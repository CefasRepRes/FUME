"""
Functions for calculating Forel Ule index of water from spectral measurements

calc_ForelUle_image - calculate Forel Ule index for multispectral image
"""


import os
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import math
import matplotlib.pyplot as plt
from matplotlib import cm
import xarray as xr
from matplotlib.colors import LinearSegmentedColormap

# Used to locate package data files
def package_path(*paths, package_directory=os.path.dirname(os.path.abspath(__file__))):
    return os.path.join(package_directory, *paths)

sensor_corr_file = package_path('data','sensor_hue_corr_WW2015.csv')

def _polynomial(coefs,x):
    order = len(coefs)-1
    return sum( [coef*x**(order-i) for i,coef in enumerate(coefs) ])


def calc_ForelUle_image(wavelength, 
                        reflectance, 
                        sensorcorr = None, 
                        cmf = 'data/FUI_CIE1931.tsv',
                        fucalibration = 'data/hue_angle_limits_WW2010.csv',
                        ):
   
    """ calculate Forel Ule index for multispectral image
    Arguments:
        wavelength:    Vector with wavlengths in nm of reflectances first dimension
        reflectance:   3D array with dimensions (wavelength, space_1, space_2) 
                       where space_* dimensions can be any aspace dimensions 
        sensor:        To do hue correction for (nocorr, olci, modis, seawifs)
        cmf:           file with tristimulus weights for a specific Color 
                       Matching Function. Default is CIE1931.
        fucalibration: File with hue angle limits for each Forel Ule class.
    Returns:
                       Array with Forel Ule class between 1-21, with 
                       dimensions (space_1, space_2)
                       Array with hue angle, with dimensions (space_1, space_2)
    """
    
    if sensorcorr:
        sensorcorrdf = pd.read_csv(sensor_corr_file,index_col=0,)
        sensor_coef = sensorcorrdf.loc[sensorcorr].values
   
    cmf = package_path(cmf)
    cmf = pd.read_csv(sep = "\t", filepath_or_buffer = cmf)

    fucalibration = package_path(fucalibration)
    fui = pd.read_csv(sep = ",", filepath_or_buffer = fucalibration, header=0)
    
    Delta = cmf['wavelength'][1] - cmf['wavelength'][0]
    
    # find overping wavlengths between cmf and multispectral image 
    start = max([cmf['wavelength'].iloc[0], wavelength.min()])
    end = min([cmf['wavelength'].iloc[-1], wavelength.max()])
    #print('wavelength overlap', start,end)
    cmfi = cmf[ (cmf['wavelength'] >= start) & (cmf['wavelength'] <= end) ]
        
    cmfi = {"wavelength": cmfi['wavelength'],
            "x": np.reshape( cmfi['x'].values, [len(cmfi['x']),1,1] ),
            "y": np.reshape( cmfi['y'].values, [len(cmfi['y']),1,1] ),
            "z": np.reshape( cmfi['z'].values, [len(cmfi['z']),1,1] )}
 
    # Filter out unphysical reflectances from Atmospheric Correction
    reflectance = reflectance.where( reflectance>0 )

    # Linearly interpolate sat reflectance at the cmf table's wavelengths
    int_r1 = interp1d( wavelength, reflectance, axis = 0, kind = 'linear')(cmfi['wavelength'].values)
    
    # Creating a dictionary to hold the reflectance values as tristimulus
    r = {"x": [], "y": [], "z": []}
    
    r["x"] = int_r1 * cmfi["x"]
    r["y"] = int_r1 * cmfi["y"]
    r["z"] = int_r1 * cmfi["z"]
    
    # ------ Sum
    s = {"x": [], "y": [], "z": []}

    s["x"] = sum(r["x"] * Delta)
    s["y"] = sum(r["y"] * Delta)
    s["z"] = sum(r["z"] * Delta)
        
    sum_xyz = s["x"] + s["y"] + s["z"]
    
    # ------ chromaticity
    chrom = {"x": [], "y": [], "z": []}
    chrom["x"] = sum(r["x"] * Delta) / sum_xyz
    chrom["y"] = sum(r["y"] * Delta) / sum_xyz
    chrom["z"] = sum(r["z"] * Delta) / sum_xyz

    sum_chrom_y = chrom["x"] + chrom["y"] + chrom["z"]

    # ------ chromaticity - whiteness
    chrom_w = {"x": [], "y": []}
    chrom_w["x"] = chrom["x"] - (1 / 3)
    chrom_w["y"] = chrom["y"] - (1 / 3)

    # ------ calculate atan2
    # we use the average atan per scale
    a_i = np.arctan2( chrom_w["y"], chrom_w["x"]) * 180 / math.pi

    a_i[a_i < 0] =  a_i[a_i < 0] + 360
    

    if sensorcorr:
        # Intrument correction for Hue Angle (Woerd and Wernand 2015)
        # polynomical correction works on Hue angle (deg)/100

        sensorcorrdf = pd.read_csv(sensor_corr_file,index_col=0,)
        sensor_coef = sensorcorrdf.loc[sensorcorr].values

        anglecorr = _polynomial(sensor_coef, a_i/100)
#        correctionOLCI = (-12.5076*pow(a_i100,5) + 
#                        91.6345*pow(a_i100,4) - 
#                        249.8480*pow(a_i100,3) + 
#                        308.6561*pow(a_i100,2) - 
#                        165.4818*a_i100 + 28.5608 )

        a_i = a_i + anglecorr
    
    # ----- fui approximation
    fu_i = np.zeros(a_i.shape)
    fu_i[ a_i >= fui["lowerlimit"][0]] = 1   #FUI = 1 its > Average
    fu_i[ np.isnan(a_i) ] = float('nan')           # FUI = NAN
    fu_i[ a_i <= fui["lowerlimit"].iloc[-1]] = 21  #FUI = 1 its > Average
    
    # find closest FU index from angle
    for c in range(1, len(fui)-1):
        fu_i[ (a_i < fui["lowerlimit"].iloc[c-1]) & (a_i >= fui["lowerlimit"].iloc[c]) ] = fui['FU'].iloc[c]            
        
    return (fu_i,a_i)


def calc_fu_WW2015(wavelength, reflec,sensor):
    # calculates Forel Ule following Woerd and Wernand 2015 supplementary material
    #
    # wavelength - 1D array with sensor wavelength in  nm
    # reflec - 1D array with measuredremote sensing reflectance (Rrs) in sr-1
    # sensor - Name of sensor used: 'olci','meris','modis','seawifs'

    nbands = {'olci':11,'meris':9,'modis':7,'seawifs':6}
    if ( len(wavelength) != nbands[sensor] ):
        raise Exception(f'ERROR: Incorrect number of bands for sensor {sensor}. Found {len(wavelength)}, expected {nbands[sensor]}.')
        
    if sensor=='olci':
        X3 = (0.154*reflec[0] + 2.957*reflec[1] + 10.861*reflec[2] + 3.744*reflec[3] + 
              3.750*reflec[4] + 34.687*reflec[5] + 41.853*reflec[6] + 7.323*reflec[7] + 
              0.591*reflec[8] + 0.549*reflec[9] + 0.189*reflec[10])
        Y3 = (0.004*reflec[0] + 0.112*reflec[1] + 1.711*reflec[2] + 5.672*reflec[3] + 
              23.263*reflec[4] + 48.791*reflec[5] + 23.949*reflec[6] + 2.836*reflec[7] + 
              0.216*reflec[8] + 0.199*reflec[9] + 0.068*reflec[10])
        Z3 = (0.731*reflec[0] + 14.354*reflec[1] + 58.356*reflec[2] + 28.227*reflec[3] + 
              4.022*reflec[4] + 0.618*reflec[5] + 0.026*reflec[6] + 0.000*reflec[7] + 
              0.000*reflec[8] + 0.000*reflec[9] + 0.000*reflec[10])
    elif sensor=='meris':
        X3 = (2.957*reflec[0] + 10.861*reflec[1] + 3.744*reflec[2] + 3.750*reflec[3] + 
              34.687*reflec[4] + 41.853*reflec[5] + 7.619*reflec[6] + 0.844*reflec[7] + 
              0.189*reflec[8])
        Y3 = (0.112*reflec[0] + 1.711*reflec[1] + 5.672*reflec[2] + 23.263*reflec[3] + 
              48.791*reflec[4] + 23.949*reflec[5] + 2.944*reflec[6] + 0.307*reflec[7] + 
              0.068*reflec[8])
        Z3 = (14.354*reflec[0] + 58.356*reflec[1] + 28.227*reflec[2] + 4.022*reflec[3] + 
              0.618*reflec[4] + 0.026*reflec[5] + 0.000*reflec[6] + 0.000*reflec[7] + 
              0.000*reflec[8])
    elif sensor == 'modis':
        X3 = (2.957*reflec[0] + 10.861*reflec[1] + 4.031*reflec[2] + 3.989*reflec[3] + 
              49.037*reflec[4] + 34.586*reflec[5] + 0.829*reflec[6])
        Y3 = (0.112*reflec[0] + 1.711*reflec[1] + 11.106*reflec[2] + 22.579*reflec[3] + 
              51.477*reflec[4] + 19.452*reflec[5] + 0.301*reflec[6])
        Z3 = (14.354*reflec[0] + 58.356*reflec[1] + 29.993*reflec[2] + 2.618*reflec[3] + 
              0.262*reflec[4] + 0.000*reflec[5] + 0.000*reflec[6])
    elif sensor == 'seawifs':
        X3 = (2.957*reflec[0] + 10.861*reflec[1] + 3.744*reflec[2] + 3.455*reflec[3] +
              52.304*reflec[4] + 32.825*reflec[5])
        Y3 = (0.112*reflec[0] + 1.711*reflec[1] + 5.672*reflec[2] + 21.929*reflec[3] + 
              59.454*reflec[4] + 17.810*reflec[5])
        Z3 = (14.354*reflec[0] + 58.356*reflec[1] + 28.227*reflec[2] + 3.967*reflec[3] + 
              0.682*reflec[4] + 0.018*reflec[5])

        
    Chrx = X3/(X3 + Y3 + Z3)
    Chry = Y3/(X3 + Y3 + Z3)

    hueangl = ((math.atan2((Chry - 0.333333), (Chrx - 0.333333)))*180/math.pi) # degrees
    if hueangl < 0:
        hueangl = hueangl+360

    hueangl_100 = hueangl/100

    # Calc sensor specific polynomial correction
    if sensor == 'olci':
        Polyhueangl = (-12.5076*pow(hueangl_100,5) + 
                        91.6345*pow(hueangl_100,4) - 
                        249.8480*pow(hueangl_100,3) + 
                        308.6561*pow(hueangl_100,2) - 
                        165.4818*hueangl_100 + 28.5608 )
    elif sensor == 'meris':
        Polyhueangl = (-12.0506*pow(hueangl_100,5) + 
                       88.9325*pow(hueangl_100,4) - 
                       244.6960*pow(hueangl_100,3) + 
                       305.2361*pow(hueangl_100,2) - 
                       164.6960* hueangl_100 + 28.5255)
    elif sensor == 'modis':
        Polyhueangl = (-48.0880*pow(hueangl_100,5) + 
                       362.6179*pow(hueangl_100,4) -
                       1011.7151*pow(hueangl_100,3) + 
                       1262.0348*pow(hueangl_100,2) - 
                       666.5981*hueangl_100 + 113.9215)
    elif sensor == 'seawifs':
        Polyhueangl = (-49.4377*pow(hueangl_100,5) + 
                       363.2770*pow(hueangl_100,4) -
                       978.1648*pow(hueangl_100,3) + 
                       1154.6030*pow(hueangl_100,2) - 
                       552.2701*hueangl_100 + 78.2940)
                       
    hueanglPcorr = hueangl + Polyhueangl

    if (hueanglPcorr > 232):
        FUSentinel3Pcorr = 0
    elif hueanglPcorr > 227.168:
        FUSentinel3Pcorr = 1
    elif (hueanglPcorr > 220.977):
        FUSentinel3Pcorr = 2
    elif (hueanglPcorr > 209.994):
        FUSentinel3Pcorr = 3
    elif (hueanglPcorr > 190.779):
        FUSentinel3Pcorr = 4
    elif (hueanglPcorr > 163.084):
        FUSentinel3Pcorr = 5
    elif (hueanglPcorr > 132.999):
        FUSentinel3Pcorr = 6
    elif (hueanglPcorr > 109.054):
        FUSentinel3Pcorr = 7
    elif (hueanglPcorr > 94.037):
        FUSentinel3Pcorr = 8
    elif (hueanglPcorr > 83.346):
        FUSentinel3Pcorr = 9
    elif (hueanglPcorr > 74.572):
        FUSentinel3Pcorr = 10
    elif (hueanglPcorr > 67.957):
        FUSentinel3Pcorr =  11
    elif (hueanglPcorr > 62.186):
        FUSentinel3Pcorr =  12
    elif (hueanglPcorr > 56.435):
        FUSentinel3Pcorr = 13
    elif (hueanglPcorr > 50.665):
        FUSentinel3Pcorr =  14
    elif (hueanglPcorr > 45.129):
        FUSentinel3Pcorr =  15
    elif (hueanglPcorr > 39.769):
        FUSentinel3Pcorr =  16
    elif (hueanglPcorr > 34.906):
        FUSentinel3Pcorr =  17
    elif (hueanglPcorr > 30.439):
        FUSentinel3Pcorr =  18
    elif (hueanglPcorr > 26.337):
        FUSentinel3Pcorr =  19
    elif (hueanglPcorr > 22.741):
        FUSentinel3Pcorr =  20
    elif (hueanglPcorr > 19):
        FUSentinel3Pcorr =  21
    elif (hueanglPcorr < 19):
        FUSentinel3Pcorr =  21
    
    return (FUSentinel3Pcorr, hueanglPcorr)


def forelulecmap():
    # Returns Forel Ule colormap
    
    fuh ={"1":"#2158bc",
          "2":"#316dc5",
          "3":"#327cbb",
          "4":"#4b80a0",
          "5":"#568f96",
          "6":"#6d9298",
          "7":"#698c86",
          "8":"#759e72",
          "9":"#7ba654",
          "10":"#7dae38",
          "11":"#95b645",
          "12":"#94b660",
          "13":"#a5bc76",
          "14":"#aab86d",
          "15":"#adb55f",
          "16":"#a8a965",
          "17":"#ae9f5c",
          "18":"#b3a053",
          "19":"#af8a44",
          "20":"#a46905",
          "21":"#a14d04"}

    def _hex_to_rgb(value):
        '''
        Converts hex to rgb colours
        value: string of 6 characters representing a hex colour.
        Returns: list length 3 of RGB values'''
        value = value.strip("#") # removes hash symbol if present
        lv = len(value)
        return tuple(int(value[i:i + lv // 3], 16)/256 for i in range(0, lv, lv // 3))

    furgb = [_hex_to_rgb(hexstr) for hexstr in fuh.values()]

    cm = LinearSegmentedColormap.from_list('ForelUle',furgb,N=21)
    cm.set_under(color='k', alpha=None)
    cm.set_over(color='k', alpha=None)
    cm.set_bad(color='k', alpha=0)
    
    return cm
