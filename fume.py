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

sensor_corr_file = 'data/sensor_hue_corr_WW2015'


def _polynomial(coefs,x):
    order = len(coefs)-1
    return sum( [coef*x**(order-i) for i,coef in enumerate(coefs) ])


def calc_ForelUle_image(wavelength, reflectance, sensorcorr=None, cmf='data/FUI_CIE1931_JV.tsv'):
    """ calculate Forel Ule index for multispectral image
    Arguments:
        wavelength: vector with wavlengths in nm of reflectances first dimension
        reflectance: 3D array with dimensions (wavelength, space_1, space_2) where space_* dimensions can be any aspace dimensions 
        sensor: to do hue correction for (nocorr, olci, modis, seawifs) 
    Returns:
        Array with Forel Ule class between 1-21, with dimensions (space_1, space_2)
    """
    
    if sensorcorr:
        sensor_corr_file = 'data/sensor_hue_corr_WW2015'
        sensorcorrdf = pd.read_csv(sensor_corr_file,index_col=0,)
        sensor_coef = sensorcorrdf.loc[sensorcorr].values
    
    cmf = pd.read_csv(sep = "\t", filepath_or_buffer = cmf)
    fui = pd.read_csv(sep = "\t", filepath_or_buffer = "data/FUI_ATAN210.tsv", names = ["value", "atan"])
    
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
    fu_i[ a_i >= fui["atan"][0]] = 1   #FUI = 1 its > Average
    fu_i[ np.isnan(a_i) ] = 0           # FUI = NAN = 0
    fu_i[ a_i <= fui["atan"].iloc[-1]] = 21  #FUI = 1 its > Average
    
    # find closest FU index from angle
    for c in range(0, len(fui)-1):
        fu_i[ (a_i <= fui["atan"].iloc[c]) & (a_i > fui["atan"].iloc[c+1]) ] = fui['value'].iloc[c]            
    
    return fu_i
