# Forel-Ule calculator
Calculates the Forel-Ule Index from surface reflectances
Works with both point data (observations) and images (satellite). 

Notebooks:

* Test_FU - compares this implementation with the published formulae in Woerd and Wernand (2015)
* Test_FU_images - Applies FU to test images
* forel_ule_tutorial.ipynb - Step by step calculation, with formulae and references.
* prep_test_images.ipynb - Utility to generate test images included.

Module fume
```
    calc_ForelUle_image(wavelength, reflectance, 
                        sensorcorr=None, 
                        cmf='data/FUI_CIE1931.tsv',
                        fucalibration='data/hue_angle_limits_WW2010.csv')
                        
    calc_fu_WW2015(wavelength, reflec,sensor)
        # calculates Forel Ule following Woerd and Wernand 2015 supplementary material. For point data only.
        
    forelulecmap()
        # Generates a colormap for plotting Forel Ule classes
        
```

Install in develop mode with:
    git clone https://github.com/jobel-openscience/FUME.git
    cd FUME
    pip install -e . 

To use the module:
```
import fume
```

To use FUME WHW2013

```
fume.calc_ForelUle_image(wavelength, 
                         satRrs,
                         sensorcorr=None,
                         fucalibration='data/hue_angle_limits_WW2010_decimalFU.csv')
```

To use FUME WW2015
```
fume.calc_ForelUle_image(wavelength, 
                         satRrs,
                         sensorcorr=sensor,
                         fucalibration='data/hue_angle_limits_NWW2013.csv')
```

## References

WW2010

Wernand, M. R., and van der Woerd, H. J. 2010. Spectral analysis of the Forel-Ule ocean colour comparator scale. Journal of the European Optical Society, 5.

WHW2013

M.R. Wernand, A. Hommersom, and H. J. van der Woerd, 2013. MERIS-based ocean colour classification with the discrete Forel–Ule scale. Ocean Science. doi:10.5194/os-9-477-2013

NWW2013

Novoa, S., Wernand, M. R., and Woerd, H. J. van der. 2013. The Forel-Ule scale revisited spectrally: preparation protocol, transmission measurements and chromaticity. Journal of the European Optical Society: Rapid Publications, 8: 13057. https://www.jeos.org/index.php/jeos_rp/article/view/13057.

WW2015

Hendrik J. van der Woerd and Marcel R. Wernand 2015. True Colour Classification of Natural Waters with Medium-Spectral Resolution Satellites: SeaWiFS, MODIS, MERIS and OLCI. Sensors, 15, 25663-25680.
