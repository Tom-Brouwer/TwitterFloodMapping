# -*- coding: utf-8 -*-
"""
File for creating a map displaying the different classess in the confusion
matrix based on a flood map used as input.

The function uses the tif file of a flood map as wel as the tif file of
validation data to determine which cells in the map belong to which section of
the confusion matrix.

Input:
    fmloc ---
    String with file location of the flood map containing water depths (.tif
    format). Contains water depth values for flooded areas and nan values for
    non-flooded and nodata areas.
    
    validmaploc ---
    String with file location of the validation data which contains ones for
    flooded areas, zeros for non-flooded areas and nan values for Nodata cells
    outside the study area. The location of these cells also determines the
    outline of the result
    
    confmatrixsaveloc ---
    String with file location to save the .tif file of the confusion matrix map
    
Output:
    Confusion matrix map ---
    This map is save to the location in confmatrixsaveloc. 1 means flooded in
    both arrays, 2 means flooded in WDarray but not in valarray
    (overestimation), 3 means flooded in valarray but not in WDarray
    (underestimation) and 4 means flooded in neither of the arrays.
    
Metadata:
    Version ---
    1.0 (12-08-16 09:53)
    
    Type ---
    Standalone script
    
    Dependencies ---
    FloodInterpolationTools (custom module)
    
    Author ---
    Tom Brouwer
"""
import FloodInterpolationTools

#%% Input:
fmloc = r'..\Output\floodmap.tif'
validmaploc = r'..\Data\york_validation_innercity.tif'
confmatrixsaveloc = r'floodmap_ConfMatrix.tif'

#%% Script:
WDarray, properties = FloodInterpolationTools.tif2array(fmloc) #Load Floodmap & validation
validarray, properties = FloodInterpolationTools.tif2array(validmaploc)

confmaparray = FloodInterpolationTools.confusionmatrixmap(WDarray, validarray) #Create map

FloodInterpolationTools.array2tif(confmaparray, properties, confmatrixsaveloc,
                                  tif_ndv = 99,dtype = 'GDT_Byte') #Save the map