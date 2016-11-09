# -*- coding: utf-8 -*-
"""
Create Flood Maps by interpolating along flood paths (Simplified)

This is only a wrapper for loading the observations, grids and saving the
results. The actual flood mapping script is included in the module
FloodInterpolationTools.py.

Input:
    JSONloc --- 
    string with the location of the JSON file with the observations, in format
    of the output of the FloodTags API.
                    
    HANDloc ---
    string with file location of the HAND map, used in flood mapping
    
    DTMloc ---
    Location of the DTM to use in grouping (In the HAND map sinks are not
    accurately represented, meaning that the filtering of sinks leads to
    bad results).
    
    power ---
    Integer specifying the value of the power parameter used for interpolation.
    
    smoothing ---
    Integer specifying the value of the smoothing parameter (number of cells)
    in interpolation.
    
    DWD ---
    Default water depth used for observations without water depths. If None,
    only observations with water depths are used.
    
    area_filtersinks ---
    Area of sinks that are filtered by the LDD creation process, used in
    grouping
    
    tifsaveloc ---
    string with location to save tif file. None if not used.

Output:
    Tif with Water depths ---
    Tif with water depths is exported in case tifsaveloc has a value.
                                        
    Figures ---
    Figures with Water depths, and water levels along the flow paths are
    created by default.

Metadata:
    Version ---
    1.21 (09-08-16 14:48)
    
    Type ---
    Standalone script
    
    Dependencies ---
    json, pcraster, numpy, matplotlib, gdal, copy, time (typical modules);
    FTMappingTools (custom module)
    
    Author ---
    Tom Brouwer
"""
import json
import matplotlib.pyplot as plt
import FloodInterpolationTools

#%% Input:
JSONloc = r'..\Data\observations_york.json'
HANDloc = r'..\Data\20mHAND.tif'
DTMloc = r'..\Data\20mDTM.tif'
power = 4
smoothing = 30
DWD = 50
area_filtersinks = 1200
tifsaveloc = None

#%% Script
with open(JSONloc,'r') as jsonfile: #Load the observations
    observations = json.load(jsonfile)

HAND, properties = FloodInterpolationTools.tif2array(HANDloc) #Load the HAND

TLx = properties['TLx']; TLy = properties['TLy']; Xres = properties['Xres'] #Derive HAND properties
Yres = properties['Yres']; Xsize = properties['Xsize']
Ysize = properties['Ysize']; extent = [TLx,TLx+Xsize*Xres,TLy+Ysize*Yres,TLy]

FloodInterpolationTools.tif2pcr(DTMloc) #Make DTM pcraster file, to use for grouping
pcrDTMloc = DTMloc[:-3] + 'map'

WDarray = FloodInterpolationTools.InterpAlongFlowPaths(observations, HAND,
                                                       properties, pcrDTMloc,
                                                       power, smoothing,
                                                       area_filtersinks,
                                                       DWD = DWD) #Do the interpolation

plt.close(1) #Plot the results
plt.figure(1)
plt.imshow(WDarray,interpolation = 'nearest',cmap = 'Blues',extent = extent)
plt.gca().set_axis_bgcolor('black')
plt.title('Water depth (m)')
plt.colorbar()
plt.show()

if tifsaveloc: #Save the results
    FloodInterpolationTools.array2tif(WDarray,properties,tifsaveloc)