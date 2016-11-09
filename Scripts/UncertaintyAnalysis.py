# -*- coding: utf-8 -*-
"""
Creates Probabilistic Flood Inundation Maps Using specified error sources

The script loads all observations and gives them random errors, based on the
error sources specified. The script creates input files for every simulation,
which contain both the input observations, references to the DTM to use, as
well as the parameters for mapping. This allows the parameters to be varied
for each simulation.

Input:
    number ---
    Integer specifying the number of random simulations to perform. This should
    be devidable by nprocesses and if random grids are used, there should be
    enough grids in the folder specified by HANDloc or DTMloc
    
    nprocesses ---
    Integer specifies the amount of processes to spawn. This allows for using
    multple cores, and is generally equal to the amount of cores of the
    processor
    
    JSONloc ---
    String with file location of json file. The file contains a dictionary
    with the observations being in the list 'tags'. Each of these observations
    in the list contains a LatLon location and optionally a water depth or
    range of water depths in its enrichments.
    
    power ---
    The power parameter for the IDW equation. This can be a single integer in
    which case a single value is used, a list containing a lower and upper
    bound in which case a uniform distribution is used or a list containing
    mean, standard deviation and a lower and upper bound, in which case a
    clipped normal distribution is used.
    
    smoothing ---
    The smooting parameter for the IDW equation. This can be a singleinteger in
    which case a single value is used, a list containing a lower and upper
    bound in which case a uniform distribution is used or a list containing
    mean, standard deviation and a lower and upper bound, in which case a
    clipped normal distribution is used.
    
    DWD ---
    Value of the DWD parameter used to calculate the WL for 'depthless'
    observations. This can be a single integer in which case a single value is
    used, a list containing a lower and upper bound in which case a uniform
    distribution is used or a list containing mean, standard deviation and a
    lower and upper bound, in which case a clipped normal distribution is used.
    
    area_filtersinks ---
    The surface area (cells) under which sinks are filtered from the elevation
    model
    
    locationerrors ---
    A list containing the mean and standard deviation in errors of points
    referenced in the dataset (intersections & POIs) as well as the mean
    and standard deviation in error of the streets referenced in the
    dataset. None if not used.
    
    deptherrors ---
    A list with mean and standard deviation in water depth errors (in cm). None
    if not used.
    
    gridtype ---
    String specifying to use a single HAND and DTM ('single') or a list of
    randomly generated HAND and DTM maps ('random').
    
    HANDloc ---
    String with the file location of the HAND map (gridtype = 'single') or the
    folder with the HAND maps with random errors (gridtype = 'random', file
    name format: #_HAND.tif). File name suffix can be changed by additionally
    specifying this suffix in the generaterandomobs and generatesimulationinput
    functions.
    
    DTMloc ---
    String with the file location of the DTM map (gridtype = 'single') or the
    folder with the DTMs with random errors (gridtype = 'random', file name
    format: #_DTM.tif). File name suffix can be changed by additionally
    specifying this suffix in the generatesimulationinput function.
    
    inputfolder ---
    String specifying the folder used to store the input files for each random
    simulation. Folder is created if it doesn't exist yet. THIS VARIABLE SHOULD
    ALSO BE TAKEN OVER IN MulticoreFuntions.generatefloodextentgrids!
    
    outputfolder ---
    String specifying the folder used to store the grids resulting from the
    random generation. Folder is created if it doesn't exist yet.

Metadata:
    Version ---
    1.01 (25-08-16)
    
    Type ----
    Standalone Script
    
    Dependencies ---
    multiprocessing, gdal, numpy, json, sys, copy, os (typical modules);
    RNDModule, FTMappingTools (Custom modules)
                        
    Author ---
    Tom Brouwer
"""
from multiprocessing import Pool
import MulticoreFunctions

if __name__ == "__main__":   
    import os
    import json
    import FloodInterpolationTools

    #%% Input:
    number = 4
    nprocesses = 4
    JSONloc = r'..\Data\observations_york.json'
    streetpointsjsonloc = r'..\Data\street_points_york.json'
    power = 4
    smoothing = 30
    DWD = 50
    area_filtersinks = 1200
    locationerrors = None#[0,50,0,200]
    deptherrors = None
    gridtype = 'single'
    HANDloc = r'..\Data\20mHAND.tif'
    DTMloc = r'..\Data\20mDTM.tif'
    inputfolder = r'D:\input\\'
    outputfolder = r'D:\output\\'
    
    #%% Script
    if not os.path.isdir(inputfolder): #Create directories
        os.makedirs(inputfolder)  
    if not os.path.isdir(outputfolder):
        os.makedirs(outputfolder)
    
    with open(JSONloc,'r') as jsonfile: #Load the original observations
        observations = json.load(jsonfile)
        
    with open(streetpointsjsonloc,'r') as streetpointsfile:
        streetpointdict = json.load(streetpointsfile)
    
    if gridtype == 'single': #Load one DEM to use the projection in determining XY locations
        HAND, properties = FloodInterpolationTools.tif2array(HANDloc)
    else:
        HAND, properties = FloodInterpolationTools.tif2array(HANDloc+'0_HAND.tif')
    
    XYWDlocmatch = FloodInterpolationTools.DeriveXYWDlocmatch(observations,properties) #Derive X,Y,WD and locational reference from observations
    
    XYZlist = FloodInterpolationTools.generaterandomobs_alongstreet(XYWDlocmatch,
                                                                    HANDloc,
                                                                    number,
                                                                    streetpointdict,
                                                                    DWD = DWD,
                                                                    deptherrors = deptherrors,
                                                                    locationerrors = locationerrors,
                                                                    gridtype = gridtype,) #List containing lists with XYZ tuples for each realization
                                                        
    FloodInterpolationTools.generatesimulationinput(XYZlist,power,smoothing,
                                                    area_filtersinks,HANDloc,
                                                    DTMloc,nprocesses,
                                                    inputfolder,
                                                    gridtype = gridtype)#Writes .inp files for all realizations, including parameter settings
                                                    
    inp_list = []
    rngstart = 0
    for nr in range(nprocesses):#Creates list specifying which process creates which maps
        rngend = rngstart+(number/nprocesses)
        settings = {'outputfolder':outputfolder,'gridtype':gridtype,
                    'inputfolder':inputfolder,'range':[rngstart,rngend]}
        inp_list.append(settings)
        rngstart = rngend
    
    p = Pool(nprocesses) #Spawns the flood interpolation processes
    results = p.map(MulticoreFunctions.generatefloodextentgrids,inp_list)
    p.close()
    
    Probabilitygrid = FloodInterpolationTools.mergefloodmaps(outputfolder,number)
    
    FloodInterpolationTools.array2tif(Probabilitygrid,properties,outputfolder + 'ProbabilisticMap.tif')