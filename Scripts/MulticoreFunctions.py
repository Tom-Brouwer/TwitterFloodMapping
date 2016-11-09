# -*- coding: utf-8 -*-
"""
Functions used via Multiprocessing.pool in the uncertainty analysis

Metadata:
    Version ---
    1.2 (6-9-2016)
    
    Type ---
    Module
    
    Dependencies ---
    json, numpy, FloodInterpolationTools
    
    Author ---
    Tom Brouwer
    
"""
import json
import numpy as np
import FloodInterpolationTools


def generatefloodextentgrids(inp_dict):
    """
    Function generates a number of flood maps, using the settings and
    observations in the input files. Flood extents are written to the output-
    folder specified in settings.inp in the inputfolder. Inputfolder is set in
    this function also.
    
    Input:
        inp_dict ---
        Dictionary containing 'range': list with lowest and last number (-1) of
        flood extent grids to generate, 'inputfolder': string specifying the
        location with the .inp files, 'gridtype': string specifying if one
        'single' grid is used or a collection of grids with 'random' errors is
        used and 'outputfolder' specifying the folder where the resulting flood
        extent maps are stored.
        
    Output:
        flood extent grids ---
        Using the input settings for each run, flood extent grids are created,
        and exported to the outputfolder specified in the .inp file in the
        inputfolder
    """
    inputfolder = inp_dict['inputfolder']
    outputfolder = inp_dict['outputfolder']
    gridtype = inp_dict['gridtype']
    inp_list = inp_dict['range']
    
    if gridtype == 'single': #In case only one grid is used, it is already loaded
        with open(inputfolder + str(inp_list[0]) + '.inp','r') as inputfile:
            runinput = json.load(inputfile)
        HANDloc = str(runinput['HANDloc'])
        DTMloc = str(runinput['DTMloc'])
        HAND, properties = FloodInterpolationTools.tif2array(HANDloc)
        FloodInterpolationTools.tif2pcr(DTMloc)
        pcrDTMloc = DTMloc[:-3] + 'map'
        
    for runnr in range(inp_list[0],inp_list[1]): #Create realizations of the flood extents
        
        with open(inputfolder + str(runnr) + '.inp','r') as inputfile:#Load the input
            runinput = json.load(inputfile)
        
        observations = runinput['observations']
        power = runinput['power']
        smoothing = runinput['smoothing']
        area_filtersinks = runinput['area_filtersinks']
        
        if gridtype == 'random': #If the grid changes per run, a new one is loaded
            HANDloc = str(runinput['HANDloc'])
            DTMloc = str(runinput['DTMloc'])
            HAND, properties = FloodInterpolationTools.tif2array(HANDloc)
            FloodInterpolationTools.tif2pcr(DTMloc)
            pcrDTMloc = DTMloc[:-3] + 'map'
        
        WD_array = FloodInterpolationTools.InterpAlongFlowPaths(observations,
                                                                HAND, properties,
                                                                pcrDTMloc, power,
                                                                smoothing,
                                                                area_filtersinks)#Interpolation
        
        floodextent = np.zeros((properties['Ysize'],properties['Xsize']))#Generate binary flood extent grid
        floodextent[np.isnan(HAND)] = np.nan
        floodextent[WD_array > 0] = 1
        
        FloodInterpolationTools.array2tif(floodextent,properties,outputfolder + 'realization_' + str(runnr) + '.tif',tif_ndv = 99, dtype = 'GDT_Byte')#Write result
        
def generategridrealizations(inp_dict):
    """
    Generates a grid with random errors from the original input grids.
    
    Input:
        inp_dict ---
        dictionary with 'orgDEM': string with location of original DEM,
        'sigma_z': float with standard deviation in elevation,'corr_dist':
        integer with correlation distance, 'outputdir': String with location to
        store output and 'range':list with the start number and end number of
        the random grids to generate.
        
    Output:
        grids with random errors ---
        These grids are written to the directory specified by the 'outputdir'
        variable. File name is realization_#.tif.
    """
    orggridloc = inp_dict['orggridloc']
    sigma_z = inp_dict['sigma_z']
    T_x = T_y = inp_dict['corr_dist']
    outputdir = inp_dict['outputdir']
    exportrange = inp_dict['range']
    
    orggridarray,properties = FloodInterpolationTools.tif2array(orggridloc) #Load original grid once
    
    Xres = properties['Xres']; Yres = properties['Yres']; #Load relevant properties
    Xsize = properties['Xsize']; Ysize = properties['Ysize']
    
    delta_x = Xres
    delta_y = -Yres
    
    chuncksize = 1000
    #Determine parameters:
    if T_y != 0:
        s = np.exp(-delta_y/float(T_y))
    else:
        s = 0
        
    if T_x != 0:
        r = np.exp(-delta_x/float(T_x))
    else:
        r = 0
    
    if T_x < 9*delta_x or T_y < 9*delta_y:
        print "Warning!: Correlation distance less than 9 times the grid size!"
    
    sigma_u = np.sqrt((1-(s**2))*(1-(r**2))*(sigma_z**2))
    
    start , end = exportrange
    
    for fl in range(start,end):
        #Create an array of random errors:
        z_array = np.zeros((Ysize,Xsize),dtype=np.float32)
        
        #Initialize first column and first row:
        cols = range(1,Xsize)
        rows = range(1,Ysize)
        
        coldata = np.random.randn(1,Xsize) * sigma_u
        rowdata = np.random.randn(Ysize,1) * sigma_u
        coldata[0,0] = rowdata[0,0] = np.random.randn()*sigma_z
        
        for cl in cols:
            coldata[0,cl] += r*coldata[0,cl-1]
        for rw in rows:
            rowdata[rw,0] += s*rowdata[rw-1,0]
        
        #Now process the chuncks:
        n_chuncks = int(np.ceil(Ysize/float(chuncksize)))
        for chnr in range(n_chuncks): #Ceates the error in chuncks of rows for speed
            if 0 < chnr < (n_chuncks - 1):
                int_array = np.random.randn(chuncksize+1,Xsize)
                int_array *= sigma_u
                int_array[0,:] = z_array[(chnr*1000)-1,:]
                int_array[1:,0] = rowdata[chnr*1000:(chnr+1)*1000,0]
                for rw in range(1,chuncksize+1):
                    rwdata = int_array[rw-1:rw+1,:]    
                    for cl in cols:
                        rwdata[1,cl] += r*rwdata[1,cl-1] + s*rwdata[0,cl] - r*s*rwdata[0,cl-1]
                z_array[(chnr*1000)-1:(chnr+1)*1000,:] = int_array
            elif chnr == 0:
                int_array = np.random.randn(chuncksize,Xsize)
                int_array *= sigma_u
                int_array[0,:] = coldata
                int_array[:,0] = rowdata[chnr*1000:(chnr+1)*1000,0]
                for rw in range(1,chuncksize):
                    rwdata = int_array[rw-1:rw+1,:]    
                    for cl in cols:
                        rwdata[1,cl] += r*rwdata[1,cl-1] + s*rwdata[0,cl] - r*s*rwdata[0,cl-1]
                z_array[chnr*1000:(chnr+1)*1000,:] = int_array
            else:
                upcut = (chnr+1)*1000 - Ysize
                int_array = np.random.randn(chuncksize-upcut+1,Xsize)
                int_array *= sigma_u
                int_array[0,:] = z_array[(chnr*1000)-1,:]
                int_array[1:,0] = rowdata[chnr*1000:(chnr+1)*1000-upcut,0]
                for rw in range(1,chuncksize-upcut+1):
                    rwdata = int_array[rw-1:rw+1,:]    
                    for cl in cols:
                        rwdata[1,cl] += r*rwdata[1,cl-1] + s*rwdata[0,cl] - r*s*rwdata[0,cl-1]
                z_array[(chnr*1000)-1:(chnr+1)*1000-upcut,:] = int_array
        
        realization = orggridarray + z_array
        
        outpfilename = '%srealization_%d.tif' %(outputdir,fl)
        FloodInterpolationTools.array2tif(realization,properties,outpfilename)
        
def generateHANDmaps(inp_dict):
    """
    Generates the HAND maps
    
    Input:
        inp_dict ---
        dictionary containing 'outputdir': Directory where the files are stored
        'range': the range of files to create. 'accuthreshold': the threshold
        used to generate the HAND map
        
    Output:
        HAND maps ---
        Created from the DTMs in the folder
    """
    outputdir = inp_dict['outputdir']
    outputrange = inp_dict['range']
    accuthreshold = inp_dict['accuthreshold']
    
    start, end = outputrange
    
    for handnr in range(start,end):
        dtmfileloc = '%slargeDTM_%d.tif' %(outputdir,handnr)
        nsdtmfileloc = '%slargensDTM_%d.tif' %(outputdir,handnr)
        handfileloc = '%slargeHAND_%d.tif' %(outputdir,handnr)
        FloodInterpolationTools.createHAND(dtmfileloc,nsdtmfileloc,
                                           accuthreshold,handfileloc)
        