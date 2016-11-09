"""
Generates the grids with random errors for the York case study

Generates both the HAND and DTM at 20m resolution. Uses the same functions as
CreateDTMHANDYork.py to ensure consistent results. Since subdirectories are
used to store temporary results, the file locations should NOT be specified
relative to the current file.

Input:
    number ---
    Integer with number of grids to generate (should be a multiple of
    batchsize)
    
    batchsize ---
    Integer with the size of batches to process seperately. Used since storing
    all realizations of the 2 m DTM takes a lot of disk space
    
    nprocesses ---
    Integer with the number of processes to use. Depends on number of
    cores/threads and the available system memory.
    
    sigma_z ---
    Float with the standard deviation in ground elevation
    
    corr_dist ---
    Integer giving the correlation distance of error in metres
    
    accuthreshold ---
    Integer with threshold of upstream cells used to determine the location of
    drainage channels in producing the HAND map.
    
    orgDEMloc ---
    String with file location of the original DEM
    
    shpfileloc ---
    To ensure drainage channels are correctly represented at the side of a map,
    a larger DTM is used to produce the HAND map. This is later clipped to the
    extent of the study area. This is a string with the file location of the
    shapefile to use (same projection as DTM).
    
    outputdir ---
    String with the directory where the final DTMs and HAND maps are stored.
    
Output:
    HAND maps and DTMs ---
    These are stored to the directory indicated by outputdir
    
Metadata:
    Version ---
    1.0 (6-9-2016)
    
    Type ---
    Script
    
    Dependencies ---
    MulticoreFunctions, multiprocessing, os, sys, shutil,
    FloodInterpolationTools, SAGA GIS: saga_cmd (2.3.1 used for testing) in
    PATH system variable.
    
    Author:
    Tom Brouwer    
"""

import MulticoreFunctions
from multiprocessing import Pool

if __name__ == "__main__":   
    import os
    import sys
    import shutil
    import FloodInterpolationTools
    #%% Input:
    number = 4
    batchsize = 4
    nprocesses = 4
    sigma_z = 0.20
    corr_dist = 100
    accuthreshold = 15000
    orgDEMloc = r'2mDTM.tif'
    shpfileloc = r'..Data\york_innercity_buffer.shp'
    outputdir = r'D:\test\\'
    
    #%% Script
    tempdir = '%stemp\\' %outputdir
    batchtempdir = '%sbatchtemp\\'%tempdir
    prjfileloc = '%s2mDTM.prj' %batchtempdir
    
    if not os.path.isdir(tempdir):#Create directory
        os.makedirs(tempdir)

    if number%batchsize != 0 or batchsize%nprocesses != 0:
        sys.exit("totalgrids, batchsize or nprocessess do not align, should be\
                multiples!")
                
    for nr in range(nprocesses):#Creates copies to be used by processes
        demloc = "%sorgDEMcopy%d.tif" %(tempdir,nr)
        shutil.copy(orgDEMloc,demloc)
    
    n_batches = number / batchsize
    for batchnr in range(n_batches):
        if not os.path.isdir(batchtempdir):#Create directory
            os.makedirs(batchtempdir)
        print "Processing batch %d of %d" %(batchnr+1,n_batches)        
        batchstart = batchnr * batchsize
        batchend = (batchnr+1) * batchsize
        #First create pool workers input:
        inp_list = []
        rngstart = batchstart
        for nr in range(nprocesses):
            rngend = rngstart+(batchsize/nprocesses)
            demloc = "%sorgDEMcopy%d.tif" %(tempdir,nr)
            settings = {'orggridloc':demloc,
                        'sigma_z':sigma_z,
                        'corr_dist':corr_dist,
                        'outputdir':batchtempdir,
                        'accuthreshold':accuthreshold,
                        'range':[rngstart,rngend],
                        }
            inp_list.append(settings)
            rngstart = rngend
        
        #Now spawn the processess:
        p = Pool(nprocesses)
        results = p.map(MulticoreFunctions.generategridrealizations,inp_list)
        p.close()
        
        #After creation of grids, resample to 20m and create DTM without sinks:
        for gridnr in range(batchstart,batchend):
            print "Creating DTM %d of %d" %(gridnr+1,number)            
            inputfileloc = '%srealization_%d.tif' %(batchtempdir,gridnr)
            dtmfileloc = '%slargeDTM_%d.tif' %(batchtempdir,gridnr)
            nsdtmfileloc = '%slargensDTM_%d.tif' %(batchtempdir,gridnr)
            
            FloodInterpolationTools.resampleDTMyork(inputfileloc,dtmfileloc,
                                                    nsdtmfileloc,batchtempdir)
            
        print "Creating HAND for batch %d of %d" %(batchnr+1,n_batches)
        #Determine the HAND:
        p = Pool(nprocesses)
        results = p.map(MulticoreFunctions.generateHANDmaps,inp_list)
        p.close()
        
        #Clip the files:
        for gridnr in range(batchstart,batchend):
            print "Clipping Grid %d of %d" %(gridnr+1,number)            
            largedtmfileloc = '%slargeDTM_%d.tif' %(batchtempdir,gridnr)
            largehandfileloc = '%slargeHAND_%d.tif' %(batchtempdir,gridnr)
            dtmfileloc = '%sDTM_%d.tif' %(outputdir,gridnr)
            handfileloc = '%sHAND_%d.tif' %(outputdir,gridnr)
            
            FloodInterpolationTools.clipDTMHANDyork(largedtmfileloc,
                                                    largehandfileloc,
                                                    dtmfileloc,handfileloc,
                                                    shpfileloc,prjfileloc,
                                                    batchtempdir)
        shutil.rmtree(batchtempdir)
    shutil.rmtree(tempdir)