"""
Generates the DTM and HAND for the York case study

This function is used to generate the 20m resolution DTM and HAND map for the
York case study

Input:
    orgDEMfileloc ---
    String with the file location of the original high resolution DTM
    
    accuthreshold ---
    The threshold used for accumulated flow, to identify drainage channels in
    the DTM. Used for creating the HAND map
    
    shpfileloc ---
    String with the file location of the shapefile used to clip the DTM and
    HAND maps
    
    savedir ---
    Directory to save the 20m DTM and 20m HAND
    
Output:
    20mDTM.tif ---
    The resampled 20m DTM, saved to the savedir
    
    20mHAND.tif ---
    The HAND at 20m resolution, saved to the savedir
    
Metadata:
    Version ---
    1.0 (5-9-2016)
    
    Type ---
    Script
    
    Dependencies ---
    os,shutil,FloodInterpolationTools, saga_cmd of SAGA GIS (tested with 2.3.1)
    in path variable.
    
    Author ---
    Tom Brouwer
"""
import os
import shutil
import FloodInterpolationTools

#%% Input:
orgDEMfileloc = r'2mDTM.tif'
accuthreshold = 15000
shpfileloc = r'..\Data\YorkInnerCity.shp'
savedir = r'D:\test\\'

#%% Script:
DTMfileloc = '%s20mDTM.tif' %savedir
HANDfileloc = '%s20mHAND.tif' %savedir
tempdir = '%ssagatemp\\' %savedir
largeDTMloc = '%slrgDTM.tif' %tempdir
largensDTMloc = '%slrgnsDTM.tif' %tempdir
largeHANDloc = '%slrgHAND.tif' %tempdir
prjfileloc = '%s2mDTM.prj' %tempdir #Created as temp file by resample script

if not os.path.isdir(tempdir):#Create directory
    os.makedirs(tempdir)

FloodInterpolationTools.resampleDTMyork(orgDEMfileloc,largeDTMloc,
                                        largensDTMloc, tempdir) #Resamples DTM to 20m res

FloodInterpolationTools.createHAND(largeDTMloc,largensDTMloc,accuthreshold,
                                   largeHANDloc) #Create HAND using resampled DTM
                                   
FloodInterpolationTools.clipDTMHANDyork(largeDTMloc,largeHANDloc,DTMfileloc,
                                        HANDfileloc,shpfileloc,prjfileloc,
                                        tempdir) #Clips the files and writes them
                                        
shutil.rmtree(tempdir)