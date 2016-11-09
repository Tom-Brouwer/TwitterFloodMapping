"""
Flood Interpolation Toolset

Contains a variety of functions used in creating the flood maps

Metadata:
    Version ---
    2.00 (6-9-2016)
    
    Type ---
    Module
    
    Dependencies ---
    sys, re, subprocess, shutil, numpy, pcraster, osgeo, json, os, saga_cmd (
    tested with 2.3.1) in PATH variable for creating DTMs.
    
    Author---
    Tom Brouwer
"""
from osgeo import gdal
from osgeo import osr
from osgeo import ogr
import sys
import re
import subprocess
import shutil
import numpy as np
import pcraster as pcr
import json
import os

def tif2array(tifloc):
    """
    Loads a tif file as an array
    
    Input:
        tifloc ---
        String with location of tif file
        
    Returns:
        array ---
        Numpy array containing the data of the Tif file. Nodata cells are nan.
        All arrays returned are dtype np.float64
        
        properties ---
        dictionary containing the properties of the tif: 'TLx','TLy','Xres',
        'Yres','Xsize','Ysize','Nodata','EPSG','projectionwkt'.
    """
    ds = gdal.Open(tifloc,0)
    Geotransform = ds.GetGeoTransform()
    projectionwkt = ds.GetProjection()
    Projectionsplit = re.split('\[|\]|\"|,',projectionwkt)
    Projectionsplit = [i for i in Projectionsplit if i != ""]
    for i in reversed(range(len(Projectionsplit))):
        if Projectionsplit[i] == "EPSG":
            EPSG = int(Projectionsplit[i+1])
            break
    else:
        EPSG = 4326
        print "WARNING: No EPSG code found in projection data, assuming WGS84 for input DEM"

    properties = {'TLx':Geotransform[0],'TLy':Geotransform[3],'Xres':Geotransform[1],
                  'Yres':Geotransform[5],'Xsize':ds.RasterXSize,'Ysize':ds.RasterYSize,
                  'EPSG':EPSG,'projectionwkt':projectionwkt}
    nodata = ds.GetRasterBand(1).GetNoDataValue()
    array = ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
    array[array == nodata] = np.nan
    ds = None
    
    return array,properties
    
    
def tif2pcr(Tiflocation, Pcrlocation = None, checkexist = True):
    """
    Transforms the input .tif raster to a PCRaster .map file. If Pcrlocation = None
    the same folder as the tif is used as an output location.
    
    Input:
        Tiflocation ---
        String with file location of .tif
        
        Pcrlocation (default:None) ---
        String with output file location
        
        checkexist (default:True) ---
        Don't create if the target file already exists
        
    Output:
        PCRaster map file ---
        Written to same folder as tif file, or location specified in
        Pcrlocation
    """
    if not Pcrlocation:
        Pcrlocation = Tiflocation[:-3] + 'map'
    
    if (not os.path.isfile(Pcrlocation)) or not checkexist:
        src_ds = gdal.Open(Tiflocation)
        driver = gdal.GetDriverByName("PCRaster")
        dst_ds = driver.CreateCopy(Pcrlocation,src_ds,0)
        src_ds = dst_ds = None
        
def array2tif(array,properties,tifloc,tif_ndv = -99999,dtype = 'GDT_Float32'):
    """
    Function saves an array as a tif file
    
    Input:
        array ---
        Numpy array to write to tif. Nodata cells should be nan.
        
        properties ---
        Dictionary with properties of the original tif (see also tif2array)
        
        tifloc ---
        String with location to save the new tif file
        
        tif_ndv (default:-99999) ---
        Nodata value for the Tif
        
        data_type (default:'GDT_Float32') ---
        String specifying the gdal data type of the output. For all possible
        types see: http://www.gdal.org/gdal_8h.html#a22e22ce0a55036a96f652765793fb7a4
        
    Output:
        Tif File ---
        Tif File gets written to location specified in tifloc
        
    """
    data_type = getattr(gdal,dtype)    
    
    array[np.isnan(array)] = tif_ndv
    wrdriver = gdal.GetDriverByName("GTiff")
    outds = wrdriver.Create(tifloc,properties['Xsize'],properties['Ysize'],1,data_type)
    outds.SetGeoTransform([properties['TLx'],properties['Xres'],0,
                           properties['TLy'],0,properties['Yres']])
    outds.SetProjection(properties['projectionwkt'])
    outds.GetRasterBand(1).SetNoDataValue(tif_ndv)
    outds.GetRasterBand(1).WriteArray(array)
    stats = outds.GetRasterBand(1).ComputeStatistics(False)
    outds.GetRasterBand(1).SetStatistics(stats[0],stats[1],stats[2],stats[3])
    outds.FlushCache()
    outds = None
    
def resampleDTMyork(inputDTMfileloc,outputDTMfileloc,outputnsDTMfileloc,tempdir):
    """
    Resamples the original DTM for York to 20m resolution and creates a DTM
    without sinks.
    
    Input:
        inputDTMfileloc ---
        String with location of the input DTM
        
        outputDTMfileloc ---
        String with the location of the resampled DTM
        
        outputnsDTMfileloc ---
        String with the location of the sink filtered resampled DTM
        
        tempdir ---
        Used to store temporary files (will not be deleted automatically)
        
    Output:
        Output DTM ---
        Saved to location specified by outputDTMfileloc
        
        Ouput DTM with no sinks ---
        Save to location specified by outputnsDTMfileloc
    """        
    sagascript = """#Load the Grid:
io_gdal 0 -GRIDS=\"%s2mDTM.sgrd\" -FILES=\"%s\"

#Resample the grid to 20m:
grid_tools 0 -INPUT=\"%s2mDTM.sgrd\" -OUTPUT=\"%s20mDTMlrg.sgrd\" -TARGET_USER_SIZE=20 -TARGET_USER_XMIN=451969.00012207031 -TARGET_USER_XMAX=467149.00012207031 -TARGET_USER_YMIN=441433.00012207031 -TARGET_USER_YMAX=463693.00012207031

#Filter the sinks:
ta_preprocessor 4 -ELEV=\"%s20mDTMlrg.sgrd\" -FILLED=\"%s20mDTMlrg_NoSinks.sgrd\" -MINSLOPE=0.0

#Save the results to be used for creating the HAND:
io_gdal 2 -GRIDS=\"%s20mDTMlrg.sgrd\" -FILE=\"%s\"
io_gdal 2 -GRIDS=\"%s20mDTMlrg_NoSinks.sgrd\" -FILE=\"%s\"""" %(tempdir,inputDTMfileloc,
                                                                tempdir,tempdir,
                                                                tempdir,tempdir,
                                                                tempdir,outputDTMfileloc,
                                                                tempdir,outputnsDTMfileloc)
                                                        
    sagascript = sagascript.replace('\\\\','\\')
    
    with open('%sresamplescript.txt' %tempdir,'w') as sagafile:
        sagafile.write(sagascript)
        
    command = 'saga_cmd %sresamplescript.txt' %tempdir
    command = command.replace('\\\\','\\')           
    state = subprocess.call(command)
    if state:
        sys.exit('saga_cmd could not execute resamplescript.txt')
    
def createHAND(inputDTMfileloc,inputnsDTMfileloc, accuthreshold, outputHANDfileloc):
    """
    Uses pcraster to create the HAND of the input DTM, and reintroduces sinks
    after the process
    
    Input:
        inputDTMfileloc ---
        String with location of the input DTM
        
        inputnsDTMfileloc ---
        String with the location of the DTM file without sinks
        
        accuthreshold ---
        Threshold of accumulated flow (nr of cells) used to identify drainage
        channels
        
        outputHANDfileloc ---
        String with the file location of the HAND map (should end with .tif)
        
    Ouput:
        HAND map ---
        Saved to the location specified in outputHANDfileloc
    """
    DTMarray,properties = tif2array(inputDTMfileloc)
    tif2pcr(inputDTMfileloc)
    tif2pcr(inputnsDTMfileloc)
    
    pcrDEM=pcr.readmap(inputDTMfileloc[:-3]+'map')
    pcrDEMnosinks=pcr.readmap(inputnsDTMfileloc[:-3]+'map')
    
    ldd = pcr.lddcreate(pcrDEMnosinks,0,0,0,0)
    
    stream = pcr.ifthenelse(pcr.accuflux(ldd, 1) >= accuthreshold,
                            pcr.boolean(1), pcr.boolean(0))
    
    height_river = pcr.ifthenelse(stream, pcr.ordinal(pcrDEMnosinks*1000), 0)
    
    up_elevation = pcr.scalar(pcr.subcatchment(ldd, height_river))

    hand = pcr.max(pcr.scalar(pcr.ordinal(pcrDEMnosinks*1000))-up_elevation, 0)/1000
    
    HANDmap = hand - (pcrDEMnosinks - pcrDEM)
    
    HANDmaparray = pcr.pcr2numpy(HANDmap,float('NaN'))
    
    array2tif(HANDmaparray,properties,outputHANDfileloc)
    
def clipDTMHANDyork(inputDTMfileloc,inputHANDfileloc,outputDTMfileloc,
                    outputHANDfileloc,shpfileloc,prjfileloc,tempdir):
    """
    Clips the DTM and HAND for the York case, to the extent specified by
    in a shapefile
    
    Input:
        inputDTMfileloc ---
        String with the location of the large DTM file
        
        inputHANDfileloc ---
        String with the location of the large HAND file
        
        outputDTMfileloc ---
        String with the location to save the clipped DTM file
        
        outputHANDfileloc ---
        String with the location to save the clipped HAND file
        
        shpfileloc ---
        String the location of the shapefile used to clip the DTM and HAND.
        
        prjfileloc ---
        String with the location of the prj file. This is copied since SAGA
        does not add projection information after clipping the DTM & HAND. Also
        the .sdat.aux.xml file should have the same name and should be in the
        same directory
        
        tempdir ---
        Directory used to store the temporary files from the SAGA script
    
    Output:
        Clipped DTM ---
        tif file saved to the location specified by outputDTMfileloc
        
        Clipped HAND ---
        tif file saved to the location specified by outputHANDfileloc
    """
    sagascript = """#Load the created HAND Grid:
io_gdal 0 -GRIDS=\"%sHANDlrg.sgrd\" -FILES=\"%s\"
io_gdal 0 -GRIDS=\"%sDTMlrg.sgrd\" -FILES=\"%s\"

#Now clip all grids using YorkInnerCity.shp:
shapes_grid 7 -OUTPUT=\"%sDTM_clip.sgrd\" -INPUT=\"%sDTMlrg.sgrd\" -POLYGONS=\"%s\"
shapes_grid 7 -OUTPUT=\"%sHAND_clip.sgrd\" -INPUT=\"%sHANDlrg.sgrd\" -POLYGONS=\"%s\"""" %(tempdir,inputHANDfileloc,
                                                                                           tempdir,inputDTMfileloc,
                                                                                           tempdir,tempdir,shpfileloc,
                                                                                           tempdir,tempdir,shpfileloc)
    sagascript = sagascript.replace('\\\\','\\')
    
    with open('%sclipscript.txt' %tempdir,'w') as sagafile:
        sagafile.write(sagascript)
    
    command = "saga_cmd %sclipscript.txt" %tempdir
    command = command.replace('\\\\','\\')
    state = subprocess.call(command)
    
    if state:
        sys.exit('saga_cmd could not execute clipscript.txt')
    
    prjfnbase = prjfileloc[:-4]
    sagaHANDfnbase = '%sHAND_clip' %tempdir
    sagaDTMfnbase = '%sDTM_clip' %tempdir
    
    shutil.copy('%s.prj' %prjfnbase,'%s.prj' %sagaHANDfnbase)
    shutil.copy('%s.sdat.aux.xml' %prjfnbase,'%s.sdat.aux.xml' %sagaHANDfnbase)
    shutil.copy('%s.prj' %prjfnbase,'%s.prj' %sagaDTMfnbase)
    shutil.copy('%s.sdat.aux.xml' %prjfnbase,'%s.sdat.aux.xml' %sagaDTMfnbase)
        
    sagascript2 = """#Save the result:
io_gdal 2 -GRIDS=\"%sDTM_clip.sgrd\" -FILE=\"%s\"
io_gdal 2 -GRIDS=\"%sHAND_clip.sgrd\" -FILE=\"%s\"""" %(tempdir,outputDTMfileloc,
                                                        tempdir,outputHANDfileloc)
    sagascript2 = sagascript2.replace('\\\\','\\')
    
    with open('%sclipscript2.txt' %tempdir,'w') as sagafile2:
        sagafile2.write(sagascript2)
        
    command = "saga_cmd %sclipscript2.txt" %tempdir
    command = command.replace('\\\\','\\')
    state = subprocess.call(command)
    
    if state:
        sys.exit('saga_cmd could not execute clipscript2.txt')
    
def DeriveXY(latlons,EPSG):
    """
    Derives the XY coordinates of the list of LatLons in the EPSG projection
    number specified.
    
    Input:
        LatLons ---
        A list containing either lists with the lat lon values, or tuples with
        the lat lon values.
        
        EPSG ---
        EPSG code of the projection to project to
        
    Ouput:
        XY ---
        A list containing tuples with XY values of each LatLon in the input.
    """
    sourceproj = osr.SpatialReference()
    sourceproj.ImportFromEPSG(4326)
    targetproj = osr.SpatialReference()
    targetproj.ImportFromEPSG(EPSG)
    transformation = osr.CoordinateTransformation(sourceproj,targetproj)
    
    XY = []
    for latlon in latlons:
        wkt = "POINT (%.6f %.6f)" %(latlon[1],latlon[0])
        point = ogr.CreateGeometryFromWkt(wkt)
        point.Transform(transformation)
        X=point.GetX()
        Y=point.GetY()
        XY.append((X,Y))
    
    return XY

def DeriveXYZ(observations, DEM, properties, DWD = None, usedepth = 'avg'):
    """
    Gives a list of tuples with XYZ (Z being water level) derived from
    observations.
    
    Input:
        observations ---
        Dictionary in format downloaded from the FloodTags API, with each
        observation in observations['tags'] containing enrichments with at 
        least a LatLon location and optionally a water depth or range of water
        depths.
        
        DEM ---
        Array containing an elevation model where the water levels values are
        based upon
        
        properties ---
        properties of the DEM, loaded using the tif2array function.
        
        DWD (default:None) ---
        Specifies the default water depth, added to observations without one
        
        usedepth (default:'avg') ---
        String specifying which depth to use, in case a range is specified.
        'avg' means average,'low' means lower estimate and 'high' means higher
        estimate.
        
    Returns:
        XYZ ---
        A list containing tuples with XYZ values. The z-value is None if no
        water depth was specified.
    """
    #Load properties of the DEM:
    TLx = properties['TLx']
    TLy = properties['TLy']
    Xres = properties['Xres']
    Yres = properties['Yres']
    
    X_size = properties['Xsize']
    Y_size = properties['Ysize']
    
    EPSG = properties['EPSG']
    
    Enrichments = [enr for ob in observations['tags'] for enr in ob['Enrichments']]
    
    latlons = [enr['LatLon'] for enr in Enrichments]
    
    if EPSG != 4326:
        XY = DeriveXY(latlons,EPSG)
    else:
        XY = latlons
    
    XYZ = []
    for enr_nr in range(len(Enrichments)):
        X = XY[enr_nr][0]
        Y = XY[enr_nr][1]
        row = int(np.floor((Y-TLy)/Yres))
        col = int(np.floor((X-TLx)/Xres))
        
        if row < 0 or row > (Y_size-1) or col < 0 or col > (X_size -1):
            sys.exit("ERROR: Observation outside DEM BBOX")
        
        GL = DEM[row,col]
        
        if np.isnan(GL):
            sys.exit("ERROR: Observation in NODATA area")
        
        if Enrichments[enr_nr]['WD']:
            depth = Enrichments[enr_nr]['WD']
        elif Enrichments[enr_nr]['WDLow']:
            if usedepth == 'avg':
                depth = (Enrichments[enr_nr]['WDLow']+Enrichments[enr_nr]['WDHigh'])/2
            elif usedepth == 'low':
                depth = Enrichments[enr_nr]['WDLow']
            elif usedepth == 'high':
                depth = Enrichments[enr_nr]['WDHigh']
            else:
                sys.exit("ERROR: Incorrect value of usedepth")
        else:
            depth = DWD
        
        if depth:
            Z = GL + (float(depth)/100)
            XYZ.append((X,Y,Z))
        else:
            XYZ.append((X,Y,None))

    return XYZ
    
    
def DeriveXYWDmethod(observations,properties,usedepth = 'avg'):
    """
    Gives a list of tuples with XY ,water depth and type of locational
    reference derived from observations, and the method used to derive the
    location. Output is used in random generation.
    
    Input:
        observations ---
        Dictionary in format downloaded from the FloodTags API, with each
        observation in observations['tags'] containing enrichments with at 
        least a LatLon location and optionally a water depth or range of water
        depths
        
        properties ---
        properties from loading the elevation model, loaded using tif2array.
        
        usedepth (default:'avg') ---
        String specifying which depth to use, in case a range is specified.
        'avg' means average,'low' means lower estimate and 'high' means higher
        estimate. Default is 'avg'.
        
    Returns:
        XYWDmethod ---
        A list containing tuples with XY, WD and method values. The WD-value is
        None if no water depth was specified.
    """
    #Load properties of the DEM:
    TLx = properties['TLx']
    TLy = properties['TLy']
    Xres = properties['Xres']
    Yres = properties['Yres']
    
    X_size = properties['Xsize']
    Y_size = properties['Ysize']
    
    EPSG = properties['EPSG']
    
    Enrichments = [enr for ob in observations['tags'] for enr in ob['Enrichments']]
    
    latlons = [enr['LatLon'] for enr in Enrichments]
    
    if EPSG != 4326:
        XY = DeriveXY(latlons,EPSG)
    else:
        XY = latlons
    
    XYWDmethod = []
    for enr_nr in range(len(Enrichments)):
        X = XY[enr_nr][0]
        Y = XY[enr_nr][1]
        method = Enrichments[enr_nr]['Method']
            
        row = int(np.floor((Y-TLy)/Yres))
        col = int(np.floor((X-TLx)/Xres))
        
        if row < 0 or row > (Y_size-1) or col < 0 or col > (X_size -1):
            sys.exit("ERROR: Observation outside DEM BBOX")
        
        if Enrichments[enr_nr]['WD']:
            depth = Enrichments[enr_nr]['WD']
        elif Enrichments[enr_nr]['WDLow']:
            if usedepth == 'avg':
                depth = (Enrichments[enr_nr]['WDLow']+Enrichments[enr_nr]['WDHigh'])/2
            elif usedepth == 'low':
                depth = Enrichments[enr_nr]['WDLow']
            elif usedepth == 'high':
                depth = Enrichments[enr_nr]['WDHigh']
            else:
                sys.exit("ERROR: Incorrect value of usedepth")
        else:
            depth = None
        
        XYWDmethod.append((X,Y,depth,method))

    return XYWDmethod
    
def DeriveXYWDlocmatch(observations,properties,usedepth = 'avg'):
    """
    Gives a list of tuples with XY ,water depth and locational reference
    derived from observations, and the method used to derive the location.
    Output is used in random generation.
    
    Input:
        observations ---
        Dictionary in format downloaded from the FloodTags API, with each
        observation in observations['tags'] containing enrichments with at 
        least a LatLon location and optionally a water depth or range of water
        depths
        
        properties ---
        properties from loading the elevation model, loaded using tif2array.
        
        usedepth (default:'avg') ---
        String specifying which depth to use, in case a range is specified.
        'avg' means average,'low' means lower estimate and 'high' means higher
        estimate. Default is 'avg'.
        
    Returns:
        XYWDlocmatch ---
        A list containing tuples with XY, WD and locmatch values. The WD-value
        is None if no water depth was specified.
    """
    #Load properties of the DEM:
    TLx = properties['TLx']
    TLy = properties['TLy']
    Xres = properties['Xres']
    Yres = properties['Yres']
    
    X_size = properties['Xsize']
    Y_size = properties['Ysize']
    
    EPSG = properties['EPSG']
    
    Enrichments = [enr for ob in observations['tags'] for enr in ob['Enrichments']]
    
    latlons = [enr['LatLon'] for enr in Enrichments]
    
    if EPSG != 4326:
        XY = DeriveXY(latlons,EPSG)
    else:
        XY = latlons
    
    XYWDlocmatch = []
    for enr_nr in range(len(Enrichments)):
        X = XY[enr_nr][0]
        Y = XY[enr_nr][1]
        locmatch = str(Enrichments[enr_nr]['LocMatch'])
            
        row = int(np.floor((Y-TLy)/Yres))
        col = int(np.floor((X-TLx)/Xres))
        
        if row < 0 or row > (Y_size-1) or col < 0 or col > (X_size -1):
            sys.exit("ERROR: Observation outside DEM BBOX")
        
        if Enrichments[enr_nr]['WD']:
            depth = Enrichments[enr_nr]['WD']
        elif Enrichments[enr_nr]['WDLow']:
            if usedepth == 'avg':
                depth = (Enrichments[enr_nr]['WDLow']+Enrichments[enr_nr]['WDHigh'])/2
            elif usedepth == 'low':
                depth = Enrichments[enr_nr]['WDLow']
            elif usedepth == 'high':
                depth = Enrichments[enr_nr]['WDHigh']
            else:
                sys.exit("ERROR: Incorrect value of usedepth")
        else:
            depth = None
        
        XYWDlocmatch.append((X,Y,depth,locmatch))

    return XYWDlocmatch
    
def groupbyflowpaths(XYZ,pcrLDD,properties):
    """
    Groups observations based on flow path intersections downstream.
    
    Input:    
        XYZ ---
        List with tuples containing X,Y,Z values, created using DeriveXYZ
        
        pcrLDD ---
        PCraster LDD map
        
        properties ---
        Properties of the DTM, generated using tif2array
        
    Returns:
        XYZa ---
        List with tuples containing X,Y,Z and the number of the area the
        observation belongs to.
        
        groupmap ---
        Numpy array specifying which flow paths belong to which areas
       
    """
    Xsize = properties['Xsize']
    Ysize = properties['Ysize']
    TLy = properties['TLy']
    TLx = properties['TLx']
    Yres = properties['Yres']
    Xres = properties['Xres']
    Nodata = -99999
    
    #Transform LDD to array:
    LDDArray = pcr.pcr2numpy(pcrLDD,Nodata)
    LDDArray = LDDArray.astype(np.int32)
    LDDArray[LDDArray == 97] = 5
    
    #Create flow paths for each observation:
    flow_dirs=np.array([[7, 8, 9], [4, 5, 6], [1, 2, 3]])
    flow_dir_y = np.array([-1, 0, 1])
    flow_dir_x = np.array([-1, 0, 1])
    
    DSPath_array = np.zeros((Ysize,Xsize),dtype=np.uint8)
    for ob in XYZ:
        row = int(np.floor((ob[1]-TLy)/Yres))
        col = int(np.floor((ob[0]-TLx)/Xres))

        while LDDArray[row][col] != 5 and DSPath_array[row][col] != 1:
            DSPath_array[row][col]=1
            flow_dir = np.where(flow_dirs==LDDArray[row][col])
            y_dir = flow_dir_y[flow_dir[0][0]]
            x_dir = flow_dir_x[flow_dir[1][0]]
            row += y_dir
            col += x_dir
        else:
            DSPath_array[row][col] = 1
  
    #Transform to pcr object:
    Areacomp = pcr.numpy2pcr(pcr.Ordinal,DSPath_array,Nodata)
    
    #Calculate clumps:
    pcr.setglobaloption("diagonal")
    Areanums = pcr.clump(Areacomp)
    
    #Tranform to a array:
    areaarray = pcr.pcr2numpy(Areanums,Nodata)

    #Determine area number for each location:
    XYZa = []
    for ob in XYZ:
        row = int(np.floor((ob[1]-TLy)/Yres))
        col = int(np.floor((ob[0]-TLx)/Xres))
        a=areaarray[row,col]
        XYZa.append((ob[0],ob[1],ob[2],a))
    
    return XYZa, areaarray
    
    
def FloodFill(floodpronearray,startpoints,neighbours = [(-1,0),(0,1),(1,0),(0,-1)]):
    """
    Performs a flood fill operation from the points specified in Startpoints
    
    Input:
        floodpronearray ---
        Boolean array specifying which cells can be flooded, and which cannot.
        
        startpoints ---
        Points to start flood fill operations from list with lists [row, col],
        or tuples (row,col)
        
        neighbours (Default:(four-side rule)) ---
        The relative cells to consider (in row, col). Can be used to set either
        a four or eight side flood fill rule.
        
    Returns:
        array with flooded cells ---
        Boolean array specifying which cells are flooded
    """
    #Initialize arrays:
    shape = floodpronearray.shape
    FParray = np.zeros((shape[0]+2,shape[1]+2),dtype=np.bool_)
    FParray[1:-1,1:-1] = floodpronearray #To prevent flood fill from crossing borders and ending up at negative rows/columns
    Flooded = np.zeros((shape[0]+2,shape[1]+2),dtype=np.bool_)
    queu = set()
    #Start from each startpoint:
    for point in startpoints:
        row = point[0]+1
        col = point[1]+1
        if FParray[row,col]:
            queu.add((row,col))
            Flooded[row,col] = True
            while queu:
                row,col = queu.pop()
                for drow,dcol in neighbours:
                    rw = row + drow
                    cl = col + dcol
                    if not Flooded[rw, cl] and FParray[rw, cl]:
                            queu.add((rw,cl))
                            Flooded[rw,cl] = True
    
    return Flooded[1:-1,1:-1]
    

def InterpAlongFlowPaths(observations, HAND, properties, pcrDTMloc, power,
                         smoothing, area_filtersinks, DWD = None,
                         neighbours = [(-1,0),(0,1),(1,0),(0,-1)]):
    """
    Interpolates the water levels along flow paths and calculates flood extents
    
    Input:
        observations (or XYZ) ---
        dictionary which contains the observations. each observation in 
        observations['tags'] has enrichments containing the location and
        water depths. Alternatively this can be a list of tuples with X,Y and
        WL values.
        
        HAND ---
        Numpy array containing the HAND values
        
        properties ---
        Properties of the HAND, loaded using the tif2array function
        
        pcrDTMloc ---
        Location of the .map file for the DTM, necessary for creating the LDD
        
        power ---
        Integer with the power parameter of the inverse distance weighting
        equation
        
        smoothing ---
        Integer with the smooting parameter for IDW (in number of cells)
        
        
        area_filtersinks ---
        area of sinks (in number of cells) to filter before creating the LDD
        used in grouping
        
        DWD (default:None) ---
        The default water depth to use for 'depthless' observations in cm. Is
        not used if XYZ file is supplied.
        
        neighbours (default:(four-side flood rule)) ---
        List containing tupples specifying the relative rows and columns to
        review in flood fill. Can either be a four or eight side flood rule.
        
    Returns:
        Water depth array ---
        Numpy array containing the water depths
    """
    #Load Dataset properties:
    TLx = properties['TLx']; TLy = properties['TLy'];
    Xres = properties['Xres']; Yres = properties['Yres'];
    Xsize = properties['Xsize']; Ysize = properties['Ysize']
    Nodata = -99999
    
    #Calculate the X/Y coordinates and water levels for each observation:
    if isinstance(observations,dict):
        XYZ = DeriveXYZ(observations,HAND,properties,DWD = DWD)
    else:
        XYZ = observations
    
    #Group the observations based on downstream flow paths:
    pcr.setglobaloption("unitcell")
    pcr.setclone(pcrDTMloc)
    pcrDTM = pcr.readmap(pcrDTMloc)
    pcrLDD = pcr.lddcreate(pcrDTM,99999999,99999999,area_filtersinks,99999999)
    XYZa , areaarray = groupbyflowpaths(XYZ,pcrLDD,properties)

    #Determine list of unique areas
    AreaNums = [ob[3] for ob in XYZa]
    UniqueAreaNums = list(set(AreaNums))
    maxarea = max(UniqueAreaNums)
    
    #Do Inverse distance weighting:
    Cum_WLtimesWeight = np.zeros((Ysize,Xsize))
    Cum_Weight = np.zeros_like(Cum_WLtimesWeight)
    for area in UniqueAreaNums:
        print "Processing area ",area," of ",maxarea
        AreaObs = [ob for ob in XYZa if (ob[2] and ob[3] == area)]#All enrichments belonging to the area
        
        #Extract the complete flowpath of the area
        areamap = np.full((Ysize,Xsize),Nodata,dtype = np.float64)#float64 to allow for any nodata value, incl. nan
        areamap[areaarray == area] = 1
        
        if len(AreaObs) > 1: #If 1 observation in the area --> no IDW necessary
            RestMap = pcr.numpy2pcr(pcr.Scalar,areamap,Nodata)#Resistance along flow path, to calculate distance using spread
            for ob in AreaObs: 
                row = int(np.floor((ob[1]-TLy)/Yres))
                col = int(np.floor((ob[0]-TLx)/Xres))        
                pointarray = np.zeros((Ysize,Xsize))
                pointarray[row,col]=1#Used to define starting point of spread
                Pointmap = pcr.numpy2pcr(pcr.Boolean,pointarray,Nodata)
                Dists = pcr.spread(Pointmap,0.1,RestMap)#0.1 used as starting value, otherwise deviding by zero in IDW
                Distsarray = pcr.pcr2numpy(Dists,np.inf)
                Weights = 1/np.power((Distsarray + smoothing),power)
                #Write results:
                Cum_WLtimesWeight += ob[2]*Weights
                Cum_Weight += Weights
        else:
            Cum_WLtimesWeight[areamap == 1] += AreaObs[0][2]
            Cum_Weight[areamap == 1] += 1
            
    InterpRivWLarray = Cum_WLtimesWeight / Cum_Weight#Calculate water levels
    InterpRivWLs = pcr.numpy2pcr(pcr.Ordinal,InterpRivWLarray*1000,np.nan)
    
    Subcatchments = pcr.subcatchment(pcrLDD,InterpRivWLs)#Spread waterlevels over upstream cells
    Subcatchmentsarray = pcr.pcr2numpy(pcr.scalar(Subcatchments),np.nan)/1000#*1000 and /1000 due to datatype of subcatchment fucntion
    
    #Subtract DEM elevation from subcatchment:
    WDarray = Subcatchmentsarray - HAND
    WDarray[WDarray <= 0] = np.nan
    
    #Now remove unconnected areas:
    Rownrs = np.arange(Ysize)#Following steps create a list of the cells on the flowpath, which is input for floodfill
    Colnrs = np.arange(Xsize)
    Colsgrid , Rowsgrid = np.meshgrid(Colnrs,Rownrs,copy =  False)
    FloodCols = Colsgrid[areaarray > 1]
    FloodRows = Rowsgrid[areaarray > 1]
    floodrowscols = zip(FloodRows,FloodCols)
    rest = WDarray > 0
    Flooded = FloodFill(rest,floodrowscols,neighbours = neighbours)
    WDarray[Flooded == False] = np.nan

    return WDarray
    

def generaterandomobs(XYWDmethod,HANDloc,number,DWD = None, deptherrors = None,
                      locationerrors = None, gridtype = 'single',
                      loctypes = None, HANDsuffix = '_HAND.tif'):
    """
    Function generates the random observations used as input for flood mapping
    
    Input:
        XYWDmethod ---
        List with X,Y,WD and type of locational reference for each of the
        observations, generated using DeviveXYWDmethod.
        
        HANDloc ---
        String with location of HAND. In case gridtype='single', a file
        location is given, if gridtypetype = 'random', the folder containg the
        random HAND maps is given. Filename format is #_HAND.tif.
        
        number ---
        integer with number of random observation sets to generate
        
        DWD (default:None) ---
        Value of the DWD parameter used to calculate the WL for 'depthless'
        observations. This can be a single integer in which case a single value
        is used, a list containing a lower and upper bound in which case a
        uniform distribution is used or a list containing mean, standard
        deviation and a lower and upper bound, in which case a clipped normal
        distribution is used.
        
        deptherrors (default:None) ---
        list containing the mean and standard deviation of water depth errors.
        
        locationerrors (default:None) ---
        Either a list with the mean and standard deviation of error (if
        loctypes=None) or a list containing a list with mean errors for each
        type of locational reference in loctypes and a list of standard
        deviations for each type of locational reference in loctypes.
        
        gridtype (default:'single') ---
        string with type of HAND, 'single' means one HAND is used (default),
        'random'means a set of random HAND maps is used.
        
        loctypes (default:None) ---
        list containing all locational types, in which order the mean errors
        and standard deviations are specified in locationalerrors.
        
        HANDsuffix (default:'_HAND.tif') ---
        The file name suffix used for loading the randomly generated HAND
        files.
        
    Returns:
        list with XYZ lists ---
        Each random realization of the XYZ collection is added to a list
    """
    nr_of_obs = len(XYWDmethod)
    
    if gridtype == 'single': #Load one HAND to have properties as reference
        HAND, properties = tif2array(HANDloc)
        TLx = properties['TLx']; TLy = properties['TLy'];
        Xres = properties['Xres']; Yres = properties['Yres'];
        Xsize = properties['Xsize']; Ysize = properties['Ysize']
    else:
        HANDlocname = HANDloc + '0_HAND.tif'
        HAND, properties = tif2array(HANDlocname)
        TLx = properties['TLx']; TLy = properties['TLy'];
        Xres = properties['Xres']; Yres = properties['Yres'];
        Xsize = properties['Xsize']; Ysize = properties['Ysize']

    if locationerrors:
        X_random = np.random.randn(number,nr_of_obs)
        Y_random = np.random.randn(number,nr_of_obs)
        
    if deptherrors:
        WD_random = np.random.randn(number,nr_of_obs)
        
    if isinstance(DWD,list):
        if len(DWD) == 2:
            DWD_random = np.random.uniform(DWD[0],DWD[1],number)
        elif len(DWD) == 4:
            DWD_random = DWD[0]+np.random.randn(number)*DWD[1]
            lb = DWD[2]
            ub = DWD[3]
            if lb > DWD[0] > ub:
                sys.exit("Mean DWD not in between lower bound and upper bound")
            for nr in range(number):#Allows for replacing values outside bounds. Simpy using bounds for these locations would give peaks at end of distribution
                while not (lb < DWD_random[nr] < ub):
                    DWD_random[nr] = DWD[0]+np.random.randn()*DWD[1]
                
    elif DWD:
        DWD_random = np.ones(number)*DWD
        
        
    XYZlist = []
    for nr in range(number):
        if gridtype == 'random':
             HANDlocname = HANDloc + str(nr) + HANDsuffix
             HAND, properties = tif2array(HANDlocname)
        
        XYZ=[]
        for obnr in range(nr_of_obs):
            X = XYWDmethod[obnr][0]
            Y = XYWDmethod[obnr][1]
            WD = XYWDmethod[obnr][2]
            
            if loctypes: #For different errors for each type of locational reference
                method = XYWDmethod[obnr][3]
                methodind = loctypes.index(method)
                X += locationerrors[0][methodind]+X_random[nr,obnr]*locationerrors[1][methodind]
                Y += locationerrors[0][methodind]+Y_random[nr,obnr]*locationerrors[1][methodind]
            elif locationerrors:
                X += locationerrors[0]+X_random[nr,obnr]*locationerrors[1]
                Y += locationerrors[0]+Y_random[nr,obnr]*locationerrors[1]
            
            if WD and deptherrors:
                WD += deptherrors[0]+WD_random[nr,obnr]*deptherrors[1]
            
            row = int(np.floor((Y-TLy)/Yres))
            col = int(np.floor((X-TLx)/Xres))    
            
            if (0 <= row < Ysize and 0 <= col <Xsize) and not np.isnan(HAND[row,col]): #Check if XY still in HAND. If not, observation is excluded.
                if WD:
                    GL = HAND[row,col]
                    WL = GL+(float(WD)/100)
                elif DWD:
                    GL = HAND[row,col]
                    WL = GL+(DWD_random[nr]/100)
                else:
                    WL = None
                XYZ.append((X,Y,WL))
        XYZlist.append(XYZ)
    
    return XYZlist
    
def generaterandomobs_alongstreet(XYWDlocmatch,HANDloc,number,streetpointdict,
                                  DWD = None, deptherrors = None,
                                  locationerrors = None, gridtype = 'single',
                                  HANDsuffix = '_HAND.tif'):
    """
    Function generates the random observations used as input for flood mapping
    
    Input:
        XYWDlocmatch ---
        List with X,Y,WD and locational reference for each of the observations,
        generated using DeviveXYWDlocmatch.
        
        HANDloc ---
        String with location of HAND. In case gridtype='single', a file
        location is given, if gridtypetype = 'random', the folder containg the
        random HAND maps is given. Filename format is #_HAND.tif.
        
        number ---
        integer with number of random observation sets to generate
        
        streetpointdict ---
        Dictionary with the locational references of streets, each containing
        a collection of points belonging to that street.
        
        DWD (default:None) ---
        Value of the DWD parameter used to calculate the WL for 'depthless'
        observations. This can be a single integer in which case a single value
        is used, a list containing a lower and upper bound in which case a
        uniform distribution is used or a list containing mean, standard
        deviation and a lower and upper bound, in which case a clipped normal
        distribution is used.
        
        deptherrors (default:None) ---
        list containing the mean and standard deviation of water depth errors.
        
        locationerrors (default:None) ---
        A list containing the mean and standard deviation in errors of points
        referenced in the dataset (intersections & POIs) as well as the mean
        and standard deviation in error of the streets referenced in the
        dataset. None if not used.
        
        gridtype (default:'single') ---
        string with type of HAND, 'single' means one HAND is used (default),
        'random'means a set of random HAND maps is used.
        
        HANDsuffix (default:'_HAND.tif') ---
        The file name suffix used for loading the randomly generated HAND
        files.
        
    Returns:
        list with XYZ lists ---
        Each random realization of the XYZ collection is added to a list
    """
    nr_of_obs = len(XYWDlocmatch)
    
    if gridtype == 'single': #Load one HAND to have properties as reference
        HAND, properties = tif2array(HANDloc)
        TLx = properties['TLx']; TLy = properties['TLy'];
        Xres = properties['Xres']; Yres = properties['Yres'];
        Xsize = properties['Xsize']; Ysize = properties['Ysize']
    else:
        HANDlocname = HANDloc + '0_HAND.tif'
        HAND, properties = tif2array(HANDlocname)
        TLx = properties['TLx']; TLy = properties['TLy'];
        Xres = properties['Xres']; Yres = properties['Yres'];
        Xsize = properties['Xsize']; Ysize = properties['Ysize']

    if locationerrors:
        X_random = np.random.randn(number,nr_of_obs) # Also used for locational errors along streets
        Y_random = np.random.randn(number,nr_of_obs)
        
    if deptherrors:
        WD_random = np.random.randn(number,nr_of_obs)
        
    if isinstance(DWD,list):
        if len(DWD) == 2:
            DWD_random = np.random.uniform(DWD[0],DWD[1],number)
        elif len(DWD) == 4:
            DWD_random = DWD[0]+np.random.randn(number)*DWD[1]
            lb = DWD[2]
            ub = DWD[3]
            if lb > DWD[0] > ub:
                sys.exit("Mean DWD not in between lower bound and upper bound")
            for nr in range(number):#Allows for replacing values outside bounds. Simpy using bounds for these locations would give peaks at end of distribution
                while not (lb < DWD_random[nr] < ub):
                    DWD_random[nr] = DWD[0]+np.random.randn()*DWD[1]
                
    elif DWD:
        DWD_random = np.ones(number)*DWD
        
        
    XYZlist = []
    for nr in range(number):
        if gridtype == 'random':
             HANDlocname = HANDloc + str(nr) + HANDsuffix
             HAND, properties = tif2array(HANDlocname)
        
        XYZ=[]
        for obnr in range(nr_of_obs):
            X = XYWDlocmatch[obnr][0]
            Y = XYWDlocmatch[obnr][1]
            WD = XYWDlocmatch[obnr][2]
            locmatch = XYWDlocmatch[obnr][3]
            
            if locationerrors and locmatch in streetpointdict: #For different errors for streets and POIs
                distance_along_street = locationerrors[2]+X_random[nr,obnr]*locationerrors[3]
                streetpoints = streetpointdict[locmatch]
                dists = []
                for point in streetpoints: #Calculate distances to each point
                    dist = ((point[0] - X)**2+(point[1] - Y)**2)**0.5
                    if point[1] < Y:
                        dist *= -1
                    dists.append(dist)
                    
                dists = np.array(dists)
                
                if np.min(np.absolute(dists)) > 50:
                    sys.exit("Point not on specified street! at nr %d and locmatch %s" %(obnr,locmatch))
                
                while np.min(dists) > distance_along_street > np.max(dists):
                    distance_along_street = locationerrors[2]+np.random.randn()*locationerrors[3]
                    
                differences = dists - distance_along_street
                
                pointind = np.argmin(np.absolute(differences))
                
                X = streetpoints[pointind][0]
                Y = streetpoints[pointind][1]
                    
            elif locationerrors:
                X += locationerrors[0]+X_random[nr,obnr]*locationerrors[1]
                Y += locationerrors[0]+Y_random[nr,obnr]*locationerrors[1]
            
            if WD and deptherrors:
                WD += deptherrors[0]+WD_random[nr,obnr]*deptherrors[1]
            
            row = int(np.floor((Y-TLy)/Yres))
            col = int(np.floor((X-TLx)/Xres))    
            
            if (0 <= row < Ysize and 0 <= col <Xsize) and not np.isnan(HAND[row,col]): #Check if XY still in HAND. If not, observation is excluded.
                if WD:
                    GL = HAND[row,col]
                    WL = GL+(float(WD)/100)
                elif DWD:
                    GL = HAND[row,col]
                    WL = GL+(DWD_random[nr]/100)
                else:
                    WL = None
                XYZ.append((X,Y,WL))
        XYZlist.append(XYZ)
    
    return XYZlist
    
def generatesimulationinput(XYZlist,power,smoothing,area_filtersinks,
                            HANDloc,DTMloc,nprocesses,inputfolder,
                            gridtype = 'single',DTMsuffix = '_DTM.tif',
                            HANDsuffix = '_HAND.tif'):
    """
    Generates the input of each random simulation. This input file for each
    simulation contains the observation and settings for the interpolation process.
    Allows for variation of all input datasets and parameters between runs.
    
    Input:
        XYZlist ---
        List lists of XYZ values for each random realization
        
        power ---
        The power parameter for the IDW equation. This can be a single integer
        in which case a single value is used, a list containing a lower and
        upper bound in which case a uniform distribution is used or a list
        containing mean, standard deviation and a lower and upper bound, in
        which case a clipped normal distribution is used. The fact that these
        values are rounded to integers using np.round means that for example
        for a uniform distribution from 2 to 5, the setting [1.5,5.5] should be
        used, to have an equal probability of occurence.
        
        smoothing ---
        The smooting parameter for the IDW equation. This can be a single
        integer in which case a single value is used, a list containing a lower
        and upper bound in which case a uniform distribution is used or a list
        containing mean, standard deviation and a lower and upper bound, in
        which case a clipped normal distribution is used.
        
        area_filtersinks ---
        integer, parameter decides to which extent to filter sinks, before
        grouping.
        
        HANDloc ---
        File location of HAND if gridtype = 'single' or folder with random
        HAND maps. (File format #_HAND.tif)
        
        DTMloc ---
        File loction of DTM is gridtype = 'single' or folder with random DTMs.
        (File format #_DTM.tif)
        
        nprocesses ---
        The number of processes the runs will be devided over. Is used to 
        control the copying and assigning of HAND/DTM files to each run, to
        ensure each process has it's own HAND/DTM.
        
        inputfolder ---
        String with folder to write the input files
        
        gridtype (default:'single') ---
        Specifies if a single grid or random grids are used, using 'single' or
        'random' respectively.
        
        DTMsuffix (default:'_DTM.tif') ---
        The file name suffix used for loading the randomly generated DTM files.
        
        HANDsuffix (default:'_HAND.tif') ---
        The file name suffix used for loading the randomly generated HAND
        files.
    
    Output:
        Input files ---
        Input files for each simulation written to the inputfolder
    """
    number = len(XYZlist)
    
    if isinstance(power,list):#Random power variables
        if len(power) == 2:
            power_random = np.round(np.random.uniform(power[0],power[1],number))
        elif len(power) == 4:
            power_random = np.round(power[0]+np.random.randn(number)*power[1])
            lb = power[2]
            ub = power[3]
            if lb > power[0] > ub:
                sys.exit("Mean power not in between lower bound and upper bound")
            for nr in range(number):#Allows for replacing values outside bounds. Simpy using bounds for these locations would give peaks at end of distribution
                while not (lb < power_random[nr] < ub):
                    power_random[nr] = np.round(power[0]+np.random.randn()*power[1])
    else:
        power_random = np.ones((number))*power
        
    if isinstance(smoothing,list):
        if len(smoothing) == 2:
            smoothing_random = np.random.uniform(smoothing[0],smoothing[1],number)
        elif len(smoothing) == 4:
            smoothing_random = smoothing[0]+np.random.randn(number)*smoothing[1]
            lb = smoothing[2]
            ub = smoothing[3]
            if lb > smoothing[0] > ub:
                sys.exit("Mean smoothing not in between lower bound and upper bound")
            for nr in range(number):#Allows for replacing values outside bounds. Simpy using bounds for these locations would give peaks at end of distribution
                while not (lb < smoothing_random[nr] < ub):
                    smoothing_random[nr] = smoothing[0]+np.random.randn()*smoothing[1]
            
    else:
        smoothing_random = np.ones((number))*smoothing
        
    inputfiles = []
    if gridtype == 'single':
        HANDcopyfolder = os.path.dirname(HANDloc) + '\\Copies\\'
        DTMcopyfolder = os.path.dirname(DTMloc) + '\\Copies\\'
        if not os.path.isdir(HANDcopyfolder): #Create directories
            os.makedirs(HANDcopyfolder)  
        if not os.path.isdir(DTMcopyfolder):
            os.makedirs(DTMcopyfolder)
        HANDfn = os.path.basename(HANDloc)
        DTMfn = os.path.basename(DTMloc)
        
        rngstart = 0
        for nr in range(nprocesses):#Copies the HAND and DTM, so that each process uses a seperate copy
            newHANDloc = HANDcopyfolder + HANDfn[:-4] + '_copy' + str(nr) + '.tif'
            newDTMloc = DTMcopyfolder + DTMfn[:-4] + '_copy' + str(nr) + '.tif'
            os.system ("copy %s %s" %(HANDloc, newHANDloc))
            os.system ("copy %s %s" %(DTMloc, newDTMloc))
            rngend = rngstart+(number/nprocesses)
            for nr in range(rngstart,rngend):
                XYZ = XYZlist[nr]
                powerparam = power_random[nr]
                smoothingparam = smoothing_random[nr]
                inputfiles.append({'observations':XYZ, 'HANDloc':newHANDloc,
                               'DTMloc':newDTMloc,'power':powerparam,
                               'smoothing':smoothingparam,
                               'area_filtersinks':area_filtersinks})
            rngstart = rngend
    else:
        for nr in range(number):#Combine all input
            XYZ = XYZlist[nr]
            powerparam = power_random[nr]
            smoothingparam = smoothing_random[nr]
            HANDlocation = HANDloc + str(nr) + HANDsuffix
            DTMlocation = DTMloc + str(nr) + DTMsuffix
            inputfiles.append({'observations':XYZ, 'HANDloc':HANDlocation,
                               'DTMloc':DTMlocation,'power':powerparam,
                               'smoothing':smoothingparam,
                               'area_filtersinks':area_filtersinks})
    
    for nr in range(number):#Write all input
        with open(inputfolder+str(nr)+'.inp','w') as inputfile:
            json.dump(inputfiles[nr],inputfile)
            
def mergefloodmaps(outputfolder,number):
    """
    Merges all realizations of flood maps till 'number' in the outputfolder.
    The floodmaps are titled realization_#.tif.
    
    Input:
        outputfolder ---
        The folder in which the output is located
        
        number ---
        The number of realizations in the folder
        
    Return:
        Probabilitygrid ---
        Array with the probability of each cell being flooded
    """
    ds = gdal.Open(outputfolder + 'realization_0.tif',0)
    Nodata=ds.GetRasterBand(1).GetNoDataValue()
    Nrtimesflooded = ds.GetRasterBand(1).ReadAsArray()
    Nrtimesflooded = Nrtimesflooded.astype(np.float64)
    Nrtimesflooded[Nrtimesflooded == Nodata] = np.nan
                  
    for nr in range(1,number):
        if nr % 100 == 0:
            print "Merging grid %d of %d" %(nr,number)
        ds = gdal.Open(outputfolder + 'realization_' + str(nr) + '.tif',0)
        Nrtimesflooded += ds.GetRasterBand(1).ReadAsArray()
    
    Probabilitygrid = Nrtimesflooded / float(number)
    
    return Probabilitygrid
    
def confusionmatrixmap(WDarray,valarray):
    """
    Function creates a map that classifies the cells in the grid of water
    depths based on the confusion matrix.
    
    Input:
        WDarray ---
        Array containing water depths, nan for non-flooded areas and nodata
        
        valarray ---
        Array containing the validation data, ones for flooded cells, 0 for not
        flooded, and nan outside of area. Location of nan values in this array
        is used to determine the outline of the confusion matrix array.
        
    Returns:
        confmaparray ---
        Array with different values for the different quandrants of the
        confusion matrix. Contains Nodata (nan) in the cells in which the
        WDarray contained Nodata values. 1 means flooded in both arrays, 2
        means flooded in WDarray but not in valarray (overestimation), 3 means
        flooded in valarray but not in WDarray (underestimation) and 4 means
        flooded in neither of the arrays.
    """
    shape = WDarray.shape
    confmaparray = np.zeros(shape)
    floodextent = np.zeros(shape)
    floodextent[WDarray > 0] = 1
    
    confmaparray[(floodextent == 1) & (valarray == 1)] = 1
    confmaparray[(floodextent == 1) & (valarray == 0)] = 2
    confmaparray[(floodextent == 0) & (valarray == 1)] = 3
    confmaparray[(floodextent == 0) & (valarray == 0)] = 4
    confmaparray[np.isnan(valarray)] = np.nan
    
    return confmaparray

def reliabilitydiagram(probabilityarray,valarray,nr_of_intvs = 10):
    """
    Function creates a reliability diagram of a probabilistic flood extent map
    using recorded flood extents. Results of the probabilistic map are devided
    in equal intervals of probability. Of these intervals it is determined in 
    how many of the cells belonging to them, there were actually floods. This
    data can be used to construct a reliability diagram.
    
    Input:
        probabilityarray ---
        Array with the probability of each cell being flooded. Contains values
        between 0 and 1 in the study area and nan values outside.
        
        valarray ---
        Array with recorded flood extents. A cell is 1 if it's flooded and 0 if
        it's not. Outside the study area there are nan values.
        
        nr_of_intvs (default = 10) ---
        Specifies the number of intervals that is used. Actually one more than
        this number is plotted, since cells having 0% of flooding are
        considered seperately.
        
    Returns:
        modprob ---
        A list with X-values (average modelled probability)
        
        obsprob ---
        A list with Y-values (percentage of cells recorded flooded in interval)
        
        nr_of_cells ---
        A list with the number of cells per bin
    """
    probabilityarray[np.isnan(valarray)] = np.nan
    modprob = []
    obsprob = []
    nr_of_cells = []
    
    modprob.append(0)
    valdata = valarray[probabilityarray == 0]
    obsprob.append(np.mean(valdata))
    nr_of_cells.append(valdata.size)
    
    for intv_start in np.linspace(0,1,nr_of_intvs,endpoint = False):
        intv_end = intv_start + (1/float(nr_of_intvs))
        modprob.append((intv_start+intv_end)/2)
        valdata = valarray[(probabilityarray > intv_start) & (probabilityarray <= intv_end)]
        obsprob.append(np.mean(valdata))
        nr_of_cells.append(valdata.size)
    
    return modprob, obsprob, nr_of_cells
    
def Fstatistic(floodmaparray,valarray):
    """
    Calculates the F(2) statistic described by Aronica et al. (2002).
    dx.doi.org/10.1002/hyp.398
    
    Input:
        floodmaparray ---
        Array containing values > 1 for areas that are flooded, and 0 for areas
        which are not flooded. np.nan for nodata only outside the area
        considered
        
        valarray ---
        Array containing values of 1 for observed flooded areas and 0 for
        observed non-flooded areas. np.nan cells will not be evaluated in
        calculating the F-statistic. The cells having either zeros or ones here
        should not contain nodata in the floodmaparray
        
    Returns:
        Fvalue ---
        Float with the f-value calculated from the two maps
    """
    area_correctly_flooded = np.sum(((floodmaparray > 0) & (valarray == 1)))
    area_flooded_both = np.sum((((floodmaparray > 0) & (valarray >= 0)) | (valarray == 1))) #Only area where valarray has values is reviewed
    
    Fvalue  = float(area_correctly_flooded)/area_flooded_both
    
    return Fvalue

def create_F2_ECDF(outputfolder,number,valarray):
    """
    Calculates the Empirical-Cumulative distribution function of the F(2). Writes the
    resulting X (F-value) and Y (percentage) to a json-file in the outputfolder
    
    Input:
        outputfolder ---
        String with the folder where the simulation results are stored. Here
        json file with the E-cdf description is written.
        
        number ---
        Integer with the number of files in the outputfolder
        
        valarray ---
        String with location of the tiff file used for validation. This
        contains 1 for flooded cells, 0 for non-flooded cells and nan for cells
        which are not considered by the analysis. 
        
    Output:
        F_ECDF.json ---
        This file contains the x and Y values to describe the E-CDF of all runs
    """
    F2_vals = []
    Y_vals = []
    for nr in range(0,number):
        if nr % 100 == 0:
            print "Evaluating grid %d of %d" %(nr,number)
        ds = gdal.Open(outputfolder + 'realization_' + str(nr) + '.tif',0)
        flooded =  ds.GetRasterBand(1).ReadAsArray()
        probability_underestimation = float(nr)/number
        F_2_statistic = Fstatistic(flooded,valarray)
        F2_vals.append(F_2_statistic)
        Y_vals.append(probability_underestimation)
        
    F2_vals.sort()
    ECDF_XY = [F2_vals,Y_vals]
        
    with open(outputfolder + 'F_ECDF.json','w') as jsonfile:
        json.dump(ECDF_XY,jsonfile)