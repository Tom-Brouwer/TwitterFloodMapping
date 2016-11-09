# -*- coding: utf-8 -*-
"""
Create reliability diagram

Function creates a reliability diagram using an uncertainty map containing values between 0 and 1 and a
validation dataset with 1/0 values for flooded and not flooded respectively.

Input:
    uncertaintyloc ---
    String with location of probabilistic flood extent map, with values between
    1 and 0. np.nan values in nodata areas
    
    validmaploc ---
    String with location of the validation map, which contains a 1 if a cells
    is flooded, and a 0 if it's not. Nan values outside study area, nan values
    in this map are indicative of the processing extent.
    
    excludeloc ---
    String with location of raster specifying which cells to exclude. This can
    for example be used to exclude a buffer of cells from the edge of the map.
    contains 0 on the part of the map to be analyzed, and ones on the part of
    the map to be ignored.
    
    figuresaveloc ---
    String with file location to save file. None if not used
    
Output:
    Figure ---
    The result is presented in figure 1, and if figuresaveloc != None, written
    to disk
    
Metadata:
    Version ---
    1.0
    
    Type ---
    Script
    
    Dependencies ---
    FloodInterpolationTools
    
    Author ---
    Tom Brouwer
"""
import numpy as np
import matplotlib.pyplot as plt
import FloodInterpolationTools

#%% Input:
uncertaintyloc = r'..\Output\probabilistic_map_all_error.tif'
validmaploc = r'..\Data\york_validation_innercity.tif'
excludeloc = None
figuresaveloc = None

#%% Script:
probabilityarray, properties = FloodInterpolationTools.tif2array(uncertaintyloc)
validationarray, properties = FloodInterpolationTools.tif2array(validmaploc)

if excludeloc:
    excludearray, pptys = FloodInterpolationTools.tif2array(excludeloc)
    probabilityarray[excludearray == 1] = np.nan
    
modprob, obsprob, nr_of_cells = FloodInterpolationTools.reliabilitydiagram(probabilityarray,validationarray)

tnrfont = {'fontname':'Times New Roman','size':10}
plt.close(1)
fig = plt.figure(1,figsize = (4.5/2.54,4.5/2.54))
ax = fig.add_subplot(111)
plt.plot([0,1],[0,1],'--k',modprob,obsprob,'.-b')
plt.xlabel('Modelled probability',fontdict = tnrfont)
plt.ylabel('Observed probability',fontdict = tnrfont)
ax.set_xticklabels(ax.get_xticks(), fontdict = tnrfont)
ax.set_yticklabels(ax.get_yticks(), fontdict = tnrfont)


tnrfont = {'fontname':'Times New Roman','size':9}
ax2 = fig.add_axes([.30,.66,0.24,0.24],axisbg='w')
ax2.bar(modprob[2:],nr_of_cells[2:],width = 0.1,color = 'white',align = 'center')
ax2.set_xlim([0.1,1])
plt.xticks([0.1,1])
upperytick = int(ax2.get_yticks()[-1])
plt.yticks([0,upperytick])
ax2.set_xticklabels(ax2.get_xticks(), fontdict = tnrfont)
ax2.set_yticklabels(ax2.get_yticks(), fontdict = tnrfont)
ax2.yaxis.tick_right()

plt.gcf().subplots_adjust(bottom = 0.22, top = 0.96, left = 0.25,right = 0.95)
plt.show()
if figuresaveloc:
    plt.savefig(figuresaveloc)