# -*- coding: utf-8 -*-
"""
Calculate X/Y coordinate errors

Script calculates the errors in X/Y coordinates of the observations

Input:
    jsonloc ---
    String with file location of the json file containing the observations
    Lat/Lon derived from text as well as the Lat/Lon derived from a photograph.
    
    EPSG ---
    Integer specifying the projection number to use
    
    figureloc ---
    String with the location to save the figure. None if not used
    
Output:
    Figure ---
    Figure is optionally saved to the location specified in figureloc
    
Metadata:
    Version ---
    1.0 (15-08-16 12:07)
    
    Type ---
    Script
    
    Dependencies ---
    FloodInterpolationTools (custom module), Numpy, Matplotlib, json
    
    Author ---
    Tom Brouwer
"""
import numpy as np
import matplotlib.pyplot as plt
import json
import FloodInterpolationTools

#%% Input:
jsonloc = r'..\Data\observations_york.json'
EPSG = 27700
y_range = [0,30]
x_range=[-2000,2001]
hist_range=[-1950,1950]
hist_intval = 100
tick_intvalxy = 1000
figure_spacing = {'bottom':0.15,
                  'top':0.95,
                  'left':0.13,
                  'right':0.95}
texts_x_align = 0.70
figureloc = None

#%% Script:
with open(jsonloc,'r') as jsonfile: # Load the json file
    observations = json.load(jsonfile)
    
Enrichments = [enr for ob in observations['tags'] for enr in ob['Enrichments'] if enr['PhotoLatLon']]

latlons = [enr['LatLon'] for enr in Enrichments] #Extract locations derived from text and photo
photolatlons = [enr['PhotoLatLon'] for enr in Enrichments]
methods = [enr['Method'] for enr in Enrichments]

XY_text = FloodInterpolationTools.DeriveXY(latlons,EPSG) #Reproject both
XY_photo = FloodInterpolationTools.DeriveXY(photolatlons,EPSG)
all_inds = range(len(XY_text))

indices_streetrefs = [ind for ind,method in zip(all_inds,methods) if 'street' in method.lower() and method.lower() != 'street near poi']
indices_pointrefs = [ind for ind,method in zip(all_inds,methods) if not ('street' in method.lower() and method.lower() != 'street near poi')]

differences_all = []
for enrnr in range(len(XY_text)):
    differences_all.append(XY_text[enrnr][0]-XY_photo[enrnr][0])
    differences_all.append(XY_text[enrnr][1]-XY_photo[enrnr][1])
    
differences_streetrefs = []
for enrnr in indices_streetrefs:
    differences_streetrefs.append(XY_text[enrnr][0]-XY_photo[enrnr][0])
    differences_streetrefs.append(XY_text[enrnr][1]-XY_photo[enrnr][1])
    
differences_pointrefs = []
for enrnr in indices_pointrefs:
    differences_pointrefs.append(XY_text[enrnr][0]-XY_photo[enrnr][0])
    differences_pointrefs.append(XY_text[enrnr][1]-XY_photo[enrnr][1])
    
  
tnrfont = {'fontname':'Times New Roman','size':10}
plt.close(1)
fig = plt.figure(1,figsize = (8.5/2.54,9/2.54))

ax1 = plt.subplot2grid((2,4),(0,1),colspan = 2)
ax1.hist(differences_all,bins=range(-1950,1951,100))
plt.xticks(range(-2000,2001,1000))
plt.ylabel("Frequency",fontdict = tnrfont)
plt.xlabel("Deviation (m)",fontdict = tnrfont)
plt.figtext(0.37,0.93,'(a)',fontdict = tnrfont)
plt.figtext(0.56,0.87,'$\sigma=%d$ $m$ \n$n=%d$' %(np.std(differences_all,ddof = 1),len(differences_all)),fontsize=9)
ax1.set_ylim([0,30])
ax1.set_xticklabels(ax1.get_xticks(), fontdict = tnrfont)
ax1.set_yticklabels(ax1.get_yticks().astype(int), fontdict = tnrfont)

ax2 = plt.subplot2grid((2,2),(1,0))
ax2.hist(differences_streetrefs,bins=range(-1950,1951,100))
plt.xticks(range(-2000,2001,1000))
plt.ylabel("Frequency",fontdict = tnrfont)
plt.xlabel("Deviation (m)",fontdict = tnrfont)
plt.figtext(0.15,0.43,'(b)',fontdict = tnrfont)
plt.figtext(0.32,0.36,'$\sigma=%d$ $m$ \n$n=%d$' %(np.std(differences_streetrefs,ddof = 1),len(differences_streetrefs)),fontsize=9)
ax2.set_ylim([0,22])
ax2.set_xticklabels(ax2.get_xticks(), fontdict = tnrfont)
ax2.set_yticklabels(ax2.get_yticks().astype(int), fontdict = tnrfont)

ax3 = plt.subplot2grid((2,2),(1,1))
ax3.hist(differences_pointrefs,bins=range(-190,191,20))
plt.xticks(range(-200,201,100))
plt.xlabel("Deviation (m)",fontdict = tnrfont)
plt.figtext(0.62,0.43,'(c)',fontdict = tnrfont)
plt.figtext(0.80,0.37,'$\sigma=%d$ $m$ \n$n=%d$' %(np.std(differences_pointrefs,ddof = 1),len(differences_pointrefs)),fontsize=9)
ax3.set_ylim([0,6])
ax3.set_xticklabels(ax3.get_xticks(), fontdict = tnrfont)
ax3.set_yticklabels(ax3.get_yticks().astype(int), fontdict = tnrfont)

if figure_spacing:
    fig.subplots_adjust(bottom = 0.11, top = 0.97, left = 0.12,
                                right = 0.96, wspace = 0.25, hspace = 0.36)
if figureloc:
    plt.savefig(figureloc)
plt.show()