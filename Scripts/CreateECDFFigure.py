# -*- coding: utf-8 -*-
"""
Function creates a graph of E-CDF of F-values of simulating the three different
error sources

Input:
    DTM_ECDFloc ---
    String with location of the json file, in which the ECDF of the DTM error is
    stored.
    
    param_ECDFloc ---
    String with the location of the json file, in which the ECDF of the error
    in parameters is stored.
    
    Loc_ECDFloc ---
    String with the location of the json file, in which the ECDF of the
    locational errors is stored.
    
    figsaveloc ---
    String with location to save the resulting figure. If not, then None
    
Output:
    Figure ---
    Graph with the E-CDFs of the F-statistic values of each of the sources of
    uncertainty. Saved to figsaveloc.
    
Metadata:
    Version ---
    1.0 (25-08-2016)
    
    Type ---
    Script
    
    Dependencies ---
    json, matplotlib
    
    Author ---
    Tom Brouwer
"""
import json
import matplotlib.pyplot as plt
import matplotlib as mpl

#%% input:
DTM_ECDFloc = r'..\Output\DTM_error_ECDF.json'
param_ECDFloc = r'..\Output\param_error_ECDF.json'
Loc_ECDFloc = r'..\Output\loc_error_ECDF.json'
All_ECDFloc = r'..\Output\all_error_ECDF.json'
figsaveloc = None

#%% script
with open(DTM_ECDFloc,'r') as DTMjson:
    DTM_ECDF = json.load(DTMjson)
    
with open(param_ECDFloc,'r') as paramjson:
    param_ECDF = json.load(paramjson)
    
with open(Loc_ECDFloc,'r') as Locjson:
    Loc_ECDF = json.load(Locjson)
    
with open(All_ECDFloc,'r') as Alljson:
    All_ECDF = json.load(Alljson)
    
tnrfont = {'family':'Times New Roman','size':10}
mpl.rc('font',family = 'Times New Roman')
plt.close(1)
plt.figure(1,figsize = (8.5/2.54,6/2.54))
ax = plt.gca()
ax.set_xlim([0,1])
ax.set_ylim([0,1])
ax.set_xticklabels(ax.get_xticks(), fontdict = tnrfont)
ax.set_yticklabels(ax.get_yticks(), fontdict = tnrfont)
plt.xlabel('F$^2$ value', fontdict = tnrfont)
plt.ylabel('Probability X <= x', fontdict = tnrfont)
plt.plot(DTM_ECDF[0],DTM_ECDF[1],'-b',
         param_ECDF[0],param_ECDF[1],'-r',
         Loc_ECDF[0],Loc_ECDF[1],'-k',
         All_ECDF[0],All_ECDF[1],'--g')
plt.legend(['Elev. errors','Param. errors','Loc. errors','All errors'],loc = 'upper left',fontsize = 9)
plt.gcf().subplots_adjust(bottom = 0.19, top = 0.96, left = 0.14,right = 0.97)
if figsaveloc:    
    plt.savefig(figsaveloc)
plt.show()

