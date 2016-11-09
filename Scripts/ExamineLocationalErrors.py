# -*- coding: utf-8 -*-
"""
Analyze XYZlist errors

Analyzes the real XY errors that were introduced to the observations, by
comparing the XYZlist with the actual observations. Can only be run after
the UncertaintyAnalysis scripts.

Metadata:
    Version ---
    1.0 (18-08-2016)
    
    Type ---
    Script, dependent on other script running in advance
    
    Dependencies ---
    Running UncertaintyAnalysis or UncertaintyAnalysis_alongstreet prior to
    execution & numpy
    
    Author ---
    Tom Brouwer
"""
import numpy as np

XYWDmethod = FloodInterpolationTools.DeriveXYWDmethod(observations,properties)

indices = range(len(XYWDmethod))
streetindices = [ind for ind,XYWDmet in zip(indices,XYWDmethod) if ('street' in XYWDmet[3].lower()) and not ('poi' in XYWDmet[3].lower())]
pointindices = [ind for ind,XYWDmet in zip(indices,XYWDmethod) if not (('street' in XYWDmet[3].lower()) and not ('poi' in XYWDmet[3].lower()))]

devsstreet = []
for nr in range(len(XYZlist)):
    for ind in streetindices:
        devsstreet.append(XYZlist[nr][ind][0] - XYWDmethod[ind][0])
        devsstreet.append(XYZlist[nr][ind][1] - XYWDmethod[ind][1])
        
print "Standard deviation for locational errors of Tweets referring to streets: %.0f" %np.std(devsstreet,ddof = 1)

devspoint = []
for nr in range(len(XYZlist)):
    for ind in pointindices:
        devspoint.append(XYZlist[nr][ind][0] - XYWDmethod[ind][0])
        devspoint.append(XYZlist[nr][ind][1] - XYWDmethod[ind][1])
        
print "Standard deviation for locational errors of Tweets referring to points: %.0f" %np.std(devspoint,ddof = 1)

alldevs = devsstreet + devspoint

print "Standard deviation for all locational errors of Tweets : %.0f" %np.std(alldevs,ddof = 1)