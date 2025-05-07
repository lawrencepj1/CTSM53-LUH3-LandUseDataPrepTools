#! /usr/bin/env python
import sys
import os.path
import string
import subprocess
import datetime as date
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import netCDF4 as netcdf4
from scipy import interpolate
from math import sqrt
from numpy.polynomial import Chebyshev as T
from numpy.polynomial import Polynomial as P

# load proper modules first, i.e.
# cgd machines
'''
module load lang/python/2.7.14

'''

def readdatavalue(datasetfilename,variable1name,variable2name):
    datasetfile = open(datasetfilename,'r')
    datasetlist = datasetfile.readlines()
    numdatasetlist = len(datasetlist)
    variablelist = datasetlist[0].strip().split(",")
    numvariablelist = len(variablelist)
    variable1index = -1
    variable2index = -1
    datasetvalue = 0.0
    for variablesearchindex in range(numvariablelist):
        if (variable1name == variablelist[variablesearchindex]):
            variable1index = variablesearchindex    
    if (variable1index == -1):
        print("Warning: " + datasetfilename + " Variable Not Found " + variable1name)
    else:
        for datasetlistindex in range(numdatasetlist):
            valuelist = datasetlist[datasetlistindex].split(",")
            if (variable2name == valuelist[0]):
                variable2index = datasetlistindex
                datasetvalue = float(valuelist[variable1index])
    if (variable2index == -1):
        print("Warning: " + datasetfilename + " Variable Not Found " + variable2name)

    return datasetvalue


#--  end of function definitions  ---------------------------------

# process the input arguments

arguments = len(sys.argv) - 1

if (arguments != 10):
    print("Error: Usage createrawglobalmonthlyareasumdata.py pft1filename pft1climateid pft2filename pft2climateid pft3filename pft3climateid output1filename output2filename output3filename climatemaskfilename")
    sys.exit()
else:
    print("Processing: " + str(sys.argv[1]) + " => " + str(sys.argv[7]))

pft1filename = str(sys.argv[1])
pft1climateid = str(sys.argv[2])
pft2filename = str(sys.argv[3])
pft2climateid = str(sys.argv[4])
pft3filename = str(sys.argv[5])
pft3climateid = str(sys.argv[6])

numprecipclasses = 100
precipclassstart = 0.0
precipclasssize = 50.0
numtempclasses = 100
tempclassstart = -15.0
tempclasssize = 0.5
full1matrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
full2matrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
full3matrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
poly1matrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
poly2matrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
poly3matrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft1matrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft2matrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft3matrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft1output = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft2output = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft3output = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
outputcsv1filename = str(sys.argv[7])
outputcsv2filename = str(sys.argv[8])
outputcsv3filename = str(sys.argv[9])
outputstr = "{0:.3f}"


for precipannindex in range(numprecipclasses):
    precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
    precipannstr = str(precipannvalue)
    for tempaverageindex in range(numtempclasses):
        tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
        tempaveragestr = str(tempaveragevalue)
        full1matrix[precipannindex,tempaverageindex] = readdatavalue(pft1filename,tempaveragestr,precipannstr)
        
for precipannindex in range(numprecipclasses):
    precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
    for tempaverageindex in range(numtempclasses):
        tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
        numtempvalues = 0
        tempnear = 0
        precipnear = 0
        for tempsearchindex in range(numtempclasses):
            if (full1matrix[precipannindex,tempsearchindex] > 0.0):
                numtempvalues += 1
        if (numtempvalues > 1):
            tempvaluearray = np.zeros(shape=(numtempvalues),dtype=np.float64)
            pftvaluearray = np.zeros(shape=(numtempvalues),dtype=np.float64)
            tempvalueindex = 0
            for tempsearchindex in range(numtempclasses):
                if (full1matrix[precipannindex,tempsearchindex] > 0.0):
                    tempvaluearray[tempvalueindex] = tempclassstart + float(tempsearchindex) * tempclasssize
                    pftvaluearray[tempvalueindex] = full1matrix[precipannindex,tempsearchindex]
                    tempdistance = abs(tempsearchindex - tempaverageindex)
                    if (tempdistance < 2):
                        tempnear = 1
                    tempvalueindex += 1
            tcoeff = P.fit(tempvaluearray, pftvaluearray, deg=3)
            tcoeffvals = tcoeff.convert().coef
            if (len(tcoeffvals) >= 4):
                pftvaluetemppoly = tcoeffvals[0] + tcoeffvals[1] * tempaveragevalue + tcoeffvals[2] * pow(tempaveragevalue,2) + tcoeffvals[3] * pow(tempaveragevalue,3) 
            if (len(tcoeffvals) == 3):
                pftvaluetemppoly = tcoeffvals[0] + tcoeffvals[1] * tempaveragevalue + tcoeffvals[2] * pow(tempaveragevalue,2)
            if (pftvaluetemppoly < 0.0):
                tempnear = 0
        else:
            pftvaluetemppoly = 0.0
        numprecipvalues = 0
        for precipsearchindex in range(numprecipclasses):
            if (full1matrix[precipsearchindex,tempaverageindex] > 0.0):
                numprecipvalues += 1
        if (numprecipvalues > 1):
            precipvaluearray = np.zeros(shape=(numprecipvalues),dtype=np.float64)
            pftvaluearray = np.zeros(shape=(numprecipvalues),dtype=np.float64)
            precipvalueindex = 0
            for precipsearchindex in range(numprecipclasses):
                if (full1matrix[precipsearchindex,tempaverageindex] > 0.0):
                    precipvaluearray[precipvalueindex] = precipclassstart + float(precipsearchindex) * precipclasssize
                    pftvaluearray[precipvalueindex] = full1matrix[precipsearchindex,tempaverageindex]
                    precipdistance = abs(precipsearchindex - precipannindex)
                    if (precipdistance < 3):
                        precipnear = 1
                    precipvalueindex += 1
            pcoeff = P.fit(precipvaluearray, pftvaluearray, deg=3)
            pcoeffvals = pcoeff.convert().coef
            if (len(pcoeffvals) >= 4):
                pftvalueprecippoly = pcoeffvals[0] + pcoeffvals[1] * precipannvalue + pcoeffvals[2] * pow(precipannvalue,2) + pcoeffvals[3] * pow(precipannvalue,3)
            if (len(pcoeffvals) == 3):
                pftvalueprecippoly = pcoeffvals[0] + pcoeffvals[1] * precipannvalue + pcoeffvals[2] * pow(precipannvalue,2)
            if (pftvalueprecippoly < 0.0):
                precipnear = 0
        else:
            pftvalueprecippoly = 0.0
        pftvaluepoly = 0.0
        if (tempnear == 1 and precipnear == 1):
            pftvaluepoly = (pftvaluetemppoly + pftvalueprecippoly) / 2.0
        if (tempnear == 0 and precipnear == 1): 
            pftvaluepoly = pftvalueprecippoly
        if (tempnear == 1 and precipnear == 0): 
            pftvaluepoly = pftvaluetemppoly
                
        poly1matrix[precipannindex,tempaverageindex] = pftvaluepoly
        if (poly1matrix[precipannindex,tempaverageindex] < 0.0):
            poly1matrix[precipannindex,tempaverageindex] = 0.0

for precipannindex in range(numprecipclasses):
    precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
    for tempaverageindex in range(numtempclasses):
        tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
        numtempvalues = 0
        tempnear = 0
        precipnear = 0
        for tempsearchindex in range(numtempclasses):
            if (poly1matrix[precipannindex,tempsearchindex] > 0.0):
                numtempvalues += 1
        if (numtempvalues > 1):
            tempvaluearray = np.zeros(shape=(numtempvalues),dtype=np.float64)
            pftvaluearray = np.zeros(shape=(numtempvalues),dtype=np.float64)
            tempvalueindex = 0
            for tempsearchindex in range(numtempclasses):
                if (poly1matrix[precipannindex,tempsearchindex] > 0.0):
                    tempvaluearray[tempvalueindex] = tempclassstart + float(tempsearchindex) * tempclasssize
                    pftvaluearray[tempvalueindex] = poly1matrix[precipannindex,tempsearchindex]
                    tempdistance = abs(tempsearchindex - tempaverageindex)
                    if (tempdistance < 6):
                        tempnear = 1
                    tempvalueindex += 1
            tcoeff = P.fit(tempvaluearray, pftvaluearray, deg=3)
            tcoeffvals = tcoeff.convert().coef
            if (len(tcoeffvals) >= 4):
                pftvaluetemppoly = tcoeffvals[0] + tcoeffvals[1] * tempaveragevalue + tcoeffvals[2] * pow(tempaveragevalue,2) + tcoeffvals[3] * pow(tempaveragevalue,3) 
            if (len(tcoeffvals) == 3):
                pftvaluetemppoly = tcoeffvals[0] + tcoeffvals[1] * tempaveragevalue + tcoeffvals[2] * pow(tempaveragevalue,2)
            if (pftvaluetemppoly < 0.0):
                tempnear = 0
        else:
            pftvaluetemppoly = 0.0
        numprecipvalues = 0
        for precipsearchindex in range(numprecipclasses):
            if (poly1matrix[precipsearchindex,tempaverageindex] > 0.0):
                numprecipvalues += 1
        if (numprecipvalues > 1):
            precipvaluearray = np.zeros(shape=(numprecipvalues),dtype=np.float64)
            pftvaluearray = np.zeros(shape=(numprecipvalues),dtype=np.float64)
            precipvalueindex = 0
            for precipsearchindex in range(numprecipclasses):
                if (poly1matrix[precipsearchindex,tempaverageindex] > 0.0):
                    precipvaluearray[precipvalueindex] = precipclassstart + float(precipsearchindex) * precipclasssize
                    pftvaluearray[precipvalueindex] = poly1matrix[precipsearchindex,tempaverageindex]
                    precipdistance = abs(precipsearchindex - precipannindex)
                    if (precipdistance < 6):
                        precipnear = 1
                    precipvalueindex += 1
            pcoeff = P.fit(precipvaluearray, pftvaluearray, deg=3)
            pcoeffvals = pcoeff.convert().coef
            if (len(pcoeffvals) >= 4):
                pftvalueprecippoly = pcoeffvals[0] + pcoeffvals[1] * precipannvalue + pcoeffvals[2] * pow(precipannvalue,2) + pcoeffvals[3] * pow(precipannvalue,3)
            if (len(pcoeffvals) == 3):
                pftvalueprecippoly = pcoeffvals[0] + pcoeffvals[1] * precipannvalue + pcoeffvals[2] * pow(precipannvalue,2)
            if (pftvalueprecippoly < 0.0):
                precipnear = 0
        else:
            pftvalueprecippoly = 0.0
        pftvaluepoly = 0.0
        if (tempnear == 1 and precipnear == 1):
            pftvaluepoly = (pftvaluetemppoly + pftvalueprecippoly) / 2.0
        if (tempnear == 0 and precipnear == 1): 
            pftvaluepoly = pftvalueprecippoly
        if (tempnear == 1 and precipnear == 0): 
            pftvaluepoly = pftvaluetemppoly
                
        pft1matrix[precipannindex,tempaverageindex] = pftvaluepoly
        if (pft1matrix[precipannindex,tempaverageindex] < 0.0):
            pft1matrix[precipannindex,tempaverageindex] = 0.0
        if (pft1matrix[precipannindex,tempaverageindex] > 100.0):
            pft1matrix[precipannindex,tempaverageindex] = 100.0

for precipannindex in range(numprecipclasses):
    for tempaverageindex in range(numtempclasses):
        landuseallfound = 0
        if (full1matrix[precipannindex,tempaverageindex] > 0.0):
            landuseallfound = 1
            
        if (landuseallfound == 0):
            if (precipannindex > 0):
                landuseleftaverage = np.average(full1matrix[0:precipannindex+1,tempaverageindex])
            else:
                landuseleftaverage = 0.0
            if (precipannindex < numprecipclasses-1):
                landuserightaverage = np.average(full1matrix[precipannindex:numprecipclasses,tempaverageindex])
            else:
                landuserightaverage = 0.0
            if (tempaverageindex > 0):
                landusebottomaverage = np.average(full1matrix[precipannindex,0:tempaverageindex+1])
            else:
                landuseleftaverage = 0.0
            if (tempaverageindex < numprecipclasses-1):
                landusetopaverage = np.average(full1matrix[precipannindex,tempaverageindex:numtempclasses])
            else:
                landuserightaverage = 0.0
            if (landuseleftaverage > 0.0 and landuserightaverage > 0.0 and landusebottomaverage > 0.0 and landusetopaverage > 0.0):
                landuseallfound = 1

        if (landuseallfound == 0):
            pft1matrix[precipannindex,tempaverageindex] = 0.0


for precipannindex in range(numprecipclasses):
    precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
    precipannstr = str(precipannvalue)
    for tempaverageindex in range(numtempclasses):
        tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
        tempaveragestr = str(tempaveragevalue)
        full2matrix[precipannindex,tempaverageindex] = readdatavalue(pft2filename,tempaveragestr,precipannstr)
        
for precipannindex in range(numprecipclasses):
    precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
    for tempaverageindex in range(numtempclasses):
        tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
        numtempvalues = 0
        tempnear = 0
        precipnear = 0
        for tempsearchindex in range(numtempclasses):
            if (full2matrix[precipannindex,tempsearchindex] > 0.0):
                numtempvalues += 1
        if (numtempvalues > 1):
            tempvaluearray = np.zeros(shape=(numtempvalues),dtype=np.float64)
            pftvaluearray = np.zeros(shape=(numtempvalues),dtype=np.float64)
            tempvalueindex = 0
            for tempsearchindex in range(numtempclasses):
                if (full2matrix[precipannindex,tempsearchindex] > 0.0):
                    tempvaluearray[tempvalueindex] = tempclassstart + float(tempsearchindex) * tempclasssize
                    pftvaluearray[tempvalueindex] = full2matrix[precipannindex,tempsearchindex]
                    tempdistance = abs(tempsearchindex - tempaverageindex)
                    if (tempdistance < 2):
                        tempnear = 1
                    tempvalueindex += 1
            tcoeff = P.fit(tempvaluearray, pftvaluearray, deg=3)
            tcoeffvals = tcoeff.convert().coef
            if (len(tcoeffvals) >= 4):
                pftvaluetemppoly = tcoeffvals[0] + tcoeffvals[1] * tempaveragevalue + tcoeffvals[2] * pow(tempaveragevalue,2) + tcoeffvals[3] * pow(tempaveragevalue,3) 
            if (len(tcoeffvals) == 3):
                pftvaluetemppoly = tcoeffvals[0] + tcoeffvals[1] * tempaveragevalue + tcoeffvals[2] * pow(tempaveragevalue,2)
            if (pftvaluetemppoly < 0.0):
                tempnear = 0
        else:
            pftvaluetemppoly = 0.0
        numprecipvalues = 0
        for precipsearchindex in range(numprecipclasses):
            if (full2matrix[precipsearchindex,tempaverageindex] > 0.0):
                numprecipvalues += 1
        if (numprecipvalues > 1):
            precipvaluearray = np.zeros(shape=(numprecipvalues),dtype=np.float64)
            pftvaluearray = np.zeros(shape=(numprecipvalues),dtype=np.float64)
            precipvalueindex = 0
            for precipsearchindex in range(numprecipclasses):
                if (full2matrix[precipsearchindex,tempaverageindex] > 0.0):
                    precipvaluearray[precipvalueindex] = precipclassstart + float(precipsearchindex) * precipclasssize
                    pftvaluearray[precipvalueindex] = full2matrix[precipsearchindex,tempaverageindex]
                    precipdistance = abs(precipsearchindex - precipannindex)
                    if (precipdistance < 3):
                        precipnear = 1
                    precipvalueindex += 1
            pcoeff = P.fit(precipvaluearray, pftvaluearray, deg=3)
            pcoeffvals = pcoeff.convert().coef
            if (len(pcoeffvals) >= 4):
                pftvalueprecippoly = pcoeffvals[0] + pcoeffvals[1] * precipannvalue + pcoeffvals[2] * pow(precipannvalue,2) + pcoeffvals[3] * pow(precipannvalue,3)
            if (len(pcoeffvals) == 3):
                pftvalueprecippoly = pcoeffvals[0] + pcoeffvals[1] * precipannvalue + pcoeffvals[2] * pow(precipannvalue,2)
            if (pftvalueprecippoly < 0.0):
                precipnear = 0
        else:
            pftvalueprecippoly = 0.0
        pftvaluepoly = 0.0
        if (tempnear == 1 and precipnear == 1):
            pftvaluepoly = (pftvaluetemppoly + pftvalueprecippoly) / 2.0
        if (tempnear == 0 and precipnear == 1): 
            pftvaluepoly = pftvalueprecippoly
        if (tempnear == 1 and precipnear == 0): 
            pftvaluepoly = pftvaluetemppoly
                
        poly2matrix[precipannindex,tempaverageindex] = pftvaluepoly
        if (poly2matrix[precipannindex,tempaverageindex] < 0.0):
            poly2matrix[precipannindex,tempaverageindex] = 0.0

for precipannindex in range(numprecipclasses):
    precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
    for tempaverageindex in range(numtempclasses):
        tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
        numtempvalues = 0
        tempnear = 0
        precipnear = 0
        for tempsearchindex in range(numtempclasses):
            if (poly2matrix[precipannindex,tempsearchindex] > 0.0):
                numtempvalues += 1
        if (numtempvalues > 1):
            tempvaluearray = np.zeros(shape=(numtempvalues),dtype=np.float64)
            pftvaluearray = np.zeros(shape=(numtempvalues),dtype=np.float64)
            tempvalueindex = 0
            for tempsearchindex in range(numtempclasses):
                if (poly2matrix[precipannindex,tempsearchindex] > 0.0):
                    tempvaluearray[tempvalueindex] = tempclassstart + float(tempsearchindex) * tempclasssize
                    pftvaluearray[tempvalueindex] = poly2matrix[precipannindex,tempsearchindex]
                    tempdistance = abs(tempsearchindex - tempaverageindex)
                    if (tempdistance < 6):
                        tempnear = 1
                    tempvalueindex += 1
            tcoeff = P.fit(tempvaluearray, pftvaluearray, deg=3)
            tcoeffvals = tcoeff.convert().coef
            if (len(tcoeffvals) >= 4):
                pftvaluetemppoly = tcoeffvals[0] + tcoeffvals[1] * tempaveragevalue + tcoeffvals[2] * pow(tempaveragevalue,2) + tcoeffvals[3] * pow(tempaveragevalue,3) 
            if (len(tcoeffvals) == 3):
                pftvaluetemppoly = tcoeffvals[0] + tcoeffvals[1] * tempaveragevalue + tcoeffvals[2] * pow(tempaveragevalue,2)
            if (pftvaluetemppoly < 0.0):
                tempnear = 0
        else:
            pftvaluetemppoly = 0.0
        numprecipvalues = 0
        for precipsearchindex in range(numprecipclasses):
            if (poly2matrix[precipsearchindex,tempaverageindex] > 0.0):
                numprecipvalues += 1
        if (numprecipvalues > 1):
            precipvaluearray = np.zeros(shape=(numprecipvalues),dtype=np.float64)
            pftvaluearray = np.zeros(shape=(numprecipvalues),dtype=np.float64)
            precipvalueindex = 0
            for precipsearchindex in range(numprecipclasses):
                if (poly2matrix[precipsearchindex,tempaverageindex] > 0.0):
                    precipvaluearray[precipvalueindex] = precipclassstart + float(precipsearchindex) * precipclasssize
                    pftvaluearray[precipvalueindex] = poly2matrix[precipsearchindex,tempaverageindex]
                    precipdistance = abs(precipsearchindex - precipannindex)
                    if (precipdistance < 6):
                        precipnear = 1
                    precipvalueindex += 1
            pcoeff = P.fit(precipvaluearray, pftvaluearray, deg=3)
            pcoeffvals = pcoeff.convert().coef
            if (len(pcoeffvals) >= 4):
                pftvalueprecippoly = pcoeffvals[0] + pcoeffvals[1] * precipannvalue + pcoeffvals[2] * pow(precipannvalue,2) + pcoeffvals[3] * pow(precipannvalue,3)
            if (len(pcoeffvals) == 3):
                pftvalueprecippoly = pcoeffvals[0] + pcoeffvals[1] * precipannvalue + pcoeffvals[2] * pow(precipannvalue,2)
            if (pftvalueprecippoly < 0.0):
                precipnear = 0
        else:
            pftvalueprecippoly = 0.0
        pftvaluepoly = 0.0
        if (tempnear == 1 and precipnear == 1):
            pftvaluepoly = (pftvaluetemppoly + pftvalueprecippoly) / 2.0
        if (tempnear == 0 and precipnear == 1): 
            pftvaluepoly = pftvalueprecippoly
        if (tempnear == 1 and precipnear == 0): 
            pftvaluepoly = pftvaluetemppoly
                
        pft2matrix[precipannindex,tempaverageindex] = pftvaluepoly
        if (pft2matrix[precipannindex,tempaverageindex] < 0.0):
            pft2matrix[precipannindex,tempaverageindex] = 0.0
        if (pft2matrix[precipannindex,tempaverageindex] > 100.0):
            pft2matrix[precipannindex,tempaverageindex] = 100.0

for precipannindex in range(numprecipclasses):
    for tempaverageindex in range(numtempclasses):
        landuseallfound = 0
        if (full2matrix[precipannindex,tempaverageindex] > 0.0):
            landuseallfound = 1
            
        if (landuseallfound == 0):
            if (precipannindex > 0):
                landuseleftaverage = np.average(full2matrix[0:precipannindex+1,tempaverageindex])
            else:
                landuseleftaverage = 0.0
            if (precipannindex < numprecipclasses-1):
                landuserightaverage = np.average(full2matrix[precipannindex:numprecipclasses,tempaverageindex])
            else:
                landuserightaverage = 0.0
            if (tempaverageindex > 0):
                landusebottomaverage = np.average(full2matrix[precipannindex,0:tempaverageindex+1])
            else:
                landuseleftaverage = 0.0
            if (tempaverageindex < numprecipclasses-1):
                landusetopaverage = np.average(full2matrix[precipannindex,tempaverageindex:numtempclasses])
            else:
                landuserightaverage = 0.0
            if (landuseleftaverage > 0.0 and landuserightaverage > 0.0 and landusebottomaverage > 0.0 and landusetopaverage > 0.0):
                landuseallfound = 1

        if (landuseallfound == 0):
            pft2matrix[precipannindex,tempaverageindex] = 0.0


for precipannindex in range(numprecipclasses):
    precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
    precipannstr = str(precipannvalue)
    for tempaverageindex in range(numtempclasses):
        tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
        tempaveragestr = str(tempaveragevalue)
        full3matrix[precipannindex,tempaverageindex] = readdatavalue(pft3filename,tempaveragestr,precipannstr)
        
for precipannindex in range(numprecipclasses):
    precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
    for tempaverageindex in range(numtempclasses):
        tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
        numtempvalues = 0
        tempnear = 0
        precipnear = 0
        for tempsearchindex in range(numtempclasses):
            if (full3matrix[precipannindex,tempsearchindex] > 0.0):
                numtempvalues += 1
        if (numtempvalues > 1):
            tempvaluearray = np.zeros(shape=(numtempvalues),dtype=np.float64)
            pftvaluearray = np.zeros(shape=(numtempvalues),dtype=np.float64)
            tempvalueindex = 0
            for tempsearchindex in range(numtempclasses):
                if (full3matrix[precipannindex,tempsearchindex] > 0.0):
                    tempvaluearray[tempvalueindex] = tempclassstart + float(tempsearchindex) * tempclasssize
                    pftvaluearray[tempvalueindex] = full3matrix[precipannindex,tempsearchindex]
                    tempdistance = abs(tempsearchindex - tempaverageindex)
                    if (tempdistance < 2):
                        tempnear = 1
                    tempvalueindex += 1
            tcoeff = P.fit(tempvaluearray, pftvaluearray, deg=3)
            tcoeffvals = tcoeff.convert().coef
            if (len(tcoeffvals) >= 4):
                pftvaluetemppoly = tcoeffvals[0] + tcoeffvals[1] * tempaveragevalue + tcoeffvals[2] * pow(tempaveragevalue,2) + tcoeffvals[3] * pow(tempaveragevalue,3) 
            if (len(tcoeffvals) == 3):
                pftvaluetemppoly = tcoeffvals[0] + tcoeffvals[1] * tempaveragevalue + tcoeffvals[2] * pow(tempaveragevalue,2)
            if (pftvaluetemppoly < 0.0):
                tempnear = 0
        else:
            pftvaluetemppoly = 0.0
        numprecipvalues = 0
        for precipsearchindex in range(numprecipclasses):
            if (full3matrix[precipsearchindex,tempaverageindex] > 0.0):
                numprecipvalues += 1
        if (numprecipvalues > 1):
            precipvaluearray = np.zeros(shape=(numprecipvalues),dtype=np.float64)
            pftvaluearray = np.zeros(shape=(numprecipvalues),dtype=np.float64)
            precipvalueindex = 0
            for precipsearchindex in range(numprecipclasses):
                if (full3matrix[precipsearchindex,tempaverageindex] > 0.0):
                    precipvaluearray[precipvalueindex] = precipclassstart + float(precipsearchindex) * precipclasssize
                    pftvaluearray[precipvalueindex] = full3matrix[precipsearchindex,tempaverageindex]
                    precipdistance = abs(precipsearchindex - precipannindex)
                    if (precipdistance < 3):
                        precipnear = 1
                    precipvalueindex += 1
            pcoeff = P.fit(precipvaluearray, pftvaluearray, deg=3)
            pcoeffvals = pcoeff.convert().coef
            if (len(pcoeffvals) >= 4):
                pftvalueprecippoly = pcoeffvals[0] + pcoeffvals[1] * precipannvalue + pcoeffvals[2] * pow(precipannvalue,2) + pcoeffvals[3] * pow(precipannvalue,3)
            if (len(pcoeffvals) == 3):
                pftvalueprecippoly = pcoeffvals[0] + pcoeffvals[1] * precipannvalue + pcoeffvals[2] * pow(precipannvalue,2)
            if (pftvalueprecippoly < 0.0):
                precipnear = 0
        else:
            pftvalueprecippoly = 0.0
        pftvaluepoly = 0.0
        if (tempnear == 1 and precipnear == 1):
            pftvaluepoly = (pftvaluetemppoly + pftvalueprecippoly) / 2.0
        if (tempnear == 0 and precipnear == 1): 
            pftvaluepoly = pftvalueprecippoly
        if (tempnear == 1 and precipnear == 0): 
            pftvaluepoly = pftvaluetemppoly
                
        poly3matrix[precipannindex,tempaverageindex] = pftvaluepoly
        if (poly3matrix[precipannindex,tempaverageindex] < 0.0):
            poly3matrix[precipannindex,tempaverageindex] = 0.0

for precipannindex in range(numprecipclasses):
    precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
    for tempaverageindex in range(numtempclasses):
        tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
        numtempvalues = 0
        tempnear = 0
        precipnear = 0
        for tempsearchindex in range(numtempclasses):
            if (poly3matrix[precipannindex,tempsearchindex] > 0.0):
                numtempvalues += 1
        if (numtempvalues > 1):
            tempvaluearray = np.zeros(shape=(numtempvalues),dtype=np.float64)
            pftvaluearray = np.zeros(shape=(numtempvalues),dtype=np.float64)
            tempvalueindex = 0
            for tempsearchindex in range(numtempclasses):
                if (poly3matrix[precipannindex,tempsearchindex] > 0.0):
                    tempvaluearray[tempvalueindex] = tempclassstart + float(tempsearchindex) * tempclasssize
                    pftvaluearray[tempvalueindex] = poly3matrix[precipannindex,tempsearchindex]
                    tempdistance = abs(tempsearchindex - tempaverageindex)
                    if (tempdistance < 6):
                        tempnear = 1
                    tempvalueindex += 1
            tcoeff = P.fit(tempvaluearray, pftvaluearray, deg=3)
            tcoeffvals = tcoeff.convert().coef
            if (len(tcoeffvals) >= 4):
                pftvaluetemppoly = tcoeffvals[0] + tcoeffvals[1] * tempaveragevalue + tcoeffvals[2] * pow(tempaveragevalue,2) + tcoeffvals[3] * pow(tempaveragevalue,3) 
            if (len(tcoeffvals) == 3):
                pftvaluetemppoly = tcoeffvals[0] + tcoeffvals[1] * tempaveragevalue + tcoeffvals[2] * pow(tempaveragevalue,2)
            if (pftvaluetemppoly < 0.0):
                tempnear = 0
        else:
            pftvaluetemppoly = 0.0
        numprecipvalues = 0
        for precipsearchindex in range(numprecipclasses):
            if (poly3matrix[precipsearchindex,tempaverageindex] > 0.0):
                numprecipvalues += 1
        if (numprecipvalues > 1):
            precipvaluearray = np.zeros(shape=(numprecipvalues),dtype=np.float64)
            pftvaluearray = np.zeros(shape=(numprecipvalues),dtype=np.float64)
            precipvalueindex = 0
            for precipsearchindex in range(numprecipclasses):
                if (poly3matrix[precipsearchindex,tempaverageindex] > 0.0):
                    precipvaluearray[precipvalueindex] = precipclassstart + float(precipsearchindex) * precipclasssize
                    pftvaluearray[precipvalueindex] = poly3matrix[precipsearchindex,tempaverageindex]
                    precipdistance = abs(precipsearchindex - precipannindex)
                    if (precipdistance < 6):
                        precipnear = 1
                    precipvalueindex += 1
            pcoeff = P.fit(precipvaluearray, pftvaluearray, deg=3)
            pcoeffvals = pcoeff.convert().coef
            if (len(pcoeffvals) >= 4):
                pftvalueprecippoly = pcoeffvals[0] + pcoeffvals[1] * precipannvalue + pcoeffvals[2] * pow(precipannvalue,2) + pcoeffvals[3] * pow(precipannvalue,3)
            if (len(pcoeffvals) == 3):
                pftvalueprecippoly = pcoeffvals[0] + pcoeffvals[1] * precipannvalue + pcoeffvals[2] * pow(precipannvalue,2)
            if (pftvalueprecippoly < 0.0):
                precipnear = 0
        else:
            pftvalueprecippoly = 0.0
        pftvaluepoly = 0.0
        if (tempnear == 1 and precipnear == 1):
            pftvaluepoly = (pftvaluetemppoly + pftvalueprecippoly) / 2.0
        if (tempnear == 0 and precipnear == 1): 
            pftvaluepoly = pftvalueprecippoly
        if (tempnear == 1 and precipnear == 0): 
            pftvaluepoly = pftvaluetemppoly
                
        pft3matrix[precipannindex,tempaverageindex] = pftvaluepoly
        if (pft3matrix[precipannindex,tempaverageindex] < 0.0):
            pft3matrix[precipannindex,tempaverageindex] = 0.0
        if (pft3matrix[precipannindex,tempaverageindex] > 100.0):
            pft3matrix[precipannindex,tempaverageindex] = 100.0

for precipannindex in range(numprecipclasses):
    for tempaverageindex in range(numtempclasses):
        landuseallfound = 0
        if (full3matrix[precipannindex,tempaverageindex] > 0.0):
            landuseallfound = 1
            
        if (landuseallfound == 0):
            if (precipannindex > 0):
                landuseleftaverage = np.average(full3matrix[0:precipannindex+1,tempaverageindex])
            else:
                landuseleftaverage = 0.0
            if (precipannindex < numprecipclasses-1):
                landuserightaverage = np.average(full3matrix[precipannindex:numprecipclasses,tempaverageindex])
            else:
                landuserightaverage = 0.0
            if (tempaverageindex > 0):
                landusebottomaverage = np.average(full3matrix[precipannindex,0:tempaverageindex+1])
            else:
                landuseleftaverage = 0.0
            if (tempaverageindex < numprecipclasses-1):
                landusetopaverage = np.average(full3matrix[precipannindex,tempaverageindex:numtempclasses])
            else:
                landuserightaverage = 0.0
            if (landuseleftaverage > 0.0 and landuserightaverage > 0.0 and landusebottomaverage > 0.0 and landusetopaverage > 0.0):
                landuseallfound = 1

        if (landuseallfound == 0):
            pft3matrix[precipannindex,tempaverageindex] = 0.0


climatemaskfilename = str(sys.argv[10])

climatemaskfile = open(climatemaskfilename,'r')
climatemasklist = climatemaskfile.readlines()
numclimatemasks = len(climatemasklist)
climatemaskids = [""] * numclimatemasks
climatemasktypes = [""] * numclimatemasks
climatemasklowtemps = [""] * numclimatemasks
climatemaskhightemps = [""] * numclimatemasks
climatemasklowprecips = [""] * numclimatemasks
climatemaskhighprecips = [""] * numclimatemasks
pft1climatemaskindex = -1
pft2climatemaskindex = -1
pft3climatemaskindex = -1
for climatemasklistindex in range(numclimatemasks):
    climatemasklistvalues = climatemasklist[climatemasklistindex].split()
    climatemaskids[climatemasklistindex] = climatemasklistvalues[0]
    climatemasktypes[climatemasklistindex] = climatemasklistvalues[1]
    climatemasklowtemps[climatemasklistindex] = climatemasklistvalues[2]
    climatemaskhightemps[climatemasklistindex] = climatemasklistvalues[3]
    climatemasklowprecips[climatemasklistindex] = climatemasklistvalues[4]
    climatemaskhighprecips[climatemasklistindex] = climatemasklistvalues[5]
    if (pft1climateid == climatemaskids[climatemasklistindex]):
        pft1climatemaskindex = climatemasklistindex
    if (pft2climateid == climatemaskids[climatemasklistindex]):
        pft2climatemaskindex = climatemasklistindex
    if (pft3climateid == climatemaskids[climatemasklistindex]):
        pft3climatemaskindex = climatemasklistindex

if (pft1climatemaskindex == -1):
    print("Error: PFT 1 Climate Mask Id not found " + pft1filename + " - " + pft1climateid)
    sys.exit()

pft1climatemasktype = climatemasktypes[pft1climatemaskindex]
pft1climatemasklowtemp = float(climatemasklowtemps[pft1climatemaskindex])
pft1climatemaskhightemp = float(climatemaskhightemps[pft1climatemaskindex])
pft1climatemasklowprecip = float(climatemasklowprecips[pft1climatemaskindex])
pft1climatemaskhighprecip = float(climatemaskhighprecips[pft1climatemaskindex])

if (pft1climatemasktype == "+"):
    for precipannindex in range(numprecipclasses):
        precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
        for tempaverageindex in range(numtempclasses):
            tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
            if (tempaveragevalue < pft1climatemasklowtemp or tempaveragevalue > pft1climatemaskhightemp or precipannvalue < pft1climatemasklowprecip or precipannvalue > pft1climatemaskhighprecip):
                pft1matrix[precipannindex,tempaverageindex] = 0.0

if (pft1climatemasktype == "-"):
    for precipannindex in range(numprecipclasses):
        precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
        for tempaverageindex in range(numtempclasses):
            tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
            if (tempaveragevalue >= pft1climatemasklowtemp and tempaveragevalue <= pft1climatemaskhightemp and precipannvalue >= pft1climatemasklowprecip and precipannvalue <= pft1climatemaskhighprecip):
                pft1matrix[precipannindex,tempaverageindex] = 0.0

if (pft2climatemaskindex == -1):
    print("Error: PFT 2 Climate Mask Id not found " + pft2filename + " - " + pft2climateid)
    sys.exit()

pft2climatemasktype = climatemasktypes[pft2climatemaskindex]
pft2climatemasklowtemp = float(climatemasklowtemps[pft2climatemaskindex])
pft2climatemaskhightemp = float(climatemaskhightemps[pft2climatemaskindex])
pft2climatemasklowprecip = float(climatemasklowprecips[pft2climatemaskindex])
pft2climatemaskhighprecip = float(climatemaskhighprecips[pft2climatemaskindex])

if (pft2climatemasktype == "+"):
    for precipannindex in range(numprecipclasses):
        precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
        for tempaverageindex in range(numtempclasses):
            tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
            if (tempaveragevalue < pft2climatemasklowtemp or tempaveragevalue > pft2climatemaskhightemp or precipannvalue < pft2climatemasklowprecip or precipannvalue > pft2climatemaskhighprecip):
                pft2matrix[precipannindex,tempaverageindex] = 0.0

if (pft2climatemasktype == "-"):
    for precipannindex in range(numprecipclasses):
        precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
        for tempaverageindex in range(numtempclasses):
            tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
            if (tempaveragevalue >= pft2climatemasklowtemp and tempaveragevalue <= pft2climatemaskhightemp and precipannvalue >= pft2climatemasklowprecip and precipannvalue <= pft2climatemaskhighprecip):
                pft2matrix[precipannindex,tempaverageindex] = 0.0

if (pft3climatemaskindex == -1):
    print("Error: PFT 3 Climate Mask Id not found " + pft3filename + " - " + pft3climateid)
    sys.exit()

pft3climatemasktype = climatemasktypes[pft3climatemaskindex]
pft3climatemasklowtemp = float(climatemasklowtemps[pft3climatemaskindex])
pft3climatemaskhightemp = float(climatemaskhightemps[pft3climatemaskindex])
pft3climatemasklowprecip = float(climatemasklowprecips[pft3climatemaskindex])
pft3climatemaskhighprecip = float(climatemaskhighprecips[pft3climatemaskindex])

if (pft3climatemasktype == "+"):
    for precipannindex in range(numprecipclasses):
        precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
        for tempaverageindex in range(numtempclasses):
            tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
            if (tempaveragevalue < pft3climatemasklowtemp or tempaveragevalue > pft3climatemaskhightemp or precipannvalue < pft3climatemasklowprecip or precipannvalue > pft3climatemaskhighprecip):
                pft3matrix[precipannindex,tempaverageindex] = 0.0

if (pft3climatemasktype == "-"):
    for precipannindex in range(numprecipclasses):
        precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
        for tempaverageindex in range(numtempclasses):
            tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
            if (tempaveragevalue >= pft3climatemasklowtemp and tempaveragevalue <= pft3climatemaskhightemp and precipannvalue >= pft3climatemasklowprecip and precipannvalue <= pft3climatemaskhighprecip):
                pft3matrix[precipannindex,tempaverageindex] = 0.0

for precipannindex in range(numprecipclasses):
    for tempaverageindex in range(numtempclasses):
        pft1val = pft1matrix[precipannindex,tempaverageindex]
        pft2val = pft2matrix[precipannindex,tempaverageindex]
        pft3val = pft3matrix[precipannindex,tempaverageindex]
        pftsum = pft1val + pft2val + pft3val
        if (pftsum <= 0.0):
            pft1normalized = 0.0
            pft2normalized = 0.0
            pft3normalized = 0.0
        else:
            if (pft1val >= 100.0):
                pft1normalized = 100.0
                pft2normalized = 0.0
                pft3normalized = 0.0
            else:
                if (pft1val + pft2val >= 100.0):
                    pft1normalized = pft1val
                    pft2normalized = 100.0 - pft1val 
                    pft3normalized = 0.0
                else:
                    if (pft1val + pft2val + pft3val >= 100.0):
                        pft1normalized = pft1val
                        pft2normalized = pft2val 
                        pft3normalized = 100.0 - (pft1val + pft2val)
                    else:
                        pft1normalized = pft1val / pftsum * 100.0
                        pft2normalized = pft2val / pftsum * 100.0
                        pft3normalized = pft3val / pftsum * 100.0
        pft1output[precipannindex,tempaverageindex] = pft1normalized
        pft2output[precipannindex,tempaverageindex] = pft2normalized
        pft3output[precipannindex,tempaverageindex] = pft3normalized
        

# write timeseries file

outputfile = open(outputcsv1filename,"w")
outputline = "precipann"
for tempaverageindex in range(numtempclasses):
    outputline = outputline + ","
    tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
    outputline = outputline + str(tempaveragevalue)
outputline = outputline + "\n"
outputfile.write(outputline)

for precipannindex in range(numprecipclasses):
    precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
    outputline = str(precipannvalue)
    for tempaverageindex in range(numtempclasses):
        outputline = outputline + ","
        outputline = outputline + outputstr.format(pft1output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

# write timeseries file

outputfile = open(outputcsv2filename,"w")
outputline = "precipann"
for tempaverageindex in range(numtempclasses):
    outputline = outputline + ","
    tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
    outputline = outputline + str(tempaveragevalue)
outputline = outputline + "\n"
outputfile.write(outputline)

for precipannindex in range(numprecipclasses):
    precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
    outputline = str(precipannvalue)
    for tempaverageindex in range(numtempclasses):
        outputline = outputline + ","
        outputline = outputline + outputstr.format(pft2output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

# write timeseries file

outputfile = open(outputcsv3filename,"w")
outputline = "precipann"
for tempaverageindex in range(numtempclasses):
    outputline = outputline + ","
    tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
    outputline = outputline + str(tempaveragevalue)
outputline = outputline + "\n"
outputfile.write(outputline)

for precipannindex in range(numprecipclasses):
    precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
    outputline = str(precipannvalue)
    for tempaverageindex in range(numtempclasses):
        outputline = outputline + ","
        outputline = outputline + outputstr.format(pft3output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

