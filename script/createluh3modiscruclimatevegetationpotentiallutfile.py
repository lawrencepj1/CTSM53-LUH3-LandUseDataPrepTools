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

if (arguments != 8):
    print("Error: Usage createrawglobalmonthlyareasumdata.py landuse1filename landuse1weight landuse2filename landuse2weight landuse3filename landuse3weight landuseallfilename outputfilename")
    sys.exit()
else:
    print("Processing: " + str(sys.argv[1]) + " => " + str(sys.argv[7]))

landuse1filename = str(sys.argv[1])
landuse1weight = str(sys.argv[2])
landuse1weightvalue = float(landuse1weight)
landuse2filename = str(sys.argv[3])
landuse2weight = str(sys.argv[4])
landuse2weightvalue = float(landuse2weight)
landuse3filename = str(sys.argv[5])
landuse3weight = str(sys.argv[6])
landuse3weightvalue = float(landuse3weight)
landuseallfilename = str(sys.argv[7])

numprecipclasses = 100
precipclassstart = 0.0
precipclasssize = 50.0
numtempclasses = 100
tempclassstart = -15.0
tempclasssize = 0.5
landuse1matrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
landuse2matrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
landuse3matrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
landuseallmatrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
fullmatrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
polymatrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
finalmatrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
outputfilename = str(sys.argv[8]) 
outputstr = "{0:.3f}"

for precipannindex in range(numprecipclasses):
    precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
    precipannstr = str(precipannvalue)
    for tempaverageindex in range(numtempclasses):
        tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
        tempaveragestr = str(tempaveragevalue)
        landuse1matrix[precipannindex,tempaverageindex] = readdatavalue(landuse1filename,tempaveragestr,precipannstr) * landuse1weightvalue
        landuse2matrix[precipannindex,tempaverageindex] = readdatavalue(landuse2filename,tempaveragestr,precipannstr) * landuse2weightvalue
        landuse3matrix[precipannindex,tempaverageindex] = readdatavalue(landuse3filename,tempaveragestr,precipannstr) * landuse3weightvalue
        landuseallmatrix[precipannindex,tempaverageindex] = readdatavalue(landuseallfilename,tempaveragestr,precipannstr)

for precipannindex in range(numprecipclasses):
    for tempaverageindex in range(numtempclasses):
        if (landuse1matrix[precipannindex,tempaverageindex] > 0.0):
            fullmatrix[precipannindex,tempaverageindex] = landuse1matrix[precipannindex,tempaverageindex]
        else:
            if (landuse2matrix[precipannindex,tempaverageindex] > 0.0):
                fullmatrix[precipannindex,tempaverageindex] = landuse2matrix[precipannindex,tempaverageindex]
            else:
                if (landuse3matrix[precipannindex,tempaverageindex] > 0.0):
                    fullmatrix[precipannindex,tempaverageindex] = landuse3matrix[precipannindex,tempaverageindex]

for precipannindex in range(numprecipclasses):
    precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
    for tempaverageindex in range(numtempclasses):
        tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
        numtempvalues = 0
        tempnear = 0
        precipnear = 0
        for tempsearchindex in range(numtempclasses):
            if (fullmatrix[precipannindex,tempsearchindex] > 0.0):
                numtempvalues += 1
        if (numtempvalues > 1):
            tempvaluearray = np.zeros(shape=(numtempvalues),dtype=np.float64)
            pftvaluearray = np.zeros(shape=(numtempvalues),dtype=np.float64)
            tempvalueindex = 0
            for tempsearchindex in range(numtempclasses):
                if (fullmatrix[precipannindex,tempsearchindex] > 0.0):
                    tempvaluearray[tempvalueindex] = tempclassstart + float(tempsearchindex) * tempclasssize
                    pftvaluearray[tempvalueindex] = fullmatrix[precipannindex,tempsearchindex]
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
            if (fullmatrix[precipsearchindex,tempaverageindex] > 0.0):
                numprecipvalues += 1
        if (numprecipvalues > 1):
            precipvaluearray = np.zeros(shape=(numprecipvalues),dtype=np.float64)
            pftvaluearray = np.zeros(shape=(numprecipvalues),dtype=np.float64)
            precipvalueindex = 0
            for precipsearchindex in range(numprecipclasses):
                if (fullmatrix[precipsearchindex,tempaverageindex] > 0.0):
                    precipvaluearray[precipvalueindex] = precipclassstart + float(precipsearchindex) * precipclasssize
                    pftvaluearray[precipvalueindex] = fullmatrix[precipsearchindex,tempaverageindex]
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
                
        polymatrix[precipannindex,tempaverageindex] = pftvaluepoly
        if (polymatrix[precipannindex,tempaverageindex] < 0.0):
            polymatrix[precipannindex,tempaverageindex] = 0.0

for precipannindex in range(numprecipclasses):
    precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
    for tempaverageindex in range(numtempclasses):
        tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
        numtempvalues = 0
        tempnear = 0
        precipnear = 0
        for tempsearchindex in range(numtempclasses):
            if (polymatrix[precipannindex,tempsearchindex] > 0.0):
                numtempvalues += 1
        if (numtempvalues > 1):
            tempvaluearray = np.zeros(shape=(numtempvalues),dtype=np.float64)
            pftvaluearray = np.zeros(shape=(numtempvalues),dtype=np.float64)
            tempvalueindex = 0
            for tempsearchindex in range(numtempclasses):
                if (polymatrix[precipannindex,tempsearchindex] > 0.0):
                    tempvaluearray[tempvalueindex] = tempclassstart + float(tempsearchindex) * tempclasssize
                    pftvaluearray[tempvalueindex] = polymatrix[precipannindex,tempsearchindex]
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
            if (polymatrix[precipsearchindex,tempaverageindex] > 0.0):
                numprecipvalues += 1
        if (numprecipvalues > 1):
            precipvaluearray = np.zeros(shape=(numprecipvalues),dtype=np.float64)
            pftvaluearray = np.zeros(shape=(numprecipvalues),dtype=np.float64)
            precipvalueindex = 0
            for precipsearchindex in range(numprecipclasses):
                if (polymatrix[precipsearchindex,tempaverageindex] > 0.0):
                    precipvaluearray[precipvalueindex] = precipclassstart + float(precipsearchindex) * precipclasssize
                    pftvaluearray[precipvalueindex] = polymatrix[precipsearchindex,tempaverageindex]
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
                
        finalmatrix[precipannindex,tempaverageindex] = pftvaluepoly
        if (finalmatrix[precipannindex,tempaverageindex] < 0.0):
            finalmatrix[precipannindex,tempaverageindex] = 0.0
        if (finalmatrix[precipannindex,tempaverageindex] > 100.0):
            finalmatrix[precipannindex,tempaverageindex] = 100.0

for precipannindex in range(numprecipclasses):
    for tempaverageindex in range(numtempclasses):
        landuseallfound = 0
        if (landuseallmatrix[precipannindex,tempaverageindex] > 0.0):
            landuseallfound = 1
            
        if (landuseallfound == 0):
            if (precipannindex > 0):
                landuseleftaverage = np.average(landuseallmatrix[0:precipannindex+1,tempaverageindex])
            else:
                landuseleftaverage = 0.0
            if (precipannindex < numprecipclasses-1):
                landuserightaverage = np.average(landuseallmatrix[precipannindex:numprecipclasses,tempaverageindex])
            else:
                landuserightaverage = 0.0
            if (tempaverageindex > 0):
                landusebottomaverage = np.average(landuseallmatrix[precipannindex,0:tempaverageindex+1])
            else:
                landuseleftaverage = 0.0
            if (tempaverageindex < numprecipclasses-1):
                landusetopaverage = np.average(landuseallmatrix[precipannindex,tempaverageindex:numtempclasses])
            else:
                landuserightaverage = 0.0
            if (landuseleftaverage > 0.0 and landuserightaverage > 0.0 and landusebottomaverage > 0.0 and landusetopaverage > 0.0):
                landuseallfound = 1

        if (landuseallfound == 0):
            finalmatrix[precipannindex,tempaverageindex] = 0.0

# write timeseries file

outputfile = open(outputfilename,"w")
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
        outputline = outputline + outputstr.format(finalmatrix[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)
