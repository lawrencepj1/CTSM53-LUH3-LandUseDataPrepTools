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
pft1matrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft2matrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft3matrix = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft1output = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft2output = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft3output = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
output1filename = str(sys.argv[7]) 
output2filename = str(sys.argv[8]) 
output3filename = str(sys.argv[9]) 
outputstr = "{0:.3f}"

for precipannindex in range(numprecipclasses):
    precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
    precipannstr = str(precipannvalue)
    for tempaverageindex in range(numtempclasses):
        tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
        tempaveragestr = str(tempaveragevalue)
        pft1matrix[precipannindex,tempaverageindex] = readdatavalue(pft1filename,tempaveragestr,precipannstr)
        pft2matrix[precipannindex,tempaverageindex] = readdatavalue(pft2filename,tempaveragestr,precipannstr)
        pft3matrix[precipannindex,tempaverageindex] = readdatavalue(pft3filename,tempaveragestr,precipannstr)
        
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

outputfile = open(output1filename,"w")
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

outputfile = open(output2filename,"w")
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

outputfile = open(output3filename,"w")
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
