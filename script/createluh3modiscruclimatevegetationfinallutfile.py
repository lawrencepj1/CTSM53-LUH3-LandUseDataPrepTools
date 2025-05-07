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

if (arguments != 9):
    print("Error: Usage createrawglobalmonthlyareasumdata.py pft1filename pft2filename pft3filename pft4filename pft5filename pft6filename pft7filename pft8filename pft9filename")
    sys.exit()
else:
    print("Processing: " + str(sys.argv[1]) + " => " + str(sys.argv[7]))

pft1filename = str(sys.argv[1]) + ".potentialrefitted.csv"
pft2filename = str(sys.argv[2]) + ".potentialrefitted.csv"
pft3filename = str(sys.argv[3]) + ".potentialrefitted.csv"
pft4filename = str(sys.argv[4]) + ".potentialrefitted.csv"
pft5filename = str(sys.argv[5]) + ".potentialrefitted.csv"
pft6filename = str(sys.argv[6]) + ".potentialrefitted.csv"
pft7filename = str(sys.argv[7]) + ".potentialrefitted.csv"
pft8filename = str(sys.argv[8]) + ".potentialrefitted.csv"
pft9filename = str(sys.argv[9]) + ".potentialrefitted.csv"

numprecipclasses = 100
precipclassstart = 0.0
precipclasssize = 50.0
numtempclasses = 100
tempclassstart = -15.0
tempclasssize = 0.5
pft1input = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft2input = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft3input = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft4input = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft5input = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft6input = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft7input = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft8input = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft9input = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft1output = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft2output = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft3output = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft4output = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft5output = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft6output = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft7output = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft8output = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
pft9output = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float)
outputcsv1filename = str(sys.argv[1]) + ".potentialfinal.csv"
outputcsv2filename = str(sys.argv[2]) + ".potentialfinal.csv"
outputcsv3filename = str(sys.argv[3]) + ".potentialfinal.csv"
outputcsv4filename = str(sys.argv[4]) + ".potentialfinal.csv"
outputcsv5filename = str(sys.argv[5]) + ".potentialfinal.csv"
outputcsv6filename = str(sys.argv[6]) + ".potentialfinal.csv"
outputcsv7filename = str(sys.argv[7]) + ".potentialfinal.csv"
outputcsv8filename = str(sys.argv[8]) + ".potentialfinal.csv"
outputcsv9filename = str(sys.argv[9]) + ".potentialfinal.csv"
outputtxt1filename = str(sys.argv[1]) + ".potentialfinal.txt"
outputtxt2filename = str(sys.argv[2]) + ".potentialfinal.txt"
outputtxt3filename = str(sys.argv[3]) + ".potentialfinal.txt"
outputtxt4filename = str(sys.argv[4]) + ".potentialfinal.txt"
outputtxt5filename = str(sys.argv[5]) + ".potentialfinal.txt"
outputtxt6filename = str(sys.argv[6]) + ".potentialfinal.txt"
outputtxt7filename = str(sys.argv[7]) + ".potentialfinal.txt"
outputtxt8filename = str(sys.argv[8]) + ".potentialfinal.txt"
outputtxt9filename = str(sys.argv[9]) + ".potentialfinal.txt"
outputstr = "{0:.3f}"


for precipannindex in range(numprecipclasses):
    precipannvalue = precipclassstart + float(precipannindex) * precipclasssize
    precipannstr = str(precipannvalue)
    for tempaverageindex in range(numtempclasses):
        tempaveragevalue = tempclassstart + float(tempaverageindex) * tempclasssize
        tempaveragestr = str(tempaveragevalue)
        pft1input[precipannindex,tempaverageindex] = readdatavalue(pft1filename,tempaveragestr,precipannstr)
        pft2input[precipannindex,tempaverageindex] = readdatavalue(pft2filename,tempaveragestr,precipannstr)
        pft3input[precipannindex,tempaverageindex] = readdatavalue(pft3filename,tempaveragestr,precipannstr)
        pft4input[precipannindex,tempaverageindex] = readdatavalue(pft4filename,tempaveragestr,precipannstr)
        pft5input[precipannindex,tempaverageindex] = readdatavalue(pft5filename,tempaveragestr,precipannstr)
        pft6input[precipannindex,tempaverageindex] = readdatavalue(pft6filename,tempaveragestr,precipannstr)
        pft7input[precipannindex,tempaverageindex] = readdatavalue(pft7filename,tempaveragestr,precipannstr)
        pft8input[precipannindex,tempaverageindex] = readdatavalue(pft8filename,tempaveragestr,precipannstr)
        pft9input[precipannindex,tempaverageindex] = readdatavalue(pft9filename,tempaveragestr,precipannstr)

for precipannindex in range(numprecipclasses):
    for tempaverageindex in range(numtempclasses):
        foresttreepct = pft1input[precipannindex,tempaverageindex]
        nonforesttreepct = pft4input[precipannindex,tempaverageindex]
        rangetreepct = pft7input[precipannindex,tempaverageindex]
        if (nonforesttreepct < rangetreepct):
            nonforesttreepct = rangetreepct
        if (foresttreepct < nonforesttreepct):
            foresttreepct = nonforesttreepct
        pft1output[precipannindex,tempaverageindex] = foresttreepct
        pft4output[precipannindex,tempaverageindex] = nonforesttreepct
        pft7output[precipannindex,tempaverageindex] = rangetreepct

for precipannindex in range(numprecipclasses):
    for tempaverageindex in range(numtempclasses):
        foresttreepct = pft1output[precipannindex,tempaverageindex]
        nonforesttreepct = pft4output[precipannindex,tempaverageindex]
        rangetreepct = pft7output[precipannindex,tempaverageindex]
        forestherbpct = pft2input[precipannindex,tempaverageindex]
        nonforestherbpct = pft5input[precipannindex,tempaverageindex]
        rangeherbpct = pft8input[precipannindex,tempaverageindex]
        if (rangeherbpct > 100.0 - rangetreepct):
            rangeherbpct = 100.0 - rangetreepct
        if (nonforestherbpct > rangeherbpct):
            nonforestherbpct = rangeherbpct
        if (nonforestherbpct > 100.0 - nonforesttreepct):
            nonforestherbpct = 100.0 - nonforesttreepct
        if (forestherbpct > nonforestherbpct):
            forestherbpct = nonforestherbpct
        if (forestherbpct > 100.0 - foresttreepct):
            forestherbpct = 100.0 - foresttreepct
        pft2output[precipannindex,tempaverageindex] = forestherbpct
        pft5output[precipannindex,tempaverageindex] = nonforestherbpct
        pft8output[precipannindex,tempaverageindex] = rangeherbpct

for precipannindex in range(numprecipclasses):
    for tempaverageindex in range(numtempclasses):
        foresttreepct = pft1output[precipannindex,tempaverageindex]
        nonforesttreepct = pft4output[precipannindex,tempaverageindex]
        rangetreepct = pft7output[precipannindex,tempaverageindex]
        forestherbpct = pft2output[precipannindex,tempaverageindex]
        nonforestherbpct = pft5output[precipannindex,tempaverageindex]
        rangeherbpct = pft8output[precipannindex,tempaverageindex]
        forestbarepct = pft3input[precipannindex,tempaverageindex]
        nonforestbarepct = pft6input[precipannindex,tempaverageindex]
        rangebarepct = pft9input[precipannindex,tempaverageindex]
        if (foresttreepct > 0.0 or forestherbpct > 0.0 or forestbarepct > 0.0):
            forestbarepct = 100.0 - foresttreepct - forestherbpct
        if (nonforesttreepct > 0.0 or nonforestherbpct > 0.0 or nonforestbarepct > 0.0):
            nonforestbarepct = 100.0 - nonforesttreepct - nonforestherbpct
        if (rangetreepct > 0.0 or rangeherbpct > 0.0 or rangebarepct > 0.0):
            rangebarepct = 100.0 - rangetreepct - rangeherbpct
        pft3output[precipannindex,tempaverageindex] = forestbarepct
        pft6output[precipannindex,tempaverageindex] = nonforestbarepct
        pft9output[precipannindex,tempaverageindex] = rangebarepct

# write timeseries file 1

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

# write timeseries file 1

outputfile = open(outputtxt1filename,"w")

for precipannindex in range(numprecipclasses):
    outputline = ""
    for tempaverageindex in range(numtempclasses):
        if (tempaverageindex != 0):
            outputline = outputline + " "
        outputline = outputline + outputstr.format(pft1output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

# write timeseries file 2

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

# write timeseries file 2

outputfile = open(outputtxt2filename,"w")

for precipannindex in range(numprecipclasses):
    outputline = ""
    for tempaverageindex in range(numtempclasses):
        if (tempaverageindex != 0):
            outputline = outputline + " "
        outputline = outputline + outputstr.format(pft2output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

# write timeseries file 3

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

# write timeseries file 3

outputfile = open(outputtxt3filename,"w")

for precipannindex in range(numprecipclasses):
    outputline = ""
    for tempaverageindex in range(numtempclasses):
        if (tempaverageindex != 0):
            outputline = outputline + " "
        outputline = outputline + outputstr.format(pft3output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

# write timeseries file 4

outputfile = open(outputcsv4filename,"w")
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
        outputline = outputline + outputstr.format(pft4output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

# write timeseries file 4

outputfile = open(outputtxt4filename,"w")

for precipannindex in range(numprecipclasses):
    outputline = ""
    for tempaverageindex in range(numtempclasses):
        if (tempaverageindex != 0):
            outputline = outputline + " "
        outputline = outputline + outputstr.format(pft4output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

# write timeseries file 5

outputfile = open(outputcsv5filename,"w")
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
        outputline = outputline + outputstr.format(pft5output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

# write timeseries file 5

outputfile = open(outputtxt5filename,"w")

for precipannindex in range(numprecipclasses):
    outputline = ""
    for tempaverageindex in range(numtempclasses):
        if (tempaverageindex != 0):
            outputline = outputline + " "
        outputline = outputline + outputstr.format(pft5output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

# write timeseries file 6

outputfile = open(outputcsv6filename,"w")
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
        outputline = outputline + outputstr.format(pft6output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

# write timeseries file 6

outputfile = open(outputtxt6filename,"w")

for precipannindex in range(numprecipclasses):
    outputline = ""
    for tempaverageindex in range(numtempclasses):
        if (tempaverageindex != 0):
            outputline = outputline + " "
        outputline = outputline + outputstr.format(pft6output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

# write timeseries file 7

outputfile = open(outputcsv7filename,"w")
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
        outputline = outputline + outputstr.format(pft7output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

# write timeseries file 7

outputfile = open(outputtxt7filename,"w")

for precipannindex in range(numprecipclasses):
    outputline = ""
    for tempaverageindex in range(numtempclasses):
        if (tempaverageindex != 0):
            outputline = outputline + " "
        outputline = outputline + outputstr.format(pft7output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

# write timeseries file 8

outputfile = open(outputcsv8filename,"w")
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
        outputline = outputline + outputstr.format(pft8output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

# write timeseries file 8

outputfile = open(outputtxt8filename,"w")

for precipannindex in range(numprecipclasses):
    outputline = ""
    for tempaverageindex in range(numtempclasses):
        if (tempaverageindex != 0):
            outputline = outputline + " "
        outputline = outputline + outputstr.format(pft8output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

# write timeseries file 9

outputfile = open(outputcsv9filename,"w")
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
        outputline = outputline + outputstr.format(pft9output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

# write timeseries file 9

outputfile = open(outputtxt9filename,"w")

for precipannindex in range(numprecipclasses):
    outputline = ""
    for tempaverageindex in range(numtempclasses):
        if (tempaverageindex != 0):
            outputline = outputline + " "
        outputline = outputline + outputstr.format(pft9output[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)

