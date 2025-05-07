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

# load proper modules first, i.e.
# cgd machines
'''
module load lang/python/2.7.14

'''

#--  end of function definitions  ---------------------------------

# process the input arguments

arguments = len(sys.argv) - 1

if (arguments != 10):
    print("Error: Usage createrawglobalmonthlyareasumdata.py experimentname yearname experimentlist luh3surfdir worksurfdir landusename landuselist pftname pftlist outputfilename")
    sys.exit()
else:
    print("Processing: " + str(sys.argv[1]) + " => " + str(sys.argv[8]))

experimentname = str(sys.argv[1])
yearname = str(sys.argv[2])
experimentlistfilename = str(sys.argv[3])
luh3surfdirname = str(sys.argv[4])
worksurfdirname = str(sys.argv[5])

filtersize = 0.8

experimentlistfile = open(experimentlistfilename,'r')
experimentlist = experimentlistfile.readlines()
numexperiments = len(experimentlist)
experimentnames = [""] * numexperiments
experimentyears = [""] * numexperiments
experimentfiletypes = [""] * numexperiments
experimentstartyears = [""] * numexperiments
experimentendyears = [""] * numexperiments
experimentdirs = [""] * numexperiments
experimentfiles = [""] * numexperiments
luh3fileindex = -1
modisfileindex = -1
climatefileindex = -1
for experimentlistindex in range(numexperiments):
    experimentlistvalues = experimentlist[experimentlistindex].split()
    experimentnames[experimentlistindex] = experimentlistvalues[0]
    experimentyears[experimentlistindex] = experimentlistvalues[1]
    experimentfiletypes[experimentlistindex] = experimentlistvalues[2]
    experimentstartyears[experimentlistindex] = experimentlistvalues[3]
    experimentendyears[experimentlistindex] = experimentlistvalues[4]
    experimentdirs[experimentlistindex] = experimentlistvalues[5]
    if (experimentdirs[experimentlistindex] == "<luh3surfdir>"):
        experimentdirs[experimentlistindex] = luh3surfdirname
    if (experimentdirs[experimentlistindex] == "<worksurfdir>"):
        experimentdirs[experimentlistindex] = worksurfdirname
    experimentfiles[experimentlistindex] = experimentlistvalues[6]
    if (experimentname == experimentnames[experimentlistindex] and yearname == experimentyears[experimentlistindex] and "LUH3" == experimentfiletypes[experimentlistindex]):
        luh3fileindex = experimentlistindex
    if (experimentname == experimentnames[experimentlistindex] and yearname == experimentyears[experimentlistindex] and "MODIS" == experimentfiletypes[experimentlistindex]):
        modisfileindex = experimentlistindex
    if (experimentname == experimentnames[experimentlistindex] and yearname == experimentyears[experimentlistindex] and "CLIMATE" == experimentfiletypes[experimentlistindex]):
        climatefileindex = experimentlistindex

if (luh3fileindex == -1):
    print("Error: Experiment and Year LUH3 file not found " + experimentname + " - " + yearname)
    sys.exit()

if (modisfileindex == -1):
    print("Error: Experiment and Year MODIS file not found " + experimentname + " - " + yearname)
    sys.exit()
    
if (climatefileindex == -1):
    print("Error: Experiment and Year CLIMATE file not found " + experimentname + " - " + yearname)
    sys.exit()
    
luh3filename = experimentdirs[luh3fileindex] + experimentfiles[luh3fileindex]
luh3year = int(experimentyears[luh3fileindex])
luh3startyear = int(experimentstartyears[luh3fileindex])
luh3endyear = int(experimentendyears[luh3fileindex])
luh3yearindex = luh3year - luh3startyear

modisfilename = experimentdirs[modisfileindex] + experimentfiles[modisfileindex]

climatefilename = experimentdirs[climatefileindex] + experimentfiles[climatefileindex]

# read landuse data from landuse file 

landusename = str(sys.argv[6])
landuselistfilename = str(sys.argv[7])

landuselistfile = open(landuselistfilename,'r')
landuselist = landuselistfile.readlines()
numlanduses = len(landuselist)
landusenames = [""] * numlanduses
landuseindex = -1
numlanduseids = -1
numlandusefilters = -1
allvegindex = -1
numallvegids = -1
numallvegfilters = -1
for landuselistindex in range(numlanduses):
    landuselistvalues = landuselist[landuselistindex].split()
    landusenames[landuselistindex] = landuselistvalues[0]
    if (landusename == landusenames[landuselistindex]):
        landuseindex = landuselistindex
        numlanduseids = int(landuselistvalues[1])
        landuseids = [""] * numlanduseids
        for landuseidindex in range(numlanduseids):
            landuseids[landuseidindex] = landuselistvalues[2+landuseidindex]
        numlandusefilters = int(landuselistvalues[2+numlanduseids])
        landusefilters = [""] * numlandusefilters
        for landuseidindex in range(numlandusefilters):
            landusefilters[landuseidindex] = landuselistvalues[3+numlanduseids+landuseidindex]
    if ("allveg" == landusenames[landuselistindex]):
        allvegindex = landuselistindex
        numallvegids = int(landuselistvalues[1])
        allvegids = [""] * numallvegids
        for allvegidindex in range(numallvegids):
            allvegids[allvegidindex] = landuselistvalues[2+allvegidindex]
        numallvegfilters = int(landuselistvalues[2+numallvegids])
        allvegfilters = [""] * numallvegfilters
        for allvegidindex in range(numallvegfilters):
            allvegfilters[allvegidindex] = landuselistvalues[3+numallvegids+allvegidindex]

if (landuseindex == -1):
    print("Error: Landuse not found " + landusename)
    sys.exit()

# read pft data from pft file 

pftname = str(sys.argv[8])
pftlistfilename = str(sys.argv[9])

pftlistfile = open(pftlistfilename,'r')
pftlist = pftlistfile.readlines()
numpfts = len(pftlist)
pftnames = [""] * numpfts
pftindex = -1
pftnumids = -1
for pftlistindex in range(numpfts):
    pftlistvalues = pftlist[pftlistindex].split()
    pftnames[pftlistindex] = pftlistvalues[0]
    if (pftname == pftnames[pftlistindex]):
        pftindex = pftlistindex
        numpftids = int(pftlistvalues[1])
        pftids = np.zeros(numpftids,dtype=int)
        for pftidindex in range(numpftids):
            pftids[pftidindex] = int(pftlistvalues[2+pftidindex])

if (pftindex == -1):
    print("Error: PFT not found " + pftname)
    sys.exit()


print("Reading File: " + luh3filename)
luh3file = netcdf4.Dataset(luh3filename,'r')

luh3landfrac=np.asfarray(luh3file.variables["LANDFRAC"][:,:],np.float32)
luh3landmask=np.asfarray(luh3file.variables["LANDMASK"][:,:],np.float32)
luh3area=np.asfarray(luh3file.variables["AREA"][:,:],np.float32)
luh3lon   = np.asfarray(luh3file.variables["LON"][:],np.float32)
luh3lat   = np.asfarray(luh3file.variables["LAT"][:],np.float32)
nluh3lat  = luh3lat.size
nluh3lon  = luh3lon.size
landusevariabletype = luh3file.variables["primf"].datatype
landusevariablefillvalue = luh3file.variables["primf"]._FillValue
landusevariablevalues = np.zeros(shape=(nluh3lat,nluh3lon),dtype=float)
for landuseidindex in range(numlanduseids):
    landuseid = landuseids[landuseidindex]
    landusevalues = np.asfarray(luh3file.variables[landuseid][luh3yearindex,:,:],landusevariabletype)
    landusevariablevalues += landusevalues
landusevariablevaluesarea = luh3landfrac * landusevariablevalues * luh3area
landusefiltervalues = np.zeros(shape=(nluh3lat,nluh3lon),dtype=float)
for landuseidindex in range(numlandusefilters):
    landuseid = landusefilters[landuseidindex]
    landusevalues = np.asfarray(luh3file.variables[landuseid][luh3yearindex,:,:],landusevariabletype)
    landusefiltervalues += landusevalues
landusefiltervaluesarea = luh3landfrac * landusevariablevalues * luh3area
luh3file.close()


print("Reading File: " + modisfilename)
modisfile = netcdf4.Dataset(modisfilename,'r')

modislandfrac=np.asfarray(modisfile.variables["LANDFRAC"][:,:],np.float32)
modislandmask=np.asfarray(modisfile.variables["LANDMASK"][:,:],np.float32)
modisarea=np.asfarray(modisfile.variables["AREA"][:,:],np.float32)
modislon   = np.asfarray(modisfile.variables["LON"][:],np.float32)
modislat   = np.asfarray(modisfile.variables["LAT"][:],np.float32)
nmodislat  = modislat.size
nmodislon  = modislon.size
modisvariabletype = modisfile.variables["PCT_NAT_PFT"].datatype
modisvariablefillvalue = modisfile.variables["PCT_NAT_PFT"]._FillValue
modisvariablevalues = np.zeros(shape=(nmodislat,nmodislon),dtype=float)
for pftidindex in range(numpftids):
    pftid = pftids[pftidindex]
    pftvalues = np.asfarray(modisfile.variables["PCT_NAT_PFT"][pftid,:,:],modisvariabletype)
    modisvariablevalues += pftvalues[0,:,:]
modisvariablefracvalues = modisvariablevalues / 100.0
modisfile.close()

print("Reading File: " + climatefilename)
climatefile = netcdf4.Dataset(climatefilename,'r')

climlandfrac=np.asfarray(climatefile.variables["LANDFRAC"][:,:],np.float32)
climlandmask=np.asfarray(climatefile.variables["LANDMASK"][:,:],np.float32)
climarea=np.asfarray(climatefile.variables["AREA"][:,:],np.float32)
climlon   = np.asfarray(climatefile.variables["LON"][:],np.float32)
climlat   = np.asfarray(climatefile.variables["LAT"][:],np.float32)
nclimlat  = climlat.size
nclimlon  = climlon.size
climprecipann = np.asfarray(climatefile.variables["precipann"][:,:],np.float32)
climtempaverage = np.asfarray(climatefile.variables["tempaverage"][:,:],np.float32)


numprecipclasses = 100
precipclassstart = 0.0
precipclasssize = 50.0
numtempclasses = 100
tempclassstart = -15.0
tempclasssize = 0.5
outresults = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float64)
outareasums = np.zeros(shape=(numprecipclasses,numtempclasses),dtype=np.float64)
outputfilename = str(sys.argv[10]) + ".csv"
outputstr = "{0:.3f}"

for latindex in range(nluh3lat):
    for lonindex in range(nluh3lon):
        landusefrac = luh3landfrac[latindex,lonindex]
        if (landusefrac > 0.0):
            landusegridarea = luh3area[latindex,lonindex]
            landuselandarea = landusefrac * landusegridarea
            landusefilterarea = landusefiltervaluesarea[latindex,lonindex]
            if (landusefilterarea > filtersize * landuselandarea):
                landusevariablearea = landusevariablevaluesarea[latindex,lonindex]
                modisvariablefrac = modisvariablefracvalues[latindex,lonindex]
                landusemodisvariablearea = modisvariablefrac * landusevariablearea
                currprecipann = climprecipann[latindex,lonindex]
                currtempaverage = climtempaverage[latindex,lonindex]
                precipannindex = int((currprecipann - precipclassstart) / precipclasssize)
                if (precipannindex < 0):
                    precipannindex = 0
                if (precipannindex >= numprecipclasses):
                    precipannindex = numprecipclasses - 1
                tempaverageindex = int((currtempaverage - tempclassstart) / tempclasssize)
                if (tempaverageindex < 0):
                    tempaverageindex = 0
                if (tempaverageindex >= numtempclasses):
                    tempaverageindex = numtempclasses - 1
                outresults[precipannindex,tempaverageindex] += landusemodisvariablearea
                outareasums[precipannindex,tempaverageindex] += landusevariablearea

for precipannindex in range(numprecipclasses):
    for tempaverageindex in range(numtempclasses):
        tempresult = outresults[precipannindex,tempaverageindex]
        temparea = outareasums[precipannindex,tempaverageindex]
        if (temparea > 0.0):
            tempresult = tempresult / temparea * 100.0
        else:
            tempresult = 0.0
        outresults[precipannindex,tempaverageindex] = tempresult

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
        outputline = outputline + outputstr.format(outresults[precipannindex,tempaverageindex])
    outputline = outputline + "\n"
    outputfile.write(outputline)
