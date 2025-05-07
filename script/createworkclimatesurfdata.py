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

# load proper modules first, i.e.
# cgd machines
'''
module load lang/python/2.7.14

'''

#--  end of function definitions  ---------------------------------

# process the input arguments

arguments = len(sys.argv) - 1

if (arguments != 6):
    print("Error: Usage createdatafile.py climatenamelistfile timeseriesrawfile timeseriesrawdir regionfile worksurffile worksurfdir")
    sys.exit()
else:
    print("Processing: " + str(sys.argv[1]))
    
pftnum = 15
cftnum = 64

climatenamelistfilename = str(sys.argv[1])
timeseriesrawfilename = str(sys.argv[2])
timeseriesrawdirname = str(sys.argv[3])
subregionfilename = str(sys.argv[4])
worksurffilename = str(sys.argv[5])
worksurfdirname = str(sys.argv[6])

climatenamelistfile = open(climatenamelistfilename)
climatenamelist = climatenamelistfile.readlines()
climateauthorparts = climatenamelist[0].split()
climateauthorpartnum = len(climateauthorparts)
climateauthor = ""
for climateauthorpartindex in range(climateauthorpartnum):
    if (climateauthorpartindex > 1):
        climateauthor = climateauthor + " "
    if (climateauthorpartindex > 0):
        climateauthor = climateauthor + climateauthorparts[climateauthorpartindex]
climatenameparts = climatenamelist[1].split()
climatename = climatenameparts[1]
worksurfoutputparts = climatenamelist[2].split()
worksurfoutput = worksurfoutputparts[1]
worksurfoutputnameparts = climatenamelist[3].split()
worksurfoutputname = worksurfoutputnameparts[1]
clm5modisrawinputparts = climatenamelist[4].split()
clm5modisrawinput = clm5modisrawinputparts[1]
clm5modisyearparts = climatenamelist[5].split()
clm5modisyearname = clm5modisyearparts[1]

timeseriesrawlistfile = open(timeseriesrawfilename,'r')
timeseriesrawlist = timeseriesrawlistfile.readlines()
numtimeseriesraw = len(timeseriesrawlist)
timeseriesrawnames = [""] * numtimeseriesraw
timeseriesrawdirs = [""] * numtimeseriesraw
timeseriesrawstartyears = [""] * numtimeseriesraw
timeseriesrawendyears = [""] * numtimeseriesraw
clm5modisrawindex = -1
for timeseriesrawlistindex in range(numtimeseriesraw):
    timeseriesrawlistvalues = timeseriesrawlist[timeseriesrawlistindex].split()
    timeseriesrawnames[timeseriesrawlistindex] = timeseriesrawlistvalues[0]
    timeseriesrawdirs[timeseriesrawlistindex] = timeseriesrawlistvalues[1]
    if (timeseriesrawdirs[timeseriesrawlistindex] == "<timeseriesrawdir>"):
        timeseriesrawdirs[timeseriesrawlistindex] = timeseriesrawdirname
    timeseriesrawstartyears[timeseriesrawlistindex] = timeseriesrawlistvalues[2]
    timeseriesrawendyears[timeseriesrawlistindex] = timeseriesrawlistvalues[3]
    if (clm5modisrawinput == timeseriesrawnames[timeseriesrawlistindex]):
        clm5modisrawindex = timeseriesrawlistindex

if (clm5modisrawindex == -1):
    print("Error: clm5modisrawinput not found " + clm5modisrawinput)
    sys.exit()
else:
    print("Climate Raw Input: " + timeseriesrawdirs[clm5modisrawindex] + timeseriesrawnames[clm5modisrawindex])


worksurflistfile = open(worksurffilename,'r')
worksurflist = worksurflistfile.readlines()
numworksurf = len(worksurflist)
worksurfnames = [""] * numworksurf
worksurfdirs = [""] * numworksurf
worksurfstartyears = [""] * numworksurf
worksurfendyears = [""] * numworksurf
worksurfoutputindex = -1
for worksurflistindex in range(numworksurf):
    worksurflistvalues = worksurflist[worksurflistindex].split()
    worksurfnames[worksurflistindex] = worksurflistvalues[0]
    worksurfdirs[worksurflistindex] = worksurflistvalues[1]
    if (worksurfdirs[worksurflistindex] == "<worksurfdir>"):
        worksurfdirs[worksurflistindex] = worksurfdirname
    worksurfstartyears[worksurflistindex] = worksurflistvalues[2]
    worksurfendyears[worksurflistindex] = worksurflistvalues[3]
    if (worksurfoutput == worksurfnames[worksurflistindex]):
        worksurfoutputindex = worksurflistindex

if (worksurfoutputindex == -1):
    print("Error: worksurfoutput not found " + worksurfoutput)
    sys.exit()
else:
    print("Climate Surface Output: " + worksurfdirs[worksurfoutputindex] + worksurfnames[worksurfoutputindex])

subregionfile = open(subregionfilename,'r')
subregionlist = subregionfile.readlines()

subregionlllon = float(subregionlist[0])
subregionlllat = float(subregionlist[1])
subregionurlon = float(subregionlist[2])
subregionurlat = float(subregionlist[3])
cellsize = float(subregionlist[4])
   
subregionloncells = int((subregionurlon - subregionlllon) / cellsize)
subregionlatcells = int((subregionurlat - subregionlllat) / cellsize)

lsmnatpftid = np.zeros(pftnum,dtype=int)
for pftindex in range(pftnum):
    lsmnatpftid[pftindex] = pftindex
    
lsmcftid = np.zeros(cftnum,dtype=int)
for cftindex in range(cftnum):
    lsmcftid[cftindex] = pftnum + cftindex
    
lsmedgen = subregionurlat
lsmedgee = subregionurlon
lsmedges = subregionlllat
lsmedgew = subregionlllon

lsmLON   = np.zeros(subregionloncells,dtype=np.float32)
lsmLONGXY = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)
lsmLAT   = np.zeros(subregionlatcells,dtype=np.float32)
lsmLATIXY = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

for lonindex in range(subregionloncells):
    loncount = float(lonindex)
    lonvalue = subregionlllon + cellsize / 2.0 + loncount * cellsize
    lsmLON[lonindex] = lonvalue

for latindex in range(subregionlatcells):
    latcount = float(latindex)
    latvalue = subregionlllat + cellsize / 2.0 + latcount * cellsize
    lsmLAT[latindex] = latvalue

for lonindex in range(subregionloncells):
    for latindex in range(subregionlatcells):
        lsmLONGXY[latindex,lonindex] = lsmLON[lonindex]
        lsmLATIXY[latindex,lonindex] = lsmLAT[latindex]

climyear = "CLIM"

currentdirname = timeseriesrawdirs[clm5modisrawindex] + timeseriesrawnames[clm5modisrawindex] + "/"
currentclm5modisname = timeseriesrawnames[clm5modisrawindex]

outputfilename = worksurfdirs[worksurfoutputindex] + worksurfnames[worksurfoutputindex] + "/" + worksurfoutputname

# process raw data files

lsmLANDMASK = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

LANDMASKfilename = currentdirname + currentclm5modisname + ".LANDMASK." + str(climyear) + ".dat"
print("Processing LANDMASK > " + LANDMASKfilename)

insurffloatdata = np.fromfile(LANDMASKfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmLANDMASK[:] = insurffloatgrid[::-1,:]

lsmLANDFRAC = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

LANDFRACfilename = currentdirname + currentclm5modisname + ".LANDFRAC." + str(climyear) + ".dat"
print("Processing LANDFRAC > " + LANDFRACfilename)

insurffloatdata = np.fromfile(LANDFRACfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmLANDFRAC[:] = insurffloatgrid[::-1,:]

lsmAREA = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

AREAfilename = currentdirname + currentclm5modisname + ".AREA." + str(climyear) + ".dat"
print("Processing AREA > " + AREAfilename)

insurffloatdata = np.fromfile(AREAfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmAREA[:] = insurffloatgrid[::-1,:]

lsmprecipann = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

precipannfilename = currentdirname + currentclm5modisname + ".precipann-" + str(climyear) + ".dat"
print("Processing precipann > " + precipannfilename)

insurffloatdata = np.fromfile(precipannfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmprecipann[:] = insurffloatgrid[::-1,:]

nullindex = np.where(lsmprecipann == 90000.0)
lsmprecipann[nullindex] = 0.0

lsmprecipmaxtg22 = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

precipmaxtg22filename = currentdirname + currentclm5modisname + ".precipmaxtg22-" + str(climyear) + ".dat"
print("Processing precipmaxtg22 > " + precipmaxtg22filename)

insurffloatdata = np.fromfile(precipmaxtg22filename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmprecipmaxtg22[:] = insurffloatgrid[::-1,:]

nullindex = np.where(lsmprecipmaxtg22 == 90000.0)
lsmprecipmaxtg22[nullindex] = 0.0

lsmprecipmin = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

precipminfilename = currentdirname + currentclm5modisname + ".precipmin-" + str(climyear) + ".dat"
print("Processing precipmin > " + precipminfilename)

insurffloatdata = np.fromfile(precipminfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmprecipmin[:] = insurffloatgrid[::-1,:]

nullindex = np.where(lsmprecipmin == 90000.0)
lsmprecipmin[nullindex] = 0.0

lsmprecipwin = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

precipwinfilename = currentdirname + currentclm5modisname + ".precipwin-" + str(climyear) + ".dat"
print("Processing precipwin > " + precipwinfilename)

insurffloatdata = np.fromfile(precipwinfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmprecipwin[:] = insurffloatgrid[::-1,:]

nullindex = np.where(lsmprecipwin == 90000.0)
lsmprecipwin[nullindex] = 0.0

lsmtempaverage = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

tempaveragefilename = currentdirname + currentclm5modisname + ".tempaverage-" + str(climyear) + ".dat"
print("Processing tempaverage > " + tempaveragefilename)

insurffloatdata = np.fromfile(tempaveragefilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmtempaverage[:] = insurffloatgrid[::-1,:]

nullindex = np.where(lsmtempaverage == 90000.0)
lsmtempaverage[nullindex] = 0.0

lsmtempcoldest = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

tempcoldestfilename = currentdirname + currentclm5modisname + ".tempcoldest-" + str(climyear) + ".dat"
print("Processing tempcoldest > " + tempcoldestfilename)

insurffloatdata = np.fromfile(tempcoldestfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmtempcoldest[:] = insurffloatgrid[::-1,:]

nullindex = np.where(lsmtempcoldest == 90000.0)
lsmtempcoldest[nullindex] = 0.0

lsmtempwarmest = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

tempwarmestfilename = currentdirname + currentclm5modisname + ".tempwarmest-" + str(climyear) + ".dat"
print("Processing tempwarmest > " + tempwarmestfilename)

insurffloatdata = np.fromfile(tempwarmestfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmtempwarmest[:] = insurffloatgrid[::-1,:]

nullindex = np.where(lsmtempwarmest == 90000.0)
lsmtempwarmest[nullindex] = 0.0


# generate timeseries data

print('Creating: ' + outputfilename)
outputfile = netcdf4.Dataset(outputfilename, 'w')

outputfile.Conventions = 'NCAR-CSM'
outputfile.Author = 'Peter Lawrence, Terrestrial Sciences Section, National Center for Atmospheric Research'
datenow = date.datetime.now()
datenowstr = datenow.strftime("%m-%d-%Y %H:%M:%S")
outputfile.History_Log = 'created on: ' + datenowstr
outputfile.Data_Log = currentdirname

outputfile.createDimension('lat',int(subregionlatcells))
outputfile.createDimension('lon',int(subregionloncells))

wEDGEN = outputfile.createVariable('EDGEN',np.float32)
wEDGEE = outputfile.createVariable('EDGEE',np.float32)
wEDGES = outputfile.createVariable('EDGES',np.float32)
wEDGEW = outputfile.createVariable('EDGEW',np.float32)
wLAT  = outputfile.createVariable('LAT',np.float32,('lat',))
wLATIXY = outputfile.createVariable('LATIXY',np.float32,('lat','lon'))
wLON  = outputfile.createVariable('LON',np.float32,('lon',))
wLONGXY  = outputfile.createVariable('LONGXY',np.float32,('lat','lon'))
wLANDMASK = outputfile.createVariable('LANDMASK',np.float32,('lat','lon'))
wLANDFRAC = outputfile.createVariable('LANDFRAC',np.float32,('lat','lon'))
wAREA = outputfile.createVariable('AREA',np.float32,('lat','lon'))
wprecipann = outputfile.createVariable('precipann',np.float32,('lat','lon'),fill_value=-9999.0)
wprecipmaxtg22 = outputfile.createVariable('precipmaxtg22',np.float32,('lat','lon'),fill_value=-9999.0)
wprecipmin = outputfile.createVariable('precipmin',np.float32,('lat','lon'),fill_value=-9999.0)
wprecipwin = outputfile.createVariable('precipwin',np.float32,('lat','lon'),fill_value=-9999.0)
wtempaverage = outputfile.createVariable('tempaverage',np.float32,('lat','lon'),fill_value=-9999.0)
wtempcoldest = outputfile.createVariable('tempcoldest',np.float32,('lat','lon'),fill_value=-9999.0)
wtempwarmest = outputfile.createVariable('tempwarmest',np.float32,('lat','lon'),fill_value=-9999.0)

wEDGEN[...] = float(lsmedgen)
wEDGEN.long_name = 'northern edge of surface grid'
wEDGEN.units = 'degrees north'

wEDGEE[...] = float(lsmedgee)
wEDGEE.long_name = 'eastern edge of surface grid'
wEDGEE.units = 'degrees east'

wEDGES[...] = float(lsmedges)
wEDGES.long_name = 'southern edge of surface grid'
wEDGES.units = 'degrees north'

wEDGEW[...] = float(lsmedgew)
wEDGEW.long_name = 'western edge of surface grid'
wEDGEW.units = 'degrees east'

wLAT[:] = lsmLAT
wLAT.long_name = 'lat'
wLAT.units = 'degrees north'

wLATIXY[:] = lsmLATIXY
wLATIXY.long_name = 'latitude-2d'
wLATIXY.units = 'degrees north'

wLON[:] = lsmLON
wLON.long_name = 'lon'
wLON.units = 'degrees east'

wLONGXY[:] = lsmLONGXY
wLONGXY.long_name = 'longitude-2d'
wLONGXY.units = 'degrees east'

wLANDMASK[:] = lsmLANDMASK
wLANDMASK.long_name = 'and mask'
wLANDMASK.units = 'unitless'

wLANDFRAC[:] = lsmLANDFRAC
wLANDFRAC.long_name = 'land fraction of gridcell'
wLANDFRAC.units = 'unitless'

wAREA[:] = lsmAREA
wAREA.long_name = 'area of gridcell'
wAREA.units = 'km^2'

wprecipann[:] = lsmprecipann
wprecipann.long_name = 'average annual precipitation'
wprecipann.units = 'mm'

wprecipmaxtg22[:] = lsmprecipmaxtg22
wprecipmaxtg22.long_name = 'average precipitation for months with maximum temperature warmer than 22C'
wprecipmaxtg22.units = 'mm'

wprecipmin[:] = lsmprecipmin
wprecipmin.long_name = 'average minimum monthly precipitation'
wprecipmin.units = 'mm'

wprecipwin[:] = lsmprecipwin
wprecipwin.long_name = 'average winter precipitation'
wprecipwin.units = 'mm'

wtempaverage[:] = lsmtempaverage
wtempaverage.long_name = 'average annual temperature'
wtempaverage.units = 'C'

wtempcoldest[:] = lsmtempcoldest
wtempcoldest.long_name = 'average coldest month temperature'
wtempcoldest.units = 'C'

wtempwarmest[:] = lsmtempwarmest
wtempwarmest.long_name = 'average warmest month temperature'
wtempwarmest.units = 'C'
