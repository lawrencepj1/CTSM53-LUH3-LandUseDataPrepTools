#! /usr/bin/env python
import sys
import os.path
import math
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
    print("Error: Usage createdatafile.py luh2descnamelistfile timeseriesrawfilename timeseriesrawdirname regionfile timeseriessurffile timeseriessurfdirname")
    sys.exit()
else:
    print("Processing: " + str(sys.argv[1]))
    
pftnum = 15
cftnum = 64
monthnum = 12
monthdays = [ 14.0, 45.0, 73.0, 104.0, 134.0, 165.0, 195.0, 226.0, 257.0, 287.0, 318.0, 348.0 ]
startyear = 2005

luh2descnamelistfilename = str(sys.argv[1])
timeseriesrawfilename = str(sys.argv[2])
timeseriesrawdirname = str(sys.argv[3])
subregionfilename = str(sys.argv[4])
timeseriessurffilename = str(sys.argv[5])
timeseriessurfdirname = str(sys.argv[6])

luh2descnamelistfile = open(luh2descnamelistfilename)
luh2descnamelist = luh2descnamelistfile.readlines()
luh2descauthorparts = luh2descnamelist[0].split()
luh2descauthorpartnum = len(luh2descauthorparts)
luh2descauthor = ""
for luh2descauthorpartindex in range(luh2descauthorpartnum):
    if (luh2descauthorpartindex > 1):
        luh2descauthor = luh2descauthor + " "
    if (luh2descauthorpartindex > 0):
        luh2descauthor = luh2descauthor + luh2descauthorparts[luh2descauthorpartindex]
luh2descnameparts = luh2descnamelist[1].split()
luh2descname = luh2descnameparts[1]
timeseriessurfoutputparts = luh2descnamelist[2].split()
timeseriessurfoutput = timeseriessurfoutputparts[1]
timeseriessurfoutputnameparts = luh2descnamelist[3].split()
timeseriessurfoutputname = timeseriessurfoutputnameparts[1]
clm5modisrawinputparts = luh2descnamelist[4].split()
clm5modisrawinput = clm5modisrawinputparts[1]
clm5modisrawinputyearnameparts = luh2descnamelist[5].split()
clm5modisrawinputyearname = clm5modisrawinputyearnameparts[1]

print("author: " + luh2descauthor)
print("luh2descname: " + luh2descname)
print("timeseriessurfoutput: " + timeseriessurfoutput)
print("timeseriessurfoutputname: " + timeseriessurfoutputname)
print("clm5modisrawinput: " + clm5modisrawinput)

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
    print("Raw Input: " + timeseriesrawdirs[clm5modisrawindex] + timeseriesrawnames[clm5modisrawindex])

timeseriessurflistfile = open(timeseriessurffilename,'r')
timeseriessurflist = timeseriessurflistfile.readlines()
numtimeseriessurf = len(timeseriessurflist)
timeseriessurfnames = [""] * numtimeseriessurf
timeseriessurfdirs = [""] * numtimeseriessurf
timeseriessurfstartyears = [""] * numtimeseriessurf
timeseriessurfendyears = [""] * numtimeseriessurf
timeseriessurfoutputindex = -1
for timeseriessurflistindex in range(numtimeseriessurf):
    timeseriessurflistvalues = timeseriessurflist[timeseriessurflistindex].split()
    timeseriessurfnames[timeseriessurflistindex] = timeseriessurflistvalues[0]
    timeseriessurfdirs[timeseriessurflistindex] = timeseriessurflistvalues[1]
    if (timeseriessurfdirs[timeseriessurflistindex] == "<timeseriessurfdir>"):
        timeseriessurfdirs[timeseriessurflistindex] = timeseriessurfdirname
    timeseriessurfstartyears[timeseriessurflistindex] = timeseriessurflistvalues[2]
    timeseriessurfendyears[timeseriessurflistindex] = timeseriessurflistvalues[3]
    if (timeseriessurfoutput == timeseriessurfnames[timeseriessurflistindex]):
        timeseriessurfoutputindex = timeseriessurflistindex

if (timeseriessurfoutputindex == -1):
    print("Error: timeseriessurfoutput not found " + timeseriessurfoutput)
    sys.exit()
else:
    print("Surface Output: " + timeseriessurfdirs[timeseriessurfoutputindex] + timeseriessurfnames[timeseriessurfoutputindex])

subregionfilename = str(sys.argv[4])
subregionfile = open(subregionfilename,'r')
subregionlist = subregionfile.readlines()

subregionlllon = float(subregionlist[0])
subregionlllat = float(subregionlist[1])
subregionurlon = float(subregionlist[2])
subregionurlat = float(subregionlist[3])
cellsize = float(subregionlist[4])

subregionloncells = int((subregionurlon - subregionlllon) / cellsize)
subregionlatcells = int((subregionurlat - subregionlllat) / cellsize)

lsmpftid = np.zeros(pftnum,dtype=int)
for pftindex in range(pftnum):
    lsmpftid[pftindex] = pftindex
    
lsmtime = np.zeros(monthnum,dtype=float)

for timeindex in range(monthnum):
    lsmtime[timeindex] = monthdays[timeindex]
      
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
        
clm5modisrawdirname = timeseriesrawdirs[clm5modisrawindex] + timeseriesrawnames[clm5modisrawindex] + "/"

outputfilename = timeseriessurfdirs[timeseriessurfoutputindex] + timeseriessurfnames[timeseriessurfoutputindex] + "/" + timeseriessurfoutputname

# process raw data files

lsmLANDMASK = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

LANDMASKfilename = clm5modisrawdirname + clm5modisrawinput + ".LANDMASK.CLIM.dat"
print("Processing LANDMASK > " + LANDMASKfilename)

insurffloatdata = np.fromfile(LANDMASKfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmLANDMASK[:] = insurffloatgrid[::-1,:]

lsmLANDFRAC = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

LANDFRACfilename = clm5modisrawdirname + clm5modisrawinput + ".LANDFRAC.CLIM.dat"
print("Processing LANDFRAC > " + LANDFRACfilename)

insurffloatdata = np.fromfile(LANDFRACfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmLANDFRAC[:] = insurffloatgrid[::-1,:]

lsmAREA = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

AREAfilename = clm5modisrawdirname + clm5modisrawinput + ".AREA.CLIM.dat"
print("Processing AREA > " + AREAfilename)

insurffloatdata = np.fromfile(AREAfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmAREA[:] = insurffloatgrid[::-1,:]

lsmMONTHLYHEIGHTBOT = np.zeros(shape=(monthnum,pftnum,subregionlatcells,subregionloncells),dtype=np.float32)

for pftindex in range(1,pftnum):

    if (pftindex < 10):
        MONTHLYHEIGHTBOTfilename = clm5modisrawdirname + clm5modisrawinput + ".PFTCANOPYBOT0" + str(pftindex) + ".CLIM.dat"
    else:
        MONTHLYHEIGHTBOTfilename = clm5modisrawdirname + clm5modisrawinput + ".PFTCANOPYBOT" + str(pftindex) + ".CLIM.dat"
    
    print("Processing MONTHLYHEIGHTBOT > " + MONTHLYHEIGHTBOTfilename)

    insurffloatdata = np.fromfile(MONTHLYHEIGHTBOTfilename, dtype="float32", count=-1)
    insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
    
    for timeindex in range(monthnum):
        lsmMONTHLYHEIGHTBOT[timeindex,pftindex,:,:] = insurffloatgrid[::-1,:]

lsmMONTHLYHEIGHTTOP = np.zeros(shape=(monthnum,pftnum,subregionlatcells,subregionloncells),dtype=np.float32)

for pftindex in range(1,pftnum):

    if (pftindex < 10):
        MONTHLYHEIGHTTOPfilename = clm5modisrawdirname + clm5modisrawinput + ".PFTCANOPYTOP0" + str(pftindex) + ".CLIM.dat"
    else:
        MONTHLYHEIGHTTOPfilename = clm5modisrawdirname + clm5modisrawinput + ".PFTCANOPYTOP" + str(pftindex) + ".CLIM.dat"
    
    print("Processing MONTHLYHEIGHTTOP > " + MONTHLYHEIGHTTOPfilename)

    insurffloatdata = np.fromfile(MONTHLYHEIGHTTOPfilename, dtype="float32", count=-1)
    insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
    
    for timeindex in range(monthnum):
        lsmMONTHLYHEIGHTTOP[timeindex,pftindex,:,:] = insurffloatgrid[::-1,:]

lsmMONTHLYLAI = np.zeros(shape=(monthnum,pftnum,subregionlatcells,subregionloncells),dtype=np.float32)

for pftindex in range(1,pftnum):

    for timeindex in range(monthnum):
        monthvalue = timeindex + 1
        if (monthvalue < 10):
            monthstr = "0" + str(monthvalue)
        else:
            monthstr = str(monthvalue)
    
        if (pftindex < 10):
            MONTHLYLAIfilename = clm5modisrawdirname + clm5modisrawinput + ".PFTLAI0" + str(pftindex) + "." + monthstr + ".dat"
        else:
            MONTHLYLAIfilename = clm5modisrawdirname + clm5modisrawinput + ".PFTLAI" + str(pftindex) + "." + monthstr + ".dat"
    
        print("Processing MONTHLYLAI > " + MONTHLYLAIfilename)

        insurffloatdata = np.fromfile(MONTHLYLAIfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
    
        lsmMONTHLYLAI[timeindex,pftindex,:,:] = insurffloatgrid[::-1,:]

lsmMONTHLYSAI = np.zeros(shape=(monthnum,pftnum,subregionlatcells,subregionloncells),dtype=np.float32)

for pftindex in range(1,pftnum):

    for timeindex in range(monthnum):
        monthvalue = timeindex + 1
        if (monthvalue < 10):
            monthstr = "0" + str(monthvalue)
        else:
            monthstr = str(monthvalue)
    
        if (pftindex < 10):
            MONTHLYSAIfilename = clm5modisrawdirname + clm5modisrawinput + ".PFTSAI0" + str(pftindex) + "." + monthstr + ".dat"
        else:
            MONTHLYSAIfilename = clm5modisrawdirname + clm5modisrawinput + ".PFTSAI" + str(pftindex) + "." + monthstr + ".dat"
    
        print("Processing MONTHLYSAI > " + MONTHLYSAIfilename)

        insurffloatdata = np.fromfile(MONTHLYSAIfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
    
        lsmMONTHLYSAI[timeindex,pftindex,:,:] = insurffloatgrid[::-1,:]

# generate timeseries data

print('Creating: ' + outputfilename)
outputfile = netcdf4.Dataset(outputfilename, 'w', format='NETCDF3_64BIT_DATA')

outputfile.Conventions = 'NCAR-CSM'
outputfile.Author = 'Peter Lawrence, Terrestrial Sciences Section, National Center for Atmospheric Research'
datenow = date.datetime.now()
datenowstr = datenow.strftime("%m-%d-%Y %H:%M:%S")
outputfile.History_Log = 'created on: ' + datenowstr
outputfile.Data_Log = clm5modisrawdirname

outputfile.createDimension('pft',int(pftnum))
outputfile.createDimension('time',None)
outputfile.createDimension('lat',int(subregionlatcells))
outputfile.createDimension('lon',int(subregionloncells))

wpft = outputfile.createVariable('pft','i',('pft',))
wtime = outputfile.createVariable('time',np.float32,('time',))
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
wMONTHLYHEIGHTBOT = outputfile.createVariable('MONTHLY_HEIGHT_BOT',np.float32,('time','pft','lat','lon'),fill_value=-9999.0)
wMONTHLYHEIGHTTOP = outputfile.createVariable('MONTHLY_HEIGHT_TOP',np.float32,('time','pft','lat','lon'),fill_value=-9999.0)
wMONTHLYLAI = outputfile.createVariable('MONTHLY_LAI',np.float32,('time','pft','lat','lon'),fill_value=-9999.0)
wMONTHLYSAI = outputfile.createVariable('MONTHLY_SAI',np.float32,('time','pft','lat','lon'),fill_value=-9999.0)

wpft[:] = lsmpftid
wpft.long_name = 'indices of natural PFTs'
wpft.units = 'index'

wtime[:] = lsmtime
wtime.long_name = 'time'
wtime.bounds = 'time_bnds'
wtime.calendar = 'noleap'
wtime.units = 'days since ' + str(startyear) + '-01-01 00:00:00'

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

wMONTHLYHEIGHTBOT[:] = lsmMONTHLYHEIGHTBOT
wMONTHLYHEIGHTBOT.long_name = 'monthly height bottom'
wMONTHLYHEIGHTBOT.units = 'meters'

wMONTHLYHEIGHTTOP[:] = lsmMONTHLYHEIGHTTOP
wMONTHLYHEIGHTTOP.long_name = 'monthly height top'
wMONTHLYHEIGHTTOP.units = 'meters'

wMONTHLYLAI[:] = lsmMONTHLYLAI
wMONTHLYLAI.long_name = 'monthly leaf area index'
wMONTHLYLAI.units = 'm2/m2'

wMONTHLYSAI[:] = lsmMONTHLYSAI
wMONTHLYSAI.long_name = 'monthly stem area index'
wMONTHLYSAI.units = 'm2/m2'
