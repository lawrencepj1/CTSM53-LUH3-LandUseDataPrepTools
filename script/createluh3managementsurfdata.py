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

if (arguments != 9):
    print("Error: Usage createdatafile.py luh3namelistfile timeseriesrawfile timeseriesrawdir luh3rawfile luh3rawdir variablefile regionfile luh3surffile luh3surfdir")
    sys.exit()
else:
    print("Processing: " + str(sys.argv[1]))
    
luh3namelistfilename = str(sys.argv[1])
timeseriesrawfilename = str(sys.argv[2])
timeseriesrawdirname = str(sys.argv[3])
luh3rawfilename = str(sys.argv[4])
luh3rawdirname = str(sys.argv[5])
variablelistfilename = str(sys.argv[6])
subregionfilename = str(sys.argv[7])
luh3surffilename = str(sys.argv[8])
luh3surfdirname = str(sys.argv[9])

luh3namelistfile = open(luh3namelistfilename)
luh3namelist = luh3namelistfile.readlines()
luh3authorparts = luh3namelist[0].split()
luh3authorpartnum = len(luh3authorparts)
luh3author = ""
for luh3authorpartindex in range(luh3authorpartnum):
    if (luh3authorpartindex > 1):
        luh3author = luh3author + " "
    if (luh3authorpartindex > 0):
        luh3author = luh3author + luh3authorparts[luh3authorpartindex]
luh3nameparts = luh3namelist[1].split()
luh3name = luh3nameparts[1]
luh3surfoutputparts = luh3namelist[2].split()
luh3surfoutput = luh3surfoutputparts[1]
luh3rawinputparts = luh3namelist[3].split()
luh3rawinput = luh3rawinputparts[1]
timeseriesrawinputparts = luh3namelist[4].split()
timeseriesrawinput = timeseriesrawinputparts[1]
luh3startyearparts = luh3namelist[5].split()
luh3startyearname = luh3startyearparts[1]
luh3startyear = int(luh3startyearname)
luh3endyearparts = luh3namelist[6].split()
luh3endyearname = luh3endyearparts[1]
luh3endyear = int(luh3endyearname)
clm5modisyearparts = luh3namelist[7].split()
clm5modisyearname = clm5modisyearparts[1]

print("LUH3 Name: " + luh3name)
print("LUH3 Author: " + luh3author)
print("LUH3 Surface Output: " + luh3surfoutput)
print("LUH3 Raw Input: " + luh3rawinput)
print("CLM5 MODIS Raw Input: " + timeseriesrawinput)
print("LUH3 Start Year: " + luh3startyearname)
print("LUH3 End Year: " + luh3endyearname)
print("CLM5 Year: " + clm5modisyearname)

timeseriesrawlistfile = open(timeseriesrawfilename,'r')
timeseriesrawlist = timeseriesrawlistfile.readlines()
numtimeseriesraw = len(timeseriesrawlist)
timeseriesrawnames = [""] * numtimeseriesraw
timeseriesrawdirs = [""] * numtimeseriesraw
timeseriesrawstartyears = [""] * numtimeseriesraw
timeseriesrawendyears = [""] * numtimeseriesraw
timeseriesrawinputindex = -1
for timeseriesrawlistindex in range(numtimeseriesraw):
    timeseriesrawlistvalues = timeseriesrawlist[timeseriesrawlistindex].split()
    timeseriesrawnames[timeseriesrawlistindex] = timeseriesrawlistvalues[0]
    timeseriesrawdirs[timeseriesrawlistindex] = timeseriesrawlistvalues[1]
    if (timeseriesrawdirs[timeseriesrawlistindex] == "<timeseriesrawdir>"):
        timeseriesrawdirs[timeseriesrawlistindex] = timeseriesrawdirname
    timeseriesrawstartyears[timeseriesrawlistindex] = timeseriesrawlistvalues[2]
    timeseriesrawendyears[timeseriesrawlistindex] = timeseriesrawlistvalues[3]
    if (timeseriesrawinput == timeseriesrawnames[timeseriesrawlistindex]):
        timeseriesrawinputindex = timeseriesrawlistindex

if (timeseriesrawinputindex == -1):
    print("Error: timeseriesrawinput not found " + timeseriesrawinput)
    sys.exit()
else:
    print("CLM5 MODIS Raw Input: " + timeseriesrawdirs[timeseriesrawinputindex] + timeseriesrawnames[timeseriesrawinputindex])

luh3rawlistfile = open(luh3rawfilename,'r')
luh3rawlist = luh3rawlistfile.readlines()
numluh3raw = len(luh3rawlist)
luh3rawnames = [""] * numluh3raw
luh3rawdirs = [""] * numluh3raw
luh3rawstartyears = [""] * numluh3raw
luh3rawendyears = [""] * numluh3raw
luh3rawinputindex = -1
for luh3rawlistindex in range(numluh3raw):
    luh3rawlistvalues = luh3rawlist[luh3rawlistindex].split()
    luh3rawnames[luh3rawlistindex] = luh3rawlistvalues[0]
    luh3rawdirs[luh3rawlistindex] = luh3rawlistvalues[1]
    if (luh3rawdirs[luh3rawlistindex] == "<luh3rawdir>"):
        luh3rawdirs[luh3rawlistindex] = luh3rawdirname
    luh3rawstartyears[luh3rawlistindex] = luh3rawlistvalues[2]
    luh3rawendyears[luh3rawlistindex] = luh3rawlistvalues[3]
    if (luh3rawinput == luh3rawnames[luh3rawlistindex]):
        luh3rawinputindex = luh3rawlistindex

if (luh3rawinputindex == -1):
    print("Error: luh3rawinput not found " + luh3rawinput)
    sys.exit()
else:
    print("LUH3 Raw Input: " + luh3rawdirs[luh3rawinputindex] + luh3rawnames[luh3rawinputindex])

luh3surflistfile = open(luh3surffilename,'r')
luh3surflist = luh3surflistfile.readlines()
numluh3surf = len(luh3surflist)
luh3surfnames = [""] * numluh3surf
luh3surfdirs = [""] * numluh3surf
luh3surfstartyears = [""] * numluh3surf
luh3surfendyears = [""] * numluh3surf
luh3surfoutputindex = -1
for luh3surflistindex in range(numluh3surf):
    luh3surflistvalues = luh3surflist[luh3surflistindex].split()
    luh3surfnames[luh3surflistindex] = luh3surflistvalues[0]
    luh3surfdirs[luh3surflistindex] = luh3surflistvalues[1]
    if (luh3surfdirs[luh3surflistindex] == "<luh3surfdir>"):
        luh3surfdirs[luh3surflistindex] = luh3surfdirname
    luh3surfstartyears[luh3surflistindex] = luh3surflistvalues[2]
    luh3surfendyears[luh3surflistindex] = luh3surflistvalues[3]
    if (luh3surfoutput == luh3surfnames[luh3surflistindex]):
        luh3surfoutputindex = luh3surflistindex

if (luh3surfoutputindex == -1):
    print("Error: luh3surfoutput not found " + luh3surfoutput)
    sys.exit()
else:
    print("LUH3 Surface Output: " + luh3surfdirs[luh3surfoutputindex] + luh3surfnames[luh3surfoutputindex])

timeseriesrawinputcommonfilename = timeseriesrawdirs[timeseriesrawinputindex] + timeseriesrawnames[timeseriesrawinputindex] + "/" + timeseriesrawnames[timeseriesrawinputindex]
luh3rawinputcommonfilename = luh3rawdirs[luh3rawinputindex] + luh3rawnames[luh3rawinputindex] + "/" + luh3rawnames[luh3rawinputindex]

outputfilename = luh3surfdirs[luh3surfoutputindex] + luh3surfnames[luh3surfoutputindex] + "/management.nc"

variablelistfile = open(variablelistfilename,'r')
variablelist = variablelistfile.readlines()
numvariables = len(variablelist)
variablencnames = [""] * numvariables
variableoutnames = [""] * numvariables
variabledatatypes = [""] * numvariables
variableshapes = np.zeros(numvariables,dtype=int)
for variablelistindex in range(numvariables):
    variablelistvalues = variablelist[variablelistindex].split()
    variablencnames[variablelistindex] = variablelistvalues[0]
    variableoutnames[variablelistindex] = variablelistvalues[1]
    variabledatatypes[variablelistindex] = variablelistvalues[2]
    variableshapes[variablelistindex] = int(variablelistvalues[3])

subregionfile = open(subregionfilename,'r')
subregionlist = subregionfile.readlines()

subregionlllon = float(subregionlist[0])
subregionlllat = float(subregionlist[1])
subregionurlon = float(subregionlist[2])
subregionurlat = float(subregionlist[3])   
cellsize = float(subregionlist[4])

subregionloncells = int((subregionurlon - subregionlllon) / cellsize)
subregionlatcells = int((subregionurlat - subregionlllat) / cellsize)

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

climyear = clm5modisyearname
startyear = luh3startyear
endyear = luh3endyear
numyears = endyear - startyear + 1

lsmtime = np.zeros(numyears,dtype=float)

currenttime = 0.0
for timeindex in range(numyears):
    lsmtime[timeindex] = currenttime
    currenttime += 365.0    
    
# process raw data files

lsmLANDMASK = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

LANDMASKfilename = timeseriesrawinputcommonfilename + ".LANDMASK." + climyear + ".dat"
print("Processing LANDMASK > " + LANDMASKfilename)

insurffloatdata = np.fromfile(LANDMASKfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmLANDMASK[:] = insurffloatgrid[::-1,:]

lsmLANDFRAC = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

LANDFRACfilename = timeseriesrawinputcommonfilename + ".LANDFRAC." + climyear + ".dat"
print("Processing LANDFRAC > " + LANDFRACfilename)

insurffloatdata = np.fromfile(LANDFRACfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmLANDFRAC[:] = insurffloatgrid[::-1,:]

lsmAREA = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

AREAfilename = timeseriesrawinputcommonfilename + ".AREA." + climyear + ".dat"
print("Processing AREA > " + AREAfilename)

insurffloatdata = np.fromfile(AREAfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmAREA[:] = insurffloatgrid[::-1,:]

lsmPCTGLACIER = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTGLACIERfilename = timeseriesrawinputcommonfilename + ".PCTGLACIER." + climyear + ".dat"
print("Processing PCTGLACIER > " + PCTGLACIERfilename)

insurffloatdata = np.fromfile(PCTGLACIERfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTGLACIER[:] = insurffloatgrid[::-1,:]

lsmPCTLAKE = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTLAKEfilename = timeseriesrawinputcommonfilename + ".PCTLAKE." + climyear + ".dat"
print("Processing PCTLAKE > " + PCTLAKEfilename)

insurffloatdata = np.fromfile(PCTLAKEfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTLAKE[:] = insurffloatgrid[::-1,:]

lsmLUH3values = np.zeros(shape=(numyears,subregionlatcells,subregionloncells),dtype=np.float32)

# generate luh3surf data

print('Creating: ' + outputfilename)
outputfile = netcdf4.Dataset(outputfilename, 'w')

outputfile.Conventions = 'NCAR-CSM'
outputfile.Author = luh3author
datenow = date.datetime.now()
datenowstr = datenow.strftime("%m-%d-%Y %H:%M:%S")
outputfile.History_Log = 'created on: ' + datenowstr
outputfile.Data_Log = luh3rawinputcommonfilename

outputfile.createDimension('time',int(numyears))
outputfile.createDimension('lat',int(subregionlatcells))
outputfile.createDimension('lon',int(subregionloncells))

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
wGLACIER = outputfile.createVariable('GLACIER',np.float32,('lat','lon'))
wLAKE = outputfile.createVariable('LAKE',np.float32,('lat','lon'))

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

wGLACIER[:] = lsmPCTGLACIER / 100.0
wGLACIER.long_name = 'glacier fraction of gridcell'
wGLACIER.units = 'unitless'

wLAKE[:] = lsmPCTLAKE / 100.0
wLAKE.long_name = 'lake fraction of gridcell'
wLAKE.units = 'unitless'

for variablelistindex in range(numvariables):
    variablencname = variablencnames[variablelistindex]
    variableoutname = variableoutnames[variablelistindex]
    variabledatatype = variabledatatypes[variablelistindex]
    variableshape = variableshapes[variablelistindex]
    if (variabledatatype == "M"):
        if (variableshape == 1):
            wLUH3VAR = outputfile.createVariable(variablencname,np.float32,('time','lat','lon'),fill_value=-9999.0)

            for yearindex in range(numyears):
                yearnumber = startyear + yearindex

                LUH3variablefilename = luh3rawinputcommonfilename + "." + variableoutname + "." + str(yearnumber) + ".dat"

                print("Processing LUH3variable > " + LUH3variablefilename)

                insurffloatdata = np.fromfile(LUH3variablefilename, dtype="float32", count=-1)
                insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
                lsmLUH3values[yearindex,:,:] = insurffloatgrid[::-1,:]
            
            wLUH3VAR[:] = lsmLUH3values
            wLUH3VAR.long_name = variablencname
            wLUH3VAR.units = 'fraction'
