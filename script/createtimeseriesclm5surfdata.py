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
    print("Error: Usage createdatafile.py modiscurrentnamelistfile workrawfile workrawdir regionfile timeseriessurffile timeseriessurfdir")
    sys.exit()
else:
    print("Processing: " + str(sys.argv[1]))
    
pftnum = 15
cftnum = 64

modiscurrentnamelistfilename = str(sys.argv[1])
workrawfilename = str(sys.argv[2])
workrawdirname = str(sys.argv[3])
subregionfilename = str(sys.argv[4])
timeseriessurffilename = str(sys.argv[5])
timeseriessurfdirname = str(sys.argv[6])

modiscurrentnamelistfile = open(modiscurrentnamelistfilename)
modiscurrentnamelist = modiscurrentnamelistfile.readlines()
modiscurrentauthorparts = modiscurrentnamelist[0].split()
modiscurrentauthorpartnum = len(modiscurrentauthorparts)
modiscurrentauthor = ""
for modiscurrentauthorpartindex in range(modiscurrentauthorpartnum):
    if (modiscurrentauthorpartindex > 1):
        modiscurrentauthor = modiscurrentauthor + " "
    if (modiscurrentauthorpartindex > 0):
        modiscurrentauthor = modiscurrentauthor + modiscurrentauthorparts[modiscurrentauthorpartindex]
modiscurrentnameparts = modiscurrentnamelist[1].split()
modiscurrentname = modiscurrentnameparts[1]
timeseriessurfoutputparts = modiscurrentnamelist[2].split()
timeseriessurfoutput = timeseriessurfoutputparts[1]
timeseriessurfoutputnameparts = modiscurrentnamelist[3].split()
timeseriessurfoutputname = timeseriessurfoutputnameparts[1]
currentrawinputparts = modiscurrentnamelist[4].split()
currentrawinput = currentrawinputparts[1]
currentrawinputnameparts = modiscurrentnamelist[5].split()
currentrawinputname = currentrawinputnameparts[1]
clm5modisyearparts = modiscurrentnamelist[6].split()
clm5modisyearname = clm5modisyearparts[1]

print("author: " + modiscurrentauthor)
print("modiscurrentname: " + modiscurrentname)
print("timeseriessurfoutput: " + timeseriessurfoutput)
print("currentrawinput: " + currentrawinput)
print("currentrawinputname: " + currentrawinputname)
print("clm5modisyearname: " + clm5modisyearname)

workrawlistfile = open(workrawfilename,'r')
workrawlist = workrawlistfile.readlines()
numworkraw = len(workrawlist)
workrawnames = [""] * numworkraw
workrawdirs = [""] * numworkraw
workrawstartyears = [""] * numworkraw
workrawendyears = [""] * numworkraw
currentrawindex = -1
for workrawlistindex in range(numworkraw):
    workrawlistvalues = workrawlist[workrawlistindex].split()
    workrawnames[workrawlistindex] = workrawlistvalues[0]
    workrawdirs[workrawlistindex] = workrawlistvalues[1]
    if (workrawdirs[workrawlistindex] == "<workrawdir>"):
        workrawdirs[workrawlistindex] = workrawdirname
    workrawstartyears[workrawlistindex] = workrawlistvalues[2]
    workrawendyears[workrawlistindex] = workrawlistvalues[3]
    if (currentrawinput == workrawnames[workrawlistindex]):
        currentrawindex = workrawlistindex

if (currentrawindex == -1):
    print("Error: currentrawinput not found " + currentrawinput)
    sys.exit()
else:
    print("CLM5MODISEarthStat Raw Input: " + workrawdirs[currentrawindex] + workrawnames[currentrawindex])


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
    print("CLM5MODISEarthStat Surface Output: " + timeseriessurfdirs[timeseriessurfoutputindex] + timeseriessurfnames[timeseriessurfoutputindex])

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

climyear = int(clm5modisyearname)
startyear = int(clm5modisyearname)
endyear = int(clm5modisyearname)
numyears = endyear - startyear + 1

lsmtime = np.zeros(numyears,dtype=float)

currenttime = 0.0
for timeindex in range(numyears):
    lsmtime[timeindex] = currenttime
    currenttime += 365.0    

currentdirname = workrawdirs[currentrawindex] + workrawnames[currentrawindex] + "/" + currentrawinputname + "/"

outputfilename = timeseriessurfdirs[timeseriessurfoutputindex] + timeseriessurfnames[timeseriessurfoutputindex] + "/" + timeseriessurfoutputname

# process raw data files

lsmLANDMASK = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

LANDMASKfilename = currentdirname + currentrawinputname + ".LANDMASK." + str(climyear) + ".dat"
print("Processing LANDMASK > " + LANDMASKfilename)

insurffloatdata = np.fromfile(LANDMASKfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmLANDMASK[:] = insurffloatgrid[::-1,:]

lsmLANDFRAC = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

LANDFRACfilename = currentdirname + currentrawinputname + ".LANDFRAC." + str(climyear) + ".dat"
print("Processing LANDFRAC > " + LANDFRACfilename)

insurffloatdata = np.fromfile(LANDFRACfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmLANDFRAC[:] = insurffloatgrid[::-1,:]

lsmAREA = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

AREAfilename = currentdirname + currentrawinputname + ".AREA." + str(climyear) + ".dat"
print("Processing AREA > " + AREAfilename)

insurffloatdata = np.fromfile(AREAfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmAREA[:] = insurffloatgrid[::-1,:]

lsmPCTGLACIER = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTGLACIERfilename = currentdirname + currentrawinputname + ".PCTGLACIER." + str(climyear) + ".dat"
print("Processing PCTGLACIER > " + PCTGLACIERfilename)

insurffloatdata = np.fromfile(PCTGLACIERfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTGLACIER[:] = insurffloatgrid[::-1,:]

lsmPCTLAKE = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTLAKEfilename = currentdirname + currentrawinputname + ".PCTLAKE." + str(climyear) + ".dat"
print("Processing PCTLAKE > " + PCTLAKEfilename)

insurffloatdata = np.fromfile(PCTLAKEfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTLAKE[:] = insurffloatgrid[::-1,:]

lsmPCTWETLAND = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTWETLANDfilename = currentdirname + currentrawinputname + ".PCTWETLAND." + str(climyear) + ".dat"
print("Processing PCTWETLAND > " + PCTWETLANDfilename)

insurffloatdata = np.fromfile(PCTWETLANDfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTWETLAND[:] = insurffloatgrid[::-1,:]

lsmPCTURBAN = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTURBANfilename = currentdirname + currentrawinputname + ".PCTURBAN." + str(climyear) + ".dat"
print("Processing PCTURBAN > " + PCTURBANfilename)

insurffloatdata = np.fromfile(PCTURBANfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTURBAN[:] = insurffloatgrid[::-1,:]

lsmPCTNATVEG = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTNATVEGfilename = currentdirname + currentrawinputname + ".PCTNATVEG." + str(climyear) + ".dat"
print("Processing PCTNATVEG > " + PCTNATVEGfilename)

insurffloatdata = np.fromfile(PCTNATVEGfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTNATVEG[:] = insurffloatgrid[::-1,:]

lsmPCTCROP = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTCROPfilename = currentdirname + currentrawinputname + ".PCTCROP." + str(climyear) + ".dat"
print("Processing PCTCROP > " + PCTCROPfilename)

insurffloatdata = np.fromfile(PCTCROPfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTCROP[:] = insurffloatgrid[::-1,:]

lsmPCTNATPFT = np.zeros(shape=(pftnum,numyears,subregionlatcells,subregionloncells),dtype=np.float32)

for pftindex in range(pftnum):

    if (pftindex < 10):
        PCTNATPFTfilename = currentdirname + currentrawinputname + ".PCTNATPFT0" + str(pftindex) + "." + str(climyear) + ".dat"
    else:
        PCTNATPFTfilename = currentdirname + currentrawinputname + ".PCTNATPFT" + str(pftindex) + "." + str(climyear) + ".dat"
    
    print("Processing PCTNATPFT > " + PCTNATPFTfilename)

    insurffloatdata = np.fromfile(PCTNATPFTfilename, dtype="float32", count=-1)
    insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
    lsmPCTNATPFT[pftindex,0,:,:] = insurffloatgrid[::-1,:]

lsmPCTCFT = np.zeros(shape=(cftnum,numyears,subregionlatcells,subregionloncells),dtype=np.float32)

for cftindex in range(cftnum):

    if (cftindex < 10):
        PCTCFTfilename = currentdirname + currentrawinputname + ".PCTCFT0" + str(cftindex) + "." + str(climyear) + ".dat"
    else:
        PCTCFTfilename = currentdirname + currentrawinputname + ".PCTCFT" + str(cftindex) + "." + str(climyear) + ".dat"
    
    print("Processing PCTCFT > " + PCTCFTfilename)

    insurffloatdata = np.fromfile(PCTCFTfilename, dtype="float32", count=-1)
    insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
    lsmPCTCFT[cftindex,0,:,:] = insurffloatgrid[::-1,:]

# generate workraw data

print('Creating: ' + outputfilename)
outputfile = netcdf4.Dataset(outputfilename, 'w')

outputfile.Conventions = 'NCAR-CSM'
outputfile.Author = 'Peter Lawrence, Terrestrial Sciences Section, National Center for Atmospheric Research'
datenow = date.datetime.now()
datenowstr = datenow.strftime("%m-%d-%Y %H:%M:%S")
outputfile.History_Log = 'created on: ' + datenowstr
outputfile.Data_Log = currentdirname

outputfile.createDimension('natpft',int(pftnum))
outputfile.createDimension('cft',int(cftnum))
outputfile.createDimension('time',int(numyears))
outputfile.createDimension('lat',int(subregionlatcells))
outputfile.createDimension('lon',int(subregionloncells))

wnatpft = outputfile.createVariable('natpft','i',('natpft',))
wcft = outputfile.createVariable('cft','i',('cft',))
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
wPCTGLACIER = outputfile.createVariable('PCT_GLACIER',np.float32,('lat','lon'),fill_value=-9999.0)
wPCTLAKE = outputfile.createVariable('PCT_LAKE',np.float32,('lat','lon'),fill_value=-9999.0)
wPCTWETLAND = outputfile.createVariable('PCT_WETLAND',np.float32,('lat','lon'),fill_value=-9999.0)
wPCTURBAN = outputfile.createVariable('PCT_URBAN',np.float32,('lat','lon'),fill_value=-9999.0)
wPCTNATVEG = outputfile.createVariable('PCT_NATVEG',np.float32,('time','lat','lon'),fill_value=-9999.0)
wPCTCROP = outputfile.createVariable('PCT_CROP',np.float32,('time','lat','lon'),fill_value=-9999.0)
wPCTNATPFT = outputfile.createVariable('PCT_NAT_PFT',np.float32,('natpft','time','lat','lon'),fill_value=-9999.0)
wPCTCFT = outputfile.createVariable('PCT_CFT',np.float32,('cft','time','lat','lon'),fill_value=-9999.0)

wnatpft[:] = lsmnatpftid
wnatpft.long_name = 'indices of natural PFTs'
wnatpft.units = 'index'

wcft[:] = lsmcftid
wcft.long_name = 'indices of CFTs'
wcft.units = 'index'

wtime[0] = 0.0
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

wPCTGLACIER[:] = lsmPCTGLACIER
wPCTGLACIER.long_name = 'total percent glacier landunit'
wPCTGLACIER.units = 'unitless'

wPCTLAKE[:] = lsmPCTLAKE
wPCTLAKE.long_name = 'total percent lake landunit'
wPCTLAKE.units = 'unitless'

wPCTWETLAND[:] = lsmPCTWETLAND
wPCTWETLAND.long_name = 'total percent wetland landunit'
wPCTWETLAND.units = 'unitless'

wPCTURBAN[:] = lsmPCTURBAN
wPCTURBAN.long_name = 'total percent urban landunit'
wPCTURBAN.units = 'unitless'

wPCTNATVEG[:] = lsmPCTNATVEG
wPCTNATVEG.long_name = 'total percent natural vegetation landunit'
wPCTNATVEG.units = 'unitless'

wPCTNATPFT[:] = lsmPCTNATPFT
wPCTNATPFT.long_name = 'percent plant functional type on the natural veg landunit (% of landunit)'
wPCTNATPFT.units = 'unitless'

wPCTCROP[:] = lsmPCTCROP
wPCTCROP.long_name = 'total percent crop landunit'
wPCTCROP.units = 'unitless'

wPCTCFT[:] = lsmPCTCFT
wPCTCFT.long_name = 'percent crop functional type on the crop landunit (% of landunit)'
wPCTCFT.units = 'unitless'
