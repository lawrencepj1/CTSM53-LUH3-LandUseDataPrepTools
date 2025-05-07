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
    print("Error: Usage createdatafile.py workvcfnamelistfile workrawfile workrawdir regionfile worksurffile worksurfdir")
    sys.exit()
else:
    print("Processing: " + str(sys.argv[1]))
    
pftnum = 15
cftnum = 64

workvcfnamelistfilename = str(sys.argv[1])
workrawfilename = str(sys.argv[2])
workrawdirname = str(sys.argv[3])
subregionfilename = str(sys.argv[4])
worksurffilename = str(sys.argv[5])
worksurfdirname = str(sys.argv[6])

workvcfnamelistfile = open(workvcfnamelistfilename)
workvcfnamelist = workvcfnamelistfile.readlines()
workvcfauthorparts = workvcfnamelist[0].split()
workvcfauthorpartnum = len(workvcfauthorparts)
workvcfauthor = ""
for workvcfauthorpartindex in range(workvcfauthorpartnum):
    if (workvcfauthorpartindex > 1):
        workvcfauthor = workvcfauthor + " "
    if (workvcfauthorpartindex > 0):
        workvcfauthor = workvcfauthor + workvcfauthorparts[workvcfauthorpartindex]
workvcfnameparts = workvcfnamelist[1].split()
workvcfname = workvcfnameparts[1]
worksurfoutputparts = workvcfnamelist[2].split()
worksurfoutput = worksurfoutputparts[1]
worksurfoutputnameparts = workvcfnamelist[3].split()
worksurfoutputname = worksurfoutputnameparts[1]
workrawinputparts = workvcfnamelist[4].split()
workrawinput = workrawinputparts[1]
workrawinputnameparts = workvcfnamelist[5].split()
workrawinputname = workrawinputnameparts[1]
workrawinputyearnameparts = workvcfnamelist[6].split()
workrawinputyearname = workrawinputyearnameparts[1]

print("author: " + workvcfauthor)
print("workvcfname: " + workvcfname)
print("worksurfoutput: " + worksurfoutput)
print("worksurfoutputname: " + worksurfoutputname)
print("workrawinput: " + workrawinput)
print("workrawinputname: " + workrawinputname)

workrawlistfile = open(workrawfilename,'r')
workrawlist = workrawlistfile.readlines()
numworkraw = len(workrawlist)
workrawnames = [""] * numworkraw
workrawdirs = [""] * numworkraw
workrawstartyears = [""] * numworkraw
workrawendyears = [""] * numworkraw
workrawindex = -1
for workrawlistindex in range(numworkraw):
    workrawlistvalues = workrawlist[workrawlistindex].split()
    workrawnames[workrawlistindex] = workrawlistvalues[0]
    workrawdirs[workrawlistindex] = workrawlistvalues[1]
    if (workrawdirs[workrawlistindex] == "<workrawdir>"):
        workrawdirs[workrawlistindex] = workrawdirname
    workrawstartyears[workrawlistindex] = workrawlistvalues[2]
    workrawendyears[workrawlistindex] = workrawlistvalues[3]
    if (workrawinput == workrawnames[workrawlistindex]):
        workrawindex = workrawlistindex

if (workrawindex == -1):
    print("Error: workrawinput not found " + workrawinput)
    sys.exit()
else:
    print("CLM5MODISEarthStat Raw Input: " + workrawdirs[workrawindex] + workrawnames[workrawindex])


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
    print("CLM5MODISEarthStat Surface Output: " + worksurfdirs[worksurfoutputindex] + worksurfnames[worksurfoutputindex])

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

climyear = int(workrawinputyearname)
startyear = int(workrawinputyearname)
endyear = int(workrawinputyearname)
numyears = endyear - startyear + 1

lsmtime = np.zeros(numyears,dtype=float)

currenttime = 0.0
for timeindex in range(numyears):
    lsmtime[timeindex] = currenttime
    currenttime += 365.0    

currentdirname = workrawdirs[workrawindex] + workrawnames[workrawindex] + "/" + workrawinputname + "/"

outputfilename = worksurfdirs[worksurfoutputindex] + worksurfnames[worksurfoutputindex] + "/" + worksurfoutputname


# process raw data files

lsmLANDMASK = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

LANDMASKfilename = currentdirname + workrawinputname + ".LANDMASK." + str(climyear) + ".dat"
print("Processing LANDMASK > " + LANDMASKfilename)

insurffloatdata = np.fromfile(LANDMASKfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmLANDMASK[:] = insurffloatgrid[::-1,:]

lsmLANDFRAC = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

LANDFRACfilename = currentdirname + workrawinputname + ".LANDFRAC." + str(climyear) + ".dat"
print("Processing LANDFRAC > " + LANDFRACfilename)

insurffloatdata = np.fromfile(LANDFRACfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmLANDFRAC[:] = insurffloatgrid[::-1,:]

lsmAREA = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

AREAfilename = currentdirname + workrawinputname + ".AREA." + str(climyear) + ".dat"
print("Processing AREA > " + AREAfilename)

insurffloatdata = np.fromfile(AREAfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmAREA[:] = insurffloatgrid[::-1,:]

lsmPCTGLACIER = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTGLACIERfilename = currentdirname + workrawinputname + ".PCTGLACIER." + str(climyear) + ".dat"
print("Processing PCTGLACIER > " + PCTGLACIERfilename)

insurffloatdata = np.fromfile(PCTGLACIERfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTGLACIER[:] = insurffloatgrid[::-1,:]

lsmPCTLAKE = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTLAKEfilename = currentdirname + workrawinputname + ".PCTLAKE." + str(climyear) + ".dat"
print("Processing PCTLAKE > " + PCTLAKEfilename)

insurffloatdata = np.fromfile(PCTLAKEfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTLAKE[:] = insurffloatgrid[::-1,:]

lsmPCTWETLAND = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTWETLANDfilename = currentdirname + workrawinputname + ".PCTWETLAND." + str(climyear) + ".dat"
print("Processing PCTWETLAND > " + PCTWETLANDfilename)

insurffloatdata = np.fromfile(PCTWETLANDfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTWETLAND[:] = insurffloatgrid[::-1,:]

lsmPCTURBAN = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTURBANfilename = currentdirname + workrawinputname + ".PCTURBAN." + str(climyear) + ".dat"
print("Processing PCTURBAN > " + PCTURBANfilename)

insurffloatdata = np.fromfile(PCTURBANfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTURBAN[:] = insurffloatgrid[::-1,:]

lsmPCTNATVEG = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTNATVEGfilename = currentdirname + workrawinputname + ".PCTNATVEG." + str(climyear) + ".dat"
print("Processing PCTNATVEG > " + PCTNATVEGfilename)

insurffloatdata = np.fromfile(PCTNATVEGfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTNATVEG[:] = insurffloatgrid[::-1,:]

lsmPCTCROP = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTCROPfilename = currentdirname + workrawinputname + ".PCTCROP." + str(climyear) + ".dat"
print("Processing PCTCROP > " + PCTCROPfilename)

insurffloatdata = np.fromfile(PCTCROPfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTCROP[:] = insurffloatgrid[::-1,:]

lsmPCTTREE = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTTREEfilename = currentdirname + workrawinputname + ".PCTTREE." + str(climyear) + ".dat"
print("Processing PCTTREE > " + PCTTREEfilename)

insurffloatdata = np.fromfile(PCTTREEfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTTREE[:] = insurffloatgrid[::-1,:]

lsmPCTHERB = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTHERBfilename = currentdirname + workrawinputname + ".PCTHERB." + str(climyear) + ".dat"
print("Processing PCTHERB > " + PCTHERBfilename)

insurffloatdata = np.fromfile(PCTHERBfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTHERB[:] = insurffloatgrid[::-1,:]

lsmPCTBARE = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTBAREfilename = currentdirname + workrawinputname + ".PCTBARE." + str(climyear) + ".dat"
print("Processing PCTBARE > " + PCTBAREfilename)

insurffloatdata = np.fromfile(PCTBAREfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTBARE[:] = insurffloatgrid[::-1,:]


# generate timeseries data

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
wPCTTREE = outputfile.createVariable('PCT_TREE',np.float32,('time','lat','lon'),fill_value=-9999.0)
wPCTHERB = outputfile.createVariable('PCT_HERB',np.float32,('time','lat','lon'),fill_value=-9999.0)
wPCTBARE = outputfile.createVariable('PCT_BARE',np.float32,('time','lat','lon'),fill_value=-9999.0)

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

wPCTTREE[:] = lsmPCTTREE
wPCTTREE.long_name = 'total tree percent of natural vegetation landunit'
wPCTTREE.units = 'unitless'

wPCTHERB[:] = lsmPCTHERB
wPCTHERB.long_name = 'total herbaceous percent of natural vegetation landunit'
wPCTHERB.units = 'unitless'

wPCTBARE[:] = lsmPCTBARE
wPCTBARE.long_name = 'total bare percent of natural vegetation landunit'
wPCTBARE.units = 'unitless'
