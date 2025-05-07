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

if (arguments != 8):
    print("Error: Usage createdatafile.py luh3descnamelistfile luh3descrawfile luh3descrawdir workrawfilename workrawfiledir regionfile luh3descsurffile luh3descsurfdir")
    sys.exit()
else:
    print("Processing: " + str(sys.argv[1]))
    
pftnum = 15
cftnum = 64

luh3descnamelistfilename = str(sys.argv[1])
luh3descrawfilename = str(sys.argv[2])
luh3descrawdirname = str(sys.argv[3])
workrawfilename = str(sys.argv[4])
workrawdirname = str(sys.argv[5])
subregionfilename = str(sys.argv[6])
luh3descsurffilename = str(sys.argv[7])
luh3descsurfdirname = str(sys.argv[8])

luh3descnamelistfile = open(luh3descnamelistfilename)
luh3descnamelist = luh3descnamelistfile.readlines()
luh3descauthorparts = luh3descnamelist[0].split()
luh3descauthorpartnum = len(luh3descauthorparts)
luh3descauthor = ""
for luh3descauthorpartindex in range(luh3descauthorpartnum):
    if (luh3descauthorpartindex > 1):
        luh3descauthor = luh3descauthor + " "
    if (luh3descauthorpartindex > 0):
        luh3descauthor = luh3descauthor + luh3descauthorparts[luh3descauthorpartindex]
luh3descnameparts = luh3descnamelist[1].split()
luh3descname = luh3descnameparts[1]
luh3descsurfoutputparts = luh3descnamelist[2].split()
luh3descsurfoutput = luh3descsurfoutputparts[1]
luh3descsurfoutputnameparts = luh3descnamelist[3].split()
luh3descsurfoutputname = luh3descsurfoutputnameparts[1]
workclm5currentrawinputparts = luh3descnamelist[4].split()
workclm5currentrawinput = workclm5currentrawinputparts[1]
workclm5currentrawinputnameparts = luh3descnamelist[5].split()
workclm5currentrawinputname = workclm5currentrawinputnameparts[1]
workclm5currentrawinputyearnameparts = luh3descnamelist[6].split()
workclm5currentrawinputyearname = workclm5currentrawinputyearnameparts[1]
luh3descrawinputparts = luh3descnamelist[7].split()
luh3descrawinput = luh3descrawinputparts[1]
luh3descrawinputnameparts = luh3descnamelist[8].split()
luh3descrawinputname = luh3descrawinputnameparts[1]
luh3descrawinputstartyearnameparts = luh3descnamelist[9].split()
luh3descrawinputstartyearname = luh3descrawinputstartyearnameparts[1]
luh3descrawinputendyearnameparts = luh3descnamelist[10].split()
luh3descrawinputendyearname = luh3descrawinputendyearnameparts[1]

print("author: " + luh3descauthor)
print("luh3descname: " + luh3descname)
print("luh3descsurfoutput: " + luh3descsurfoutput)
print("luh3descsurfoutputname: " + luh3descsurfoutputname)
print("workclm5currentrawinput: " + workclm5currentrawinput)
print("workclm5currentrawinputname: " + workclm5currentrawinputname)
print("luh3descrawinput: " + luh3descrawinput)
print("luh3descrawinputname: " + luh3descrawinputname)

luh3descrawlistfile = open(luh3descrawfilename,'r')
luh3descrawlist = luh3descrawlistfile.readlines()
numluh3descraw = len(luh3descrawlist)
luh3descrawnames = [""] * numluh3descraw
luh3descrawdirs = [""] * numluh3descraw
luh3descrawstartyears = [""] * numluh3descraw
luh3descrawendyears = [""] * numluh3descraw
luh3descrawindex = -1
for luh3descrawlistindex in range(numluh3descraw):
    luh3descrawlistvalues = luh3descrawlist[luh3descrawlistindex].split()
    luh3descrawnames[luh3descrawlistindex] = luh3descrawlistvalues[0]
    luh3descrawdirs[luh3descrawlistindex] = luh3descrawlistvalues[1]
    if (luh3descrawdirs[luh3descrawlistindex] == "<luh3descrawdir>"):
        luh3descrawdirs[luh3descrawlistindex] = luh3descrawdirname
    luh3descrawstartyears[luh3descrawlistindex] = luh3descrawlistvalues[2]
    luh3descrawendyears[luh3descrawlistindex] = luh3descrawlistvalues[3]
    if (luh3descrawinput == luh3descrawnames[luh3descrawlistindex]):
        luh3descrawindex = luh3descrawlistindex

if (luh3descrawindex == -1):
    print("Error: luh3descrawinput not found " + luh3descrawinput)
    sys.exit()
else:
    print("Raw Input: " + luh3descrawdirs[luh3descrawindex] + luh3descrawnames[luh3descrawindex])

workrawlistfile = open(workrawfilename,'r')
workrawlist = workrawlistfile.readlines()
numworkraw = len(workrawlist)
workrawnames = [""] * numworkraw
workrawdirs = [""] * numworkraw
workrawstartyears = [""] * numworkraw
workrawendyears = [""] * numworkraw
workclm5currentrawindex = -1
for workrawlistindex in range(numworkraw):
    workrawlistvalues = workrawlist[workrawlistindex].split()
    workrawnames[workrawlistindex] = workrawlistvalues[0]
    workrawdirs[workrawlistindex] = workrawlistvalues[1]
    if (workrawdirs[workrawlistindex] == "<workrawdir>"):
        workrawdirs[workrawlistindex] = workrawdirname
    workrawstartyears[workrawlistindex] = workrawlistvalues[2]
    workrawendyears[workrawlistindex] = workrawlistvalues[3]
    if (workclm5currentrawinput == workrawnames[workrawlistindex]):
        workclm5currentrawindex = workrawlistindex

if (workclm5currentrawindex == -1):
    print("Error: workclm5currentrawinput not found " + workclm5currentrawinput)
    sys.exit()
else:
    print("Raw Input: " + workrawdirs[workclm5currentrawindex] + workrawnames[workclm5currentrawindex])

luh3descsurflistfile = open(luh3descsurffilename,'r')
luh3descsurflist = luh3descsurflistfile.readlines()
numluh3descsurf = len(luh3descsurflist)
luh3descsurfnames = [""] * numluh3descsurf
luh3descsurfdirs = [""] * numluh3descsurf
luh3descsurfstartyears = [""] * numluh3descsurf
luh3descsurfendyears = [""] * numluh3descsurf
luh3descsurfoutputindex = -1
for luh3descsurflistindex in range(numluh3descsurf):
    luh3descsurflistvalues = luh3descsurflist[luh3descsurflistindex].split()
    luh3descsurfnames[luh3descsurflistindex] = luh3descsurflistvalues[0]
    luh3descsurfdirs[luh3descsurflistindex] = luh3descsurflistvalues[1]
    if (luh3descsurfdirs[luh3descsurflistindex] == "<luh3descsurfdir>"):
        luh3descsurfdirs[luh3descsurflistindex] = luh3descsurfdirname
    luh3descsurfstartyears[luh3descsurflistindex] = luh3descsurflistvalues[2]
    luh3descsurfendyears[luh3descsurflistindex] = luh3descsurflistvalues[3]
    if (luh3descsurfoutput == luh3descsurfnames[luh3descsurflistindex]):
        luh3descsurfoutputindex = luh3descsurflistindex

if (luh3descsurfoutputindex == -1):
    print("Error: luh3descsurfoutput not found " + luh3descsurfoutput)
    sys.exit()
else:
    print("Surface Output: " + luh3descsurfdirs[luh3descsurfoutputindex] + luh3descsurfnames[luh3descsurfoutputindex])

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

climyear = int(workclm5currentrawinputyearname)
startyear = int(luh3descrawinputstartyearname)
endyear = int(luh3descrawinputendyearname)
numyears = endyear - startyear + 1

lsmtime = np.zeros(numyears,dtype=float)

currenttime = 0.0
for timeindex in range(numyears):
    lsmtime[timeindex] = currenttime
    currenttime += 365.0    

clm5currentdirname = workrawdirs[workclm5currentrawindex] + workrawnames[workclm5currentrawindex] + "/" + workclm5currentrawinputname + "/"
luh3descdirname = luh3descrawdirs[luh3descrawindex] + luh3descrawnames[luh3descrawindex] + "/" + luh3descrawinputname + "/"

outputfilename = luh3descsurfdirs[luh3descsurfoutputindex] + luh3descsurfnames[luh3descsurfoutputindex] + "/" + luh3descsurfoutputname


# process raw data files

lsmLANDMASK = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

LANDMASKfilename = clm5currentdirname + workclm5currentrawinputname + ".LANDMASK." + str(climyear) + ".dat"
print("Processing LANDMASK > " + LANDMASKfilename)

insurffloatdata = np.fromfile(LANDMASKfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmLANDMASK[:] = insurffloatgrid[::-1,:]

lsmLANDFRAC = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

LANDFRACfilename = clm5currentdirname + workclm5currentrawinputname + ".LANDFRAC." + str(climyear) + ".dat"
print("Processing LANDFRAC > " + LANDFRACfilename)

insurffloatdata = np.fromfile(LANDFRACfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmLANDFRAC[:] = insurffloatgrid[::-1,:]

lsmAREA = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

AREAfilename = clm5currentdirname + workclm5currentrawinputname + ".AREA." + str(climyear) + ".dat"
print("Processing AREA > " + AREAfilename)

insurffloatdata = np.fromfile(AREAfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmAREA[:] = insurffloatgrid[::-1,:]

lsmPCTGLACIER = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTGLACIERfilename = clm5currentdirname + workclm5currentrawinputname + ".PCTGLACIER." + str(climyear) + ".dat"
print("Processing PCTGLACIER > " + PCTGLACIERfilename)

insurffloatdata = np.fromfile(PCTGLACIERfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTGLACIER[:] = insurffloatgrid[::-1,:]

lsmPCTLAKE = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTLAKEfilename = clm5currentdirname + workclm5currentrawinputname + ".PCTLAKE." + str(climyear) + ".dat"
print("Processing PCTLAKE > " + PCTLAKEfilename)

insurffloatdata = np.fromfile(PCTLAKEfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTLAKE[:] = insurffloatgrid[::-1,:]

lsmPCTWETLAND = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTWETLANDfilename = clm5currentdirname + workclm5currentrawinputname + ".PCTWETLAND." + str(climyear) + ".dat"
print("Processing PCTWETLAND > " + PCTWETLANDfilename)

insurffloatdata = np.fromfile(PCTWETLANDfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTWETLAND[:] = insurffloatgrid[::-1,:]

lsmPCTURBAN = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTURBANfilename = clm5currentdirname + workclm5currentrawinputname + ".PCTURBAN." + str(climyear) + ".dat"
print("Processing PCTURBAN > " + PCTURBANfilename)

insurffloatdata = np.fromfile(PCTURBANfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTURBAN[:] = insurffloatgrid[::-1,:]

lsmPCTNATVEG = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTNATVEGfilename = clm5currentdirname + workclm5currentrawinputname + ".PCTNATVEG." + str(climyear) + ".dat"
print("Processing PCTNATVEG > " + PCTNATVEGfilename)

insurffloatdata = np.fromfile(PCTNATVEGfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTNATVEG[:] = insurffloatgrid[::-1,:]

lsmPCTCROP = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTCROPfilename = clm5currentdirname + workclm5currentrawinputname + ".PCTCROP." + str(climyear) + ".dat"
print("Processing PCTCROP > " + PCTCROPfilename)

insurffloatdata = np.fromfile(PCTCROPfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTCROP[:] = insurffloatgrid[::-1,:]

lsmPCTNATPFT = np.zeros(shape=(pftnum,numyears,subregionlatcells,subregionloncells),dtype=np.float32)

for timeindex in range(numyears):
    timeyear = timeindex + startyear
    
    for pftindex in range(pftnum):

        if (pftindex < 10):
            PCTNATPFTfilename = luh3descdirname + luh3descrawinputname + ".PCTNATPFT0" + str(pftindex) + "." + str(timeyear) + ".dat"
        else:
            PCTNATPFTfilename = luh3descdirname + luh3descrawinputname + ".PCTNATPFT" + str(pftindex) + "." + str(timeyear) + ".dat"
    
        print("Processing PCTNATPFT > " + PCTNATPFTfilename)

        insurffloatdata = np.fromfile(PCTNATPFTfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
        lsmPCTNATPFT[pftindex,timeindex,:,:] = insurffloatgrid[::-1,:]

# generate timeseries data

print('Creating: ' + outputfilename)
outputfile = netcdf4.Dataset(outputfilename, 'w')

outputfile.Conventions = 'NCAR-CSM'
outputfile.Author = 'Peter Lawrence, Terrestrial Sciences Section, National Center for Atmospheric Research'
datenow = date.datetime.now()
datenowstr = datenow.strftime("%m-%d-%Y %H:%M:%S")
outputfile.History_Log = 'created on: ' + datenowstr
outputfile.Data_Log = luh3descdirname

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
wPCTNATVEG = outputfile.createVariable('PCT_NATVEG',np.float32,('lat','lon'),fill_value=-9999.0)
wPCTCROP = outputfile.createVariable('PCT_CROP',np.float32,('lat','lon'),fill_value=-9999.0)
wPCTNATPFT = outputfile.createVariable('PCT_NAT_PFT',np.float32,('natpft','time','lat','lon'),fill_value=-9999.0)

wnatpft[:] = lsmnatpftid
wnatpft.long_name = 'indices of natural PFTs'
wnatpft.units = 'index'

wcft[:] = lsmcftid
wcft.long_name = 'indices of CFTs'
wcft.units = 'index'

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
