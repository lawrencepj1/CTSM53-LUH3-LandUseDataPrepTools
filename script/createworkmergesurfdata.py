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
    print("Error: Usage createdatafile.py workmergenamelistfile workrawfile workrawdir regionfile worksurffile worksurfdir")
    sys.exit()
else:
    print("Processing: " + str(sys.argv[1]))
    
pftnum = 15
cftnum = 64

workmergenamelistfilename = str(sys.argv[1])
workrawfilename = str(sys.argv[2])
workrawdirname = str(sys.argv[3])
subregionfilename = str(sys.argv[4])
worksurffilename = str(sys.argv[5])
worksurfdirname = str(sys.argv[6])

workmergenamelistfile = open(workmergenamelistfilename)
workmergenamelist = workmergenamelistfile.readlines()
workmergeauthorparts = workmergenamelist[0].split()
workmergeauthorpartnum = len(workmergeauthorparts)
workmergeauthor = ""
for workmergeauthorpartindex in range(workmergeauthorpartnum):
    if (workmergeauthorpartindex > 1):
        workmergeauthor = workmergeauthor + " "
    if (workmergeauthorpartindex > 0):
        workmergeauthor = workmergeauthor + workmergeauthorparts[workmergeauthorpartindex]
workmergenameparts = workmergenamelist[1].split()
workmergename = workmergenameparts[1]
workmergeluh3typesparts = workmergenamelist[2].split()
workmergeluh3types = workmergeluh3typesparts[1]
worksurfoutputparts = workmergenamelist[3].split()
worksurfoutput = worksurfoutputparts[1]
worksurfoutputnameparts = workmergenamelist[4].split()
worksurfoutputname = worksurfoutputnameparts[1]
clm5currentrawinputparts = workmergenamelist[5].split()
clm5currentrawinput = clm5currentrawinputparts[1]
clm5currentrawinputnameparts = workmergenamelist[6].split()
clm5currentrawinputname = clm5currentrawinputnameparts[1]
clm5currentrawinputyearnameparts = workmergenamelist[7].split()
clm5currentrawinputyearname = clm5currentrawinputyearnameparts[1]
workrawinputparts = workmergenamelist[8].split()
workrawinput = workrawinputparts[1]
workrawinputnameparts = workmergenamelist[9].split()
workrawinputname = workrawinputnameparts[1]
workrawinputstartyearnameparts = workmergenamelist[10].split()
workrawinputstartyearname = workrawinputstartyearnameparts[1]
workrawinputendyearnameparts = workmergenamelist[11].split()
workrawinputendyearname = workrawinputendyearnameparts[1]

print("author: " + workmergeauthor)
print("workmergename: " + workmergename)
print("worksurfoutput: " + worksurfoutput)
print("worksurfoutputname: " + worksurfoutputname)
print("clm5currentrawinput: " + clm5currentrawinput)
print("clm5currentrawinputname: " + clm5currentrawinputname)
print("workrawinput: " + workrawinput)
print("workrawinputname: " + workrawinputname)

workrawlistfile = open(workrawfilename,'r')
workrawlist = workrawlistfile.readlines()
numworkraw = len(workrawlist)
workrawnames = [""] * numworkraw
workrawdirs = [""] * numworkraw
workrawstartyears = [""] * numworkraw
workrawendyears = [""] * numworkraw
clm5currentrawindex = -1
workrawindex = -1
for workrawlistindex in range(numworkraw):
    workrawlistvalues = workrawlist[workrawlistindex].split()
    workrawnames[workrawlistindex] = workrawlistvalues[0]
    workrawdirs[workrawlistindex] = workrawlistvalues[1]
    if (workrawdirs[workrawlistindex] == "<workrawdir>"):
        workrawdirs[workrawlistindex] = workrawdirname
    workrawstartyears[workrawlistindex] = workrawlistvalues[2]
    workrawendyears[workrawlistindex] = workrawlistvalues[3]
    if (clm5currentrawinput == workrawnames[workrawlistindex]):
        clm5currentrawindex = workrawlistindex
    if (workrawinput == workrawnames[workrawlistindex]):
        workrawindex = workrawlistindex

if (clm5currentrawindex == -1):
    print("Error: clm5currentrawinput not found " + clm5currentrawinput)
    sys.exit()
else:
    print("Raw Input: " + workrawdirs[clm5currentrawindex] + workrawnames[clm5currentrawindex])

if (workrawindex == -1):
    print("Error: workrawinput not found " + workrawinput)
    sys.exit()
else:
    print("Raw Input: " + workrawdirs[workrawindex] + workrawnames[workrawindex])


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
    print("Surface Output: " + worksurfdirs[worksurfoutputindex] + worksurfnames[worksurfoutputindex])

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

climyear = int(clm5currentrawinputyearname)
startyear = int(workrawinputstartyearname)
endyear = int(workrawinputendyearname)
numyears = endyear - startyear + 1

lsmtime = np.zeros(numyears,dtype=float)

currenttime = 0.0
for timeindex in range(numyears):
    lsmtime[timeindex] = currenttime
    currenttime += 365.0    

clm5currentdirname = workrawdirs[clm5currentrawindex] + workrawnames[clm5currentrawindex] + "/" + clm5currentrawinputname + "/"
mergedirname = workrawdirs[workrawindex] + workrawnames[workrawindex] + "/" + workrawinputname + "/"

outputfilename = worksurfdirs[worksurfoutputindex] + worksurfnames[worksurfoutputindex] + "/" + worksurfoutputname


# process raw data files

lsmLANDMASK = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

LANDMASKfilename = clm5currentdirname + clm5currentrawinputname + ".LANDMASK." + str(climyear) + ".dat"
print("Processing LANDMASK > " + LANDMASKfilename)

insurffloatdata = np.fromfile(LANDMASKfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmLANDMASK[:] = insurffloatgrid[::-1,:]

lsmLANDFRAC = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

LANDFRACfilename = clm5currentdirname + clm5currentrawinputname + ".LANDFRAC." + str(climyear) + ".dat"
print("Processing LANDFRAC > " + LANDFRACfilename)

insurffloatdata = np.fromfile(LANDFRACfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmLANDFRAC[:] = insurffloatgrid[::-1,:]

lsmAREA = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

AREAfilename = clm5currentdirname + clm5currentrawinputname + ".AREA." + str(climyear) + ".dat"
print("Processing AREA > " + AREAfilename)

insurffloatdata = np.fromfile(AREAfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmAREA[:] = insurffloatgrid[::-1,:]

lsmPCTGLACIER = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTGLACIERfilename = clm5currentdirname + clm5currentrawinputname + ".PCTGLACIER." + str(climyear) + ".dat"
print("Processing PCTGLACIER > " + PCTGLACIERfilename)

insurffloatdata = np.fromfile(PCTGLACIERfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTGLACIER[:] = insurffloatgrid[::-1,:]

lsmPCTLAKE = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTLAKEfilename = clm5currentdirname + clm5currentrawinputname + ".PCTLAKE." + str(climyear) + ".dat"
print("Processing PCTLAKE > " + PCTLAKEfilename)

insurffloatdata = np.fromfile(PCTLAKEfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTLAKE[:] = insurffloatgrid[::-1,:]

lsmPCTWETLAND = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTWETLANDfilename = clm5currentdirname + clm5currentrawinputname + ".PCTWETLAND." + str(climyear) + ".dat"
print("Processing PCTWETLAND > " + PCTWETLANDfilename)

insurffloatdata = np.fromfile(PCTWETLANDfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTWETLAND[:] = insurffloatgrid[::-1,:]

lsmPCTURBAN = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTURBANfilename = clm5currentdirname + clm5currentrawinputname + ".PCTURBAN." + str(climyear) + ".dat"
print("Processing PCTURBAN > " + PCTURBANfilename)

insurffloatdata = np.fromfile(PCTURBANfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTURBAN[:] = insurffloatgrid[::-1,:]

lsmPCTNATVEG = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTNATVEGfilename = clm5currentdirname + clm5currentrawinputname + ".PCTNATVEG." + str(climyear) + ".dat"
print("Processing PCTNATVEG > " + PCTNATVEGfilename)

insurffloatdata = np.fromfile(PCTNATVEGfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTNATVEG[:] = insurffloatgrid[::-1,:]

lsmPCTCROP = np.zeros(shape=(subregionlatcells,subregionloncells),dtype=np.float32)

PCTCROPfilename = clm5currentdirname + clm5currentrawinputname + ".PCTCROP." + str(climyear) + ".dat"
print("Processing PCTCROP > " + PCTCROPfilename)

insurffloatdata = np.fromfile(PCTCROPfilename, dtype="float32", count=-1)
insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
lsmPCTCROP[:] = insurffloatgrid[::-1,:]

if (workmergeluh3types == "PRIMF" or workmergeluh3types == "ALL"):

    lsmPRIMFFRACTREE = np.zeros(shape=(numyears,subregionlatcells,subregionloncells),dtype=np.float32)

    for yearindex in range(numyears):
        currentyear = startyear + yearindex
        PRIMFFRACTREEfilename = mergedirname + workrawinputname + ".PRIMFFRACTREE." + str(currentyear) + ".dat"
        print("Processing PRIMFFRACTREE > " + PRIMFFRACTREEfilename)

        insurffloatdata = np.fromfile(PRIMFFRACTREEfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
        lsmPRIMFFRACTREE[yearindex,:,:] = insurffloatgrid[::-1,:]

    lsmPRIMFFRACHERB = np.zeros(shape=(numyears,subregionlatcells,subregionloncells),dtype=np.float32)

    for yearindex in range(numyears):
        currentyear = startyear + yearindex
        PRIMFFRACHERBfilename = mergedirname + workrawinputname + ".PRIMFFRACHERB." + str(currentyear) + ".dat"
        print("Processing PRIMFFRACHERB > " + PRIMFFRACHERBfilename)

        insurffloatdata = np.fromfile(PRIMFFRACHERBfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
        lsmPRIMFFRACHERB[yearindex,:,:] = insurffloatgrid[::-1,:]

    lsmPRIMFFRACBARE = np.zeros(shape=(numyears,subregionlatcells,subregionloncells),dtype=np.float32)

    for yearindex in range(numyears):
        currentyear = startyear + yearindex
        PRIMFFRACBAREfilename = mergedirname + workrawinputname + ".PRIMFFRACBARE." + str(currentyear) + ".dat"
        print("Processing PRIMFFRACBARE > " + PRIMFFRACBAREfilename)

        insurffloatdata = np.fromfile(PRIMFFRACBAREfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
        lsmPRIMFFRACBARE[yearindex,:,:] = insurffloatgrid[::-1,:]

if (workmergeluh3types == "SECDF" or workmergeluh3types == "ALL"):

    lsmSECDFFRACTREE = np.zeros(shape=(numyears,subregionlatcells,subregionloncells),dtype=np.float32)

    for yearindex in range(numyears):
        currentyear = startyear + yearindex
        SECDFFRACTREEfilename = mergedirname + workrawinputname + ".SECDFFRACTREE." + str(currentyear) + ".dat"
        print("Processing SECDFFRACTREE > " + SECDFFRACTREEfilename)

        insurffloatdata = np.fromfile(SECDFFRACTREEfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
        lsmSECDFFRACTREE[yearindex,:,:] = insurffloatgrid[::-1,:]

    lsmSECDFFRACHERB = np.zeros(shape=(numyears,subregionlatcells,subregionloncells),dtype=np.float32)

    for yearindex in range(numyears):
        currentyear = startyear + yearindex
        SECDFFRACHERBfilename = mergedirname + workrawinputname + ".SECDFFRACHERB." + str(currentyear) + ".dat"
        print("Processing SECDFFRACHERB > " + SECDFFRACHERBfilename)

        insurffloatdata = np.fromfile(SECDFFRACHERBfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
        lsmSECDFFRACHERB[yearindex,:,:] = insurffloatgrid[::-1,:]

    lsmSECDFFRACBARE = np.zeros(shape=(numyears,subregionlatcells,subregionloncells),dtype=np.float32)

    for yearindex in range(numyears):
        currentyear = startyear + yearindex
        SECDFFRACBAREfilename = mergedirname + workrawinputname + ".SECDFFRACBARE." + str(currentyear) + ".dat"
        print("Processing SECDFFRACBARE > " + SECDFFRACBAREfilename)

        insurffloatdata = np.fromfile(SECDFFRACBAREfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
        lsmSECDFFRACBARE[yearindex,:,:] = insurffloatgrid[::-1,:]

if (workmergeluh3types == "PRIMN" or workmergeluh3types == "ALL"):

    lsmPRIMNFRACTREE = np.zeros(shape=(numyears,subregionlatcells,subregionloncells),dtype=np.float32)

    for yearindex in range(numyears):
        currentyear = startyear + yearindex
        PRIMNFRACTREEfilename = mergedirname + workrawinputname + ".PRIMNFRACTREE." + str(currentyear) + ".dat"
        print("Processing PRIMNFRACTREE > " + PRIMNFRACTREEfilename)

        insurffloatdata = np.fromfile(PRIMNFRACTREEfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
        lsmPRIMNFRACTREE[yearindex,:,:] = insurffloatgrid[::-1,:]

    lsmPRIMNFRACHERB = np.zeros(shape=(numyears,subregionlatcells,subregionloncells),dtype=np.float32)

    for yearindex in range(numyears):
        currentyear = startyear + yearindex
        PRIMNFRACHERBfilename = mergedirname + workrawinputname + ".PRIMNFRACHERB." + str(currentyear) + ".dat"
        print("Processing PRIMNFRACHERB > " + PRIMNFRACHERBfilename)

        insurffloatdata = np.fromfile(PRIMNFRACHERBfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
        lsmPRIMNFRACHERB[yearindex,:,:] = insurffloatgrid[::-1,:]

    lsmPRIMNFRACBARE = np.zeros(shape=(numyears,subregionlatcells,subregionloncells),dtype=np.float32)

    for yearindex in range(numyears):
        currentyear = startyear + yearindex
        PRIMNFRACBAREfilename = mergedirname + workrawinputname + ".PRIMNFRACBARE." + str(currentyear) + ".dat"
        print("Processing PRIMNFRACBARE > " + PRIMNFRACBAREfilename)

        insurffloatdata = np.fromfile(PRIMNFRACBAREfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
        lsmPRIMNFRACBARE[yearindex,:,:] = insurffloatgrid[::-1,:]

if (workmergeluh3types == "SECDN" or workmergeluh3types == "ALL"):

    lsmSECDNFRACTREE = np.zeros(shape=(numyears,subregionlatcells,subregionloncells),dtype=np.float32)

    for yearindex in range(numyears):
        currentyear = startyear + yearindex
        SECDNFRACTREEfilename = mergedirname + workrawinputname + ".SECDNFRACTREE." + str(currentyear) + ".dat"
        print("Processing SECDNFRACTREE > " + SECDNFRACTREEfilename)

        insurffloatdata = np.fromfile(SECDNFRACTREEfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
        lsmSECDNFRACTREE[yearindex,:,:] = insurffloatgrid[::-1,:]

    lsmSECDNFRACHERB = np.zeros(shape=(numyears,subregionlatcells,subregionloncells),dtype=np.float32)

    for yearindex in range(numyears):
        currentyear = startyear + yearindex
        SECDNFRACHERBfilename = mergedirname + workrawinputname + ".SECDNFRACHERB." + str(currentyear) + ".dat"
        print("Processing SECDNFRACHERB > " + SECDNFRACHERBfilename)

        insurffloatdata = np.fromfile(SECDNFRACHERBfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
        lsmSECDNFRACHERB[yearindex,:,:] = insurffloatgrid[::-1,:]

    lsmSECDNFRACBARE = np.zeros(shape=(numyears,subregionlatcells,subregionloncells),dtype=np.float32)

    for yearindex in range(numyears):
        currentyear = startyear + yearindex
        SECDNFRACBAREfilename = mergedirname + workrawinputname + ".SECDNFRACBARE." + str(currentyear) + ".dat"
        print("Processing SECDNFRACBARE > " + SECDNFRACBAREfilename)

        insurffloatdata = np.fromfile(SECDNFRACBAREfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
        lsmSECDNFRACBARE[yearindex,:,:] = insurffloatgrid[::-1,:]

if (workmergeluh3types == "RANGE" or workmergeluh3types == "ALL"):

    lsmRANGEFRACTREE = np.zeros(shape=(numyears,subregionlatcells,subregionloncells),dtype=np.float32)

    for yearindex in range(numyears):
        currentyear = startyear + yearindex
        RANGEFRACTREEfilename = mergedirname + workrawinputname + ".RANGEFRACTREE." + str(currentyear) + ".dat"
        print("Processing RANGEFRACTREE > " + RANGEFRACTREEfilename)

        insurffloatdata = np.fromfile(RANGEFRACTREEfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
        lsmRANGEFRACTREE[yearindex,:,:] = insurffloatgrid[::-1,:]

    lsmRANGEFRACHERB = np.zeros(shape=(numyears,subregionlatcells,subregionloncells),dtype=np.float32)

    for yearindex in range(numyears):
        currentyear = startyear + yearindex
        RANGEFRACHERBfilename = mergedirname + workrawinputname + ".RANGEFRACHERB." + str(currentyear) + ".dat"
        print("Processing RANGEFRACHERB > " + RANGEFRACHERBfilename)

        insurffloatdata = np.fromfile(RANGEFRACHERBfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
        lsmRANGEFRACHERB[yearindex,:,:] = insurffloatgrid[::-1,:]

    lsmRANGEFRACBARE = np.zeros(shape=(numyears,subregionlatcells,subregionloncells),dtype=np.float32)

    for yearindex in range(numyears):
        currentyear = startyear + yearindex
        RANGEFRACBAREfilename = mergedirname + workrawinputname + ".RANGEFRACBARE." + str(currentyear) + ".dat"
        print("Processing RANGEFRACBARE > " + RANGEFRACBAREfilename)

        insurffloatdata = np.fromfile(RANGEFRACBAREfilename, dtype="float32", count=-1)
        insurffloatgrid = np.reshape(insurffloatdata,(subregionlatcells,subregionloncells))
        lsmRANGEFRACBARE[yearindex,:,:] = insurffloatgrid[::-1,:]


# generate timeseries data

print('Creating: ' + outputfilename)
outputfile = netcdf4.Dataset(outputfilename, 'w')

outputfile.Conventions = 'NCAR-CSM'
outputfile.Author = 'Peter Lawrence, Terrestrial Sciences Section, National Center for Atmospheric Research'
datenow = date.datetime.now()
datenowstr = datenow.strftime("%m-%d-%Y %H:%M:%S")
outputfile.History_Log = 'created on: ' + datenowstr
outputfile.Data_Log = mergedirname

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
if (workmergeluh3types == "PRIMF" or workmergeluh3types == "ALL"):
    wPRIMFFRACTREE = outputfile.createVariable('PRIMF_FRAC_TREE',np.float32,('time','lat','lon'),fill_value=-9999.0)
    wPRIMFFRACHERB = outputfile.createVariable('PRIMF_FRAC_HERB',np.float32,('time','lat','lon'),fill_value=-9999.0)
    wPRIMFFRACBARE = outputfile.createVariable('PRIMF_FRAC_BARE',np.float32,('time','lat','lon'),fill_value=-9999.0)
if (workmergeluh3types == "SECDF" or workmergeluh3types == "ALL"):
    wSECDFFRACTREE = outputfile.createVariable('SECDF_FRAC_TREE',np.float32,('time','lat','lon'),fill_value=-9999.0)
    wSECDFFRACHERB = outputfile.createVariable('SECDF_FRAC_HERB',np.float32,('time','lat','lon'),fill_value=-9999.0)
    wSECDFFRACBARE = outputfile.createVariable('SECDF_FRAC_BARE',np.float32,('time','lat','lon'),fill_value=-9999.0)
if (workmergeluh3types == "PRIMN" or workmergeluh3types == "ALL"):
    wPRIMNFRACTREE = outputfile.createVariable('PRIMN_FRAC_TREE',np.float32,('time','lat','lon'),fill_value=-9999.0)
    wPRIMNFRACHERB = outputfile.createVariable('PRIMN_FRAC_HERB',np.float32,('time','lat','lon'),fill_value=-9999.0)
    wPRIMNFRACBARE = outputfile.createVariable('PRIMN_FRAC_BARE',np.float32,('time','lat','lon'),fill_value=-9999.0)
if (workmergeluh3types == "SECDN" or workmergeluh3types == "ALL"):
    wSECDNFRACTREE = outputfile.createVariable('SECDN_FRAC_TREE',np.float32,('time','lat','lon'),fill_value=-9999.0)
    wSECDNFRACHERB = outputfile.createVariable('SECDN_FRAC_HERB',np.float32,('time','lat','lon'),fill_value=-9999.0)
    wSECDNFRACBARE = outputfile.createVariable('SECDN_FRAC_BARE',np.float32,('time','lat','lon'),fill_value=-9999.0)
if (workmergeluh3types == "RANGE" or workmergeluh3types == "ALL"):
    wRANGEFRACTREE = outputfile.createVariable('RANGE_FRAC_TREE',np.float32,('time','lat','lon'),fill_value=-9999.0)
    wRANGEFRACHERB = outputfile.createVariable('RANGE_FRAC_HERB',np.float32,('time','lat','lon'),fill_value=-9999.0)
    wRANGEFRACBARE = outputfile.createVariable('RANGE_FRAC_BARE',np.float32,('time','lat','lon'),fill_value=-9999.0)

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

if (workmergeluh3types == "PRIMF" or workmergeluh3types == "ALL"):

    wPRIMFFRACTREE[:] = lsmPRIMFFRACTREE
    wPRIMFFRACTREE.long_name = 'total tree fraction of primary forest landunit'
    wPRIMFFRACTREE.units = 'unitless'

    wPRIMFFRACHERB[:] = lsmPRIMFFRACHERB
    wPRIMFFRACHERB.long_name = 'total herb fraction of primary forest landunit'
    wPRIMFFRACHERB.units = 'unitless'

    wPRIMFFRACBARE[:] = lsmPRIMFFRACBARE
    wPRIMFFRACBARE.long_name = 'total bare fraction of primary forest landunit'
    wPRIMFFRACBARE.units = 'unitless'

if (workmergeluh3types == "SECDF" or workmergeluh3types == "ALL"):

    wSECDFFRACTREE[:] = lsmSECDFFRACTREE
    wSECDFFRACTREE.long_name = 'total tree fraction of secondary forest landunit'
    wSECDFFRACTREE.units = 'unitless'

    wSECDFFRACHERB[:] = lsmSECDFFRACHERB
    wSECDFFRACHERB.long_name = 'total herb fraction of secondary forest landunit'
    wSECDFFRACHERB.units = 'unitless'

    wSECDFFRACBARE[:] = lsmSECDFFRACBARE
    wSECDFFRACBARE.long_name = 'total bare fraction of secondary forest landunit'
    wSECDFFRACBARE.units = 'unitless'

if (workmergeluh3types == "PRIMN" or workmergeluh3types == "ALL"):

    wPRIMNFRACTREE[:] = lsmPRIMNFRACTREE
    wPRIMNFRACTREE.long_name = 'total tree fraction of primary non forest landunit'
    wPRIMNFRACTREE.units = 'unitless'

    wPRIMNFRACHERB[:] = lsmPRIMNFRACHERB
    wPRIMNFRACHERB.long_name = 'total herb fraction of primary non forest landunit'
    wPRIMNFRACHERB.units = 'unitless'

    wPRIMNFRACBARE[:] = lsmPRIMNFRACBARE
    wPRIMNFRACBARE.long_name = 'total bare fraction of primary non forest landunit'
    wPRIMNFRACBARE.units = 'unitless'

if (workmergeluh3types == "RANGE" or workmergeluh3types == "ALL"):

    wRANGEFRACTREE[:] = lsmRANGEFRACTREE
    wRANGEFRACTREE.long_name = 'total tree fraction of rangeland landunit'
    wRANGEFRACTREE.units = 'unitless'

    wRANGEFRACHERB[:] = lsmRANGEFRACHERB
    wRANGEFRACHERB.long_name = 'total herb fraction of rangeland landunit'
    wRANGEFRACHERB.units = 'unitless'

    wRANGEFRACBARE[:] = lsmRANGEFRACBARE
    wRANGEFRACBARE.long_name = 'total bare fraction of rangeland landunit'
    wRANGEFRACBARE.units = 'unitless'
