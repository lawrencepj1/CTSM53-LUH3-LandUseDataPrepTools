/* Generate Global CLM Surface Data from LUH format time series and MODIS and MIRCA2000 current day reference data */
/* Author Peter Lawrence - Terrestrial Sciences Section - National Center for Atmospheric Research */
/* Contact lawrence@ucar.edu 303 - 497 1727 */

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#define MAXCLMPIX 1440
#define MAXCLMLIN 720

#define CLMPIXSIZE 0.25
#define CLMLLX -180.0
#define CLMLLY -90.0

#define MAXLUHPIX 1440
#define MAXLUHLIN 720

#define LUHPIXSIZE 0.25
#define LUHLLX -180.0
#define LUHLLY -90.0

#define PI 4.0*atan(1.0)
#define EarthCir 40075.017

#define MAXPFT 15
#define MAXCFTRAW 31
#define MAXCFT 64

#define firsttreepft 1
#define lasttreepft 8

#define PRECIPCLASSES 100
#define PRECIPCLASSSTART 0.0
#define PRECIPCLASSSIZE 50.0

#define TEMPCLASSES 100
#define TEMPCLASSSTART -15.0
#define TEMPCLASSSIZE 0.5

#define CLMPCTThreshold 0.1
#define CLMPCTTruncate 10.0

#define LUH3StateThreshold 0.01
#define LUH3StateTruncate 1000.0

long MAXOUTPIX = MAXCLMPIX;
long MAXOUTLIN = MAXCLMLIN;
long OUTLONOFFSET = 0;
long OUTLATOFFSET = 0;
float OUTPIXSIZE = 0.25;
float OUTLLX = -180.0;
float OUTLLY = -90.0;

float *tempoutGrid;

long OUTDATASIZE = sizeof(float) * MAXCLMPIX * MAXCLMLIN;

char namelistauthor[1024];
char namelistforestallmaxvcfname[1024];
char namelistworkrawvcfoutput[1024];
char namelistworkrawvcfoutputforestname[1024];
char namelistworkrawvcfoutputnonforestname[1024];
char namelistworkrawvcfoutputrangename[1024];
char namelistworkrawclm5input[1024];
char namelistworkrawclm5inputname[1024];
char namelistworkrawclm5inputreferenceyear[1024];
char namelistworkrawmaxvcfinput[1024];
char namelistworkrawmaxvcfinputforestname[1024];
char namelistworkrawmaxvcfinputnonforestname[1024];
char namelistworkrawmaxvcfinputrangename[1024];
char namelistworkrawmaxvcfinputreferenceyear[1024];
char namelistworkrawmaxvcfinputmaxtreescale[1024];
char namelistclm5modisrawinput[1024];
char namelistclm5modisrawinputreferenceyear[1024];

int clm5modisrawinputreferenceyear;
int workrawclm5inputreferenceyear;
int workrawmaxvcfinputreferenceyear;
float workrawmaxvcfinputmaxtreescale;

char clm5modisrawinputname[1024];
char clm5modisrawinputdb[1024];
int clm5modisrawinputstartyear;
int clm5modisrawinputendyear;

char workrawclm5inputname[1024];
char workrawclm5inputdb[1024];
int workrawclm5inputstartyear;
int workrawclm5inputendyear;

char workrawmaxvcfinputname[1024];
char workrawmaxvcfinputdb[1024];
int workrawmaxvcfinputstartyear;
int workrawmaxvcfinputendyear;

char workrawvcfoutputname[1024];
char workrawvcfoutputdb[1024];
int workrawvcfoutputstartyear;
int workrawvcfoutputendyear;

char PFTluhtype[MAXPFT][256];
char CFTRAWluhtype[MAXCFTRAW][256];
char CFTluhtype[MAXCFT][256];

float tempLUHGrid[MAXLUHPIX * MAXLUHLIN];
float inLUHLANDMASKGrid[MAXLUHPIX * MAXLUHLIN];

float tempGrid[MAXCLMPIX * MAXCLMLIN];

float inLANDMASKGrid[MAXCLMPIX * MAXCLMLIN];
float inLANDFRACGrid[MAXCLMPIX * MAXCLMLIN];
float inAREAGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMPCTGLACIERGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMPCTLAKEGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMPCTWETLANDGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMPCTURBANGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMPCTNATVEGGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMPCTCROPGrid[MAXCLMPIX * MAXCLMLIN];

float inCTSMPCTTREEGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMPCTHERBGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMPCTBAREGrid[MAXCLMPIX * MAXCLMLIN];

float inLUH3PRIMFGrid[MAXCLMPIX * MAXCLMLIN];
float inLUH3PRIMNGrid[MAXCLMPIX * MAXCLMLIN];
float inLUH3SECDFGrid[MAXCLMPIX * MAXCLMLIN];
float inLUH3SECDNGrid[MAXCLMPIX * MAXCLMLIN];
float inLUH3RANGEGrid[MAXCLMPIX * MAXCLMLIN];
float inLUH3PASTRGrid[MAXCLMPIX * MAXCLMLIN];

float inPOTVEGFORESTPCTTREEGrid[MAXCLMPIX * MAXCLMLIN];
float inPOTVEGFORESTPCTHERBGrid[MAXCLMPIX * MAXCLMLIN];
float inPOTVEGFORESTPCTBAREGrid[MAXCLMPIX * MAXCLMLIN];

float inPOTVEGNONFORESTPCTTREEGrid[MAXCLMPIX * MAXCLMLIN];
float inPOTVEGNONFORESTPCTHERBGrid[MAXCLMPIX * MAXCLMLIN];
float inPOTVEGNONFORESTPCTBAREGrid[MAXCLMPIX * MAXCLMLIN];

float inPOTVEGRANGEPCTTREEGrid[MAXCLMPIX * MAXCLMLIN];
float inPOTVEGRANGEPCTHERBGrid[MAXCLMPIX * MAXCLMLIN];
float inPOTVEGRANGEPCTBAREGrid[MAXCLMPIX * MAXCLMLIN];

float inTEMPAVGGrid[MAXCLMPIX * MAXCLMLIN];
float inTEMPWARMGrid[MAXCLMPIX * MAXCLMLIN];
float inTEMPCOLDGrid[MAXCLMPIX * MAXCLMLIN];
float inTEMPGDDGrid[MAXCLMPIX * MAXCLMLIN];
float inPRECIPANNGrid[MAXCLMPIX * MAXCLMLIN];
float inPRECIPMINGrid[MAXCLMPIX * MAXCLMLIN];
float inPRECIPWINGrid[MAXCLMPIX * MAXCLMLIN];
float inPRECIPMAXTG22Grid[MAXCLMPIX * MAXCLMLIN];

float inC3LAIGrid[MAXCLMPIX * MAXCLMLIN];
float inC4LAIGrid[MAXCLMPIX * MAXCLMLIN];
float inMINLAIGrid[MAXCLMPIX * MAXCLMLIN];
float inMAXLAIGrid[MAXCLMPIX * MAXCLMLIN];

float outCTSMLANDMASKGrid[MAXCLMPIX * MAXCLMLIN];
float outCTSMLANDFRACGrid[MAXCLMPIX * MAXCLMLIN];
float outCTSMAREAGrid[MAXCLMPIX * MAXCLMLIN];
float outCTSMPCTGLACIERGrid[MAXCLMPIX * MAXCLMLIN];
float outCTSMPCTLAKEGrid[MAXCLMPIX * MAXCLMLIN];
float outCTSMPCTWETLANDGrid[MAXCLMPIX * MAXCLMLIN];
float outCTSMPCTURBANGrid[MAXCLMPIX * MAXCLMLIN];
float outCTSMPCTNATVEGGrid[MAXCLMPIX * MAXCLMLIN];
float outCTSMPCTCROPGrid[MAXCLMPIX * MAXCLMLIN];
float outCTSMPCTPFTGrid[MAXPFT][MAXCLMPIX * MAXCLMLIN];
float outCTSMPCTCFTGrid[MAXCFT][MAXCLMPIX * MAXCLMLIN];

float outFORESTPCTTREEGrid[MAXCLMPIX * MAXCLMLIN];
float outFORESTPCTHERBGrid[MAXCLMPIX * MAXCLMLIN];
float outFORESTPCTBAREGrid[MAXCLMPIX * MAXCLMLIN];

float outNONFORESTPCTTREEGrid[MAXCLMPIX * MAXCLMLIN];
float outNONFORESTPCTHERBGrid[MAXCLMPIX * MAXCLMLIN];
float outNONFORESTPCTBAREGrid[MAXCLMPIX * MAXCLMLIN];

float outRANGEPCTTREEGrid[MAXCLMPIX * MAXCLMLIN];
float outRANGEPCTHERBGrid[MAXCLMPIX * MAXCLMLIN];
float outRANGEPCTBAREGrid[MAXCLMPIX * MAXCLMLIN];

int readnamelistfile(char *namelist) {

  FILE *namelistfile;
  char templine[1024];
  char fieldname[256];
  char *token;
  int tokencount;

  printf("Reading Namelist: %s\n",namelist);
  namelistfile = fopen(namelist,"r");
  
  sprintf(namelistauthor,"");
  fgets(templine,sizeof(templine),namelistfile);
  token = strtok(templine," ");
  tokencount = 1;
  while( token != NULL ) {
     if (tokencount == 2) {
         sprintf(namelistauthor,"%s",token);
     }
     if (tokencount > 2 && strcmp(token,"\n") != 0) {
        sprintf(namelistauthor,"%s %s",namelistauthor,token);
     }
     token = strtok(NULL," ");
     tokencount++;
  }
  fscanf(namelistfile,"%s %s",fieldname,namelistforestallmaxvcfname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfoutput);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfoutputforestname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfoutputnonforestname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfoutputrangename);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5input);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5inputname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5inputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawmaxvcfinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawmaxvcfinputforestname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawmaxvcfinputnonforestname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawmaxvcfinputrangename);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawmaxvcfinputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawmaxvcfinputmaxtreescale);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5modisrawinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5modisrawinputreferenceyear);
  
  return 0;

}


int readclm5modisrawinputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *clm5modisrawinputvaluefile;
  int foundclm5modisrawinput, itemsfound;
  char fullfilenamestr[1024], inclm5modisrawinputname[1024], inclm5modisrawinputdb[1024];
  char inclm5modisrawinputstartyear[1024], inclm5modisrawinputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  clm5modisrawinputvaluefile = fopen(fullfilenamestr,"r");
  
  foundclm5modisrawinput = 0;

  while (foundclm5modisrawinput == 0) {
      itemsfound = fscanf(clm5modisrawinputvaluefile,"%s %s %s %s",inclm5modisrawinputname, inclm5modisrawinputdb, inclm5modisrawinputstartyear, inclm5modisrawinputendyear);
      if (itemsfound != 4) {
          printf("Error: clm5modisrawinput %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inclm5modisrawinputname) == 0) {
          printf("Processing clm5modisrawinput: %s\n",inclm5modisrawinputname);
	  sprintf(clm5modisrawinputname,"%s",inclm5modisrawinputname);
          if (strcmp(inclm5modisrawinputdb,"<timeseriesrawdir>") == 0) {
              sprintf(clm5modisrawinputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(clm5modisrawinputdb,"%s",inclm5modisrawinputdb);
	  }
          clm5modisrawinputstartyear = atoi(inclm5modisrawinputstartyear);
          clm5modisrawinputendyear = atoi(inclm5modisrawinputendyear);
	  foundclm5modisrawinput = 1;
      }
  }  
  
  return 0;

}


int readworkrawclm5inputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *workrawclm5inputvaluefile;
  int foundworkrawclm5input, itemsfound;
  char fullfilenamestr[1024], inworkrawclm5inputname[1024], inworkrawclm5inputdb[1024];
  char inworkrawclm5inputstartyear[1024], inworkrawclm5inputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  workrawclm5inputvaluefile = fopen(fullfilenamestr,"r");
  
  foundworkrawclm5input = 0;

  while (foundworkrawclm5input == 0) {
      itemsfound = fscanf(workrawclm5inputvaluefile,"%s %s %s %s",inworkrawclm5inputname, inworkrawclm5inputdb, inworkrawclm5inputstartyear, inworkrawclm5inputendyear);
      if (itemsfound != 4) {
          printf("Error: workrawclm5input %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inworkrawclm5inputname) == 0) {
          printf("Processing workrawclm5input: %s\n",inworkrawclm5inputname);
	  sprintf(workrawclm5inputname,"%s",inworkrawclm5inputname);
          if (strcmp(inworkrawclm5inputdb,"<workrawdir>") == 0) {
              sprintf(workrawclm5inputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(workrawclm5inputdb,"%s",inworkrawclm5inputdb);
	  }
          workrawclm5inputstartyear = atoi(inworkrawclm5inputstartyear);
          workrawclm5inputendyear = atoi(inworkrawclm5inputendyear);
	  foundworkrawclm5input = 1;
      }
  }  
  
  return 0;

}


int readworkrawmaxvcfinputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *workrawmaxvcfinputvaluefile;
  int foundworkrawmaxvcfinput, itemsfound;
  char fullfilenamestr[1024], inworkrawmaxvcfinputname[1024], inworkrawmaxvcfinputdb[1024];
  char inworkrawmaxvcfinputstartyear[1024], inworkrawmaxvcfinputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  workrawmaxvcfinputvaluefile = fopen(fullfilenamestr,"r");
  
  foundworkrawmaxvcfinput = 0;

  while (foundworkrawmaxvcfinput == 0) {
      itemsfound = fscanf(workrawmaxvcfinputvaluefile,"%s %s %s %s",inworkrawmaxvcfinputname, inworkrawmaxvcfinputdb, inworkrawmaxvcfinputstartyear, inworkrawmaxvcfinputendyear);
      if (itemsfound != 4) {
          printf("Error: workrawmaxvcfinput %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inworkrawmaxvcfinputname) == 0) {
          printf("Processing workrawmaxvcfinput: %s\n",inworkrawmaxvcfinputname);
	  sprintf(workrawmaxvcfinputname,"%s",inworkrawmaxvcfinputname);
          if (strcmp(inworkrawmaxvcfinputdb,"<workrawdir>") == 0) {
              sprintf(workrawmaxvcfinputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(workrawmaxvcfinputdb,"%s",inworkrawmaxvcfinputdb);
	  }
          workrawmaxvcfinputstartyear = atoi(inworkrawmaxvcfinputstartyear);
          workrawmaxvcfinputendyear = atoi(inworkrawmaxvcfinputendyear);
	  foundworkrawmaxvcfinput = 1;
      }
  }  
  
  return 0;

}


int readworkrawvcfoutputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *workrawvcfoutputvaluefile;
  int foundworkrawvcfoutput, itemsfound;
  char fullfilenamestr[1024], inworkrawvcfoutputname[1024], inworkrawvcfoutputdb[1024];
  char inworkrawvcfoutputstartyear[1024], inworkrawvcfoutputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  workrawvcfoutputvaluefile = fopen(fullfilenamestr,"r");
  
  foundworkrawvcfoutput = 0;

  while (foundworkrawvcfoutput == 0) {
      itemsfound = fscanf(workrawvcfoutputvaluefile,"%s %s %s %s",inworkrawvcfoutputname, inworkrawvcfoutputdb, inworkrawvcfoutputstartyear, inworkrawvcfoutputendyear);
      if (itemsfound != 4) {
          printf("Error: workrawvcfoutput %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inworkrawvcfoutputname) == 0) {
          printf("Processing workrawvcfoutput: %s\n",inworkrawvcfoutputname);
	  sprintf(workrawvcfoutputname,"%s",inworkrawvcfoutputname);
          if (strcmp(inworkrawvcfoutputdb,"<workrawdir>") == 0) {
              sprintf(workrawvcfoutputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(workrawvcfoutputdb,"%s",inworkrawvcfoutputdb);
	  }
          workrawvcfoutputstartyear = atoi(inworkrawvcfoutputstartyear);
          workrawvcfoutputendyear = atoi(inworkrawvcfoutputendyear);
	  foundworkrawvcfoutput = 1;
      }
  }  
  
  return 0;

}


int initializeGrids() {

  long ctsmlin, ctsmpix;
  long pftid, cftid;
  
  for (ctsmlin = 0; ctsmlin < MAXCLMLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXCLMPIX; ctsmpix++) {
          outCTSMLANDMASKGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
          outCTSMLANDFRACGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
          outCTSMAREAGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
          outCTSMPCTGLACIERGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
          outCTSMPCTLAKEGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
          outCTSMPCTWETLANDGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
          outCTSMPCTURBANGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
          outCTSMPCTNATVEGGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
          outCTSMPCTCROPGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;

          outFORESTPCTTREEGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
          outFORESTPCTHERBGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
          outFORESTPCTBAREGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;

          outNONFORESTPCTTREEGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
          outNONFORESTPCTHERBGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
          outNONFORESTPCTBAREGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;

          outRANGEPCTTREEGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
          outRANGEPCTHERBGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
          outRANGEPCTBAREGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
      }
  }

  return 0;
  	      
}


int setoutputregion(char *databasestr, char *filenamestr) {

  FILE *pftparamfile;
  int inpft, inpftid;
  char fullfilenamestr[256], inPFTluhtype[256];
  float lllon, lllat, urlon, urlat;
  
  sprintf(fullfilenamestr,"%s/%s",databasestr,filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  pftparamfile = fopen(fullfilenamestr,"r");

  fscanf(pftparamfile,"%f",&lllon);
  fscanf(pftparamfile,"%f",&lllat);
  fscanf(pftparamfile,"%f",&urlon);
  fscanf(pftparamfile,"%f",&urlat);
  
  MAXOUTPIX = (long) (urlon - lllon) / OUTPIXSIZE;
  MAXOUTLIN = (long) (urlat - lllat) / OUTPIXSIZE;
  
  OUTLLX = lllon;
  OUTLLY = lllat;
  
  OUTLATOFFSET = (long) (90.0 - urlat) / OUTPIXSIZE;
  OUTLONOFFSET = (long) (lllon + 180.0) / OUTPIXSIZE;

  OUTDATASIZE = MAXOUTPIX * MAXOUTLIN * sizeof(float);
  
  return 0;

}


int readlandGrids(char *databasestr, char *seriesname, char *yearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  long outlin, outpix;
  double tempvalue, tempavgvalue, tempcount, tempnewavgvalue;
  
 
  sprintf(globalbinfilename,"%s/%s/%s.LANDMASK.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inLANDMASKGrid,sizeof(inLANDMASKGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.LANDFRAC.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inLANDFRACGrid,sizeof(inLANDFRACGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.AREA.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inAREAGrid,sizeof(inAREAGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTGLACIER.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inCTSMPCTGLACIERGrid,sizeof(inCTSMPCTGLACIERGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTLAKE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inCTSMPCTLAKEGrid,sizeof(inCTSMPCTLAKEGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTWETLAND.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inCTSMPCTWETLANDGrid,sizeof(inCTSMPCTWETLANDGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTURBAN.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inCTSMPCTURBANGrid,sizeof(inCTSMPCTURBANGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTNATVEG.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inCTSMPCTNATVEGGrid,sizeof(inCTSMPCTNATVEGGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTCROP.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inCTSMPCTCROPGrid,sizeof(inCTSMPCTCROPGrid),1,globalbinfile);  
  fclose(globalbinfile);

  return 0;
  
}


int readcurrentluh3Grids(char *databasestr,char *seriesname,char *referenceyearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long ctsmlin, ctsmpix;

  sprintf(globalbinfilename,"%s/%s/%s.primf.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inLUH3PRIMFGrid,sizeof(inLUH3PRIMFGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.primn.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inLUH3PRIMNGrid,sizeof(inLUH3PRIMNGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.secdf.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inLUH3SECDFGrid,sizeof(inLUH3SECDFGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.secdn.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inLUH3SECDNGrid,sizeof(inLUH3SECDNGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.range.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inLUH3RANGEGrid,sizeof(inLUH3RANGEGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.pastr.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inLUH3PASTRGrid,sizeof(inLUH3PASTRGrid),1,globalbinfile);  
  fclose(globalbinfile);

  return 0;
  
}


int readcurrentvcfGrids(char *databasestr,char *seriesname,char *referenceyearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long ctsmlin, ctsmpix;

  sprintf(globalbinfilename,"%s/%s/%s.PCTTREE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inCTSMPCTTREEGrid,sizeof(inCTSMPCTTREEGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTHERB.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inCTSMPCTHERBGrid,sizeof(inCTSMPCTHERBGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTBARE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inCTSMPCTBAREGrid,sizeof(inCTSMPCTBAREGrid),1,globalbinfile);  
  fclose(globalbinfile);

  return 0;
  
}


int readmaxforestvcfGrids(char *databasestr,char *seriesname,char *referenceyearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long ctsmlin, ctsmpix;

  sprintf(globalbinfilename,"%s/%s/%s.PCTTREE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPOTVEGFORESTPCTTREEGrid,sizeof(inPOTVEGFORESTPCTTREEGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTHERB.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPOTVEGFORESTPCTHERBGrid,sizeof(inPOTVEGFORESTPCTHERBGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTBARE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPOTVEGFORESTPCTBAREGrid,sizeof(inPOTVEGFORESTPCTBAREGrid),1,globalbinfile);  
  fclose(globalbinfile);

  return 0;
  
}


int readmaxnonforestvcfGrids(char *databasestr,char *seriesname,char *referenceyearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long ctsmlin, ctsmpix;

  sprintf(globalbinfilename,"%s/%s/%s.PCTTREE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPOTVEGNONFORESTPCTTREEGrid,sizeof(inPOTVEGNONFORESTPCTTREEGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTHERB.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPOTVEGNONFORESTPCTHERBGrid,sizeof(inPOTVEGNONFORESTPCTHERBGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTBARE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPOTVEGNONFORESTPCTBAREGrid,sizeof(inPOTVEGNONFORESTPCTBAREGrid),1,globalbinfile);  
  fclose(globalbinfile);

  return 0;
  
}


int readmaxrangevcfGrids(char *databasestr,char *seriesname,char *referenceyearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long ctsmlin, ctsmpix;

  sprintf(globalbinfilename,"%s/%s/%s.PCTTREE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPOTVEGRANGEPCTTREEGrid,sizeof(inPOTVEGRANGEPCTTREEGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTHERB.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPOTVEGRANGEPCTHERBGrid,sizeof(inPOTVEGRANGEPCTHERBGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTBARE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPOTVEGRANGEPCTBAREGrid,sizeof(inPOTVEGRANGEPCTBAREGrid),1,globalbinfile);  
  fclose(globalbinfile);

  return 0;
  
}


int readclimvarseries(char *databasestr, char *seriesname, char *fileyear) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  
  sprintf(globalbinfilename,"%s/%s/%s.tempaverage-%s.dat",databasestr,seriesname,seriesname,fileyear);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inTEMPAVGGrid,sizeof(inTEMPAVGGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.tempwarmest-%s.dat",databasestr,seriesname,seriesname,fileyear);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inTEMPWARMGrid,sizeof(inTEMPWARMGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.tempcoldest-%s.dat",databasestr,seriesname,seriesname,fileyear);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inTEMPCOLDGrid,sizeof(inTEMPCOLDGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.growdegdays-%s.dat",databasestr,seriesname,seriesname,fileyear);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inTEMPGDDGrid,sizeof(inTEMPGDDGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.precipann-%s.dat",databasestr,seriesname,seriesname,fileyear);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPRECIPANNGrid,sizeof(inPRECIPANNGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.precipmin-%s.dat",databasestr,seriesname,seriesname,fileyear);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPRECIPMINGrid,sizeof(inPRECIPMINGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.precipwin-%s.dat",databasestr,seriesname,seriesname,fileyear);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPRECIPWINGrid,sizeof(inPRECIPWINGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.precipmaxtg22-%s.dat",databasestr,seriesname,seriesname,fileyear);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPRECIPMAXTG22Grid,sizeof(inPRECIPMAXTG22Grid),1,globalbinfile);  
  fclose(globalbinfile);

  return 0;

}


int readlaivarseries(char *databasestr, char *seriesname, char *yearname) {

  FILE *globalbinlaifile;
  char filedate[256];
  char globalbinlaifilename[256];

  sprintf(globalbinlaifilename,"%s/%s/%s.MinLai_1km.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinlaifilename);
  globalbinlaifile = fopen(globalbinlaifilename,"r");
  fread(inMINLAIGrid,sizeof(inMINLAIGrid),1,globalbinlaifile);  
  fclose(globalbinlaifile);

  sprintf(globalbinlaifilename,"%s/%s/%s.MaxLai_1km.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinlaifilename);
  globalbinlaifile = fopen(globalbinlaifilename,"r");
  fread(inMAXLAIGrid,sizeof(inMAXLAIGrid),1,globalbinlaifile);  
  fclose(globalbinlaifile);

  sprintf(globalbinlaifilename,"%s/%s/%s.C3Lai_1km.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinlaifilename);
  globalbinlaifile = fopen(globalbinlaifilename,"r");
  fread(inC3LAIGrid,sizeof(inC3LAIGrid),1,globalbinlaifile);  
  fclose(globalbinlaifile);

  sprintf(globalbinlaifilename,"%s/%s/%s.C4Lai_1km.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinlaifilename);
  globalbinlaifile = fopen(globalbinlaifilename,"r");
  fread(inC4LAIGrid,sizeof(inC4LAIGrid),1,globalbinlaifile);  
  fclose(globalbinlaifile);

  return 0;

}


int genLUHLandGrids() {

  long ctsmlin, ctsmpix;
  long pftid, cftid;
  float landmask, landfrac, area, glacierpct, lakepct, wetlandpct, natvegpct;
  for (ctsmlin = 0; ctsmlin < MAXCLMLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXCLMPIX; ctsmpix++) {
          landmask = inLANDMASKGrid[ctsmlin * MAXCLMPIX + ctsmpix];
          landfrac = inLANDFRACGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	  if (landfrac > 0.0) {
              outCTSMLANDMASKGrid[ctsmlin * MAXCLMPIX + ctsmpix] = landmask;
              outCTSMLANDFRACGrid[ctsmlin * MAXCLMPIX + ctsmpix] = landfrac;
	      outCTSMAREAGrid[ctsmlin * MAXCLMPIX + ctsmpix] = inAREAGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      glacierpct = inCTSMPCTGLACIERGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      lakepct = inCTSMPCTLAKEGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      wetlandpct = inCTSMPCTWETLANDGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      natvegpct = 100.0 - glacierpct - lakepct - wetlandpct;
              outCTSMPCTGLACIERGrid[ctsmlin * MAXCLMPIX + ctsmpix] = glacierpct;
              outCTSMPCTLAKEGrid[ctsmlin * MAXCLMPIX + ctsmpix] = lakepct;
              outCTSMPCTWETLANDGrid[ctsmlin * MAXCLMPIX + ctsmpix] = wetlandpct;
              outCTSMPCTNATVEGGrid[ctsmlin * MAXCLMPIX + ctsmpix] = natvegpct;
	  }
      }
  }

  return 0;
  	      
}


int genLUHMaxVegGrids() {

  long ctsmlin, ctsmpix;
  long pftid, cftid;
  float landmask, landfrac, area, precipann, tempaverage;
  float luh3forest, luh3nonforest, luh3range, luh3pasture, luh3total;
  float luh3treepct, luh3herbpct, luh3barepct;
  float luh3scaledfrac;
  float currentnatvegfrac, currenttreepct, currentherbpct, currentbarepct;
  float currentforesttreepct, currentnonforesttreepct, currentrangetreepct;
  float currentforestherbpct, currentnonforestherbpct, currentrangeherbpct;
  float currentforestbarepct, currentnonforestbarepct, currentrangebarepct;
  float foresttreepct, forestherbpct, forestbarepct;
  float nonforesttreepct, nonforestherbpct, nonforestbarepct;
  float rangetreepct, rangeherbpct, rangebarepct;
  float scaledforesttreepct, scaledforestherbpct, scaledforestbarepct;
  float scaledforesttreefrac, scaledforestherbfrac, scaledforestbarefrac;
  float scalednonforesttreepct, scalednonforestherbpct, scalednonforestbarepct;
  float scalednonforesttreefrac, scalednonforestherbfrac, scalednonforestbarefrac;
  float scaledrangetreepct, scaledrangeherbpct, scaledrangebarepct;
  float scaledrangetreefrac, scaledrangeherbfrac, scaledrangebarefrac;
  float potentialtreepct, potentialforesttreepct, potentialnonforesttreepct, potentialrangetreepct;
  float potentialherbpct, potentialforestherbpct, potentialnonforestherbpct, potentialrangeherbpct;
  float potentialbarepct, potentialforestbarepct, potentialnonforestbarepct, potentialrangebarepct;
  float remainingforestpct, otherforestpct;
  float remainingnonforestpct, othernonforestpct;
  float remainingrangepct, otherrangepct;
  float removetreepct, removeluh3treepct, addtreepct, addluh3treepct;
  float removeherbpct, removeluh3herbpct, addherbpct, addluh3herbpct;
  float removebarepct, removeluh3barepct, addbarepct, addluh3barepct;
  int precipannindex, tempaverageindex;

  for (ctsmlin = 0; ctsmlin < MAXCLMLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXCLMPIX; ctsmpix++) {
          landmask = inLANDMASKGrid[ctsmlin * MAXCLMPIX + ctsmpix];
          landfrac = inLANDFRACGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	  if (landfrac > 0.0) {
	  
	      foresttreepct = inPOTVEGFORESTPCTTREEGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      forestherbpct = inPOTVEGFORESTPCTHERBGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      forestbarepct = inPOTVEGFORESTPCTBAREGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      
	      nonforesttreepct = inPOTVEGNONFORESTPCTTREEGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      nonforestherbpct = inPOTVEGNONFORESTPCTHERBGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      nonforestbarepct = inPOTVEGNONFORESTPCTBAREGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      
	      rangetreepct = inPOTVEGRANGEPCTTREEGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      rangeherbpct = inPOTVEGRANGEPCTHERBGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      rangebarepct = inPOTVEGRANGEPCTBAREGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      
	      scaledforesttreepct = workrawmaxvcfinputmaxtreescale * foresttreepct;
	      if (scaledforesttreepct > 100.0) {
	          scaledforesttreepct = 100.0;
              }
	      remainingforestpct = 100.0 - scaledforesttreepct;
	      otherforestpct = forestherbpct + forestbarepct;
	      if (otherforestpct > 0.0) {
	          scaledforestherbpct = remainingforestpct * forestherbpct / otherforestpct;
	          scaledforestbarepct = remainingforestpct * forestbarepct / otherforestpct;
              }
	      else {
	          scaledforestherbpct = remainingforestpct;
		  scaledforestbarepct = 0.0;
              }
	      
	      scalednonforesttreepct = workrawmaxvcfinputmaxtreescale * nonforesttreepct;
	      if (scalednonforesttreepct > 100.0) {
	          scalednonforesttreepct = 100.0;
              }
	      remainingnonforestpct = 100.0 - scalednonforesttreepct;
	      othernonforestpct = nonforestherbpct + nonforestbarepct;
	      if (othernonforestpct > 0.0) {
	          scalednonforestherbpct = remainingnonforestpct * nonforestherbpct / othernonforestpct;
	          scalednonforestbarepct = remainingnonforestpct * nonforestbarepct / othernonforestpct;
              }
	      else {
	          scalednonforestherbpct = remainingnonforestpct;
		  scalednonforestbarepct = 0.0;
              }
	      
	      scaledrangetreepct = workrawmaxvcfinputmaxtreescale * rangetreepct;
	      if (scaledrangetreepct > 100.0) {
	          scaledrangetreepct = 100.0;
              }
	      remainingrangepct = 100.0 - scaledrangetreepct;
	      otherrangepct = rangeherbpct + rangebarepct;
	      if (otherrangepct > 0.0) {
	          scaledrangeherbpct = remainingrangepct * rangeherbpct / otherrangepct;
	          scaledrangebarepct = remainingrangepct * rangebarepct / otherrangepct;
              }
	      else {
	          scaledrangeherbpct = remainingrangepct;
		  scaledrangebarepct = 0.0;
              }
	      
              if (scaledforesttreepct > foresttreepct) {
	         foresttreepct = scaledforesttreepct;
		 forestherbpct = scaledforestherbpct;
		 forestbarepct = scaledforestbarepct;
              }
	      
              if (scalednonforesttreepct > nonforesttreepct) {
	         nonforesttreepct = scalednonforesttreepct;
		 nonforestherbpct = scalednonforestherbpct;
		 nonforestbarepct = scalednonforestbarepct;
              }

              if (scaledrangetreepct > rangetreepct) {
	         rangetreepct = scaledrangetreepct;
		 rangeherbpct = scaledrangeherbpct;
		 rangebarepct = scaledrangebarepct;
              }
	      
              outFORESTPCTTREEGrid[ctsmlin * MAXCLMPIX + ctsmpix] = foresttreepct;
              outFORESTPCTHERBGrid[ctsmlin * MAXCLMPIX + ctsmpix] = forestherbpct;
	      outFORESTPCTBAREGrid[ctsmlin * MAXCLMPIX + ctsmpix] = forestbarepct;
	      
              outNONFORESTPCTTREEGrid[ctsmlin * MAXCLMPIX + ctsmpix] = nonforesttreepct;
              outNONFORESTPCTHERBGrid[ctsmlin * MAXCLMPIX + ctsmpix] = nonforestherbpct;
	      outNONFORESTPCTBAREGrid[ctsmlin * MAXCLMPIX + ctsmpix] = nonforestbarepct;
              
	      outRANGEPCTTREEGrid[ctsmlin * MAXCLMPIX + ctsmpix] = rangetreepct;
              outRANGEPCTHERBGrid[ctsmlin * MAXCLMPIX + ctsmpix] = rangeherbpct;
	      outRANGEPCTBAREGrid[ctsmlin * MAXCLMPIX + ctsmpix] = rangebarepct;
	  }
      }
  }

  return 0;
  	      
}


int rebalanceLUHMaxVegGrids() {

  long ctsmlin, ctsmpix;
  long pftid, cftid;
  float landmask, landfrac, area, precipann, tempaverage;
  float forestallpct, foresttreepct, forestherbpct, forestbarepct;
  float nonforestallpct, nonforesttreepct, nonforestherbpct, nonforestbarepct;
  float rangeallpct, rangetreepct, rangeherbpct, rangebarepct;
  float removepct, addpct;
  int precipannindex, tempaverageindex;

  for (ctsmlin = 0; ctsmlin < MAXCLMLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXCLMPIX; ctsmpix++) {
          landmask = inLANDMASKGrid[ctsmlin * MAXCLMPIX + ctsmpix];
          landfrac = inLANDFRACGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	  if (landfrac > 0.0) {
	  
	      foresttreepct = outFORESTPCTTREEGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      forestherbpct = outFORESTPCTHERBGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      forestbarepct = outFORESTPCTBAREGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      
	      nonforesttreepct = outNONFORESTPCTTREEGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      nonforestherbpct = outNONFORESTPCTHERBGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      nonforestbarepct = outNONFORESTPCTBAREGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      
	      rangetreepct = outRANGEPCTTREEGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      rangeherbpct = outRANGEPCTHERBGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      rangebarepct = outRANGEPCTBAREGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      
	      if (foresttreepct < nonforesttreepct) {
	          foresttreepct = nonforesttreepct;
              }
	      
	      if (foresttreepct < rangetreepct) {
	          foresttreepct = rangetreepct;
              }
	      
	      if (forestherbpct > nonforestherbpct) {
	          forestherbpct = nonforestherbpct;
              }
	      
	      if (forestherbpct > rangeherbpct) {
	          forestherbpct = rangeherbpct;
              }
	      
              foresttreepct = trunc(foresttreepct * CLMPCTTruncate) / CLMPCTTruncate;
              forestherbpct = trunc(forestherbpct * CLMPCTTruncate) / CLMPCTTruncate;
              forestbarepct = trunc(forestbarepct * CLMPCTTruncate) / CLMPCTTruncate;

	      forestallpct = foresttreepct + forestherbpct + forestbarepct;	      
	      if (forestallpct > 100.0) {
	          removepct = forestallpct - 100.0;
		  if (forestbarepct > removepct) {
		      forestbarepct -= removepct;
		      removepct = 0.0;
		  }
		  else {
		      removepct -= forestbarepct;
		      forestbarepct = 0.0;
		  }
		  if (removepct > 0.0) {
		      if (forestherbpct > removepct) {
		          forestherbpct -= removepct;
			  removepct = 0.0;
		      }
		      else {
		          removepct -= forestherbpct;
	                  forestherbpct = 0.0;
                      }
		      if (removepct > 0.0) {
		          if (foresttreepct > removepct) {
			      foresttreepct -= removepct;
			      removepct = 0.0;
			  }
			  else {
			      removepct -= foresttreepct;
			      foresttreepct = 0.0;
			  }
		      }
		  }
              }
	      
	      forestallpct = foresttreepct + forestherbpct + forestbarepct;
	      if (forestallpct < 100.0) {
	          addpct = 100.0 - forestallpct;
		  if (foresttreepct > forestherbpct) {
		      foresttreepct += addpct;
		      addpct = 0.0;
		  }
		  else {
		      if (forestherbpct > forestbarepct) {
		          forestherbpct += addpct;
		          addpct = 0.0;
		      }
		      else {
		          forestbarepct += addpct;
	                  addpct = 0.0;
                      }
		  }
              }

	      
	      if (nonforesttreepct < rangetreepct) {
	          nonforesttreepct = rangetreepct;
              }
	      
	      if (nonforestherbpct > rangeherbpct) {
	          nonforestherbpct = rangeherbpct;
              }
	      
              nonforesttreepct = trunc(nonforesttreepct * CLMPCTTruncate) / CLMPCTTruncate;
              nonforestherbpct = trunc(nonforestherbpct * CLMPCTTruncate) / CLMPCTTruncate;
              nonforestbarepct = trunc(nonforestbarepct * CLMPCTTruncate) / CLMPCTTruncate;

	      nonforestallpct = nonforesttreepct + nonforestherbpct + nonforestbarepct;	      
	      if (nonforestallpct > 100.0) {
	          removepct = nonforestallpct - 100.0;
		  if (nonforestbarepct > removepct) {
		      nonforestbarepct -= removepct;
		      removepct = 0.0;
		  }
		  else {
		      removepct -= nonforestbarepct;
		      nonforestbarepct = 0.0;
		  }
		  if (removepct > 0.0) {
		      if (nonforestherbpct > removepct) {
		          nonforestherbpct -= removepct;
			  removepct = 0.0;
		      }
		      else {
		          removepct -= nonforestherbpct;
	                  nonforestherbpct = 0.0;
                      }
		      if (removepct > 0.0) {
		          if (nonforesttreepct > removepct) {
			      nonforesttreepct -= removepct;
			      removepct = 0.0;
			  }
			  else {
			      removepct -= nonforesttreepct;
			      nonforesttreepct = 0.0;
			  }
		      }
		  }
              }
	      
	      nonforestallpct = nonforesttreepct + nonforestherbpct + nonforestbarepct;
	      if (nonforestallpct < 100.0) {
	          addpct = 100.0 - nonforestallpct;
		  if (nonforesttreepct > nonforestherbpct) {
		      nonforesttreepct += addpct;
		      addpct = 0.0;
		  }
		  else {
		      if (nonforestherbpct > nonforestbarepct) {
		          nonforestherbpct += addpct;
		          addpct = 0.0;
		      }
		      else {
		          nonforestbarepct += addpct;
	                  addpct = 0.0;
                      }
		  }
              }

              rangetreepct = trunc(rangetreepct * CLMPCTTruncate) / CLMPCTTruncate;
              rangeherbpct = trunc(rangeherbpct * CLMPCTTruncate) / CLMPCTTruncate;
              rangebarepct = trunc(rangebarepct * CLMPCTTruncate) / CLMPCTTruncate;

	      rangeallpct = rangetreepct + rangeherbpct + rangebarepct;	      
	      if (rangeallpct > 100.0) {
	          removepct = rangeallpct - 100.0;
		  if (rangetreepct > removepct) {
		      rangetreepct -= removepct;
		      removepct = 0.0;
		  }
		  else {
		      removepct -= rangetreepct;
		      rangetreepct = 0.0;
		  }
		  if (removepct > 0.0) {
		      if (rangeherbpct > removepct) {
		          rangeherbpct -= removepct;
			  removepct = 0.0;
		      }
		      else {
		          removepct -= rangeherbpct;
	                  rangeherbpct = 0.0;
                      }
		      if (removepct > 0.0) {
		          if (rangebarepct > removepct) {
			      rangebarepct -= removepct;
			      removepct = 0.0;
			  }
			  else {
			      removepct -= rangebarepct;
			      rangebarepct = 0.0;
			  }
		      }
		  }
              }
	      
	      rangeallpct = rangetreepct + rangeherbpct + rangebarepct;
	      if (rangeallpct < 100.0) {
	          addpct = 100.0 - rangeallpct;
		  if (rangeherbpct > rangetreepct) {
		      rangeherbpct += addpct;
		      addpct = 0.0;
		  }
		  else {
		      if (rangetreepct > rangebarepct) {
		          rangetreepct += addpct;
		          addpct = 0.0;
		      }
		      else {
		          rangebarepct += addpct;
	                  addpct = 0.0;
                      }
		  }
              }
	      
              outFORESTPCTTREEGrid[ctsmlin * MAXCLMPIX + ctsmpix] = foresttreepct;
              outFORESTPCTHERBGrid[ctsmlin * MAXCLMPIX + ctsmpix] = forestherbpct;
	      outFORESTPCTBAREGrid[ctsmlin * MAXCLMPIX + ctsmpix] = forestbarepct;
	      
              outNONFORESTPCTTREEGrid[ctsmlin * MAXCLMPIX + ctsmpix] = nonforesttreepct;
              outNONFORESTPCTHERBGrid[ctsmlin * MAXCLMPIX + ctsmpix] = nonforestherbpct;
	      outNONFORESTPCTBAREGrid[ctsmlin * MAXCLMPIX + ctsmpix] = nonforestbarepct;
              
	      outRANGEPCTTREEGrid[ctsmlin * MAXCLMPIX + ctsmpix] = rangetreepct;
              outRANGEPCTHERBGrid[ctsmlin * MAXCLMPIX + ctsmpix] = rangeherbpct;
	      outRANGEPCTBAREGrid[ctsmlin * MAXCLMPIX + ctsmpix] = rangebarepct;
	  }
      }
  }

  return 0;
  	      
}


int writemaxforestgrids(char *databasestr, char *seriesname, char *yearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  long outlin, outpix;
  long clmlin, clmpix;
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];

  tempoutGrid = (float *) malloc(OUTDATASIZE);
  
  sprintf(globalbinfilename,"%s/%s/%s.LANDMASK.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMLANDMASKGrid,sizeof(outCTSMLANDMASKGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.LANDFRAC.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMLANDFRACGrid,sizeof(outCTSMLANDFRACGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.AREA.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMAREAGrid,sizeof(outCTSMAREAGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTGLACIER.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTGLACIERGrid,sizeof(outCTSMPCTGLACIERGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTLAKE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTLAKEGrid,sizeof(outCTSMPCTLAKEGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTWETLAND.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTWETLANDGrid,sizeof(outCTSMPCTWETLANDGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTURBAN.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTURBANGrid,sizeof(outCTSMPCTURBANGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTNATVEG.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTNATVEGGrid,sizeof(outCTSMPCTNATVEGGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTCROP.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTCROPGrid,sizeof(outCTSMPCTCROPGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTTREE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outFORESTPCTTREEGrid,sizeof(outFORESTPCTTREEGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTHERB.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outFORESTPCTHERBGrid,sizeof(outFORESTPCTHERBGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTBARE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outFORESTPCTBAREGrid,sizeof(outFORESTPCTBAREGrid),1,globalbinfile);  
  fclose(globalbinfile);
      
  return 0;

}


int writemaxnonforestgrids(char *databasestr, char *seriesname, char *yearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  long outlin, outpix;
  long clmlin, clmpix;
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];

  tempoutGrid = (float *) malloc(OUTDATASIZE);
  
  sprintf(globalbinfilename,"%s/%s/%s.LANDMASK.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMLANDMASKGrid,sizeof(outCTSMLANDMASKGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.LANDFRAC.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMLANDFRACGrid,sizeof(outCTSMLANDFRACGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.AREA.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMAREAGrid,sizeof(outCTSMAREAGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTGLACIER.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTGLACIERGrid,sizeof(outCTSMPCTGLACIERGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTLAKE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTLAKEGrid,sizeof(outCTSMPCTLAKEGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTWETLAND.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTWETLANDGrid,sizeof(outCTSMPCTWETLANDGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTURBAN.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTURBANGrid,sizeof(outCTSMPCTURBANGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTNATVEG.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTNATVEGGrid,sizeof(outCTSMPCTNATVEGGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTCROP.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTCROPGrid,sizeof(outCTSMPCTCROPGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTTREE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outNONFORESTPCTTREEGrid,sizeof(outNONFORESTPCTTREEGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTHERB.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outNONFORESTPCTHERBGrid,sizeof(outNONFORESTPCTHERBGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTBARE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outNONFORESTPCTBAREGrid,sizeof(outNONFORESTPCTBAREGrid),1,globalbinfile);  
  fclose(globalbinfile);
      
  return 0;

}


int writemaxrangegrids(char *databasestr, char *seriesname, char *yearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  long outlin, outpix;
  long clmlin, clmpix;
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];

  tempoutGrid = (float *) malloc(OUTDATASIZE);
  
  sprintf(globalbinfilename,"%s/%s/%s.LANDMASK.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMLANDMASKGrid,sizeof(outCTSMLANDMASKGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.LANDFRAC.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMLANDFRACGrid,sizeof(outCTSMLANDFRACGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.AREA.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMAREAGrid,sizeof(outCTSMAREAGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTGLACIER.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTGLACIERGrid,sizeof(outCTSMPCTGLACIERGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTLAKE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTLAKEGrid,sizeof(outCTSMPCTLAKEGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTWETLAND.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTWETLANDGrid,sizeof(outCTSMPCTWETLANDGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTURBAN.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTURBANGrid,sizeof(outCTSMPCTURBANGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTNATVEG.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTNATVEGGrid,sizeof(outCTSMPCTNATVEGGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTCROP.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCTSMPCTCROPGrid,sizeof(outCTSMPCTCROPGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTTREE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outRANGEPCTTREEGrid,sizeof(outRANGEPCTTREEGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTHERB.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outRANGEPCTHERBGrid,sizeof(outRANGEPCTHERBGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTBARE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outRANGEPCTBAREGrid,sizeof(outRANGEPCTBAREGrid),1,globalbinfile);  
  fclose(globalbinfile);
      
  return 0;

}


int main(long narg, char **argv) {

  char workrawclm5inputdirname[1024];
  char workrawmaxvcfinputdirname[1024];
  char workrawvcfoutputdirname[1024];
  char clm5modisrawinputreferenceyearstr[256];
  char workrawclm5inputreferenceyearstr[256];
  char workrawmaxvcfinputreferenceyearstr[256];
  char workrawvcfoutputreferenceyearstr[256];
  
  if(narg != 6){
        printf("Usage createcurrentpftsCTSM52Deg025bin timeseriesrawfile timeseriesdir workrawfile workrawdir forestallmaxvcfnamelistfile\n");
        return 0;
  }
  
  readnamelistfile(argv[5]);

  readclm5modisrawinputlutfile(argv[1],argv[2],namelistclm5modisrawinput);
  readworkrawclm5inputlutfile(argv[3],argv[4],namelistworkrawclm5input);
  readworkrawmaxvcfinputlutfile(argv[3],argv[4],namelistworkrawmaxvcfinput);
  readworkrawvcfoutputlutfile(argv[3],argv[4],namelistworkrawvcfoutput);
  clm5modisrawinputreferenceyear = atoi(namelistclm5modisrawinputreferenceyear);
  sprintf(clm5modisrawinputreferenceyearstr,"%04d",clm5modisrawinputreferenceyear);
  workrawclm5inputreferenceyear = atoi(namelistworkrawclm5inputreferenceyear);
  sprintf(workrawclm5inputreferenceyearstr,"%04d",workrawclm5inputreferenceyear);
  workrawmaxvcfinputreferenceyear = atoi(namelistworkrawmaxvcfinputreferenceyear);
  sprintf(workrawmaxvcfinputreferenceyearstr,"%04d",workrawmaxvcfinputreferenceyear);
  workrawmaxvcfinputmaxtreescale = atof(namelistworkrawmaxvcfinputmaxtreescale);
  
  initializeGrids();
  
  sprintf(workrawclm5inputdirname,"%s/%s",workrawclm5inputdb,workrawclm5inputname);
  sprintf(workrawmaxvcfinputdirname,"%s/%s",workrawmaxvcfinputdb,workrawmaxvcfinputname);
  readlandGrids(workrawclm5inputdirname,namelistworkrawclm5inputname,workrawclm5inputreferenceyearstr);
  readmaxforestvcfGrids(workrawmaxvcfinputdirname,namelistworkrawmaxvcfinputforestname,workrawmaxvcfinputreferenceyearstr);    
  readmaxnonforestvcfGrids(workrawmaxvcfinputdirname,namelistworkrawmaxvcfinputnonforestname,workrawmaxvcfinputreferenceyearstr);    
  readmaxrangevcfGrids(workrawmaxvcfinputdirname,namelistworkrawmaxvcfinputrangename,workrawmaxvcfinputreferenceyearstr);    
  readclimvarseries(clm5modisrawinputdb,clm5modisrawinputname,"CLIM");
  readlaivarseries(clm5modisrawinputdb,clm5modisrawinputname,"FILLCLIM");
  
  genLUHLandGrids();
  genLUHMaxVegGrids();
  rebalanceLUHMaxVegGrids();
  
  sprintf(workrawvcfoutputreferenceyearstr,"%s",workrawclm5inputreferenceyearstr);
  sprintf(workrawvcfoutputdirname,"%s/%s",workrawvcfoutputdb,workrawvcfoutputname);
  writemaxforestgrids(workrawvcfoutputdirname,namelistworkrawvcfoutputforestname,workrawvcfoutputreferenceyearstr);
  writemaxnonforestgrids(workrawvcfoutputdirname,namelistworkrawvcfoutputnonforestname,workrawvcfoutputreferenceyearstr);
  writemaxrangegrids(workrawvcfoutputdirname,namelistworkrawvcfoutputrangename,workrawvcfoutputreferenceyearstr);

  return 1;
  
}
