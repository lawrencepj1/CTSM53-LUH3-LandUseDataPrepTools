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
char namelistluh3descallsecdfallpftname[1024];
char namelistluh3descrawpftoutput[1024];
char namelistluh3descrawpftoutputsecdfname[1024];
char namelistworkrawclm5input[1024];
char namelistworkrawclm5inputname[1024];
char namelistworkrawclm5inputreferenceyear[1024];
char namelistworkrawforestmaxvcfinput[1024];
char namelistworkrawforestmaxvcfinputforestname[1024];
char namelistworkrawforestmaxvcfinputreferenceyear[1024];
char namelistworkrawscaledvcfinput[1024];
char namelistworkrawscaledvcfinputforestname[1024];
char namelistworkrawscaledvcfinputreferenceyear[1024];
char namelistworkrawextrappftinput[1024];
char namelistworkrawextrappftinputtreename[1024];
char namelistworkrawextrappftinputherbname[1024];
char namelistworkrawextrappftinputreferenceyear[1024];
char namelistworkrawmergeinput[1024];
char namelistworkrawmergeinputname[1024];
char namelistworkrawmergeinputstartyear[1024];
char namelistworkrawmergeinputendyear[1024];
char namelistclm5modisrawinput[1024];
char namelistclm5modisrawinputreferenceyear[1024];

int clm5modisrawinputreferenceyear;
int workrawclm5inputreferenceyear;
int workrawforestmaxvcfinputreferenceyear;
int workrawscaledvcfinputreferenceyear;
int workrawextrappftinputreferenceyear;
int workrawmergeinputstartyear;
int workrawmergeinputendyear;

char clm5modisrawinputname[1024];
char clm5modisrawinputdb[1024];
int clm5modisrawinputstartyear;
int clm5modisrawinputendyear;

char workrawclm5inputname[1024];
char workrawclm5inputdb[1024];
int workrawclm5inputstartyear;
int workrawclm5inputendyear;

char workrawforestmaxvcfinputname[1024];
char workrawforestmaxvcfinputdb[1024];
int workrawforestmaxvcfinputstartyear;
int workrawforestmaxvcfinputendyear;

char workrawscaledvcfinputname[1024];
char workrawscaledvcfinputdb[1024];
int workrawscaledvcfinputstartyear;
int workrawscaledvcfinputendyear;

char workrawextrappftinputname[1024];
char workrawextrappftinputdb[1024];
int workrawextrappftinputstartyear;
int workrawextrappftinputendyear;

char workrawmergeinputname[1024];
char workrawmergeinputdb[1024];
int workrawmergeinputstartyear;
int workrawmergeinputendyear;

char luh3descrawpftoutputname[1024];
char luh3descrawpftoutputdb[1024];
int luh3descrawpftoutputstartyear;
int luh3descrawpftoutputendyear;

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

float inCTSMTREEEXTRAPPCTPFTGrid[MAXPFT][MAXCLMPIX * MAXCLMLIN];
float inCTSMTREEEXTRAPPCTNDLEVGTREEGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMTREEEXTRAPPCTNDLDECTREEGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMTREEEXTRAPPCTBRDEVGTREEGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMTREEEXTRAPPCTBRDDECTREEGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMHERBEXTRAPPCTPFTGrid[MAXPFT][MAXCLMPIX * MAXCLMLIN];
float inCTSMHERBEXTRAPPCTSHRGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMHERBEXTRAPPCTGRSC3Grid[MAXCLMPIX * MAXCLMLIN];
float inCTSMHERBEXTRAPPCTGRSC4Grid[MAXCLMPIX * MAXCLMLIN];

float inPOTVEGMAXFORESTPCTTREEGrid[MAXCLMPIX * MAXCLMLIN];
float inPOTVEGMAXFORESTPCTHERBGrid[MAXCLMPIX * MAXCLMLIN];
float inPOTVEGMAXFORESTPCTBAREGrid[MAXCLMPIX * MAXCLMLIN];

float inPOTVEGSCALEDFORESTPCTTREEGrid[MAXCLMPIX * MAXCLMLIN];
float inPOTVEGSCALEDFORESTPCTHERBGrid[MAXCLMPIX * MAXCLMLIN];
float inPOTVEGSCALEDFORESTPCTBAREGrid[MAXCLMPIX * MAXCLMLIN];

float inTREEEXTRAPPCTPFTGrid[MAXPFT][MAXCLMPIX * MAXCLMLIN];
float inHERBEXTRAPPCTPFTGrid[MAXPFT][MAXCLMPIX * MAXCLMLIN];

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

float inMERGEDSECDFFRACTREEGrid[MAXCLMPIX * MAXCLMLIN];
float inMERGEDSECDFFRACHERBGrid[MAXCLMPIX * MAXCLMLIN];
float inMERGEDSECDFFRACBAREGrid[MAXCLMPIX * MAXCLMLIN];

float outDESCMERGEDLANDMASKGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCMERGEDLANDFRACGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCMERGEDAREAGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCMERGEDPCTGLACIERGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCMERGEDPCTLAKEGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCMERGEDPCTWETLANDGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCMERGEDPCTURBANGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCMERGEDPCTNATVEGGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCMERGEDPCTCROPGrid[MAXCLMPIX * MAXCLMLIN];

float outDESCMERGEDPCTPFTGrid[MAXPFT][MAXCLMPIX * MAXCLMLIN];

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
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3descallsecdfallpftname);
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3descrawpftoutput);
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3descrawpftoutputsecdfname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5input);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5inputname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5inputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawforestmaxvcfinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawforestmaxvcfinputforestname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawforestmaxvcfinputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawscaledvcfinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawscaledvcfinputforestname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawscaledvcfinputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawextrappftinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawextrappftinputtreename);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawextrappftinputherbname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawextrappftinputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawmergeinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawmergeinputname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawmergeinputstartyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawmergeinputendyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5modisrawinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5modisrawinputreferenceyear);
  
  return 0;

}


int readluh3descrawpftoutputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *luh3descrawpftoutputvaluefile;
  int foundluh3descrawpftoutput, itemsfound;
  char fullfilenamestr[1024], inluh3descrawpftoutputname[1024], inluh3descrawpftoutputdb[1024];
  char inluh3descrawpftoutputstartyear[1024], inluh3descrawpftoutputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  luh3descrawpftoutputvaluefile = fopen(fullfilenamestr,"r");
  
  foundluh3descrawpftoutput = 0;

  while (foundluh3descrawpftoutput == 0) {
      itemsfound = fscanf(luh3descrawpftoutputvaluefile,"%s %s %s %s",inluh3descrawpftoutputname, inluh3descrawpftoutputdb, inluh3descrawpftoutputstartyear, inluh3descrawpftoutputendyear);
      if (itemsfound != 4) {
          printf("Error: luh3descrawpftoutput %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inluh3descrawpftoutputname) == 0) {
          printf("Processing luh3descrawpftoutput: %s\n",inluh3descrawpftoutputname);
	  sprintf(luh3descrawpftoutputname,"%s",inluh3descrawpftoutputname);
          if (strcmp(inluh3descrawpftoutputdb,"<luh3descrawdir>") == 0) {
              sprintf(luh3descrawpftoutputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(luh3descrawpftoutputdb,"%s",inluh3descrawpftoutputdb);
	  }
          luh3descrawpftoutputstartyear = atoi(inluh3descrawpftoutputstartyear);
          luh3descrawpftoutputendyear = atoi(inluh3descrawpftoutputendyear);
	  foundluh3descrawpftoutput = 1;
      }
  }  
  
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


int readworkrawforestmaxvcfinputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *workrawforestmaxvcfinputvaluefile;
  int foundworkrawforestmaxvcfinput, itemsfound;
  char fullfilenamestr[1024], inworkrawforestmaxvcfinputname[1024], inworkrawforestmaxvcfinputdb[1024];
  char inworkrawforestmaxvcfinputstartyear[1024], inworkrawforestmaxvcfinputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  workrawforestmaxvcfinputvaluefile = fopen(fullfilenamestr,"r");
  
  foundworkrawforestmaxvcfinput = 0;

  while (foundworkrawforestmaxvcfinput == 0) {
      itemsfound = fscanf(workrawforestmaxvcfinputvaluefile,"%s %s %s %s",inworkrawforestmaxvcfinputname, inworkrawforestmaxvcfinputdb, inworkrawforestmaxvcfinputstartyear, inworkrawforestmaxvcfinputendyear);
      if (itemsfound != 4) {
          printf("Error: workrawforestmaxvcfinput %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inworkrawforestmaxvcfinputname) == 0) {
          printf("Processing workrawforestmaxvcfinput: %s\n",inworkrawforestmaxvcfinputname);
	  sprintf(workrawforestmaxvcfinputname,"%s",inworkrawforestmaxvcfinputname);
          if (strcmp(inworkrawforestmaxvcfinputdb,"<workrawdir>") == 0) {
              sprintf(workrawforestmaxvcfinputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(workrawforestmaxvcfinputdb,"%s",inworkrawforestmaxvcfinputdb);
	  }
          workrawforestmaxvcfinputstartyear = atoi(inworkrawforestmaxvcfinputstartyear);
          workrawforestmaxvcfinputendyear = atoi(inworkrawforestmaxvcfinputendyear);
	  foundworkrawforestmaxvcfinput = 1;
      }
  }  
  
  return 0;

}


int readworkrawscaledvcfinputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *workrawscaledvcfinputvaluefile;
  int foundworkrawscaledvcfinput, itemsfound;
  char fullfilenamestr[1024], inworkrawscaledvcfinputname[1024], inworkrawscaledvcfinputdb[1024];
  char inworkrawscaledvcfinputstartyear[1024], inworkrawscaledvcfinputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  workrawscaledvcfinputvaluefile = fopen(fullfilenamestr,"r");
  
  foundworkrawscaledvcfinput = 0;

  while (foundworkrawscaledvcfinput == 0) {
      itemsfound = fscanf(workrawscaledvcfinputvaluefile,"%s %s %s %s",inworkrawscaledvcfinputname, inworkrawscaledvcfinputdb, inworkrawscaledvcfinputstartyear, inworkrawscaledvcfinputendyear);
      if (itemsfound != 4) {
          printf("Error: workrawscaledvcfinput %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inworkrawscaledvcfinputname) == 0) {
          printf("Processing workrawscaledvcfinput: %s\n",inworkrawscaledvcfinputname);
	  sprintf(workrawscaledvcfinputname,"%s",inworkrawscaledvcfinputname);
          if (strcmp(inworkrawscaledvcfinputdb,"<workrawdir>") == 0) {
              sprintf(workrawscaledvcfinputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(workrawscaledvcfinputdb,"%s",inworkrawscaledvcfinputdb);
	  }
          workrawscaledvcfinputstartyear = atoi(inworkrawscaledvcfinputstartyear);
          workrawscaledvcfinputendyear = atoi(inworkrawscaledvcfinputendyear);
	  foundworkrawscaledvcfinput = 1;
      }
  }  
  
  return 0;

}


int readworkrawextrappftinputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *workrawextrappftinputvaluefile;
  int foundworkrawextrappftinput, itemsfound;
  char fullfilenamestr[1024], inworkrawextrappftinputname[1024], inworkrawextrappftinputdb[1024];
  char inworkrawextrappftinputstartyear[1024], inworkrawextrappftinputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  workrawextrappftinputvaluefile = fopen(fullfilenamestr,"r");
  
  foundworkrawextrappftinput = 0;

  while (foundworkrawextrappftinput == 0) {
      itemsfound = fscanf(workrawextrappftinputvaluefile,"%s %s %s %s",inworkrawextrappftinputname, inworkrawextrappftinputdb, inworkrawextrappftinputstartyear, inworkrawextrappftinputendyear);
      if (itemsfound != 4) {
          printf("Error: workrawextrappftinput %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inworkrawextrappftinputname) == 0) {
          printf("Processing workrawextrappftinput: %s\n",inworkrawextrappftinputname);
	  sprintf(workrawextrappftinputname,"%s",inworkrawextrappftinputname);
          if (strcmp(inworkrawextrappftinputdb,"<workrawdir>") == 0) {
              sprintf(workrawextrappftinputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(workrawextrappftinputdb,"%s",inworkrawextrappftinputdb);
	  }
          workrawextrappftinputstartyear = atoi(inworkrawextrappftinputstartyear);
          workrawextrappftinputendyear = atoi(inworkrawextrappftinputendyear);
	  foundworkrawextrappftinput = 1;
      }
  }  
  
  return 0;

}


int readworkrawmergeinputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *workrawmergeinputvaluefile;
  int foundworkrawmergeinput, itemsfound;
  char fullfilenamestr[1024], inworkrawmergeinputname[1024], inworkrawmergeinputdb[1024];
  char inworkrawmergeinputstartyear[1024], inworkrawmergeinputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  workrawmergeinputvaluefile = fopen(fullfilenamestr,"r");
  
  foundworkrawmergeinput = 0;

  while (foundworkrawmergeinput == 0) {
      itemsfound = fscanf(workrawmergeinputvaluefile,"%s %s %s %s",inworkrawmergeinputname, inworkrawmergeinputdb, inworkrawmergeinputstartyear, inworkrawmergeinputendyear);
      if (itemsfound != 4) {
          printf("Error: workrawmergeinput %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inworkrawmergeinputname) == 0) {
          printf("Processing workrawmergeinput: %s\n",inworkrawmergeinputname);
	  sprintf(workrawmergeinputname,"%s",inworkrawmergeinputname);
          if (strcmp(inworkrawmergeinputdb,"<workrawdir>") == 0) {
              sprintf(workrawmergeinputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(workrawmergeinputdb,"%s",inworkrawmergeinputdb);
	  }
          workrawmergeinputstartyear = atoi(inworkrawmergeinputstartyear);
          workrawmergeinputendyear = atoi(inworkrawmergeinputendyear);
	  foundworkrawmergeinput = 1;
      }
  }  
  
  return 0;

}


int initializeGrids() {

  long clm5lin, clm5pix;
  long pftid, cftid;
  
  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
	  
          inCTSMTREEEXTRAPPCTNDLEVGTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          inCTSMTREEEXTRAPPCTNDLDECTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          inCTSMTREEEXTRAPPCTBRDEVGTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          inCTSMTREEEXTRAPPCTBRDDECTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          inCTSMHERBEXTRAPPCTSHRGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          inCTSMHERBEXTRAPPCTGRSC3Grid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          inCTSMHERBEXTRAPPCTGRSC4Grid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;

          outDESCMERGEDLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCMERGEDLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCMERGEDAREAGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCMERGEDPCTGLACIERGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCMERGEDPCTLAKEGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCMERGEDPCTWETLANDGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCMERGEDPCTURBANGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCMERGEDPCTNATVEGGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCMERGEDPCTCROPGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          for (pftid = 0; pftid < MAXPFT; pftid++) {
              outDESCMERGEDPCTPFTGrid[pftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
	  }

      }
  }

  return 0;
  	      
}


int clearGrids() {

  long clm5lin, clm5pix;
  long pftid, cftid;
  
  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
          for (pftid = 0; pftid < MAXPFT; pftid++) {
              outDESCMERGEDPCTPFTGrid[pftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
	  }
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

  return 0;
  
}

  
int readclm5referenceGrids(char *databasestr,char *seriesname,char *referenceyearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long clm5lin, clm5pix;

  sprintf(globalbinfilename,"%s/%s/%s.PCTURBAN.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inCTSMPCTURBANGrid,sizeof(inCTSMPCTURBANGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTNATVEG.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inCTSMPCTNATVEGGrid,sizeof(inCTSMPCTNATVEGGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTCROP.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inCTSMPCTCROPGrid,sizeof(inCTSMPCTCROPGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  return 0;
  
}


int readpotvegmaxforestvcfGrids(char *databasestr,char *seriesname,char *referenceyearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long clm5lin, clm5pix;

  sprintf(globalbinfilename,"%s/%s/%s.PCTTREE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPOTVEGMAXFORESTPCTTREEGrid,sizeof(inPOTVEGMAXFORESTPCTTREEGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTHERB.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPOTVEGMAXFORESTPCTHERBGrid,sizeof(inPOTVEGMAXFORESTPCTHERBGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTBARE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPOTVEGMAXFORESTPCTBAREGrid,sizeof(inPOTVEGMAXFORESTPCTBAREGrid),1,globalbinfile);  
  fclose(globalbinfile);

  return 0;
  
}


int readpotvegscaledforestvcfGrids(char *databasestr,char *seriesname,char *referenceyearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long clm5lin, clm5pix;

  sprintf(globalbinfilename,"%s/%s/%s.PCTTREE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPOTVEGSCALEDFORESTPCTTREEGrid,sizeof(inPOTVEGSCALEDFORESTPCTTREEGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTHERB.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPOTVEGSCALEDFORESTPCTHERBGrid,sizeof(inPOTVEGSCALEDFORESTPCTHERBGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTBARE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPOTVEGSCALEDFORESTPCTBAREGrid,sizeof(inPOTVEGSCALEDFORESTPCTBAREGrid),1,globalbinfile);  
  fclose(globalbinfile);

  return 0;
  
}


int readcurrworktreeextrapGrids(char *databasestr,char *seriesname,char *referenceyearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long clm5lin, clm5pix;

  for (pftid = 0; pftid < MAXPFT; pftid++) {
      sprintf(pftidstr,"%02d",pftid);
      sprintf(globalbinfilename,"%s/%s/%s.PCTNATPFT%s.%s.dat",databasestr,seriesname,seriesname,pftidstr,referenceyearname);
      printf("Reading: %s\n",globalbinfilename);
      globalbinfile = fopen(globalbinfilename,"r");
      fread(tempGrid,sizeof(tempGrid),1,globalbinfile);  
      fclose(globalbinfile);
      for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
          for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
	      inTREEEXTRAPPCTPFTGrid[pftid][clm5lin * MAXCLMPIX + clm5pix] = tempGrid[clm5lin * MAXCLMPIX + clm5pix];
	      if (pftid == 1 || pftid == 2) {
	          inCTSMTREEEXTRAPPCTNDLEVGTREEGrid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
              }
	      if (pftid == 3) {
	          inCTSMTREEEXTRAPPCTNDLDECTREEGrid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
              }
	      if (pftid == 4 || pftid == 5) {
	          inCTSMTREEEXTRAPPCTBRDEVGTREEGrid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
              }
	      if (pftid == 6 || pftid == 7 || pftid == 8) {
	          inCTSMTREEEXTRAPPCTBRDDECTREEGrid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
              }
	  }
      }
  }
  
  return 0;
  
}


int readcurrworkherbextrapGrids(char *databasestr,char *seriesname,char *referenceyearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long clm5lin, clm5pix;

  for (pftid = 0; pftid < MAXPFT; pftid++) {
      sprintf(pftidstr,"%02d",pftid);
      sprintf(globalbinfilename,"%s/%s/%s.PCTNATPFT%s.%s.dat",databasestr,seriesname,seriesname,pftidstr,referenceyearname);
      printf("Reading: %s\n",globalbinfilename);
      globalbinfile = fopen(globalbinfilename,"r");
      fread(tempGrid,sizeof(tempGrid),1,globalbinfile);  
      fclose(globalbinfile);
      for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
          for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
	      inHERBEXTRAPPCTPFTGrid[pftid][clm5lin * MAXCLMPIX + clm5pix] = tempGrid[clm5lin * MAXCLMPIX + clm5pix];
	      if (pftid == 9 || pftid == 10 || pftid == 11) {
	          inCTSMHERBEXTRAPPCTSHRGrid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
              }
	      if (pftid == 12 || pftid == 13) {
	          inCTSMHERBEXTRAPPCTGRSC3Grid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
              }
	      if (pftid == 14) {
	          inCTSMHERBEXTRAPPCTGRSC4Grid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
              }
	  }
      }
  }
  
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


int readmworkergedSECDFGrids(char *databasestr,char *seriesname,char *referenceyearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long clm5lin, clm5pix;

  sprintf(globalbinfilename,"%s/%s/%s.SECDFFRACTREE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inMERGEDSECDFFRACTREEGrid,sizeof(inMERGEDSECDFFRACTREEGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.SECDFFRACHERB.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inMERGEDSECDFFRACHERBGrid,sizeof(inMERGEDSECDFFRACHERBGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.SECDFFRACBARE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inMERGEDSECDFFRACBAREGrid,sizeof(inMERGEDSECDFFRACBAREGrid),1,globalbinfile);  
  fclose(globalbinfile);

  
  return 0;
  
}


int genLUHLandGrids() {

  long clm5lin, clm5pix;
  long pftid, cftid;
  float landmask, landfrac, area, glacierpct, lakepct, wetlandpct, natvegpct;
  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
          landmask = inLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix];
          landfrac = inLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix];
	  if (landfrac > 0.0) {
              outDESCMERGEDLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix] = landmask;
              outDESCMERGEDLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix] = landfrac;
	      outDESCMERGEDAREAGrid[clm5lin * MAXCLMPIX + clm5pix] = inAREAGrid[clm5lin * MAXCLMPIX + clm5pix];
	      glacierpct = inCTSMPCTGLACIERGrid[clm5lin * MAXCLMPIX + clm5pix];
	      lakepct = inCTSMPCTLAKEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      wetlandpct = inCTSMPCTWETLANDGrid[clm5lin * MAXCLMPIX + clm5pix];
	      natvegpct = 100.0 - glacierpct - lakepct - wetlandpct;
              outDESCMERGEDPCTGLACIERGrid[clm5lin * MAXCLMPIX + clm5pix] = glacierpct;
              outDESCMERGEDPCTLAKEGrid[clm5lin * MAXCLMPIX + clm5pix] = lakepct;
              outDESCMERGEDPCTWETLANDGrid[clm5lin * MAXCLMPIX + clm5pix] = wetlandpct;
              outDESCMERGEDPCTNATVEGGrid[clm5lin * MAXCLMPIX + clm5pix] = natvegpct;
	  }
      }
  }

  return 0;
  	      
}


int genMergedDescriptorVegGrids() {

  long clm5lin, clm5pix;
  long pftid, cftid;
  float landmask, landfrac, area, precipann, tempaverage;
  float tempcoldestvalue, tempwarmestvalue, growdegdaysvalue, precipannvalue, precipminvalue, precipwinvalue, precipmaxtg22value;
  float foresttreepct, forestherbpct, forestbarepct, forestallpct;
  float potvegmaxtreepct, potvegscaledtreepct, mergedsecdffractree;
  float potvegmaxherbpct, potvegscaledherbpct, mergedsecdffracherb;
  float potvegmaxbarepct, potvegscaledbarepct, mergedsecdffracbare;
  float addpct, removepct;
  float forestndlevgtreepct, forestndldectreepct, forestbrdevgtreepct, forestbrddectreepct, forestshrpct, forestgrsc3pct, forestgrsc4pct;
  float forestalltreepct, forestallherbpct;
  float forestndlevgtemptreepct, forestndlevgborltreepct, forestndldecborltreepct;
  float forestbrdevgtroptreepct, forestbrddectroptreepct, forestbrdevgtemptreepct, forestbrddectemptreepct, forestbrddecborltreepct;
  float forestshrevgtemppct, forestshrdectemppct, forestshrdecborlpct, forestgrsc3arcpct, forestgrsc3nonpct, forestgrsc4nonpct;

  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
          landmask = inLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix];
          landfrac = inLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix];
	  if (landfrac > 0.0) {

              tempcoldestvalue = inTEMPCOLDGrid[clm5lin * MAXCLMPIX + clm5pix];
              tempwarmestvalue = inTEMPWARMGrid[clm5lin * MAXCLMPIX + clm5pix];
              growdegdaysvalue = inTEMPGDDGrid[clm5lin * MAXCLMPIX + clm5pix];
              precipannvalue = inPRECIPANNGrid[clm5lin * MAXCLMPIX + clm5pix];
              precipminvalue = inPRECIPMINGrid[clm5lin * MAXCLMPIX + clm5pix];
              precipwinvalue = inPRECIPWINGrid[clm5lin * MAXCLMPIX + clm5pix];
              precipmaxtg22value = inPRECIPMAXTG22Grid[clm5lin * MAXCLMPIX + clm5pix];

	      potvegmaxtreepct = inPOTVEGMAXFORESTPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      potvegscaledtreepct = inPOTVEGSCALEDFORESTPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      mergedsecdffractree = inMERGEDSECDFFRACTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      foresttreepct = mergedsecdffractree * potvegmaxtreepct + (1.0 - mergedsecdffractree) * potvegscaledtreepct;
	      
	      potvegmaxherbpct = inPOTVEGMAXFORESTPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
	      potvegscaledherbpct = inPOTVEGSCALEDFORESTPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
	      mergedsecdffracherb = inMERGEDSECDFFRACHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
	      forestherbpct = mergedsecdffracherb * potvegmaxherbpct + (1.0 - mergedsecdffracherb) * potvegscaledherbpct;
	      
	      potvegmaxbarepct = inPOTVEGMAXFORESTPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix];
	      potvegscaledbarepct = inPOTVEGSCALEDFORESTPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix];
	      mergedsecdffracbare = inMERGEDSECDFFRACBAREGrid[clm5lin * MAXCLMPIX + clm5pix];
	      forestbarepct = mergedsecdffracbare * potvegmaxbarepct + (1.0 - mergedsecdffracbare) * potvegscaledbarepct;
	      
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

	      forestndlevgtreepct = inCTSMTREEEXTRAPPCTNDLEVGTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      forestndldectreepct = inCTSMTREEEXTRAPPCTNDLDECTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      forestbrdevgtreepct = inCTSMTREEEXTRAPPCTBRDEVGTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      forestbrddectreepct = inCTSMTREEEXTRAPPCTBRDDECTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      forestalltreepct = forestndlevgtreepct + forestndldectreepct + forestbrdevgtreepct + forestbrddectreepct;
	      
	      if (forestalltreepct > 0.0) {
	          forestndlevgtreepct = forestndlevgtreepct / forestalltreepct * foresttreepct;
                  forestndlevgtreepct = trunc(forestndlevgtreepct * CLMPCTTruncate) / CLMPCTTruncate;
	          forestndldectreepct = forestndldectreepct / forestalltreepct * foresttreepct;
                  forestndldectreepct = trunc(forestndldectreepct * CLMPCTTruncate) / CLMPCTTruncate;
	          forestbrdevgtreepct = forestbrdevgtreepct / forestalltreepct * foresttreepct;
                  forestbrdevgtreepct = trunc(forestbrdevgtreepct * CLMPCTTruncate) / CLMPCTTruncate;
	          forestbrddectreepct = forestbrddectreepct / forestalltreepct * foresttreepct;
                  forestbrddectreepct = trunc(forestbrddectreepct * CLMPCTTruncate) / CLMPCTTruncate;
              }
	      
              if (tempcoldestvalue > -19.0 && growdegdaysvalue > 1200.0) {
                  forestndlevgtemptreepct = forestndlevgtreepct + forestndldectreepct;
                  forestndlevgborltreepct = 0.0;
                  forestndldecborltreepct = 0.0;
              }
	      else {
	          forestndlevgtemptreepct = 0.0;
		  forestndlevgborltreepct = forestndlevgtreepct;
		  forestndldecborltreepct = forestndldectreepct;
              }

              if ((tempcoldestvalue < 10.0 && forestbrdevgtreepct < 1.0 && clm5lin < MAXCLMLIN / 2) || (tempcoldestvalue < 0.0 && forestbrdevgtreepct < 1.0)) {                   forestbrddectreepct = forestbrddectreepct + forestbrdevgtreepct;
		  forestbrdevgtreepct = 0.0;
              } 
	      
              if (tempcoldestvalue > 15.0) {
	          forestbrdevgtroptreepct = forestbrdevgtreepct;
		  forestbrdevgtemptreepct = 0.0;
              }
	      else {
		  forestbrdevgtroptreepct = 0.0;
	          forestbrdevgtemptreepct = forestbrdevgtreepct;
	      }
	      
              if (tempcoldestvalue > 15.0) {
                  forestbrddectroptreepct = forestbrddectreepct;
                  forestbrddectemptreepct = 0.0;
                  forestbrddecborltreepct = 0.0;
	      }
	      else {
                  if (tempcoldestvalue > -15.0 && growdegdaysvalue > 1200.0) {
                      forestbrddectroptreepct = 0.0;
                      forestbrddectemptreepct = forestbrddectreepct;
                      forestbrddecborltreepct = 0.0;
		  }
		  else {
                      forestbrddectroptreepct = 0.0;
                      forestbrddectemptreepct = 0.0;
                      forestbrddecborltreepct = forestbrddectreepct;
		  }
              }
	      
	      forestshrpct = inCTSMHERBEXTRAPPCTSHRGrid[clm5lin * MAXCLMPIX + clm5pix];
	      forestgrsc3pct = inCTSMHERBEXTRAPPCTGRSC3Grid[clm5lin * MAXCLMPIX + clm5pix];
	      forestgrsc4pct = inCTSMHERBEXTRAPPCTGRSC4Grid[clm5lin * MAXCLMPIX + clm5pix];
              forestallherbpct = forestshrpct + forestgrsc3pct + forestgrsc4pct;
	      
	      if (forestallherbpct > 0.0) {
	          forestshrpct = forestshrpct / forestallherbpct * forestherbpct;
                  forestshrpct = trunc(forestshrpct * CLMPCTTruncate) / CLMPCTTruncate;
	          forestgrsc3pct = forestgrsc3pct / forestallherbpct * forestherbpct;
                  forestgrsc3pct = trunc(forestgrsc3pct * CLMPCTTruncate) / CLMPCTTruncate;
	          forestgrsc4pct = forestgrsc4pct / forestallherbpct * forestherbpct;
                  forestgrsc4pct = trunc(forestgrsc4pct * CLMPCTTruncate) / CLMPCTTruncate;
              }
	      
              if ((tempcoldestvalue > -19.0 && growdegdaysvalue > 1200.0) && (precipannvalue > 520.0 && precipwinvalue > (precipannvalue * 2.0 / 3.0))) {
	          forestshrevgtemppct = forestshrpct;
		  forestshrdectemppct = 0.0;
		  forestshrdecborlpct = 0.0;
	      }
	      else {
		  if (tempcoldestvalue > -19.0 && growdegdaysvalue > 1200.0) {
                      forestshrevgtemppct = 0.0;
                      forestshrdectemppct = forestshrpct;
                      forestshrdecborlpct = 0.0;
                  }
                  else {
                      forestshrevgtemppct = 0.0;
                      forestshrdectemppct = 0.0;
                      forestshrdecborlpct = forestshrpct;
                  }
	      }
	      
              if (growdegdaysvalue <= 1000.0) {
                  forestgrsc3arcpct = forestgrsc3pct + forestgrsc4pct;
                  forestgrsc3nonpct = 0.0;
                  forestgrsc4nonpct = 0.0;
              }
              else {
                  if (tempcoldestvalue > 22.0 && precipminvalue > 25.0) {
                      forestgrsc3arcpct = 0.0;
                      forestgrsc3nonpct = 0.0;
                      forestgrsc4nonpct = forestgrsc3pct + forestgrsc4pct;
		  }
		  else {
                      forestgrsc3arcpct = 0.0;
                      forestgrsc3nonpct = forestgrsc3pct;
                      forestgrsc4nonpct = forestgrsc4pct;
		  }
              }
	      
	      outDESCMERGEDPCTPFTGrid[0][clm5lin * MAXCLMPIX + clm5pix] = forestbarepct;
	      outDESCMERGEDPCTPFTGrid[1][clm5lin * MAXCLMPIX + clm5pix] = forestndlevgtemptreepct;		  
	      outDESCMERGEDPCTPFTGrid[2][clm5lin * MAXCLMPIX + clm5pix] = forestndlevgborltreepct;		  
	      outDESCMERGEDPCTPFTGrid[3][clm5lin * MAXCLMPIX + clm5pix] = forestndldecborltreepct;		  
	      outDESCMERGEDPCTPFTGrid[4][clm5lin * MAXCLMPIX + clm5pix] = forestbrdevgtroptreepct;		  
	      outDESCMERGEDPCTPFTGrid[5][clm5lin * MAXCLMPIX + clm5pix] = forestbrdevgtemptreepct;		  
	      outDESCMERGEDPCTPFTGrid[6][clm5lin * MAXCLMPIX + clm5pix] = forestbrddectroptreepct;		  
	      outDESCMERGEDPCTPFTGrid[7][clm5lin * MAXCLMPIX + clm5pix] = forestbrddectemptreepct;		  
	      outDESCMERGEDPCTPFTGrid[8][clm5lin * MAXCLMPIX + clm5pix] = forestbrddecborltreepct;		  
	      outDESCMERGEDPCTPFTGrid[9][clm5lin * MAXCLMPIX + clm5pix] = forestshrevgtemppct;
	      outDESCMERGEDPCTPFTGrid[10][clm5lin * MAXCLMPIX + clm5pix] = forestshrdectemppct;
	      outDESCMERGEDPCTPFTGrid[11][clm5lin * MAXCLMPIX + clm5pix] = forestshrdecborlpct;
	      outDESCMERGEDPCTPFTGrid[12][clm5lin * MAXCLMPIX + clm5pix] = forestgrsc3arcpct;
	      outDESCMERGEDPCTPFTGrid[13][clm5lin * MAXCLMPIX + clm5pix] = forestgrsc3nonpct;
	      outDESCMERGEDPCTPFTGrid[14][clm5lin * MAXCLMPIX + clm5pix] = forestgrsc4nonpct; 
	  }
      }
  }

  return 0;
  	      
}


int balanceMergedDescriptorVegGrids() {

  long clm5lin, clm5pix;
  long pftid, forestmaxpftid;
  float landmask, landfrac, area, precipann, tempaverage;
  float forestpftpct, forestmaxpftpct, forestallpftpct;
  float forestremovepct, forestaddpct;

  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
          landmask = inLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix];
          landfrac = inLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix];
	  if (landfrac > 0.0) {
	  
	      forestmaxpftid = -1;
	      forestmaxpftpct = 0.0;
	      forestallpftpct = 0.0;
	      
	      for (pftid = 0; pftid < MAXPFT; pftid++) {
	          forestpftpct = outDESCMERGEDPCTPFTGrid[pftid][clm5lin * MAXCLMPIX + clm5pix];
		  forestallpftpct += forestpftpct;
		  if (forestpftpct > forestmaxpftpct) {
		      forestmaxpftpct = forestpftpct;
		      forestmaxpftid = pftid;
		  }
              }
	      
	      if (forestallpftpct > 100.0) {
	          forestremovepct = forestallpftpct - 100.0;
	          outDESCMERGEDPCTPFTGrid[forestmaxpftid][clm5lin * MAXCLMPIX + clm5pix] -= forestremovepct;
              }
	      
	      if (forestallpftpct < 100.0) {
	          forestaddpct = 100.0 - forestallpftpct;
	          outDESCMERGEDPCTPFTGrid[forestmaxpftid][clm5lin * MAXCLMPIX + clm5pix] += forestaddpct;
              }
	      
	      
	  }
      }
  }

  return 0;
  	      
}


int writedescmergedgrids(char *databasestr, char *seriesname, char *yearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  long outlin, outpix;
  long clm5lin, clm5pix;
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];

  for (pftid = 0; pftid < MAXPFT; pftid++) {

      for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
          for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
	      tempGrid[clm5lin * MAXCLMPIX + clm5pix] = outDESCMERGEDPCTPFTGrid[pftid][clm5lin * MAXCLMPIX + clm5pix];
	  }
      }
      
      sprintf(pftidstr,"%02d",pftid);
      sprintf(globalbinfilename,"%s/%s/%s.PCTNATPFT%s.%s.dat",databasestr,seriesname,seriesname,pftidstr,yearname);
      printf("Writing: %s\n",globalbinfilename);
      globalbinfile = fopen(globalbinfilename,"w+");
      fwrite(tempGrid,sizeof(tempGrid),1,globalbinfile);  
      fclose(globalbinfile);
  }

  return 0;

}


int main(long narg, char **argv) {

  char workrawclm5inputdirname[1024];
  char workrawforestmaxvcfinputdirname[1024];
  char workrawscaledvcfinputdirname[1024];
  char workrawextrappftinputdirname[1024];
  char workrawmergeinputdirname[1024];
  char luh3descrawpftoutputdirname[1024];
  char clm5modisrawinputreferenceyearstr[256];
  char workrawclm5inputreferenceyearstr[256];
  char workrawforestmaxvcfinputreferenceyearstr[256];
  char workrawscaledvcfinputreferenceyearstr[256];
  char workrawextrappftinputreferenceyearstr[256];
  char workrawmergeinputreferenceyearstr[256];
  char luh3descrawpftoutputreferenceyearstr[256];
  int currentyear;
    
  if(narg != 8){
        printf("Usage createalldescsecdfCTSM53Deg025bin luh3descrawfile luh3descrawdir timeseriesrawfile timeseriesrawdir workrawfile workrawdir luh3descsecdfpftnamelistfile\n");
        return 0;
  }
  
  readnamelistfile(argv[7]);

  readluh3descrawpftoutputlutfile(argv[1],argv[2],namelistluh3descrawpftoutput);
  readclm5modisrawinputlutfile(argv[3],argv[4],namelistclm5modisrawinput);
  readworkrawclm5inputlutfile(argv[5],argv[6],namelistworkrawclm5input);
  readworkrawforestmaxvcfinputlutfile(argv[5],argv[6],namelistworkrawforestmaxvcfinput);
  readworkrawscaledvcfinputlutfile(argv[5],argv[6],namelistworkrawscaledvcfinput);
  readworkrawextrappftinputlutfile(argv[5],argv[6],namelistworkrawextrappftinput);
  readworkrawmergeinputlutfile(argv[5],argv[6],namelistworkrawmergeinput);
  clm5modisrawinputreferenceyear = atoi(namelistclm5modisrawinputreferenceyear);
  sprintf(clm5modisrawinputreferenceyearstr,"%04d",clm5modisrawinputreferenceyear);
  workrawclm5inputreferenceyear = atoi(namelistworkrawclm5inputreferenceyear);
  sprintf(workrawclm5inputreferenceyearstr,"%04d",workrawclm5inputreferenceyear);
  workrawforestmaxvcfinputreferenceyear = atoi(namelistworkrawforestmaxvcfinputreferenceyear);
  sprintf(workrawforestmaxvcfinputreferenceyearstr,"%04d",workrawforestmaxvcfinputreferenceyear);
  workrawscaledvcfinputreferenceyear = atoi(namelistworkrawscaledvcfinputreferenceyear);
  sprintf(workrawscaledvcfinputreferenceyearstr,"%04d",workrawscaledvcfinputreferenceyear);
  workrawextrappftinputreferenceyear = atoi(namelistworkrawextrappftinputreferenceyear);
  sprintf(workrawextrappftinputreferenceyearstr,"%04d",workrawextrappftinputreferenceyear);  
  workrawmergeinputstartyear = atoi(namelistworkrawmergeinputstartyear);
  workrawmergeinputendyear = atoi(namelistworkrawmergeinputendyear);
  
  initializeGrids();

  sprintf(workrawclm5inputdirname,"%s/%s",workrawclm5inputdb,workrawclm5inputname);
  sprintf(workrawforestmaxvcfinputdirname,"%s/%s",workrawforestmaxvcfinputdb,workrawforestmaxvcfinputname);
  sprintf(workrawscaledvcfinputdirname,"%s/%s",workrawscaledvcfinputdb,workrawscaledvcfinputname);
  sprintf(workrawextrappftinputdirname,"%s/%s",workrawextrappftinputdb,workrawextrappftinputname);
  readlandGrids(workrawclm5inputdirname,namelistworkrawclm5inputname,workrawclm5inputreferenceyearstr);
  readclm5referenceGrids(workrawclm5inputdirname,namelistworkrawclm5inputname,workrawclm5inputreferenceyearstr);    
  readpotvegmaxforestvcfGrids(workrawforestmaxvcfinputdirname,namelistworkrawforestmaxvcfinputforestname,workrawforestmaxvcfinputreferenceyearstr);    
  readpotvegscaledforestvcfGrids(workrawscaledvcfinputdirname,namelistworkrawscaledvcfinputforestname,workrawscaledvcfinputreferenceyearstr);    
  readcurrworktreeextrapGrids(workrawextrappftinputdirname,namelistworkrawextrappftinputtreename,workrawextrappftinputreferenceyearstr);
  readcurrworkherbextrapGrids(workrawextrappftinputdirname,namelistworkrawextrappftinputherbname,workrawextrappftinputreferenceyearstr);
  readclimvarseries(clm5modisrawinputdb,clm5modisrawinputname,"CLIM");
  readlaivarseries(clm5modisrawinputdb,clm5modisrawinputname,"FILLCLIM");
  
  genLUHLandGrids();
  
  sprintf(workrawmergeinputdirname,"%s/%s",workrawmergeinputdb,workrawmergeinputname);
  sprintf(luh3descrawpftoutputdirname,"%s/%s",luh3descrawpftoutputdb,luh3descrawpftoutputname);
  for (currentyear = workrawmergeinputstartyear; currentyear <= workrawmergeinputendyear; currentyear++) {
      sprintf(workrawmergeinputreferenceyearstr,"%04d",currentyear);  
      sprintf(luh3descrawpftoutputreferenceyearstr,"%04d",currentyear);  
      clearGrids();
      readmworkergedSECDFGrids(workrawmergeinputdirname,namelistworkrawmergeinputname,workrawmergeinputreferenceyearstr);
      genMergedDescriptorVegGrids();
      balanceMergedDescriptorVegGrids();
      writedescmergedgrids(luh3descrawpftoutputdirname,namelistluh3descrawpftoutputsecdfname,luh3descrawpftoutputreferenceyearstr);
  }
  

  return 1; 
  
}
