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
char namelistpotvegallscaledvcfname[1024];
char namelistworkrawvcfoutput[1024];
char namelistworkrawvcfoutputforestname[1024];
char namelistworkrawvcfoutputnonforestname[1024];
char namelistworkrawvcfoutputrangename[1024];
char namelistluh3rawinput[1024];
char namelistluh3rawinputreferenceyear[1024];
char namelistworkrawclm5input[1024];
char namelistworkrawclm5inputname[1024];
char namelistworkrawclm5inputvcfname[1024];
char namelistworkrawclm5inputreferenceyear[1024];
char namelistworkrawcombvcfinput[1024];
char namelistworkrawcombvcfinputforestname[1024];
char namelistworkrawcombvcfinputnonforestname[1024];
char namelistworkrawcombvcfinputrangename[1024];
char namelistworkrawcombvcfinputreferenceyear[1024];
char namelistclm5modisrawinput[1024];
char namelistclm5modisrawinputreferenceyear[1024];

int luh3rawinputreferenceyear;
int clm5modisrawinputreferenceyear;
int workrawclm5inputreferenceyear;
int workrawcombvcfinputreferenceyear;

char luh3rawinputname[1024];
char luh3rawinputdb[1024];
int luh3rawinputstartyear;
int luh3rawinputendyear;

char clm5modisrawinputname[1024];
char clm5modisrawinputdb[1024];
int clm5modisrawinputstartyear;
int clm5modisrawinputendyear;

char workrawclm5inputname[1024];
char workrawclm5inputdb[1024];
int workrawclm5inputstartyear;
int workrawclm5inputendyear;

char workrawcombvcfinputname[1024];
char workrawcombvcfinputdb[1024];
int workrawcombvcfinputstartyear;
int workrawcombvcfinputendyear;

char workrawscaledvcfinputname[1024];
char workrawscaledvcfinputdb[1024];
int workrawscaledvcfinputstartyear;
int workrawscaledvcfinputendyear;

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

float inFORESTPCTTREEGrid[MAXCLMPIX * MAXCLMLIN];
float inFORESTPCTHERBGrid[MAXCLMPIX * MAXCLMLIN];
float inFORESTPCTBAREGrid[MAXCLMPIX * MAXCLMLIN];

float inNONFORESTPCTTREEGrid[MAXCLMPIX * MAXCLMLIN];
float inNONFORESTPCTHERBGrid[MAXCLMPIX * MAXCLMLIN];
float inNONFORESTPCTBAREGrid[MAXCLMPIX * MAXCLMLIN];

float inRANGEPCTTREEGrid[MAXCLMPIX * MAXCLMLIN];
float inRANGEPCTHERBGrid[MAXCLMPIX * MAXCLMLIN];
float inRANGEPCTBAREGrid[MAXCLMPIX * MAXCLMLIN];

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
  fscanf(namelistfile,"%s %s",fieldname,namelistpotvegallscaledvcfname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfoutput);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfoutputforestname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfoutputnonforestname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfoutputrangename);
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3rawinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3rawinputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5input);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5inputname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5inputvcfname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5inputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawcombvcfinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawcombvcfinputforestname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawcombvcfinputnonforestname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawcombvcfinputrangename);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawcombvcfinputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5modisrawinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5modisrawinputreferenceyear);
  
  return 0;

}


int readluh3rawinputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *luh3rawinputvaluefile;
  int foundluh3rawinput, itemsfound;
  char fullfilenamestr[1024], inluh3rawinputname[1024], inluh3rawinputdb[1024];
  char inluh3rawinputstartyear[1024], inluh3rawinputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  luh3rawinputvaluefile = fopen(fullfilenamestr,"r");
  
  foundluh3rawinput = 0;

  while (foundluh3rawinput == 0) {
      itemsfound = fscanf(luh3rawinputvaluefile,"%s %s %s %s",inluh3rawinputname, inluh3rawinputdb, inluh3rawinputstartyear, inluh3rawinputendyear);
      if (itemsfound != 4) {
          printf("Error: luh3rawinput %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inluh3rawinputname) == 0) {
          printf("Processing luh3rawinput: %s\n",inluh3rawinputname);
	  sprintf(luh3rawinputname,"%s",inluh3rawinputname);
          if (strcmp(inluh3rawinputdb,"<luh3rawdir>") == 0) {
              sprintf(luh3rawinputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(luh3rawinputdb,"%s",inluh3rawinputdb);
	  }
          luh3rawinputstartyear = atoi(inluh3rawinputstartyear);
          luh3rawinputendyear = atoi(inluh3rawinputendyear);
	  foundluh3rawinput = 1;
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


int readworkrawcombvcfinputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *workrawcombvcfinputvaluefile;
  int foundworkrawcombvcfinput, itemsfound;
  char fullfilenamestr[1024], inworkrawcombvcfinputname[1024], inworkrawcombvcfinputdb[1024];
  char inworkrawcombvcfinputstartyear[1024], inworkrawcombvcfinputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  workrawcombvcfinputvaluefile = fopen(fullfilenamestr,"r");
  
  foundworkrawcombvcfinput = 0;

  while (foundworkrawcombvcfinput == 0) {
      itemsfound = fscanf(workrawcombvcfinputvaluefile,"%s %s %s %s",inworkrawcombvcfinputname, inworkrawcombvcfinputdb, inworkrawcombvcfinputstartyear, inworkrawcombvcfinputendyear);
      if (itemsfound != 4) {
          printf("Error: workrawcombvcfinput %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inworkrawcombvcfinputname) == 0) {
          printf("Processing workrawcombvcfinput: %s\n",inworkrawcombvcfinputname);
	  sprintf(workrawcombvcfinputname,"%s",inworkrawcombvcfinputname);
          if (strcmp(inworkrawcombvcfinputdb,"<workrawdir>") == 0) {
              sprintf(workrawcombvcfinputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(workrawcombvcfinputdb,"%s",inworkrawcombvcfinputdb);
	  }
          workrawcombvcfinputstartyear = atoi(inworkrawcombvcfinputstartyear);
          workrawcombvcfinputendyear = atoi(inworkrawcombvcfinputendyear);
	  foundworkrawcombvcfinput = 1;
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

  long clm5lin, clm5pix;
  long pftid, cftid;
  
  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
          outCTSMLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outCTSMLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outCTSMAREAGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outCTSMPCTGLACIERGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outCTSMPCTLAKEGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outCTSMPCTWETLANDGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outCTSMPCTURBANGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outCTSMPCTNATVEGGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outCTSMPCTCROPGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;

          outFORESTPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outFORESTPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outFORESTPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;

          outNONFORESTPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outNONFORESTPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outNONFORESTPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;

          outRANGEPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outRANGEPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outRANGEPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
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
  long clm5lin, clm5pix;

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
  long clm5lin, clm5pix;

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


int readcombinedforestvcfGrids(char *databasestr,char *seriesname,char *referenceyearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long clm5lin, clm5pix;

  sprintf(globalbinfilename,"%s/%s/%s.PCTTREE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inFORESTPCTTREEGrid,sizeof(inFORESTPCTTREEGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTHERB.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inFORESTPCTHERBGrid,sizeof(inFORESTPCTHERBGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTBARE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inFORESTPCTBAREGrid,sizeof(inFORESTPCTBAREGrid),1,globalbinfile);  
  fclose(globalbinfile);

  return 0;
  
}


int readcombinednonforestvcfGrids(char *databasestr,char *seriesname,char *referenceyearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long clm5lin, clm5pix;

  sprintf(globalbinfilename,"%s/%s/%s.PCTTREE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inNONFORESTPCTTREEGrid,sizeof(inNONFORESTPCTTREEGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTHERB.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inNONFORESTPCTHERBGrid,sizeof(inNONFORESTPCTHERBGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTBARE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inNONFORESTPCTBAREGrid,sizeof(inNONFORESTPCTBAREGrid),1,globalbinfile);  
  fclose(globalbinfile);

  return 0;
  
}


int readcombinedrangevcfGrids(char *databasestr,char *seriesname,char *referenceyearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long clm5lin, clm5pix;

  sprintf(globalbinfilename,"%s/%s/%s.PCTTREE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inRANGEPCTTREEGrid,sizeof(inRANGEPCTTREEGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTHERB.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inRANGEPCTHERBGrid,sizeof(inRANGEPCTHERBGrid),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.PCTBARE.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inRANGEPCTBAREGrid,sizeof(inRANGEPCTBAREGrid),1,globalbinfile);  
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

  long clm5lin, clm5pix;
  long pftid, cftid;
  float landmask, landfrac, area, glacierpct, lakepct, wetlandpct, natvegpct;
  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
          landmask = inLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix];
          landfrac = inLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix];
	  if (landfrac > 0.0) {
              outCTSMLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix] = landmask;
              outCTSMLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix] = landfrac;
	      outCTSMAREAGrid[clm5lin * MAXCLMPIX + clm5pix] = inAREAGrid[clm5lin * MAXCLMPIX + clm5pix];
	      glacierpct = inCTSMPCTGLACIERGrid[clm5lin * MAXCLMPIX + clm5pix];
	      lakepct = inCTSMPCTLAKEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      wetlandpct = inCTSMPCTWETLANDGrid[clm5lin * MAXCLMPIX + clm5pix];
	      natvegpct = 100.0 - glacierpct - lakepct - wetlandpct;
              outCTSMPCTGLACIERGrid[clm5lin * MAXCLMPIX + clm5pix] = glacierpct;
              outCTSMPCTLAKEGrid[clm5lin * MAXCLMPIX + clm5pix] = lakepct;
              outCTSMPCTWETLANDGrid[clm5lin * MAXCLMPIX + clm5pix] = wetlandpct;
              outCTSMPCTNATVEGGrid[clm5lin * MAXCLMPIX + clm5pix] = natvegpct;
	  }
      }
  }

  return 0;
  	      
}


int genLUHScaledVegGrids() {

  long clm5lin, clm5pix;
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
  float removetreepct, removeluh3treepct, addtreepct, addluh3treepct;
  float removeherbpct, removeluh3herbpct, addherbpct, addluh3herbpct;
  float removebarepct, removeluh3barepct, addbarepct, addluh3barepct;
  int precipannindex, tempaverageindex;

  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
          landmask = inLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix];
          landfrac = inLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix];
	  if (landfrac > 0.0) {
	  
	      luh3forest = inLUH3PRIMFGrid[clm5lin * MAXCLMPIX + clm5pix] + inLUH3SECDFGrid[clm5lin * MAXCLMPIX + clm5pix];
	      luh3nonforest = inLUH3PRIMNGrid[clm5lin * MAXCLMPIX + clm5pix] + inLUH3SECDNGrid[clm5lin * MAXCLMPIX + clm5pix];
	      luh3range = inLUH3RANGEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      luh3pasture = inLUH3PASTRGrid[clm5lin * MAXCLMPIX + clm5pix];
	      luh3total = luh3forest + luh3nonforest + luh3range;

              if (luh3total > 0.0) {
	          luh3forest = luh3forest / luh3total;
		  luh3nonforest = luh3nonforest / luh3total;
		  luh3range = luh3range / luh3total;
              }
	      	  
	      luh3forest = trunc(luh3forest * LUH3StateTruncate) / LUH3StateTruncate;
	      luh3nonforest = trunc(luh3nonforest * LUH3StateTruncate) / LUH3StateTruncate;
	      luh3range = trunc(luh3range * LUH3StateTruncate) / LUH3StateTruncate;
	      luh3pasture = trunc(luh3pasture * LUH3StateTruncate) / LUH3StateTruncate;
	      luh3total = luh3forest + luh3nonforest + luh3range;
	      
	      currentnatvegfrac = inCTSMPCTNATVEGGrid[clm5lin * MAXCLMPIX + clm5pix] / 100.0;
	      currenttreepct = inCTSMPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
              currenttreepct = trunc(currenttreepct * CLMPCTTruncate) / CLMPCTTruncate;
	      currentherbpct = inCTSMPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
              currentherbpct = trunc(currentherbpct * CLMPCTTruncate) / CLMPCTTruncate;
	      currentbarepct = inCTSMPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix];
              currentbarepct = trunc(currentbarepct * CLMPCTTruncate) / CLMPCTTruncate;
	      
	      foresttreepct = inFORESTPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      forestherbpct = inFORESTPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
	      forestbarepct = inFORESTPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix];
	      
	      nonforesttreepct = inNONFORESTPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      nonforestherbpct = inNONFORESTPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
	      nonforestbarepct = inNONFORESTPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix];
	      
	      rangetreepct = inRANGEPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      rangeherbpct = inRANGEPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
	      rangebarepct = inRANGEPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix];
	      
	      scaledforesttreepct = foresttreepct;
	      scaledforestherbpct = forestherbpct;
	      scaledforestbarepct = forestbarepct;
	      
	      scalednonforesttreepct = nonforesttreepct;
	      scalednonforestherbpct = nonforestherbpct;
	      scalednonforestbarepct = nonforestbarepct;
	      
	      scaledrangetreepct = rangetreepct;
	      scaledrangeherbpct = rangeherbpct;
	      scaledrangebarepct = rangebarepct;
	      
	      if (luh3total > LUH3StateThreshold) {
	      
	          luh3treepct = luh3forest * scaledforesttreepct + luh3nonforest * scalednonforesttreepct + luh3range * scaledrangetreepct;
                  if (luh3treepct > currenttreepct) { /* remove tree pct values */
                      removetreepct = luh3treepct - currenttreepct;
                      if (luh3forest > LUH3StateThreshold && removetreepct > CLMPCTThreshold) {
	                  removeluh3treepct = removetreepct * luh3forest * scaledforesttreepct / luh3treepct;
			  if (removeluh3treepct > 0.5 * luh3forest * scaledforesttreepct) {
			      removeluh3treepct = 0.5 * luh3forest * scaledforesttreepct;
			  }
                          scaledforesttreepct = scaledforesttreepct - removeluh3treepct / luh3forest;
                          removetreepct -= removeluh3treepct;
                          if (scaledforesttreepct < 0.0) {
                              scaledforesttreepct = 0.0;
		          }
		      }
                      if (luh3nonforest > LUH3StateThreshold && removetreepct > CLMPCTThreshold) {
                          removeluh3treepct = removetreepct * luh3nonforest * scalednonforesttreepct / luh3treepct;
			  if (removeluh3treepct > luh3nonforest * scalednonforesttreepct) {
			      removeluh3treepct = luh3nonforest * scalednonforesttreepct;
			  }
                          scalednonforesttreepct = scalednonforesttreepct - removeluh3treepct / luh3nonforest;
                          removetreepct -= removeluh3treepct;
                          if (scalednonforesttreepct < 0.0) {
                              scalednonforesttreepct = 0.0;
		          }
                      }
                      if (luh3range > LUH3StateThreshold && removetreepct > CLMPCTThreshold) {
                          removeluh3treepct = removetreepct * luh3range * scaledrangetreepct / luh3treepct;
			  if (removeluh3treepct > luh3range * scaledrangetreepct) {
			      removeluh3treepct = luh3range * scaledrangetreepct;
			  }
                          scaledrangetreepct = scaledrangetreepct - removeluh3treepct / luh3range;
                          removetreepct -= removeluh3treepct;
                          if (scaledrangetreepct < 0.0) {
                              scaledrangetreepct = 0.0;
		          }
                      }
                  }
		  else {
		      if (currenttreepct > luh3treepct && luh3treepct > CLMPCTThreshold) { /* increase scaled tree pct values */
                          addtreepct = currenttreepct - luh3treepct;
                          if (luh3forest > LUH3StateThreshold && addtreepct > CLMPCTThreshold) {
		              addluh3treepct = addtreepct * luh3forest / luh3total;
                              scaledforesttreepct = scaledforesttreepct + addluh3treepct / luh3forest;
                              addtreepct -= addluh3treepct;
		              if (scaledforesttreepct > 100.0) {
			          scaledforesttreepct = 100.0;
			      }
                          }
                          if (luh3nonforest > LUH3StateThreshold && addtreepct > CLMPCTThreshold) {
			      addluh3treepct = addtreepct * luh3nonforest / luh3total;
                              scalednonforesttreepct = scalednonforesttreepct + addluh3treepct / luh3nonforest;
                              addtreepct -= addluh3treepct;
		              if (scalednonforesttreepct > 100.0) {
			          scalednonforesttreepct = 100.0;
			      }
                          }
                          if (luh3range > LUH3StateThreshold && addtreepct > CLMPCTThreshold) {
			      addluh3treepct = addtreepct * luh3range / luh3total;
                              scaledrangetreepct = scaledrangetreepct + addluh3treepct / luh3range;
                              addtreepct -= addluh3treepct;
		              if (scaledrangetreepct > 100.0) {
			          scaledrangetreepct = 100.0;
			      }
                          }
                      }
                  }
		  
	  
	          luh3herbpct = luh3forest * scaledforestherbpct + luh3nonforest * scalednonforestherbpct + luh3range * scaledrangeherbpct;
                  if (luh3herbpct > currentherbpct) { /* remove herb pct values */
                      removeherbpct = luh3herbpct - currentherbpct;
                      if (luh3forest > LUH3StateThreshold && removeherbpct > CLMPCTThreshold) {
                          removeluh3herbpct = removeherbpct * luh3forest * scaledforestherbpct / luh3herbpct;
  			  if (removeluh3herbpct > 0.2 * luh3forest * scaledforestherbpct) {
			      removeluh3herbpct = 0.2 * luh3forest * scaledforestherbpct;
			  }
                          scaledforestherbpct = scaledforestherbpct - removeluh3herbpct / luh3forest;
                          removeherbpct -= removeluh3herbpct;
                          if (scaledforestherbpct < 0.0) {
                              scaledforestherbpct = 0.0;
		          }
		      }
                      if (luh3nonforest > LUH3StateThreshold && removeherbpct > CLMPCTThreshold) {
		          removeluh3herbpct = removeherbpct * luh3nonforest * scalednonforestherbpct / luh3herbpct;
  			  if (removeluh3herbpct > 0.6 * luh3nonforest * scalednonforestherbpct) {
			      removeluh3herbpct = 0.6 * luh3nonforest * scalednonforestherbpct;
			  }
		          scalednonforestherbpct = scalednonforestherbpct - removeluh3herbpct / luh3nonforest;
		          removeherbpct -= removeluh3herbpct;
		          if (scalednonforestherbpct < 0.0) {
		              scalednonforestherbpct = 0.0;
		          }
                      }
                      if (luh3range > LUH3StateThreshold && removeherbpct > CLMPCTThreshold) {
		          removeluh3herbpct = removeherbpct * luh3range * scaledrangeherbpct / luh3herbpct;
  			  if (removeluh3herbpct > 0.8 * luh3range * scaledrangeherbpct) {
			      removeluh3herbpct = 0.8 * luh3range * scaledrangeherbpct;
			  }
		          scaledrangeherbpct = scaledrangeherbpct - removeluh3herbpct / luh3range;
		          removeherbpct -= removeluh3herbpct;
		          if (scaledrangeherbpct < 0.0) {
		              scaledrangeherbpct = 0.0;
		          }
		      }
                  }
		  else {
		      if (currentherbpct > luh3herbpct && luh3herbpct > CLMPCTThreshold) { /* increase scaled herb pct values */
                          addherbpct = currentherbpct - luh3herbpct;
                          if (luh3forest > LUH3StateThreshold && addherbpct > CLMPCTThreshold) {
                              addluh3herbpct = addherbpct * luh3forest / luh3total;
                              scaledforestherbpct = scaledforestherbpct + addluh3herbpct / luh3forest;
                              addherbpct -= addluh3herbpct;
		              if (scaledforestherbpct > 100.0) {
			          scaledforestherbpct = 100.0;
			      }
                          }
                          if (luh3nonforest > LUH3StateThreshold && addherbpct > CLMPCTThreshold) {
			      addluh3herbpct = addherbpct * luh3nonforest / luh3total;
                              scalednonforestherbpct = scalednonforestherbpct + addluh3herbpct / luh3nonforest;
                              addherbpct -= addluh3herbpct;
		              if (scalednonforestherbpct > 100.0) {
			          scalednonforestherbpct = 100.0;
			      }
                          }
                          if (luh3range > LUH3StateThreshold && addherbpct > CLMPCTThreshold) {
			      addluh3herbpct = addherbpct * luh3range / luh3total;
                              scaledrangeherbpct = scaledrangeherbpct + addluh3herbpct / luh3range;
                              addherbpct -= addluh3herbpct;
		              if (scaledrangeherbpct > 100.0) {
			          scaledrangeherbpct = 100.0;
			      }
                          }
                      }
                  }
		  		  
	          luh3barepct = luh3forest * scaledforestbarepct + luh3nonforest * scalednonforestbarepct + luh3range * scaledrangebarepct;
                  if (luh3barepct > currentbarepct) { /* remove bare pct values */
                      removebarepct = luh3barepct - currentbarepct;
                      if (luh3forest > LUH3StateThreshold && removebarepct > CLMPCTThreshold) {
                          removeluh3barepct = removebarepct * luh3forest * scaledforestbarepct / luh3barepct;
  			  if (removeluh3barepct > luh3forest * scaledforestbarepct) {
			      removeluh3barepct = luh3forest * scaledforestbarepct;
			  }
                          scaledforestbarepct = scaledforestbarepct - removeluh3barepct / luh3forest;
			  removebarepct -= removeluh3barepct;
			  if (scaledforestbarepct < 0.0) {
			      scaledforestbarepct = 0.0;
		          }
                      }
                      if (luh3nonforest > LUH3StateThreshold && removebarepct > CLMPCTThreshold) {
			  removeluh3barepct = removebarepct * luh3nonforest * scalednonforestbarepct / luh3barepct;
  			  if (removeluh3barepct > luh3nonforest * scalednonforestbarepct) {
			      removeluh3barepct = luh3nonforest * scalednonforestbarepct;
			  }
			  scalednonforestbarepct = scalednonforestbarepct - removeluh3barepct / luh3nonforest;
			  removebarepct -= removeluh3barepct;
			  if (scalednonforestbarepct < 0.0) {
			      scalednonforestbarepct = 0.0;
		          }
                      }
                      if (luh3range > LUH3StateThreshold && removebarepct > CLMPCTThreshold) {
			  removeluh3barepct = removebarepct * luh3range * scaledrangebarepct / luh3barepct;
  			  if (removeluh3barepct > luh3range * scaledrangebarepct) {
			      removeluh3barepct = luh3range * scaledrangebarepct;
			  }
			  scaledrangebarepct = scaledrangebarepct - removeluh3barepct / luh3range;
			  removebarepct -= removeluh3barepct;
			  if (scaledrangebarepct < 0.0) {
			      scaledrangebarepct = 0.0;
		          }
	              }
                  }
		  else {
		      if (currentbarepct > luh3barepct && luh3barepct > CLMPCTThreshold) { /* increase scaled bare pct values */
                          addbarepct = currentbarepct - luh3barepct;
                          if (luh3forest > LUH3StateThreshold && addbarepct > CLMPCTThreshold) {
		              addluh3barepct = addbarepct * luh3forest / luh3total;
                              scaledforestbarepct = scaledforestbarepct + addluh3barepct / luh3forest;
                              addbarepct -= addluh3barepct;
		              if (scaledforestbarepct > 100.0) {
			          scaledforestbarepct = 100.0;
			      }
                          }
                          if (luh3nonforest > LUH3StateThreshold && addbarepct > CLMPCTThreshold) {
			      addluh3barepct = addbarepct * luh3nonforest / luh3total;
                              scalednonforestbarepct = scalednonforestbarepct + addluh3barepct / luh3nonforest;
                              addbarepct -= addluh3barepct;
		              if (scalednonforestbarepct > 100.0) {
			          scalednonforestbarepct = 100.0;
			      }
                          }
                          if (luh3range > LUH3StateThreshold && addbarepct > CLMPCTThreshold) {
			      addluh3barepct = addbarepct * luh3range / luh3total;
                              scaledrangebarepct = scaledrangebarepct + addluh3barepct / luh3range;
                              addbarepct -= addluh3barepct;
		              if (scaledrangebarepct > 100.0) {
			          scaledrangebarepct = 100.0;
			      }
                          }
                      }
                  }

                  if (luh3forest <= LUH3StateThreshold && luh3nonforest > LUH3StateThreshold) {
		      scaledforesttreepct = nonforesttreepct;
                      scaledforestherbpct = nonforestherbpct;
		      scaledforestbarepct = nonforestbarepct;
		  }

		  if (luh3nonforest <= LUH3StateThreshold && luh3range > LUH3StateThreshold) {
		      scalednonforesttreepct = rangetreepct;
                      scalednonforestherbpct = rangeherbpct;
		      scalednonforestbarepct = rangebarepct;
		  }

              }


              outFORESTPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = scaledforesttreepct;
              outFORESTPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] = scaledforestherbpct;
	      outFORESTPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix] = scaledforestbarepct;
	      
              outNONFORESTPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = scalednonforesttreepct;
              outNONFORESTPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] = scalednonforestherbpct;
	      outNONFORESTPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix] = scalednonforestbarepct;
              
	      outRANGEPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = scaledrangetreepct;
              outRANGEPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] = scaledrangeherbpct;
	      outRANGEPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix] = scaledrangebarepct;
	  }
      }
  }

  return 0;
  	      
}


int rebalanceLUHScaledVegGrids() {

  long clm5lin, clm5pix;
  long pftid, cftid;
  float landmask, landfrac, area, precipann, tempaverage;
  float forestallpct, foresttreepct, forestherbpct, forestbarepct;
  float nonforestallpct, nonforesttreepct, nonforestherbpct, nonforestbarepct;
  float rangeallpct, rangetreepct, rangeherbpct, rangebarepct;
  float removepct, addpct;
  int precipannindex, tempaverageindex;

  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
          landmask = inLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix];
          landfrac = inLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix];
	  if (landfrac > 0.0) {
	  
	      foresttreepct = outFORESTPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      forestherbpct = outFORESTPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
	      forestbarepct = outFORESTPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix];
	      
	      nonforesttreepct = outNONFORESTPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      nonforestherbpct = outNONFORESTPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
	      nonforestbarepct = outNONFORESTPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix];
	      
	      rangetreepct = outRANGEPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      rangeherbpct = outRANGEPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
	      rangebarepct = outRANGEPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix];
	      
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
	      
              outFORESTPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = foresttreepct;
              outFORESTPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] = forestherbpct;
	      outFORESTPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix] = forestbarepct;
	      
              outNONFORESTPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = nonforesttreepct;
              outNONFORESTPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] = nonforestherbpct;
	      outNONFORESTPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix] = nonforestbarepct;
              
	      outRANGEPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = rangetreepct;
              outRANGEPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] = rangeherbpct;
	      outRANGEPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix] = rangebarepct;
	  }
      }
  }

  return 0;
  	      
}


int writescaledforestgrids(char *databasestr, char *seriesname, char *yearname) {

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


int writescalednonforestgrids(char *databasestr, char *seriesname, char *yearname) {

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


int writescaledrangegrids(char *databasestr, char *seriesname, char *yearname) {

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
  char workrawcombvcfinputdirname[1024];
  char workrawvcfoutputdirname[1024];
  char luh3rawinputreferenceyearstr[256];
  char clm5modisrawinputreferenceyearstr[256];
  char workrawclm5inputreferenceyearstr[256];
  char workrawcombvcfinputreferenceyearstr[256];
  char workrawvcfoutputreferenceyearstr[256];
  
  if(narg != 8){
        printf("Usage createcurrentpftsCTSM52Deg025bin luh3rawfile luh3rawdir timeseriesrawfile timeseriesrawdir workrawfile workrawdir potvegallscaledvcfnamelistfile\n");
        return 0;
  }
  
  readnamelistfile(argv[7]);

  readluh3rawinputlutfile(argv[1],argv[2],namelistluh3rawinput);
  readclm5modisrawinputlutfile(argv[3],argv[4],namelistclm5modisrawinput);
  readworkrawclm5inputlutfile(argv[5],argv[6],namelistworkrawclm5input);
  readworkrawcombvcfinputlutfile(argv[5],argv[6],namelistworkrawcombvcfinput);
  readworkrawvcfoutputlutfile(argv[5],argv[6],namelistworkrawvcfoutput);
  luh3rawinputreferenceyear = atoi(namelistluh3rawinputreferenceyear);
  sprintf(luh3rawinputreferenceyearstr,"%04d",luh3rawinputreferenceyear);
  clm5modisrawinputreferenceyear = atoi(namelistclm5modisrawinputreferenceyear);
  sprintf(clm5modisrawinputreferenceyearstr,"%04d",clm5modisrawinputreferenceyear);
  workrawclm5inputreferenceyear = atoi(namelistworkrawclm5inputreferenceyear);
  sprintf(workrawclm5inputreferenceyearstr,"%04d",workrawclm5inputreferenceyear);
  workrawcombvcfinputreferenceyear = atoi(namelistworkrawcombvcfinputreferenceyear);
  sprintf(workrawcombvcfinputreferenceyearstr,"%04d",workrawcombvcfinputreferenceyear);
  
  initializeGrids();

  sprintf(workrawclm5inputdirname,"%s/%s",workrawclm5inputdb,workrawclm5inputname);
  sprintf(workrawcombvcfinputdirname,"%s/%s",workrawcombvcfinputdb,workrawcombvcfinputname);
  readlandGrids(workrawclm5inputdirname,namelistworkrawclm5inputname,workrawclm5inputreferenceyearstr);
  readcurrentluh3Grids(luh3rawinputdb,luh3rawinputname,luh3rawinputreferenceyearstr);
  readcurrentvcfGrids(workrawclm5inputdirname,namelistworkrawclm5inputvcfname,workrawclm5inputreferenceyearstr);      
  readcombinedforestvcfGrids(workrawcombvcfinputdirname,namelistworkrawcombvcfinputforestname,workrawcombvcfinputreferenceyearstr);    
  readcombinednonforestvcfGrids(workrawcombvcfinputdirname,namelistworkrawcombvcfinputnonforestname,workrawcombvcfinputreferenceyearstr);    
  readcombinedrangevcfGrids(workrawcombvcfinputdirname,namelistworkrawcombvcfinputrangename,workrawcombvcfinputreferenceyearstr);    
  readclimvarseries(clm5modisrawinputdb,clm5modisrawinputname,"CLIM");
  readlaivarseries(clm5modisrawinputdb,clm5modisrawinputname,"FILLCLIM");
  
  genLUHLandGrids();
  genLUHScaledVegGrids();
  rebalanceLUHScaledVegGrids();
  
  sprintf(workrawvcfoutputreferenceyearstr,"%s",workrawclm5inputreferenceyearstr);
  sprintf(workrawvcfoutputdirname,"%s/%s",workrawvcfoutputdb,workrawvcfoutputname);
  writescaledforestgrids(workrawvcfoutputdirname,namelistworkrawvcfoutputforestname,workrawvcfoutputreferenceyearstr);
  writescalednonforestgrids(workrawvcfoutputdirname,namelistworkrawvcfoutputnonforestname,workrawvcfoutputreferenceyearstr);
  writescaledrangegrids(workrawvcfoutputdirname,namelistworkrawvcfoutputrangename,workrawvcfoutputreferenceyearstr);

  return 1;
  
}
