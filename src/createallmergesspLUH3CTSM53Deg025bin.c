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
char namelistluh3sspallmergename[1024];
char namelistworkrawmergeoutput[1024];
char namelistworkrawmergeoutputname[1024];
char namelistworkrawmergeoutputstartyear[1024];
char namelistworkrawmergeoutputendyear[1024];
char namelistluh3rawinput[1024];
char namelistluh3rawinputstartyear[1024];
char namelistluh3rawinputendyear[1024];
char namelistluh3rawinputreferenceyear[1024];
char namelistworkrawclm5input[1024];
char namelistworkrawclm5inputname[1024];
char namelistworkrawclm5inputreferenceyear[1024];
char namelistclm5modisrawinput[1024];
char namelistclm5modisrawinputreferenceyear[1024];

int workrawmergeoutputstartyear;
int workrawmergeoutputendyear;
int luh3rawinputstartyear;
int luh3rawinputendyear;
int luh3rawinputreferenceyear;
int clm5modisrawinputreferenceyear;
int workrawclm5inputreferenceyear;

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

char workrawmergeoutputname[1024];
char workrawmergeoutputdb[1024];
int workrawmergeoutputstartyear;
int workrawmergeoutputendyear;

char workmergename[1024];
char workmergecomponent[1024];
int workmergestartyear;
int workmergeendyear;
float workmergereforestscalar;

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

float inLUH3BASEPRIMFGrid[MAXCLMPIX * MAXCLMLIN];
float inLUH3BASESECDFGrid[MAXCLMPIX * MAXCLMLIN];

float inLUH3CURRENTPRIMFGrid[MAXCLMPIX * MAXCLMLIN];
float inLUH3CURRENTSECDFGrid[MAXCLMPIX * MAXCLMLIN];

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

float outMERGEFRACGrid[MAXCLMPIX * MAXCLMLIN];
float maxMERGEFRACGrid[MAXCLMPIX * MAXCLMLIN];

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
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3sspallmergename);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawmergeoutput);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawmergeoutputname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawmergeoutputstartyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawmergeoutputendyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3rawinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3rawinputstartyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3rawinputendyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3rawinputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5input);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5inputname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5inputreferenceyear);
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


int readworkrawmergeoutputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *workrawmergeoutputvaluefile;
  int foundworkrawmergeoutput, itemsfound;
  char fullfilenamestr[1024], inworkrawmergeoutputname[1024], inworkrawmergeoutputdb[1024];
  char inworkrawmergeoutputstartyear[1024], inworkrawmergeoutputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  workrawmergeoutputvaluefile = fopen(fullfilenamestr,"r");
  
  foundworkrawmergeoutput = 0;

  while (foundworkrawmergeoutput == 0) {
      itemsfound = fscanf(workrawmergeoutputvaluefile,"%s %s %s %s",inworkrawmergeoutputname, inworkrawmergeoutputdb, inworkrawmergeoutputstartyear, inworkrawmergeoutputendyear);
      if (itemsfound != 4) {
          printf("Error: workrawmergeoutput %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inworkrawmergeoutputname) == 0) {
          printf("Processing workrawmergeoutput: %s\n",inworkrawmergeoutputname);
	  sprintf(workrawmergeoutputname,"%s",inworkrawmergeoutputname);
          if (strcmp(inworkrawmergeoutputdb,"<workrawdir>") == 0) {
              sprintf(workrawmergeoutputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(workrawmergeoutputdb,"%s",inworkrawmergeoutputdb);
	  }
          workrawmergeoutputstartyear = atoi(inworkrawmergeoutputstartyear);
          workrawmergeoutputendyear = atoi(inworkrawmergeoutputendyear);
	  foundworkrawmergeoutput = 1;
      }
  }  
  
  return 0;

}


int readworkmergelutfile(char *filenamestr, char *namestr, char *componentstr) {

  FILE *workmergevaluefile;
  int foundworkmerge, itemsfound;
  char fullfilenamestr[1024], inworkmergename[1024], inworkmergecomponent[1024];
  char inworkmergestartyear[1024], inworkmergeendyear[1024];
  char inworkmergereforestscalar[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  workmergevaluefile = fopen(fullfilenamestr,"r");
  
  foundworkmerge = 0;

  while (foundworkmerge == 0) {
      itemsfound = fscanf(workmergevaluefile,"%s %s %s %s %s",inworkmergename, inworkmergecomponent, inworkmergestartyear, inworkmergeendyear, inworkmergereforestscalar);
      if (itemsfound != 5) {
          printf("Error: merge %s component %s Not Found\n",namestr,componentstr);
	  exit(0);
      }
      if (strcmp(namestr,inworkmergename) == 0 &&
	  strcmp(componentstr,inworkmergecomponent) == 0) {
          printf("Processing Name: %s Component: %s\n",inworkmergename,inworkmergecomponent);
	  sprintf(workmergename,"%s",inworkmergename);
          sprintf(workmergecomponent,"%s",inworkmergecomponent);
          workmergestartyear = atoi(inworkmergestartyear);
          workmergeendyear = atoi(inworkmergeendyear);
	  workmergereforestscalar = atof(inworkmergereforestscalar);
	  foundworkmerge = 1;
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

          outMERGEFRACGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
          maxMERGEFRACGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
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


int readbaseyearLUH3Grids(char *databasestr, char *seriesname, char *yearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  long outlin, outpix;
  double tempvalue, tempavgvalue, tempcount, tempnewavgvalue;
  
 
  sprintf(globalbinfilename,"%s/%s/%s.primf.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inLUH3BASEPRIMFGrid,sizeof(inLUH3BASEPRIMFGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.secdf.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inLUH3BASESECDFGrid,sizeof(inLUH3BASESECDFGrid),1,globalbinfile);  
  fclose(globalbinfile);

  return 0;
  
}


int readcurrentyearLUH3Grids(char *databasestr, char *seriesname, char *yearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  long outlin, outpix;
  double tempvalue, tempavgvalue, tempcount, tempnewavgvalue;
  
 
  sprintf(globalbinfilename,"%s/%s/%s.primf.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inLUH3CURRENTPRIMFGrid,sizeof(inLUH3CURRENTPRIMFGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.secdf.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inLUH3CURRENTSECDFGrid,sizeof(inLUH3CURRENTSECDFGrid),1,globalbinfile);  
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


int genLUHCurrentVegGrids(float startyearnumber, float endyearnumber, float yearnumber) {

  long ctsmlin, ctsmpix;
  long pftid, cftid;
  float landmask, landfrac, area, precipann, tempaverage, climatemask;
  float baseluh3primf, baseluh3secdf, currentluh3primf, currentluh3secdf;
  float changeprimf, changesecdf, changeforest;
  float yearmergefrac, mergefrac, startweight, endweight;
  
  
  for (ctsmlin = 0; ctsmlin < MAXCLMLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXCLMPIX; ctsmpix++) {
          landmask = inLANDMASKGrid[ctsmlin * MAXCLMPIX + ctsmpix];
          landfrac = inLANDFRACGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	  if (landfrac > 0.0) {
	  	      
              baseluh3primf = inLUH3BASEPRIMFGrid[ctsmlin * MAXCLMPIX + ctsmpix];
              baseluh3secdf = inLUH3BASESECDFGrid[ctsmlin * MAXCLMPIX + ctsmpix];
              currentluh3primf = inLUH3CURRENTPRIMFGrid[ctsmlin * MAXCLMPIX + ctsmpix];
              currentluh3secdf = inLUH3CURRENTSECDFGrid[ctsmlin * MAXCLMPIX + ctsmpix];
	      changeprimf = currentluh3primf - baseluh3primf;
	      changesecdf = currentluh3secdf - baseluh3secdf;
	      changeforest = changeprimf + changesecdf;
	      
	      if (changeforest > 0.0) {
	          mergefrac = workmergereforestscalar * changeforest / (baseluh3secdf + changeforest);
              }
	      else {
	          mergefrac = 0.0;
              }
	      
	      if (mergefrac < 0.0) {
	          mergefrac = 0.0;
              }
	      
	      if (mergefrac > 1.0) {
	          mergefrac = 1.0;
	      }
	      
	      if (mergefrac > maxMERGEFRACGrid[ctsmlin * MAXCLMPIX + ctsmpix]) {
	          maxMERGEFRACGrid[ctsmlin * MAXCLMPIX + ctsmpix] = mergefrac;
              }
	      if (maxMERGEFRACGrid[ctsmlin * MAXCLMPIX + ctsmpix] > mergefrac) {
	          mergefrac = maxMERGEFRACGrid[ctsmlin * MAXCLMPIX + ctsmpix];
              }
	      
              outMERGEFRACGrid[ctsmlin * MAXCLMPIX + ctsmpix] = mergefrac;
	  }
          else {
              outMERGEFRACGrid[ctsmlin * MAXCLMPIX + ctsmpix] = 0.0;
	  }
      }      
  }

  return 0;
  	      
}


int writemergegrids(char *databasestr, char *seriesname, char *componentname, char *yearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  long outlin, outpix;
  long clmlin, clmpix;
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];

  tempoutGrid = (float *) malloc(OUTDATASIZE);
  
  sprintf(globalbinfilename,"%s/%s/%s.%s.%s.dat",databasestr,seriesname,seriesname,componentname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outMERGEFRACGrid,sizeof(outMERGEFRACGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  return 0;

}


int main(long narg, char **argv) {

  char workrawclm5inputdirname[1024];
  char workrawmergeoutputdirname[1024];
  char luh3rawinputreferenceyearstr[256];
  char luh3rawinputcurrentyearstr[256];
  char clm5modisrawinputreferenceyearstr[256];
  char workrawclm5inputreferenceyearstr[256];
  char workrawmergeoutputcurrentyearstr[256];
  int startyear, endyear, currentyear;
    
  if(narg != 10){
        printf("Usage createfuturepftsCTSM53Deg025bin luh3rawfile luh3rawdir timeseriesrawfile timeseriesrawdir workrawfile workrawdir workmergefile mergecomponent luh3allvcfmergenamelistfile\n");
        return 0;
  }
  
  readnamelistfile(argv[9]);

  readluh3rawinputlutfile(argv[1],argv[2],namelistluh3rawinput);
  readclm5modisrawinputlutfile(argv[3],argv[4],namelistclm5modisrawinput);
  readworkrawclm5inputlutfile(argv[5],argv[6],namelistworkrawclm5input);
  readworkrawmergeoutputlutfile(argv[5],argv[6],namelistworkrawmergeoutput);
  readworkmergelutfile(argv[7],namelistworkrawmergeoutput,argv[8]);
  luh3rawinputreferenceyear = atoi(namelistluh3rawinputreferenceyear);
  sprintf(luh3rawinputreferenceyearstr,"%04d",luh3rawinputreferenceyear);
  clm5modisrawinputreferenceyear = atoi(namelistclm5modisrawinputreferenceyear);
  sprintf(clm5modisrawinputreferenceyearstr,"%04d",clm5modisrawinputreferenceyear);
  workrawclm5inputreferenceyear = atoi(namelistworkrawclm5inputreferenceyear);
  sprintf(workrawclm5inputreferenceyearstr,"%04d",workrawclm5inputreferenceyear);

  startyear = atoi(namelistworkrawmergeoutputstartyear);
  endyear = atoi(namelistworkrawmergeoutputendyear);
  
  initializeGrids();
  
  sprintf(workrawclm5inputdirname,"%s/%s",workrawclm5inputdb,workrawclm5inputname);
  readlandGrids(workrawclm5inputdirname,namelistworkrawclm5inputname,workrawclm5inputreferenceyearstr);
  readbaseyearLUH3Grids(luh3rawinputdb,luh3rawinputname,luh3rawinputreferenceyearstr);
  readclimvarseries(clm5modisrawinputdb,clm5modisrawinputname,"CLIM");
  readlaivarseries(clm5modisrawinputdb,clm5modisrawinputname,"FILLCLIM");

  genLUHLandGrids();
  
  sprintf(workrawmergeoutputdirname,"%s/%s",workrawmergeoutputdb,workrawmergeoutputname);
  for (currentyear = startyear; currentyear <= endyear; currentyear++) {
  
      sprintf(luh3rawinputcurrentyearstr,"%04d",currentyear);
      sprintf(workrawmergeoutputcurrentyearstr,"%04d",currentyear);
      readcurrentyearLUH3Grids(luh3rawinputdb,luh3rawinputname,luh3rawinputcurrentyearstr);
      genLUHCurrentVegGrids(workmergestartyear,workmergeendyear,currentyear);
      writemergegrids(workrawmergeoutputdirname,namelistworkrawmergeoutputname,workmergecomponent,workrawmergeoutputcurrentyearstr);

  }
  
  return 1;
  
}
