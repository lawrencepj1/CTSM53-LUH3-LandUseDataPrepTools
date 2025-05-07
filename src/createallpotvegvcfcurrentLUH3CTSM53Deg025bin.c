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
char namelistcurrentvcfname[1024];
char namelistworkrawoutput[1024];
char namelistworkrawinputname[1024];
char namelistworkrawinput[1024];
char namelistworkrawoutputname[1024];
char namelistclm5currentrawinput[1024];
char namelistclm5referenceyear[1024];

int clm5referenceyear;
int luhreferenceyear;
int luhprocessyear;

char clm5currentrawinputname[1024];
char clm5currentrawinputdb[1024];
int clm5currentrawinputstartyear;
int clm5currentrawinputendyear;

char workrawinputname[1024];
char workrawinputdb[1024];
int workrawinputstartyear;
int workrawinputendyear;

char workrawoutputname[1024];
char workrawoutputdb[1024];
int workrawoutputstartyear;
int workrawoutputendyear;

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
float inCTSMPCTPFTGrid[MAXPFT][MAXCLMPIX * MAXCLMLIN];
float inCTSMPCTCFTGrid[MAXCFT][MAXCLMPIX * MAXCLMLIN];

float inCTSMPCTTREEGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMPCTHERBGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMPCTGRASSGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMPCTBAREGrid[MAXCLMPIX * MAXCLMLIN];

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

float outPOTVEGLANDMASKGrid[MAXCLMPIX * MAXCLMLIN];
float outPOTVEGLANDFRACGrid[MAXCLMPIX * MAXCLMLIN];
float outPOTVEGAREAGrid[MAXCLMPIX * MAXCLMLIN];
float outPOTVEGPCTGLACIERGrid[MAXCLMPIX * MAXCLMLIN];
float outPOTVEGPCTLAKEGrid[MAXCLMPIX * MAXCLMLIN];
float outPOTVEGPCTWETLANDGrid[MAXCLMPIX * MAXCLMLIN];
float outPOTVEGPCTURBANGrid[MAXCLMPIX * MAXCLMLIN];
float outPOTVEGPCTNATVEGGrid[MAXCLMPIX * MAXCLMLIN];
float outPOTVEGPCTCROPGrid[MAXCLMPIX * MAXCLMLIN];
float outPOTVEGPCTPFTGrid[MAXPFT][MAXCLMPIX * MAXCLMLIN];
float outPOTVEGPCTCFTGrid[MAXCFT][MAXCLMPIX * MAXCLMLIN];

float outPOTVEGPCTTREEGrid[MAXCLMPIX * MAXCLMLIN];
float outPOTVEGPCTHERBGrid[MAXCLMPIX * MAXCLMLIN];
float outPOTVEGPCTGRASSGrid[MAXCLMPIX * MAXCLMLIN];
float outPOTVEGPCTBAREGrid[MAXCLMPIX * MAXCLMLIN];

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
  fscanf(namelistfile,"%s %s",fieldname,namelistcurrentvcfname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawoutput);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawoutputname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawinputname);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5currentrawinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5referenceyear);
  
  return 0;

}


int readclm5currentrawinputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *clm5currentrawinputvaluefile;
  int foundclm5currentrawinput, itemsfound;
  char fullfilenamestr[1024], inclm5currentrawinputname[1024], inclm5currentrawinputdb[1024];
  char inclm5currentrawinputstartyear[1024], inclm5currentrawinputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  clm5currentrawinputvaluefile = fopen(fullfilenamestr,"r");
  
  foundclm5currentrawinput = 0;

  while (foundclm5currentrawinput == 0) {
      itemsfound = fscanf(clm5currentrawinputvaluefile,"%s %s %s %s",inclm5currentrawinputname, inclm5currentrawinputdb, inclm5currentrawinputstartyear, inclm5currentrawinputendyear);
      if (itemsfound != 4) {
          printf("Error: clm5currentrawinput %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inclm5currentrawinputname) == 0) {
          printf("Processing clm5currentrawinput: %s\n",inclm5currentrawinputname);
	  sprintf(clm5currentrawinputname,"%s",inclm5currentrawinputname);
          if (strcmp(inclm5currentrawinputdb,"<timeseriesrawdir>") == 0) {
              sprintf(clm5currentrawinputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(clm5currentrawinputdb,"%s",inclm5currentrawinputdb);
	  }
          clm5currentrawinputstartyear = atoi(inclm5currentrawinputstartyear);
          clm5currentrawinputendyear = atoi(inclm5currentrawinputendyear);
	  foundclm5currentrawinput = 1;
      }
  }  
  
  return 0;

}


int readworkrawoutputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *workrawoutputvaluefile;
  int foundworkrawoutput, itemsfound;
  char fullfilenamestr[1024], inworkrawoutputname[1024], inworkrawoutputdb[1024];
  char inworkrawoutputstartyear[1024], inworkrawoutputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  workrawoutputvaluefile = fopen(fullfilenamestr,"r");
  
  foundworkrawoutput = 0;

  while (foundworkrawoutput == 0) {
      itemsfound = fscanf(workrawoutputvaluefile,"%s %s %s %s",inworkrawoutputname, inworkrawoutputdb, inworkrawoutputstartyear, inworkrawoutputendyear);
      if (itemsfound != 4) {
          printf("Error: workrawoutput %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inworkrawoutputname) == 0) {
          printf("Processing workrawoutput: %s\n",inworkrawoutputname);
	  sprintf(workrawoutputname,"%s",inworkrawoutputname);
          if (strcmp(inworkrawoutputdb,"<workrawdir>") == 0) {
              sprintf(workrawoutputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(workrawoutputdb,"%s",inworkrawoutputdb);
	  }
          workrawoutputstartyear = atoi(inworkrawoutputstartyear);
          workrawoutputendyear = atoi(inworkrawoutputendyear);
	  foundworkrawoutput = 1;
      }
  }  
  
  return 0;

}


int readworkrawinputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *workrawinputvaluefile;
  int foundworkrawinput, itemsfound;
  char fullfilenamestr[1024], inworkrawinputname[1024], inworkrawinputdb[1024];
  char inworkrawinputstartyear[1024], inworkrawinputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  workrawinputvaluefile = fopen(fullfilenamestr,"r");
  
  foundworkrawinput = 0;

  while (foundworkrawinput == 0) {
      itemsfound = fscanf(workrawinputvaluefile,"%s %s %s %s",inworkrawinputname, inworkrawinputdb, inworkrawinputstartyear, inworkrawinputendyear);
      if (itemsfound != 4) {
          printf("Error: workrawinput %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inworkrawinputname) == 0) {
          printf("Processing workrawinput: %s\n",inworkrawinputname);
	  sprintf(workrawinputname,"%s",inworkrawinputname);
          if (strcmp(inworkrawinputdb,"<workrawdir>") == 0) {
              sprintf(workrawinputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(workrawinputdb,"%s",inworkrawinputdb);
	  }
          workrawinputstartyear = atoi(inworkrawinputstartyear);
          workrawinputendyear = atoi(inworkrawinputendyear);
	  foundworkrawinput = 1;
      }
  }  
  
  return 0;

}


int initializeGrids() {

  long clm5lin, clm5pix;
  long pftid, cftid;
  
  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
          inCTSMPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          inCTSMPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
	  inCTSMPCTGRASSGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          inCTSMPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
	  
          outPOTVEGLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGAREAGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGPCTGLACIERGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGPCTLAKEGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGPCTWETLANDGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGPCTURBANGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGPCTNATVEGGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGPCTCROPGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          for (pftid = 0; pftid < MAXPFT; pftid++) {
              outPOTVEGPCTPFTGrid[pftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
	  }
          for (cftid = 0; cftid < MAXCFT; cftid++) {
              outPOTVEGPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
	  }

          outPOTVEGPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGPCTGRASSGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;

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
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      sprintf(pftidstr,"%02d",pftid);
      sprintf(globalbinfilename,"%s/%s/%s.PCTNATPFT%s.%s.dat",databasestr,seriesname,seriesname,pftidstr,referenceyearname);
      printf("Reading: %s\n",globalbinfilename);
      globalbinfile = fopen(globalbinfilename,"r");
      fread(tempGrid,sizeof(tempGrid),1,globalbinfile);  
      fclose(globalbinfile);
      for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
          for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
	      inCTSMPCTPFTGrid[pftid][clm5lin * MAXCLMPIX + clm5pix] = tempGrid[clm5lin * MAXCLMPIX + clm5pix];
	      if (pftid >= 1 && pftid <= 8) {
	          inCTSMPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
              }
	      if (pftid >= 9 && pftid <= 14) {
	          inCTSMPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
              }
	      if (pftid >= 12 && pftid <= 14) {
	          inCTSMPCTGRASSGrid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
              }
	      if (pftid == 0) {
	          inCTSMPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
              }
	  }
      }
  }
  
  for (cftid = 0; cftid < MAXCFT; cftid++) {
      sprintf(cftidstr,"%02d",cftid);
      sprintf(globalbinfilename,"%s/%s/%s.PCTCFT%s.%s.dat",databasestr,seriesname,seriesname,cftidstr,referenceyearname);
      printf("Reading: %s\n",globalbinfilename);
      globalbinfile = fopen(globalbinfilename,"r");
      fread(tempGrid,sizeof(tempGrid),1,globalbinfile);  
      fclose(globalbinfile);
      for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
          for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
	      inCTSMPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = tempGrid[clm5lin * MAXCLMPIX + clm5pix];
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


int genLUHLandGrids() {

  long clm5lin, clm5pix;
  long pftid, cftid;
  float landmask, landfrac, area, glacierpct, lakepct, wetlandpct, natvegpct;
  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
          landmask = inLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix];
          landfrac = inLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix];
	  if (landfrac > 0.0) {
              outPOTVEGLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix] = landmask;
              outPOTVEGLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix] = landfrac;
	      outPOTVEGAREAGrid[clm5lin * MAXCLMPIX + clm5pix] = inAREAGrid[clm5lin * MAXCLMPIX + clm5pix];
	      glacierpct = inCTSMPCTGLACIERGrid[clm5lin * MAXCLMPIX + clm5pix];
	      lakepct = inCTSMPCTLAKEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      wetlandpct = inCTSMPCTWETLANDGrid[clm5lin * MAXCLMPIX + clm5pix];
	      natvegpct = 100.0 - glacierpct - lakepct - wetlandpct;
              outPOTVEGPCTGLACIERGrid[clm5lin * MAXCLMPIX + clm5pix] = glacierpct;
              outPOTVEGPCTLAKEGrid[clm5lin * MAXCLMPIX + clm5pix] = lakepct;
              outPOTVEGPCTWETLANDGrid[clm5lin * MAXCLMPIX + clm5pix] = wetlandpct;
              outPOTVEGPCTNATVEGGrid[clm5lin * MAXCLMPIX + clm5pix] = natvegpct;
	  }
      }
  }

  return 0;
  	      
}


int genLUHPotentialVegGrids() {

  long clm5lin, clm5pix;
  long pftid, cftid;
  float landmask, landfrac, area, precipann, tempaverage;
  float potvegtreepct, potvegherbpct, potvegbarepct;
  float currenttreepct, currentherbpct, currentgrasspct, currentbarepct;
  float climatemask;
  float newpotveg, availpotveg, removepotveg, addpotveg;
  int precipannindex, tempaverageindex;

  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
          landmask = inLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix];
          landfrac = inLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix];
	  if (landfrac > 0.0) {

	      currenttreepct = inCTSMPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      currentherbpct = inCTSMPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
	      currentgrasspct = inCTSMPCTGRASSGrid[clm5lin * MAXCLMPIX + clm5pix];
	      currentbarepct = inCTSMPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix];
	      
	      
              outPOTVEGPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = currenttreepct;
              outPOTVEGPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] = currentherbpct;
              outPOTVEGPCTGRASSGrid[clm5lin * MAXCLMPIX + clm5pix] = currentgrasspct;
	      outPOTVEGPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix] = currentbarepct;
	  }
      }
  }

  return 0;
  	      
}


int writepotveggrids(char *databasestr, char *seriesname, char *yearname) {

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
  fwrite(outPOTVEGLANDMASKGrid,sizeof(outPOTVEGLANDMASKGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.LANDFRAC.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outPOTVEGLANDFRACGrid,sizeof(outPOTVEGLANDFRACGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.AREA.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outPOTVEGAREAGrid,sizeof(outPOTVEGAREAGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTGLACIER.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outPOTVEGPCTGLACIERGrid,sizeof(outPOTVEGPCTGLACIERGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTLAKE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outPOTVEGPCTLAKEGrid,sizeof(outPOTVEGPCTLAKEGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTWETLAND.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outPOTVEGPCTWETLANDGrid,sizeof(outPOTVEGPCTWETLANDGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTURBAN.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outPOTVEGPCTURBANGrid,sizeof(outPOTVEGPCTURBANGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTNATVEG.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outPOTVEGPCTNATVEGGrid,sizeof(outPOTVEGPCTNATVEGGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTCROP.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outPOTVEGPCTCROPGrid,sizeof(outPOTVEGPCTCROPGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTTREE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outPOTVEGPCTTREEGrid,sizeof(outPOTVEGPCTTREEGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTHERB.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outPOTVEGPCTHERBGrid,sizeof(outPOTVEGPCTHERBGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTGRASS.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outPOTVEGPCTGRASSGrid,sizeof(outPOTVEGPCTGRASSGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTBARE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outPOTVEGPCTBAREGrid,sizeof(outPOTVEGPCTBAREGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  return 0;

}


int
main(long narg, char **argv) {

  char clm5referenceyearstr[256];
  char workrawinputreferenceyearstr[256];
  char workrawoutputreferenceyearstr[256];
  char workrawinputdirname[1024];
  char workrawoutputdirname[1024];
  
  if(narg != 6){
        printf("Usage createpotvegpftsCTSM53Deg025bin timeseriesrawfile timeseriesrawdir workrawfile workrawdir generatecurrentvcfnamelist\n");
        return 0;
  }
  
  readnamelistfile(argv[5]);

  readclm5currentrawinputlutfile(argv[1],argv[2],namelistclm5currentrawinput);
  readworkrawinputlutfile(argv[3],argv[4],namelistworkrawinput);
  readworkrawoutputlutfile(argv[3],argv[4],namelistworkrawoutput);
  clm5referenceyear = atoi(namelistclm5referenceyear);
  sprintf(clm5referenceyearstr,"%04d",clm5referenceyear);
  
  initializeGrids();
  
  sprintf(workrawinputreferenceyearstr,"%s",clm5referenceyearstr);
  sprintf(workrawinputdirname,"%s/%s",workrawinputdb,workrawinputname);
  readlandGrids(workrawinputdirname,namelistworkrawinputname,workrawinputreferenceyearstr);
  readclm5referenceGrids(workrawinputdirname,namelistworkrawinputname,workrawinputreferenceyearstr);    
  readclimvarseries(clm5currentrawinputdb,clm5currentrawinputname,"CLIM");
  readlaivarseries(clm5currentrawinputdb,clm5currentrawinputname,"FILLCLIM");
  
  genLUHLandGrids();
  genLUHPotentialVegGrids();
  
  sprintf(workrawoutputreferenceyearstr,"%s",clm5referenceyearstr);
  sprintf(workrawoutputdirname,"%s/%s",workrawoutputdb,workrawoutputname);
  writepotveggrids(workrawoutputdirname,namelistworkrawoutputname,workrawoutputreferenceyearstr);

  return 1;
  
}
