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
char namelistpotvegallcombinedvcfname[1024];
char namelistworkrawvcfoutput[1024];
char namelistworkrawvcfoutputforestname[1024];
char namelistworkrawvcfoutputnonforestname[1024];
char namelistworkrawvcfoutputrangename[1024];
char namelistworkrawclm5input[1024];
char namelistworkrawclm5inputname[1024];
char namelistworkrawclm5inputreferenceyear[1024];
char namelistworkrawvcfinput[1024];
char namelistworkrawvcfinputforestname[1024];
char namelistworkrawvcfinputnonforestname[1024];
char namelistworkrawvcfinputrangename[1024];
char namelistworkrawvcfinputreferenceyear[1024];
char namelistclm5currentrawinput[1024];
char namelistclm5currentrawinputreferenceyear[1024];

int workrawclm5inputreferenceyear;
int workrawvcfinputreferenceyear;

char clm5currentrawinputname[1024];
char clm5currentrawinputdb[1024];
int clm5currentrawinputstartyear;
int clm5currentrawinputendyear;

char workrawclm5inputname[1024];
char workrawclm5inputdb[1024];
int workrawclm5inputstartyear;
int workrawclm5inputendyear;

char workrawvcfinputname[1024];
char workrawvcfinputdb[1024];
int workrawvcfinputstartyear;
int workrawvcfinputendyear;

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
float inCLM5PCTGLACIERGrid[MAXCLMPIX * MAXCLMLIN];
float inCLM5PCTLAKEGrid[MAXCLMPIX * MAXCLMLIN];
float inCLM5PCTWETLANDGrid[MAXCLMPIX * MAXCLMLIN];
float inCLM5PCTURBANGrid[MAXCLMPIX * MAXCLMLIN];
float inCLM5PCTNATVEGGrid[MAXCLMPIX * MAXCLMLIN];
float inCLM5PCTCROPGrid[MAXCLMPIX * MAXCLMLIN];

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

float outCLM5LANDMASKGrid[MAXCLMPIX * MAXCLMLIN];
float outCLM5LANDFRACGrid[MAXCLMPIX * MAXCLMLIN];
float outCLM5AREAGrid[MAXCLMPIX * MAXCLMLIN];
float outCLM5PCTGLACIERGrid[MAXCLMPIX * MAXCLMLIN];
float outCLM5PCTLAKEGrid[MAXCLMPIX * MAXCLMLIN];
float outCLM5PCTWETLANDGrid[MAXCLMPIX * MAXCLMLIN];
float outCLM5PCTURBANGrid[MAXCLMPIX * MAXCLMLIN];
float outCLM5PCTNATVEGGrid[MAXCLMPIX * MAXCLMLIN];
float outCLM5PCTCROPGrid[MAXCLMPIX * MAXCLMLIN];
float outCLM5PCTPFTGrid[MAXPFT][MAXCLMPIX * MAXCLMLIN];
float outCLM5PCTCFTGrid[MAXCFT][MAXCLMPIX * MAXCLMLIN];

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
  fscanf(namelistfile,"%s %s",fieldname,namelistpotvegallcombinedvcfname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfoutput);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfoutputforestname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfoutputnonforestname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfoutputrangename);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5input);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5inputname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5inputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfinputforestname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfinputnonforestname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfinputrangename);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfinputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5currentrawinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5currentrawinputreferenceyear);
  
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


int readworkrawvcfinputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *workrawvcfinputvaluefile;
  int foundworkrawvcfinput, itemsfound;
  char fullfilenamestr[1024], inworkrawvcfinputname[1024], inworkrawvcfinputdb[1024];
  char inworkrawvcfinputstartyear[1024], inworkrawvcfinputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  workrawvcfinputvaluefile = fopen(fullfilenamestr,"r");
  
  foundworkrawvcfinput = 0;

  while (foundworkrawvcfinput == 0) {
      itemsfound = fscanf(workrawvcfinputvaluefile,"%s %s %s %s",inworkrawvcfinputname, inworkrawvcfinputdb, inworkrawvcfinputstartyear, inworkrawvcfinputendyear);
      if (itemsfound != 4) {
          printf("Error: workrawvcfinput %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inworkrawvcfinputname) == 0) {
          printf("Processing workrawvcfinput: %s\n",inworkrawvcfinputname);
	  sprintf(workrawvcfinputname,"%s",inworkrawvcfinputname);
          if (strcmp(inworkrawvcfinputdb,"<workrawdir>") == 0) {
              sprintf(workrawvcfinputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(workrawvcfinputdb,"%s",inworkrawvcfinputdb);
	  }
          workrawvcfinputstartyear = atoi(inworkrawvcfinputstartyear);
          workrawvcfinputendyear = atoi(inworkrawvcfinputendyear);
	  foundworkrawvcfinput = 1;
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
          outCLM5LANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outCLM5LANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outCLM5AREAGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outCLM5PCTGLACIERGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outCLM5PCTLAKEGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outCLM5PCTWETLANDGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outCLM5PCTURBANGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outCLM5PCTNATVEGGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outCLM5PCTCROPGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;

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
  fread(inCLM5PCTGLACIERGrid,sizeof(inCLM5PCTGLACIERGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTLAKE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inCLM5PCTLAKEGrid,sizeof(inCLM5PCTLAKEGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTWETLAND.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inCLM5PCTWETLANDGrid,sizeof(inCLM5PCTWETLANDGrid),1,globalbinfile);  
  fclose(globalbinfile);

  return 0;
  
}


int readinitialforestvcfGrids(char *databasestr,char *seriesname,char *referenceyearname) {

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


int readinitialnonforestvcfGrids(char *databasestr,char *seriesname,char *referenceyearname) {

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


int readinitialrangevcfGrids(char *databasestr,char *seriesname,char *referenceyearname) {

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
              outCLM5LANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix] = landmask;
              outCLM5LANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix] = landfrac;
	      outCLM5AREAGrid[clm5lin * MAXCLMPIX + clm5pix] = inAREAGrid[clm5lin * MAXCLMPIX + clm5pix];
	      glacierpct = inCLM5PCTGLACIERGrid[clm5lin * MAXCLMPIX + clm5pix];
	      lakepct = inCLM5PCTLAKEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      wetlandpct = inCLM5PCTWETLANDGrid[clm5lin * MAXCLMPIX + clm5pix];
	      natvegpct = 100.0 - glacierpct - lakepct - wetlandpct;
              outCLM5PCTGLACIERGrid[clm5lin * MAXCLMPIX + clm5pix] = glacierpct;
              outCLM5PCTLAKEGrid[clm5lin * MAXCLMPIX + clm5pix] = lakepct;
              outCLM5PCTWETLANDGrid[clm5lin * MAXCLMPIX + clm5pix] = wetlandpct;
              outCLM5PCTNATVEGGrid[clm5lin * MAXCLMPIX + clm5pix] = natvegpct;
	  }
      }
  }

  return 0;
  	      
}


int genLUHCurrentVegGrids() {

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
	  
	      foresttreepct = inFORESTPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      forestherbpct = inFORESTPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
	      forestbarepct = inFORESTPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix];
	      
	      nonforesttreepct = inNONFORESTPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      nonforestherbpct = inNONFORESTPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
	      nonforestbarepct = inNONFORESTPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix];
	      
	      rangetreepct = inRANGEPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      rangeherbpct = inRANGEPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
	      rangebarepct = inRANGEPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix];
	      
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


int writecombinedforestgrids(char *databasestr, char *seriesname, char *yearname) {

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
  fwrite(outCLM5LANDMASKGrid,sizeof(outCLM5LANDMASKGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.LANDFRAC.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5LANDFRACGrid,sizeof(outCLM5LANDFRACGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.AREA.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5AREAGrid,sizeof(outCLM5AREAGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTGLACIER.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTGLACIERGrid,sizeof(outCLM5PCTGLACIERGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTLAKE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTLAKEGrid,sizeof(outCLM5PCTLAKEGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTWETLAND.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTWETLANDGrid,sizeof(outCLM5PCTWETLANDGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTURBAN.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTURBANGrid,sizeof(outCLM5PCTURBANGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTNATVEG.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTNATVEGGrid,sizeof(outCLM5PCTNATVEGGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTCROP.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTCROPGrid,sizeof(outCLM5PCTCROPGrid),1,globalbinfile);  
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


int writecombinednonforestgrids(char *databasestr, char *seriesname, char *yearname) {

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
  fwrite(outCLM5LANDMASKGrid,sizeof(outCLM5LANDMASKGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.LANDFRAC.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5LANDFRACGrid,sizeof(outCLM5LANDFRACGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.AREA.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5AREAGrid,sizeof(outCLM5AREAGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTGLACIER.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTGLACIERGrid,sizeof(outCLM5PCTGLACIERGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTLAKE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTLAKEGrid,sizeof(outCLM5PCTLAKEGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTWETLAND.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTWETLANDGrid,sizeof(outCLM5PCTWETLANDGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTURBAN.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTURBANGrid,sizeof(outCLM5PCTURBANGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTNATVEG.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTNATVEGGrid,sizeof(outCLM5PCTNATVEGGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTCROP.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTCROPGrid,sizeof(outCLM5PCTCROPGrid),1,globalbinfile);  
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


int writecombinedrangegrids(char *databasestr, char *seriesname, char *yearname) {

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
  fwrite(outCLM5LANDMASKGrid,sizeof(outCLM5LANDMASKGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.LANDFRAC.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5LANDFRACGrid,sizeof(outCLM5LANDFRACGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.AREA.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5AREAGrid,sizeof(outCLM5AREAGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTGLACIER.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTGLACIERGrid,sizeof(outCLM5PCTGLACIERGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTLAKE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTLAKEGrid,sizeof(outCLM5PCTLAKEGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTWETLAND.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTWETLANDGrid,sizeof(outCLM5PCTWETLANDGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTURBAN.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTURBANGrid,sizeof(outCLM5PCTURBANGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTNATVEG.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTNATVEGGrid,sizeof(outCLM5PCTNATVEGGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTCROP.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outCLM5PCTCROPGrid,sizeof(outCLM5PCTCROPGrid),1,globalbinfile);  
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
  char workrawvcfinputdirname[1024];
  char workrawvcfoutputdirname[1024];
  char workrawclm5inputreferenceyearstr[256];
  char workrawvcfinputreferenceyearstr[256];
  char workrawvcfoutputreferenceyearstr[256];
  
  if(narg != 6){
        printf("Usage createpotvegpftsclm552Deg025bin timeseriesrawfile timeseriesrawdir workrawfile workrawdir generateallpotvegcombinedvcfnamelist\n");
        return 0;
  }
  
  readnamelistfile(argv[5]);

  readclm5currentrawinputlutfile(argv[1],argv[2],namelistclm5currentrawinput);
  readworkrawclm5inputlutfile(argv[3],argv[4],namelistworkrawclm5input);
  readworkrawvcfinputlutfile(argv[3],argv[4],namelistworkrawvcfinput);
  readworkrawvcfoutputlutfile(argv[3],argv[4],namelistworkrawvcfoutput);
  workrawclm5inputreferenceyear = atoi(namelistworkrawclm5inputreferenceyear);
  sprintf(workrawclm5inputreferenceyearstr,"%04d",workrawclm5inputreferenceyear);
  workrawvcfinputreferenceyear = atoi(namelistworkrawvcfinputreferenceyear);
  sprintf(workrawvcfinputreferenceyearstr,"%04d",workrawvcfinputreferenceyear);

  initializeGrids();

  sprintf(workrawclm5inputdirname,"%s/%s",workrawclm5inputdb,workrawclm5inputname);
  sprintf(workrawvcfinputdirname,"%s/%s",workrawvcfinputdb,workrawvcfinputname);
  readlandGrids(workrawclm5inputdirname,namelistworkrawclm5inputname,workrawclm5inputreferenceyearstr);
  readinitialforestvcfGrids(workrawvcfinputdirname,namelistworkrawvcfinputforestname,workrawvcfinputreferenceyearstr);    
  readinitialnonforestvcfGrids(workrawvcfinputdirname,namelistworkrawvcfinputnonforestname,workrawvcfinputreferenceyearstr);    
  readinitialrangevcfGrids(workrawvcfinputdirname,namelistworkrawvcfinputrangename,workrawvcfinputreferenceyearstr);    
  readclimvarseries(clm5currentrawinputdb,clm5currentrawinputname,"CLIM");
  readlaivarseries(clm5currentrawinputdb,clm5currentrawinputname,"FILLCLIM");
  
  genLUHLandGrids();
  genLUHCurrentVegGrids();
  
  sprintf(workrawvcfoutputreferenceyearstr,"%s",workrawvcfinputreferenceyearstr);
  sprintf(workrawvcfoutputdirname,"%s/%s",workrawvcfoutputdb,workrawvcfoutputname);
  writecombinedforestgrids(workrawvcfoutputdirname,namelistworkrawvcfoutputforestname,workrawvcfoutputreferenceyearstr);
  writecombinednonforestgrids(workrawvcfoutputdirname,namelistworkrawvcfoutputnonforestname,workrawvcfoutputreferenceyearstr);
  writecombinedrangegrids(workrawvcfoutputdirname,namelistworkrawvcfoutputrangename,workrawvcfoutputreferenceyearstr);

  return 1;
  
}
