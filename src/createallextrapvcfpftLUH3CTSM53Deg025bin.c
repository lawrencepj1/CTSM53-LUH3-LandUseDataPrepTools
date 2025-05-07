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

#define TREETYPE 1
#define HERBTYPE 2
#define GRASSTYPE 3

#define BAREPFT 0

#define FIRSTTREEPFT 1
#define LASTTREEPFT 8

#define FIRSTHERBPFT 9
#define LASTHERBPFT 14

#define FIRSTGRASSPFT 12
#define LASTGRASSPFT 14

#define PRECIPCLASSES 100
#define PRECIPCLASSSTART 0.0
#define PRECIPCLASSSIZE 50.0

#define TEMPCLASSES 100
#define TEMPCLASSSTART -15.0
#define TEMPCLASSSIZE 0.5

#define CLMPCTThreshold 1.0
#define LUHPCTThreshold 0.1

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
char namelistcurrentextrapvcfpftname[1024];
char namelistworkrawextrapoutput[1024];
char namelistworkrawextrapoutputtreename[1024];
char namelistworkrawextrapoutputherbname[1024];
char namelistworkrawextrapoutputgrassname[1024];
char namelistworkrawclm5input[1024];
char namelistworkrawclm5inputname[1024];
char namelistworkrawclm5inputreferenceyear[1024];
char namelistworkrawvcfinput[1024];
char namelistworkrawvcfinputforestname[1024];
char namelistworkrawvcfinputnonforestname[1024];
char namelistworkrawvcfinputrangename[1024];
char namelistworkrawvcfinputreferenceyear[1024];
char namelistclm5modisrawinput[1024];
char namelistclm5modisrawinputreferenceyear[1024];

int clm5modisrawinputreferenceyear;
int workrawclm5inputreferenceyear;
int workrawvcfinputreferenceyear;

char clm5modisrawinputname[1024];
char clm5modisrawinputdb[1024];
int clm5modisrawinputstartyear;
int clm5modisrawinputendyear;

char workrawclm5inputname[1024];
char workrawclm5inputdb[1024];
int workrawclm5inputstartyear;
int workrawclm5inputendyear;

char workrawvcfinputname[1024];
char workrawvcfinputdb[1024];
int workrawvcfinputstartyear;
int workrawvcfinputendyear;

char workrawextrapoutputname[1024];
char workrawextrapoutputdb[1024];
int workrawextrapoutputstartyear;
int workrawextrapoutputendyear;

char PFTluhtype[MAXPFT][256];
char CFTRAWluhtype[MAXCFTRAW][256];
char CFTluhtype[MAXCFT][256];

float tempLUHGrid[MAXLUHPIX * MAXLUHLIN];
float inLUHLANDMASKGrid[MAXLUHPIX * MAXLUHLIN];

float tempGrid[MAXCLMPIX * MAXCLMLIN];

float inLANDMASKGrid[MAXCLMPIX * MAXCLMLIN];
float inLANDFRACGrid[MAXCLMPIX * MAXCLMLIN];
float inAREAGrid[MAXCLMPIX * MAXCLMLIN];
float inclm5PCTGLACIERGrid[MAXCLMPIX * MAXCLMLIN];
float inclm5PCTLAKEGrid[MAXCLMPIX * MAXCLMLIN];
float inclm5PCTWETLANDGrid[MAXCLMPIX * MAXCLMLIN];
float inclm5PCTURBANGrid[MAXCLMPIX * MAXCLMLIN];
float inclm5PCTNATVEGGrid[MAXCLMPIX * MAXCLMLIN];
float inclm5PCTCROPGrid[MAXCLMPIX * MAXCLMLIN];
float inclm5PCTPFTGrid[MAXPFT][MAXCLMPIX * MAXCLMLIN];
float inclm5PCTCFTGrid[MAXCFT][MAXCLMPIX * MAXCLMLIN];

float inclm5PCTTREEGrid[MAXCLMPIX * MAXCLMLIN];
float inclm5PCTHERBGrid[MAXCLMPIX * MAXCLMLIN];
float inclm5PCTGRASSGrid[MAXCLMPIX * MAXCLMLIN];
float inclm5PCTBAREGrid[MAXCLMPIX * MAXCLMLIN];

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

float searchPOTVEGPCTPFTValueSum[MAXPFT];
float searchPOTVEGPCTPFTAreaSum[MAXPFT];

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
  fscanf(namelistfile,"%s %s",fieldname,namelistcurrentextrapvcfpftname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawextrapoutput);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawextrapoutputtreename);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawextrapoutputherbname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawextrapoutputgrassname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5input);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5inputname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5inputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfinputforestname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfinputnonforestname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfinputrangename);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawvcfinputreferenceyear);
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


int readworkrawextrapoutputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *workrawextrapoutputvaluefile;
  int foundworkrawextrapoutput, itemsfound;
  char fullfilenamestr[1024], inworkrawextrapoutputname[1024], inworkrawextrapoutputdb[1024];
  char inworkrawextrapoutputstartyear[1024], inworkrawextrapoutputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  workrawextrapoutputvaluefile = fopen(fullfilenamestr,"r");
  
  foundworkrawextrapoutput = 0;

  while (foundworkrawextrapoutput == 0) {
      itemsfound = fscanf(workrawextrapoutputvaluefile,"%s %s %s %s",inworkrawextrapoutputname, inworkrawextrapoutputdb, inworkrawextrapoutputstartyear, inworkrawextrapoutputendyear);
      if (itemsfound != 4) {
          printf("Error: workrawextrapoutput %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inworkrawextrapoutputname) == 0) {
          printf("Processing workrawextrapoutput: %s\n",inworkrawextrapoutputname);
	  sprintf(workrawextrapoutputname,"%s",inworkrawextrapoutputname);
          if (strcmp(inworkrawextrapoutputdb,"<workrawdir>") == 0) {
              sprintf(workrawextrapoutputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(workrawextrapoutputdb,"%s",inworkrawextrapoutputdb);
	  }
          workrawextrapoutputstartyear = atoi(inworkrawextrapoutputstartyear);
          workrawextrapoutputendyear = atoi(inworkrawextrapoutputendyear);
	  foundworkrawextrapoutput = 1;
      }
  }  
  
  return 0;

}


int initializeGrids() {

  long clm5lin, clm5pix;
  long pftid, cftid;
  
  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
          inclm5PCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          inclm5PCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
	  inclm5PCTGRASSGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          inclm5PCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
	  
          outPOTVEGLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGAREAGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGPCTGLACIERGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGPCTLAKEGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGPCTWETLANDGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGPCTURBANGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGPCTNATVEGGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outPOTVEGPCTCROPGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;

      }
  }

  return 0;
  	      
}


int initializeextrapveggrids() {

  long clm5lin, clm5pix;
  long pftid, cftid;
  
  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
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

  for (pftid = 0; pftid <= MAXPFT; pftid++) {  
      searchPOTVEGPCTPFTValueSum[pftid] = 0.0;
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
  fread(inclm5PCTGLACIERGrid,sizeof(inclm5PCTGLACIERGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTLAKE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inclm5PCTLAKEGrid,sizeof(inclm5PCTLAKEGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTWETLAND.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inclm5PCTWETLANDGrid,sizeof(inclm5PCTWETLANDGrid),1,globalbinfile);  
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
  fread(inclm5PCTURBANGrid,sizeof(inclm5PCTURBANGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTNATVEG.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inclm5PCTNATVEGGrid,sizeof(inclm5PCTNATVEGGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTCROP.%s.dat",databasestr,seriesname,seriesname,referenceyearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inclm5PCTCROPGrid,sizeof(inclm5PCTCROPGrid),1,globalbinfile);  
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
	      inclm5PCTPFTGrid[pftid][clm5lin * MAXCLMPIX + clm5pix] = tempGrid[clm5lin * MAXCLMPIX + clm5pix];
	      if (pftid >= FIRSTTREEPFT && pftid <= LASTTREEPFT) {
	          inclm5PCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
              }
	      if (pftid >= FIRSTHERBPFT && pftid <= LASTHERBPFT) {
	          inclm5PCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
              }
	      if (pftid >= FIRSTGRASSPFT && pftid <= LASTGRASSPFT) {
	          inclm5PCTGRASSGrid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
              }
	      if (pftid == BAREPFT) {
	          inclm5PCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
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
	      inclm5PCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = tempGrid[clm5lin * MAXCLMPIX + clm5pix];
	  }
      }
  }

  return 0;
  
}


int readpotvegforestvcfGrids(char *databasestr,char *seriesname,char *referenceyearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long clm5lin, clm5pix;

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


int readpotvegnonforestvcfGrids(char *databasestr,char *seriesname,char *referenceyearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long clm5lin, clm5pix;

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


int readpotvegrangevcfGrids(char *databasestr,char *seriesname,char *referenceyearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long clm5lin, clm5pix;

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


int genCLM5LandGrids() {

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
	      glacierpct = inclm5PCTGLACIERGrid[clm5lin * MAXCLMPIX + clm5pix];
	      lakepct = inclm5PCTLAKEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      wetlandpct = inclm5PCTWETLANDGrid[clm5lin * MAXCLMPIX + clm5pix];
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


float calcradiusavg(int searchtype, long centerlin, long centerpix, float searchinradius, float searchoutradius, float searchminval, float searchpotvegval) {

  long searchlin, searchpix;
  long searchoutradiuslength;
  float lindistance, pixdistance, searchdistance;
  float searchareasum;
  float landfracvalue, areavalue, searchvalue, pftvalue;
  int searchstartpft, searchendpft, pftindex;
  
  switch(searchtype) {
      case TREETYPE:
          searchstartpft = FIRSTTREEPFT;
          searchendpft = LASTTREEPFT;
          break;
      case HERBTYPE: 
          searchstartpft = FIRSTHERBPFT;
          searchendpft = LASTHERBPFT;
          break;
      case GRASSTYPE: 
          searchstartpft = FIRSTGRASSPFT;
          searchendpft = LASTGRASSPFT;
          break;
  }
  
  searchareasum = 0.0;

/* This an error and the 8 should be replaced with PFTMAX but it needs to be consistent for existing model runs */

  for (pftindex = 0; pftindex <= 8; pftindex++) {  
      searchPOTVEGPCTPFTValueSum[pftindex] = 0.0;
  }

  searchoutradiuslength = (long) (searchoutradius / CLMPIXSIZE);
  for (searchlin = centerlin - searchoutradiuslength; searchlin <= centerlin + searchoutradiuslength; searchlin++) {
      if (searchlin >= 0 && searchlin < MAXCLMLIN) {
          for (searchpix = centerpix - searchoutradiuslength; searchpix <= centerpix + searchoutradiuslength; searchpix++) {
	      if (searchpix >= 0.0 && searchpix < MAXCLMPIX) {
	          lindistance = ((float) (centerlin - searchlin)) * CLMPIXSIZE;
		  pixdistance = ((float) (centerpix - searchpix)) * CLMPIXSIZE;
                  searchdistance = pow(pow(lindistance,2) + pow(pixdistance,2),0.5);
	          if (searchdistance >= searchinradius && searchdistance <= searchoutradius) {
	              landfracvalue = inLANDFRACGrid[searchlin * MAXCLMPIX + searchpix];
		      if (landfracvalue > 0.0) {
                          switch(searchtype) {
			      case TREETYPE: 
                                  searchvalue = inclm5PCTTREEGrid[searchlin * MAXCLMPIX + searchpix];
                                  break;
			      case HERBTYPE: 
                                  searchvalue = inclm5PCTHERBGrid[searchlin * MAXCLMPIX + searchpix];
                                  break;
			      case GRASSTYPE: 
                                  searchvalue = inclm5PCTGRASSGrid[searchlin * MAXCLMPIX + searchpix];
                                  break;
                          }
			  if (searchvalue >= searchminval && searchvalue <= searchpotvegval) {
                              areavalue = landfracvalue * inAREAGrid[searchlin * MAXCLMPIX + searchpix];
			      for (pftindex = searchstartpft; pftindex <= searchendpft; pftindex++) {
			          pftvalue = inclm5PCTPFTGrid[pftindex][searchlin * MAXCLMPIX + searchpix];
			          searchPOTVEGPCTPFTValueSum[pftindex] += areavalue * pftvalue;
                              }
	                      searchareasum += areavalue;
			  }
                      }
		  }
              }
          }
      }
  }

  if (searchareasum > 0.0) {
      for (pftindex = searchstartpft; pftindex <= searchendpft; pftindex++) {
          searchPOTVEGPCTPFTValueSum[pftindex] = searchPOTVEGPCTPFTValueSum[pftindex] / searchareasum;
      }
  }

  return searchareasum;
			  
}


int genCLM5ExtrapVegGrids(int potvegtype) {

  long clm5lin, clm5pix;
  long searchstartpft, searchendpft, pftid, cftid;
  float landmask, landfrac, area, precipann, tempaverage;
  float potvegtreepct, potvegherbpct, potvegbarepct;
  float currentsearchpct, potvegpct;
  float currenttreepct, currentherbpct, currentgrasspct, currentbarepct;
  float climatemask;
  float newpotveg, availpotveg, removepotveg, addpotveg;
  float areafound, searchinradius, searchoutradius;
  int precipannindex, tempaverageindex;

  switch(potvegtype) {
      case TREETYPE:
          searchstartpft = FIRSTTREEPFT;
          searchendpft = LASTTREEPFT;
          break;
      case HERBTYPE: 
          searchstartpft = FIRSTHERBPFT;
          searchendpft = LASTHERBPFT;
          break;
      case GRASSTYPE: 
          searchstartpft = FIRSTGRASSPFT;
          searchendpft = LASTGRASSPFT;
          break;
  }
  
  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      if ((clm5lin - (clm5lin / 10) * 10) == 0) {
          printf("Calculating Tree Extrapolation line %d\n",clm5lin);
      }
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
          landmask = inLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix];
          landfrac = inLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix];
	  if (landfrac > 0.0) {

	      currenttreepct = inclm5PCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      currentherbpct = inclm5PCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
	      currentgrasspct = inclm5PCTGRASSGrid[clm5lin * MAXCLMPIX + clm5pix];
	      potvegtreepct = inPOTVEGFORESTPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      potvegtreepct += inPOTVEGNONFORESTPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      potvegtreepct += inPOTVEGRANGEPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      potvegherbpct = inPOTVEGFORESTPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
	      potvegherbpct += inPOTVEGNONFORESTPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
	      potvegherbpct += inPOTVEGRANGEPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];

              switch(potvegtype) {
                  case TREETYPE:
                      currentsearchpct = currenttreepct;
                      potvegpct = potvegtreepct;
                      break;
                  case HERBTYPE: 
                      currentsearchpct = currentherbpct;
                      potvegpct = potvegtreepct + potvegherbpct;
                      break;
                  case GRASSTYPE: 
                      currentsearchpct = currentherbpct;
                      potvegpct = potvegtreepct + potvegherbpct;
                      break;
              }
	      
	      if (currentsearchpct >= CLMPCTThreshold || potvegpct == 0.0) {
	          for (pftid = searchstartpft; pftid <= searchendpft; pftid++) {
		      outPOTVEGPCTPFTGrid[pftid][clm5lin * MAXCLMPIX + clm5pix] = inclm5PCTPFTGrid[pftid][clm5lin * MAXCLMPIX + clm5pix];
		  }
              }
	      else {
	          areafound = 0.0;
		  searchinradius = 0.0;
		  searchoutradius = 0.25;
		  while (areafound == 0.0) {
		      areafound = calcradiusavg(potvegtype,clm5lin,clm5pix,searchinradius,searchoutradius,CLMPCTThreshold,100.0);
		      if (areafound > 0.0) {
	                  for (pftid = searchstartpft; pftid <= searchendpft; pftid++) {
		              outPOTVEGPCTPFTGrid[pftid][clm5lin * MAXCLMPIX + clm5pix] = searchPOTVEGPCTPFTValueSum[pftid];
			  }
		      }
		      else {
		          searchinradius = searchoutradius;
			  searchoutradius = searchoutradius + 0.25;
		      }
		  }
             }          
		      
              outPOTVEGPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = currenttreepct;
              outPOTVEGPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] = currentherbpct;
              outPOTVEGPCTGRASSGrid[clm5lin * MAXCLMPIX + clm5pix] = currentgrasspct;
	      outPOTVEGPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix] = currentbarepct;
	      
	  }
      }
  }

  return 0;
  	      
}


int writeextrapveggrids(char *databasestr, char *seriesname, char *yearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  long outlin, outpix;
  long clm5lin, clm5pix;
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
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {

      for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
          for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
	      tempGrid[clm5lin * MAXCLMPIX + clm5pix] = outPOTVEGPCTPFTGrid[pftid][clm5lin * MAXCLMPIX + clm5pix];
	  }
      }
      
      sprintf(pftidstr,"%02d",pftid);
      sprintf(globalbinfilename,"%s/%s/%s.PCTNATPFT%s.%s.dat",databasestr,seriesname,seriesname,pftidstr,yearname);
      printf("Writing: %s\n",globalbinfilename);
      globalbinfile = fopen(globalbinfilename,"w+");
      fwrite(tempGrid,sizeof(tempGrid),1,globalbinfile);  
      fclose(globalbinfile);
  }

  for (cftid = 0; cftid < MAXCFT; cftid++) {

      for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
          for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
	      tempGrid[clm5lin * MAXCLMPIX + clm5pix] = outPOTVEGPCTPFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix];
	  }
      }
      
      sprintf(cftidstr,"%02d",cftid);
      sprintf(globalbinfilename,"%s/%s/%s.PCTCFT%s.%s.dat",databasestr,seriesname,seriesname,cftidstr,yearname);
      printf("Writing: %s\n",globalbinfilename);
      globalbinfile = fopen(globalbinfilename,"w+");
      fwrite(tempGrid,sizeof(tempGrid),1,globalbinfile);  
      fclose(globalbinfile);
  }

  return 0;

}


int main(long narg, char **argv) {

  char workrawclm5inputdirname[1024];
  char workrawvcfinputdirname[1024];
  char workrawextrapoutputdirname[1024];
  char clm5modisrawinputreferenceyearstr[256];
  char workrawclm5inputreferenceyearstr[256];
  char workrawvcfinputreferenceyearstr[256];
  char workrawextrapoutputreferenceyearstr[256];
  
  if(narg != 6){
        printf("Usage createpotvegpftsclm552Deg025bin timeseriesrawfile timeseriesrawdir workrawfile workrawdir generatecurrentvcfnamelist\n");
        return 0;
  }
  
  readnamelistfile(argv[5]);

  readclm5modisrawinputlutfile(argv[1],argv[2],namelistclm5modisrawinput);
  readworkrawclm5inputlutfile(argv[3],argv[4],namelistworkrawclm5input);
  readworkrawvcfinputlutfile(argv[3],argv[4],namelistworkrawvcfinput);
  readworkrawextrapoutputlutfile(argv[3],argv[4],namelistworkrawextrapoutput);
  clm5modisrawinputreferenceyear = atoi(namelistclm5modisrawinputreferenceyear);
  sprintf(clm5modisrawinputreferenceyearstr,"%04d",clm5modisrawinputreferenceyear);
  workrawclm5inputreferenceyear = atoi(namelistworkrawclm5inputreferenceyear);
  sprintf(workrawclm5inputreferenceyearstr,"%04d",workrawclm5inputreferenceyear);
  workrawvcfinputreferenceyear = atoi(namelistworkrawvcfinputreferenceyear);
  sprintf(workrawvcfinputreferenceyearstr,"%04d",workrawvcfinputreferenceyear);
    
  initializeGrids();
  
  sprintf(workrawclm5inputdirname,"%s/%s",workrawclm5inputdb,workrawclm5inputname);
  sprintf(workrawvcfinputdirname,"%s/%s",workrawvcfinputdb,workrawvcfinputname);
  readlandGrids(workrawclm5inputdirname,namelistworkrawclm5inputname,workrawclm5inputreferenceyearstr);
  readclm5referenceGrids(workrawclm5inputdirname,namelistworkrawclm5inputname,workrawclm5inputreferenceyearstr);    
  readpotvegforestvcfGrids(workrawvcfinputdirname,namelistworkrawvcfinputforestname,workrawvcfinputreferenceyearstr);    
  readpotvegnonforestvcfGrids(workrawvcfinputdirname,namelistworkrawvcfinputnonforestname,workrawvcfinputreferenceyearstr);    
  readpotvegrangevcfGrids(workrawvcfinputdirname,namelistworkrawvcfinputrangename,workrawvcfinputreferenceyearstr);    
  readclimvarseries(clm5modisrawinputdb,clm5modisrawinputname,"CLIM");
  readlaivarseries(clm5modisrawinputdb,clm5modisrawinputname,"FILLCLIM");
  
  genCLM5LandGrids();

  sprintf(workrawextrapoutputreferenceyearstr,"%s",workrawvcfinputreferenceyearstr);
  sprintf(workrawextrapoutputdirname,"%s/%s",workrawextrapoutputdb,workrawextrapoutputname);

  initializeextrapveggrids();
  genCLM5ExtrapVegGrids(TREETYPE); 
  writeextrapveggrids(workrawextrapoutputdirname,namelistworkrawextrapoutputtreename,workrawextrapoutputreferenceyearstr);

  initializeextrapveggrids();
  genCLM5ExtrapVegGrids(HERBTYPE); 
  writeextrapveggrids(workrawextrapoutputdirname,namelistworkrawextrapoutputherbname,workrawextrapoutputreferenceyearstr);

  initializeextrapveggrids();
  genCLM5ExtrapVegGrids(GRASSTYPE); 
  writeextrapveggrids(workrawextrapoutputdirname,namelistworkrawextrapoutputgrassname,workrawextrapoutputreferenceyearstr);

  return 1;
  
}
