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
#define MAXCFTLUH3 32
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
char namelistluh3descallc4annallcftname[1024];
char namelistluh3descrawcftoutput[1024];
char namelistluh3descrawcftoutputc4annname[1024];
char namelistworkrawclm5input[1024];
char namelistworkrawclm5inputname[1024];
char namelistworkrawclm5inputreferenceyear[1024];
char namelistclm5modisrawinput[1024];
char namelistclm5modisrawinputreferenceyear[1024];
char namelistclm5modisrawinputstartyear[1024];
char namelistclm5modisrawinputendyear[1024];

int clm5modisrawinputreferenceyear;
int clm5modisrawinputstartyear;
int clm5modisrawinputendyear;
int workrawclm5inputreferenceyear;

char clm5modisrawinputname[1024];
char clm5modisrawinputdb[1024];
int clm5modisrawinputstartyear;
int clm5modisrawinputendyear;

char workrawclm5inputname[1024];
char workrawclm5inputdb[1024];
int workrawclm5inputstartyear;
int workrawclm5inputendyear;

char luh3descrawcftoutputname[1024];
char luh3descrawcftoutputdb[1024];
int luh3descrawcftoutputstartyear;
int luh3descrawcftoutputendyear;

int CFTluhsearchid[MAXCFTLUH3];
char CFTluhtype[MAXCFTLUH3][256];
float CFTextrap0search[MAXCFTLUH3];
float CFTextrap1search[MAXCFTLUH3];
float CFTextrap2search[MAXCFTLUH3];
float CFTextrap3search[MAXCFTLUH3];
float CFTextrap0scaling[MAXCFTLUH3];
float CFTextrap1scaling[MAXCFTLUH3];
float CFTextrap2scaling[MAXCFTLUH3];
float CFTextrap3scaling[MAXCFTLUH3];

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
float inCTSMPCTPFTGrid[MAXCFTLUH3][MAXCLMPIX * MAXCLMLIN];

float inCTSMPCTCROPGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMPCTC4ANNCROPGrid[MAXCLMPIX * MAXCLMLIN];
float inCTSMPCTCFTGrid[MAXCFTLUH3][MAXCLMPIX * MAXCLMLIN];

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

float searchPOTVEGPCTCFTValueSum[MAXCFTLUH3];
float searchPOTVEGPCTCFTAreaSum[MAXCFTLUH3];

float outDESCLANDMASKGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCLANDFRACGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCAREAGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCPCTGLACIERGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCPCTLAKEGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCPCTWETLANDGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCPCTURBANGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCPCTNATVEGGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCPCTCROPGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCPCTC4ANNCROPGrid[MAXCLMPIX * MAXCLMLIN];
float outDESCPCTCFTGrid[MAXCFTLUH3][MAXCLMPIX * MAXCLMLIN];
float outDESCLVLCFTGrid[MAXCFTLUH3][MAXCLMPIX * MAXCLMLIN];

float outFINALPCTCFTGrid[MAXCFTLUH3][MAXCLMPIX * MAXCLMLIN];

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
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3descallc4annallcftname);
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3descrawcftoutput);
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3descrawcftoutputc4annname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5input);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5inputname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5inputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5modisrawinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5modisrawinputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5modisrawinputstartyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5modisrawinputendyear);
  
  return 0;

}


int readluh3descrawcftoutputlutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *luh3descrawcftoutputvaluefile;
  int foundluh3descrawcftoutput, itemsfound;
  char fullfilenamestr[1024], inluh3descrawcftoutputname[1024], inluh3descrawcftoutputdb[1024];
  char inluh3descrawcftoutputstartyear[1024], inluh3descrawcftoutputendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  luh3descrawcftoutputvaluefile = fopen(fullfilenamestr,"r");
  
  foundluh3descrawcftoutput = 0;

  while (foundluh3descrawcftoutput == 0) {
      itemsfound = fscanf(luh3descrawcftoutputvaluefile,"%s %s %s %s",inluh3descrawcftoutputname, inluh3descrawcftoutputdb, inluh3descrawcftoutputstartyear, inluh3descrawcftoutputendyear);
      if (itemsfound != 4) {
          printf("Error: luh3descrawcftoutput %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,inluh3descrawcftoutputname) == 0) {
          printf("Processing luh3descrawcftoutput: %s\n",inluh3descrawcftoutputname);
	  sprintf(luh3descrawcftoutputname,"%s",inluh3descrawcftoutputname);
          if (strcmp(inluh3descrawcftoutputdb,"<luh3descrawdir>") == 0) {
              sprintf(luh3descrawcftoutputdb,"%s",dirnamestr);
	  }
	  else {
              sprintf(luh3descrawcftoutputdb,"%s",inluh3descrawcftoutputdb);
	  }
          luh3descrawcftoutputstartyear = atoi(inluh3descrawcftoutputstartyear);
          luh3descrawcftoutputendyear = atoi(inluh3descrawcftoutputendyear);
	  foundluh3descrawcftoutput = 1;
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


int
readcftparamfile(char *filenamestr) {

  FILE *cftparamfile;
  int incft, incftid;
  int inCFTsearchid;
  float inCFTextrap0search,inCFTextrap1search,inCFTextrap2search,inCFTextrap3search;
  float inCFTextrap0scaling,inCFTextrap1scaling,inCFTextrap2scaling,inCFTextrap3scaling;
  char fullfilenamestr[256], inCFTluhtype[256], inCFTname[256];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  cftparamfile = fopen(fullfilenamestr,"r");

  for (incft = 0; incft < MAXCFTLUH3; incft++) {
      fscanf(cftparamfile,"%d %s %f %f %f %f %f %f %f %f  %s",&incftid,inCFTluhtype,&inCFTextrap0search,&inCFTextrap1search,&inCFTextrap2search,&inCFTextrap3search,&inCFTextrap0scaling,&inCFTextrap1scaling,&inCFTextrap2scaling,&inCFTextrap3scaling,inCFTname);
      CFTluhsearchid[incft] = incftid;
      sprintf(CFTluhtype[incft],"%s",inCFTluhtype);
      printf("Reading %d %s %s\n",incft,inCFTname,CFTluhtype[incft]);
      CFTextrap0search[incft] = inCFTextrap0search;
      CFTextrap1search[incft] = inCFTextrap1search;
      CFTextrap2search[incft] = inCFTextrap2search;
      CFTextrap3search[incft] = inCFTextrap3search;
      CFTextrap0scaling[incft] = inCFTextrap0scaling;
      CFTextrap1scaling[incft] = inCFTextrap1scaling;
      CFTextrap2scaling[incft] = inCFTextrap2scaling;
      CFTextrap3scaling[incft] = inCFTextrap3scaling;
  }  
  
  return 0;

}


int initializeGrids() {

  long clm5lin, clm5pix;
  long pftid, cftid;
  
  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {

          inCTSMPCTCROPGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          inCTSMPCTC4ANNCROPGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
              inCTSMPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          }

          outDESCLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCAREAGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCPCTGLACIERGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCPCTLAKEGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCPCTWETLANDGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCPCTURBANGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCPCTNATVEGGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCPCTCROPGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCPCTC4ANNCROPGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
              outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
              outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
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

          inCTSMPCTCROPGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          inCTSMPCTC4ANNCROPGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
              inCTSMPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          }

          outDESCPCTCROPGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          outDESCPCTC4ANNCROPGrid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
              outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
              outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
              outFINALPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
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
	  }
      }
  }
  
  return 0;
  
}


int readclm5cropGrids(char *databasestr,char *seriesname,char *referenceyearname,char *croptype) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid, importcftid;
  char pftidstr[256], cftidstr[256];
  long clm5lin, clm5pix;

  for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
      sprintf(cftidstr,"%02d",2*(cftid));
      sprintf(globalbinfilename,"%s/%s/%s.PCTCFT%s.%s.dat",databasestr,seriesname,seriesname,cftidstr,referenceyearname);
      printf("Reading: %s\n",globalbinfilename);
      importcftid = cftid;
      globalbinfile = fopen(globalbinfilename,"r");
      fread(tempGrid,sizeof(tempGrid),1,globalbinfile);  
      fclose(globalbinfile);
      for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
          for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
              inCTSMPCTCFTGrid[importcftid][clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
	      inCTSMPCTCROPGrid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
          }
      }
      sprintf(cftidstr,"%02d",2*(cftid)+1);
      sprintf(globalbinfilename,"%s/%s/%s.PCTCFT%s.%s.dat",databasestr,seriesname,seriesname,cftidstr,referenceyearname);
      printf("Reading: %s\n",globalbinfilename);
      globalbinfile = fopen(globalbinfilename,"r");
      fread(tempGrid,sizeof(tempGrid),1,globalbinfile);  
      fclose(globalbinfile);
      for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
          for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
              inCTSMPCTCFTGrid[importcftid][clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
	      inCTSMPCTCROPGrid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
	  }
      }
  }

  for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
      if (strcmp(CFTluhtype[cftid],croptype) == 0) {
          for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
              for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
	          inCTSMPCTC4ANNCROPGrid[clm5lin * MAXCLMPIX + clm5pix] += inCTSMPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix];
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


float calcradiusavg(long centerlin, long centerpix, float searchinradius, float searchoutradius, float searchminval, float searchpotvegval, char *croptype) {

  long searchlin, searchpix;
  long searchoutradiuslength;
  float lindistance, pixdistance, searchdistance;
  float searchareasum;
  float landfracvalue, areavalue, cropvalue, cftvalue;
  int cftindex;
  
  for (cftindex = 0; cftindex < MAXCFTLUH3; cftindex++) {
      searchPOTVEGPCTCFTValueSum[cftindex] = 0.0;
      searchPOTVEGPCTCFTAreaSum[cftindex] = 0.0;
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
                          cropvalue = inCTSMPCTC4ANNCROPGrid[searchlin * MAXCLMPIX + searchpix];
			  if (cropvalue >= searchminval && cropvalue <= searchpotvegval) {
                              areavalue = landfracvalue * inAREAGrid[searchlin * MAXCLMPIX + searchpix] * cropvalue / 100.0;
			      for (cftindex = 0; cftindex <= MAXCFTLUH3; cftindex++) {
                                  if (strcmp(CFTluhtype[cftindex],croptype) == 0) {
			               cftvalue = inCTSMPCTCFTGrid[cftindex][searchlin * MAXCLMPIX + searchpix];
			               searchPOTVEGPCTCFTValueSum[cftindex] += areavalue * cftvalue;
				       searchPOTVEGPCTCFTAreaSum[cftindex] += areavalue;
				  }
                              }
			  }
                      }
		  }
              }
          }
      }
  }

  searchareasum = 0.0;
  for (cftindex = 0; cftindex <= MAXCFTLUH3; cftindex++) {
      if (strcmp(CFTluhtype[cftindex],croptype) == 0) {
          if (searchPOTVEGPCTCFTAreaSum[cftindex] > 0.0) {
              searchPOTVEGPCTCFTValueSum[cftindex] = searchPOTVEGPCTCFTValueSum[cftindex] / searchPOTVEGPCTCFTAreaSum[cftindex];
	      searchareasum += searchPOTVEGPCTCFTAreaSum[cftindex];
	  }
      }
  }

  return searchareasum;
			  
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
              outDESCLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix] = landmask;
              outDESCLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix] = landfrac;
	      outDESCAREAGrid[clm5lin * MAXCLMPIX + clm5pix] = inAREAGrid[clm5lin * MAXCLMPIX + clm5pix];
	      glacierpct = inCTSMPCTGLACIERGrid[clm5lin * MAXCLMPIX + clm5pix];
	      lakepct = inCTSMPCTLAKEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      wetlandpct = inCTSMPCTWETLANDGrid[clm5lin * MAXCLMPIX + clm5pix];
	      natvegpct = 100.0 - glacierpct - lakepct - wetlandpct;
              outDESCPCTGLACIERGrid[clm5lin * MAXCLMPIX + clm5pix] = glacierpct;
              outDESCPCTLAKEGrid[clm5lin * MAXCLMPIX + clm5pix] = lakepct;
              outDESCPCTWETLANDGrid[clm5lin * MAXCLMPIX + clm5pix] = wetlandpct;
              outDESCPCTNATVEGGrid[clm5lin * MAXCLMPIX + clm5pix] = natvegpct;
	  }
      }
  }

  return 0;
  	      
}


int genLUHDescriptorCropGrids(char *croptype) {

  long clm5lin, clm5pix;
  long pftid, cftid;
  float landmask, landfrac, area, precipann, tempaverage;
  float currentnatvegpct, currentbarepct, currentcroppct, currentcroplvl;
  float allcftsearch0, allcftsearch1, allcftsearch2, allcftsearch3;
  float allcftscaling0, allcftscaling1, allcftscaling2, allcftscaling3;
  float cftsearch, cftscaling, searchcroppct;
  float scaledcroppct, scaledcroppctdiff;
  float areafound, searchinradius, searchoutradius;

  allcftsearch0 = 0.0;
  allcftsearch1 = 0.0;
  allcftsearch2 = 0.0;
  allcftsearch3 = 0.0;
  allcftscaling0 = 0.0;
  allcftscaling1 = 0.0;
  allcftscaling2 = 0.0;
  allcftscaling3 = 0.0;
  for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
      if (strcmp(CFTluhtype[cftid],croptype) == 0) {
          allcftsearch0 += CFTextrap0search[cftid];
          allcftsearch1 += CFTextrap1search[cftid];
          allcftsearch2 += CFTextrap2search[cftid];
          allcftsearch3 += CFTextrap3search[cftid];
          allcftscaling0 += CFTextrap0scaling[cftid];
          allcftscaling1 += CFTextrap1scaling[cftid];
          allcftscaling2 += CFTextrap2scaling[cftid];
          allcftscaling3 += CFTextrap3scaling[cftid];
      }
  }
  
  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      if ((clm5lin - (clm5lin / 10) * 10) == 0) {
          printf("Calculating Crop Extrapolation line %d\n",clm5lin);
      }
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
          landmask = inLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix];
          landfrac = inLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix];
	  if (landfrac > 0.0) {

              currentnatvegpct = inCTSMPCTNATVEGGrid[clm5lin * MAXCLMPIX + clm5pix];
	      currentbarepct = inCTSMPCTPFTGrid[0][clm5lin * MAXCLMPIX + clm5pix];
	      currentcroppct = inCTSMPCTC4ANNCROPGrid[clm5lin * MAXCLMPIX + clm5pix];
	      
	      if (clm5lin < 2 * MAXCLMLIN / 18 || clm5lin > 17 * MAXCLMLIN / 18) {
	          currentbarepct = 100.0;
              }
	      
	      if (currentcroppct >= CLMPCTThreshold || currentbarepct == 100.0) {
		  if (currentcroppct >= CLMPCTThreshold) {
		      currentcroplvl = 0.0;
		      searchcroppct = 0.0;
                      for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
                          if (strcmp(CFTluhtype[cftid],croptype) == 0) {
                              cftsearch = CFTextrap0search[cftid];
			      searchPOTVEGPCTCFTValueSum[cftid] = cftsearch * inCTSMPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix];
                              searchcroppct += searchPOTVEGPCTCFTValueSum[cftid];
                          }
		          else {
                              searchPOTVEGPCTCFTValueSum[cftid] = 0.0;
                          }
                      }
                      scaledcroppct = 0.0;
		      if (searchcroppct > 0.0) {
                          for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
                              if (strcmp(CFTluhtype[cftid],croptype) == 0) {
                                  cftscaling = CFTextrap0scaling[cftid];
                                  searchPOTVEGPCTCFTValueSum[cftid] = cftscaling * searchPOTVEGPCTCFTValueSum[cftid] / searchcroppct * 100.0;
				  if (searchPOTVEGPCTCFTValueSum[cftid] > 0.0) {
                                      outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = searchPOTVEGPCTCFTValueSum[cftid];
                                      outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = currentcroplvl;	
				  }
				  else {
                                      outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
                                      outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = -1.0;	
                                  }
                                  scaledcroppct += searchPOTVEGPCTCFTValueSum[cftid];
                              }
			      else {
                                  outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
                                  outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = -1.0;	
                              }
                          }
		          scaledcroppctdiff = 100.0 - scaledcroppct;
	                  if (scaledcroppctdiff < 0.0) {
                              scaledcroppctdiff = 0.0;
                          }
		          if (scaledcroppctdiff > 100.0 - searchPOTVEGPCTCFTValueSum[15]) {
			      scaledcroppctdiff = 100.0 - searchPOTVEGPCTCFTValueSum[15];
                          }
			  if (scaledcroppctdiff > 0.0) {
		              outDESCPCTCFTGrid[15][clm5lin * MAXCLMPIX + clm5pix] += scaledcroppctdiff; 
                              outDESCLVLCFTGrid[15][clm5lin * MAXCLMPIX + clm5pix] = currentcroplvl;
			  }	
		      }
		  }
		  else {
		      currentcroplvl = 0.0;
                      searchcroppct = 100.0;
                      scaledcroppct = 100.0;
                      for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
                          if (cftid == 15) {
                              outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 100.0;
                              outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = currentcroplvl;	
                          }
                          else {
                              outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
                              outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = -1.0;	
			  }
                      }
		  }
	      }
	      
	      else {
	          areafound = 0.0;
		  currentcroplvl = 1.0;
		  searchinradius = 0.0;
		  searchoutradius = 0.25;
		  while (areafound == 0.0) {
                      for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
                          searchPOTVEGPCTCFTValueSum[cftid] = 0.0;
                      }
		      if (searchoutradius > 0.0 && searchoutradius < 1.0) {
		          currentcroplvl = 1.0;
			  if (allcftsearch1 == 0.0) {
			      areafound = -1.0;
			  }
                      }
		      if (searchoutradius >= 1.0 && searchoutradius < 2.0) {
		          currentcroplvl = 2.0;
			  if (allcftsearch2 == 0.0) {
			      areafound = -1.0;
			  }
                      }
		      if (searchoutradius >= 2.0 && searchoutradius < 4.0) {
		          currentcroplvl = 3.0;
			  if (allcftsearch3 == 0.0) {
			      areafound = -1.0;
			  }
                      }
		      if (searchoutradius < 4.0 && areafound != -1.0) {
		          areafound = calcradiusavg(clm5lin,clm5pix,searchinradius,searchoutradius,CLMPCTThreshold,100.0,croptype);
                      }
		      else {
		          areafound = -1.0;
                      }
		      if (areafound >= CLMPCTThreshold) {
		          searchcroppct = 0.0;
	                  for (cftid = 0; cftid <= MAXCFTLUH3; cftid++) {
			      if (strcmp(CFTluhtype[cftid],croptype) == 0) {
			          cftsearch = 1.0;
			          if (currentcroplvl == 1.0) {
                                      cftsearch = CFTextrap1search[cftid];
				  }
			          if (currentcroplvl == 2.0) {
                                      cftsearch = CFTextrap2search[cftid];
				  }
				  if (currentcroplvl > 2.0) {
				      cftsearch = CFTextrap3search[cftid];
				  }
				  searchPOTVEGPCTCFTValueSum[cftid] = cftsearch * searchPOTVEGPCTCFTValueSum[cftid];
                                  searchcroppct += searchPOTVEGPCTCFTValueSum[cftid];
                              }
			  }
			  if (searchcroppct > 0.0) {
			      scaledcroppct = 0.0;
                              for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
                                  if (strcmp(CFTluhtype[cftid],croptype) == 0) {
                                      cftscaling = 1.0;
			              if (currentcroplvl == 1.0) {
                                          cftscaling = CFTextrap1scaling[cftid];
				      }
			              if (currentcroplvl == 2.0) {
                                          cftscaling = CFTextrap2scaling[cftid];
				      }
				      if (currentcroplvl > 2.0) {
				          cftscaling = CFTextrap3scaling[cftid];
				      }
				      searchPOTVEGPCTCFTValueSum[cftid] = cftscaling * searchPOTVEGPCTCFTValueSum[cftid] / searchcroppct * 100.0;
				      if (searchPOTVEGPCTCFTValueSum[cftid] > 0.0) {
                                          outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = searchPOTVEGPCTCFTValueSum[cftid];
                                          outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = currentcroplvl;
				      }
				      else {	
                                          outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
                                          outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = -1.0;	
			              }
                                      scaledcroppct += searchPOTVEGPCTCFTValueSum[cftid];
				  }
				  else {
                                      outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
                                      outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = -1.0;	
				  }
                              }
                              if (scaledcroppct > 0.0) {
			          scaledcroppctdiff = 100.0 - scaledcroppct;
			          if (scaledcroppctdiff < 0.0) {
			              scaledcroppctdiff = 0.0;
                                  }
			          if (scaledcroppctdiff > 100.0 - outDESCPCTCFTGrid[15][clm5lin * MAXCLMPIX + clm5pix]) {
			              scaledcroppctdiff = 100.0 - outDESCPCTCFTGrid[15][clm5lin * MAXCLMPIX + clm5pix];
                                  }
			          if (scaledcroppctdiff > 0.0) {
			              outDESCPCTCFTGrid[15][clm5lin * MAXCLMPIX + clm5pix] += scaledcroppctdiff; 
                                      outDESCLVLCFTGrid[15][clm5lin * MAXCLMPIX + clm5pix] = currentcroplvl;
				  }
			      }	
			      else {
			          areafound == 0.0;
			          searchcroppct = 0.0;
			      }
			  }
			  else {
			      areafound == 0.0;
			      scaledcroppct = 0.0;
			  }
		      }
		      if (areafound == 0.0) {
		          searchinradius = searchoutradius;
			  searchoutradius = searchoutradius + 0.25;
		      }
		  }
		  if (areafound == -1.0) {
                      for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
                          if (cftid == 15) {
                              outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 100.0;
                              outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = currentcroplvl;	
                          }
                          else {
                              outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
                              outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = -1.0;	
			  }
                      }
		  }   
	      }
	      scaledcroppct = 0.0;
              for (cftid = 0; cftid <= MAXCFTLUH3; cftid++) {
                  if (outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] != -1.0) {
                      scaledcroppct += outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix];
		  }
              }
	      if (scaledcroppct == 0.0) {
                  for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
                      if (cftid == 15) {
                          outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 100.0;
                          outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = currentcroplvl;	
                      }
                      else {
                          outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
                          outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = -1.0;	
                      }
                  }
                  scaledcroppct = 100.0;  
              }
	      outDESCPCTCROPGrid[clm5lin * MAXCLMPIX + clm5pix] = scaledcroppct;
	  }

      }
  }


  return 0;
  	      
}


int regroupDescriptorCropGrids(char *croptype) {

  long clm5lin, clm5pix;
  int cftid;
  char extraptype[256];
  float extraplevel;
  int resetextrap;
  float landfrac, vegfrac, cropfrac;
  float totalcftpct, cftpct, cftlevel, minextraplevel;
  float luhc4croppct, luhcftpct, luhcftlevel, luhcftpctsum, foddergrasspad;
  float firstcftpct, secondcftpct, thirdcftpct, fourthcftpct, largestcftpct, rescaledcftpct;
  int cftcount, firstcftid, secondcftid, thirdcftid, fourthcftid, largestcftid;
  float c3gen, wheat, barley, rye, cotton, rice, sugarcane, temperatecorn, temperatesoy, cassava, citrus, cocoa, coffee;
  float datepalm, millet, oilpalm, canola, sunflower, sugarbeet, groundnuts, pulses, grapes, foddergrass, sorghum, tropicalcorn, tropicalsoy;
  
  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
          landfrac = inLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix];
	  if (landfrac > 0.0) {
	  
              minextraplevel = 6.0;
              sprintf(extraptype,"%s",croptype);
	      for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
                  cftpct = outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix];
                  cftlevel = outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix];
                  if (cftpct > 0.0 && cftlevel > -1.0 && cftlevel < minextraplevel) {
                      minextraplevel = cftlevel;
                  }
              }
	      extraplevel = minextraplevel;

              if (extraplevel == 0.0) {
	          luhc4croppct = outDESCPCTCROPGrid[clm5lin * MAXCLMPIX + clm5pix];
                  totalcftpct = 0.0;
	          cftcount = 0;
	          for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
                      cftpct = outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix];
                      cftlevel = outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix];
                      if (cftpct > 0.0 && cftlevel == extraplevel) {
		          totalcftpct += cftpct;
		          cftcount++;
		      }
                  }
                  if (totalcftpct > 0.0) {
                      for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
		          cftpct = 0.0;
                          cftlevel = outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix];
                          if (cftlevel == extraplevel) {
                              cftpct = outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix];
			  }
                          outFINALPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = cftpct / totalcftpct * 100.0;
                      }
	          }
		  else {
		      printf("Error Empty Crop %d %d\n",clm5lin,clm5pix);
		  }
	      }
	      else {
	          cftcount = 0;
		  firstcftid = -1;
		  firstcftpct = 0.0;
		  secondcftid = -1;
		  secondcftpct = 0.0;
		  thirdcftid = -1;
		  thirdcftpct = 0.0;
		  fourthcftid = -1;
		  fourthcftpct = 0.0;
		  for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
                      cftpct = outDESCPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix];
                      cftlevel = outDESCLVLCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix];
	              if (cftpct > 0.0 && cftlevel == extraplevel) {
                          cftcount++;
		          if (cftpct > firstcftpct) {
		              fourthcftid = thirdcftid;
		              fourthcftpct = thirdcftpct;
		              thirdcftid = secondcftid;
		              thirdcftpct = secondcftpct;
		              secondcftid = firstcftid;
		              secondcftpct = firstcftpct;
		              firstcftid = cftid;
		              firstcftpct = cftpct;
		          }
		          else {
                              if (cftpct > secondcftpct) {
                                  fourthcftid = thirdcftid;
                                  fourthcftpct = thirdcftpct;
                                  thirdcftid = secondcftid;
                                  thirdcftpct = secondcftpct;
                                  secondcftid = cftid;
                                  secondcftpct = cftpct;
                              }
                              else {
		                  if (cftpct > thirdcftpct) {
		                      fourthcftid = thirdcftid;
                                      fourthcftpct = thirdcftpct;
                                      thirdcftid = cftid;
                                      thirdcftpct = cftpct;
                                  }
                                  else {
                                      if (cftpct > fourthcftpct) {
		                          fourthcftid = cftid;
		                          fourthcftpct = cftpct;
                                      }
                                  }
			      }
                          }
		      }
                  }
              }
              if (extraplevel == 1.0 || extraplevel == 2.0) {
                  totalcftpct = firstcftpct + secondcftpct + thirdcftpct + fourthcftpct;
                  firstcftpct = firstcftpct / totalcftpct * 100.0;
                  secondcftpct = secondcftpct / totalcftpct * 100.0;
                  thirdcftpct = thirdcftpct / totalcftpct * 100.0;
                  fourthcftpct = fourthcftpct / totalcftpct * 100.0;
                  for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
	              if (cftid == firstcftid) {
                          outFINALPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = firstcftpct;
                      }
		      if (cftid == secondcftid) {
		          outFINALPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = secondcftpct;
                      }
		      if (cftid == thirdcftid) {
                          outFINALPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = thirdcftpct;
                      }
		      if (cftid == fourthcftid) {
                          outFINALPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = fourthcftpct;
                      }
                  }
              }
              if (extraplevel == 3.0 || extraplevel == 4.0) {
                  totalcftpct = firstcftpct + secondcftpct + thirdcftpct;
                  firstcftpct = firstcftpct / totalcftpct * 100.0;
                  secondcftpct = secondcftpct / totalcftpct * 100.0;
                  thirdcftpct = thirdcftpct / totalcftpct * 100.0;
                  for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
	              if (cftid == firstcftid) {
                          outFINALPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = firstcftpct;
                      }
		      if (cftid == secondcftid) {
                          outFINALPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = secondcftpct;
                      }
		      if (cftid == thirdcftid) {
                          outFINALPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = thirdcftpct;
                      }
                  }
              }
              if (extraplevel == 5.0) {
                  totalcftpct = firstcftpct + secondcftpct;
                  firstcftpct = firstcftpct / totalcftpct * 100.0;
                  secondcftpct = secondcftpct / totalcftpct * 100.0;
                  for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
	              if (cftid == firstcftid) {
                          outFINALPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = firstcftpct;
                      }
		      if (cftid == secondcftid) {
                          outFINALPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = secondcftpct;
                      }
                  }
              }
              if (extraplevel == -1.0 || extraplevel == 6.0) {
                  for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
                      outFINALPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
                  }
                  outFINALPCTCFTGrid[15][clm5lin * MAXCLMPIX + clm5pix] = 100.0;
              }
		      
              if (clm5lin < 60 * MAXCLMLIN / 180 || clm5lin >= 120 * MAXCLMLIN / 180) {
	          temperatecorn = outFINALPCTCFTGrid[1][clm5lin * MAXCLMPIX + clm5pix];
	          temperatesoy = outFINALPCTCFTGrid[4][clm5lin * MAXCLMPIX + clm5pix];
	          tropicalcorn = outFINALPCTCFTGrid[30][clm5lin * MAXCLMPIX + clm5pix];
	          tropicalsoy = outFINALPCTCFTGrid[31][clm5lin * MAXCLMPIX + clm5pix];
		  outFINALPCTCFTGrid[1][clm5lin * MAXCLMPIX + clm5pix] = temperatecorn + tropicalcorn;
		  outFINALPCTCFTGrid[4][clm5lin * MAXCLMPIX + clm5pix] = temperatesoy + tropicalsoy;
		  outFINALPCTCFTGrid[30][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
		  outFINALPCTCFTGrid[31][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
              }
	      if (clm5lin >= 60 * MAXCLMLIN / 180 && clm5lin < 120 * MAXCLMLIN / 180) {
		  temperatecorn = outFINALPCTCFTGrid[1][clm5lin * MAXCLMPIX + clm5pix];
		  temperatesoy = outFINALPCTCFTGrid[4][clm5lin * MAXCLMPIX + clm5pix];
		  tropicalcorn = outFINALPCTCFTGrid[30][clm5lin * MAXCLMPIX + clm5pix];
		  tropicalsoy = outFINALPCTCFTGrid[31][clm5lin * MAXCLMPIX + clm5pix];
		  outFINALPCTCFTGrid[1][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
		  outFINALPCTCFTGrid[4][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
		  outFINALPCTCFTGrid[30][clm5lin * MAXCLMPIX + clm5pix] = temperatecorn + tropicalcorn;
		  outFINALPCTCFTGrid[31][clm5lin * MAXCLMPIX + clm5pix] = temperatesoy + tropicalsoy;
	      }

	      if (clm5lin < 35 * MAXCLMLIN / 180 || clm5lin >= 145 * MAXCLMLIN / 180) {
		  temperatecorn = outFINALPCTCFTGrid[1][clm5lin * MAXCLMPIX + clm5pix];
		  wheat = outFINALPCTCFTGrid[2][clm5lin * MAXCLMPIX + clm5pix];
		  cassava = outFINALPCTCFTGrid[9][clm5lin * MAXCLMPIX + clm5pix];
		  citrus = outFINALPCTCFTGrid[10][clm5lin * MAXCLMPIX + clm5pix];
		  cocoa = outFINALPCTCFTGrid[11][clm5lin * MAXCLMPIX + clm5pix];
		  coffee = outFINALPCTCFTGrid[12][clm5lin * MAXCLMPIX + clm5pix];
                  datepalm = outFINALPCTCFTGrid[14][clm5lin * MAXCLMPIX + clm5pix];
                  foddergrass = outFINALPCTCFTGrid[15][clm5lin * MAXCLMPIX + clm5pix];
		  millet = outFINALPCTCFTGrid[18][clm5lin * MAXCLMPIX + clm5pix];
		  oilpalm = outFINALPCTCFTGrid[19][clm5lin * MAXCLMPIX + clm5pix];
		  rice = outFINALPCTCFTGrid[23][clm5lin * MAXCLMPIX + clm5pix];
		  sorghum = outFINALPCTCFTGrid[24][clm5lin * MAXCLMPIX + clm5pix];
		  sugarcane = outFINALPCTCFTGrid[26][clm5lin * MAXCLMPIX + clm5pix];
		  outFINALPCTCFTGrid[1][clm5lin * MAXCLMPIX + clm5pix] = temperatecorn + millet + sorghum + sugarcane;
		  outFINALPCTCFTGrid[2][clm5lin * MAXCLMPIX + clm5pix] = wheat + rice;
		  outFINALPCTCFTGrid[9][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
		  outFINALPCTCFTGrid[10][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
		  outFINALPCTCFTGrid[11][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
		  outFINALPCTCFTGrid[12][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
		  outFINALPCTCFTGrid[14][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
		  outFINALPCTCFTGrid[15][clm5lin * MAXCLMPIX + clm5pix] = foddergrass + citrus + cassava + cocoa + coffee + oilpalm + datepalm;
		  outFINALPCTCFTGrid[18][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
		  outFINALPCTCFTGrid[19][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
		  outFINALPCTCFTGrid[23][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
	          outFINALPCTCFTGrid[24][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
	          outFINALPCTCFTGrid[26][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
              }
	      
	  }

	  if (landfrac > 0.0) {
              totalcftpct = 0.0;
              for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
                  totalcftpct += outFINALPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix];
              }
	      if (totalcftpct > 0.0) {
	          largestcftpct = 0.0;
	          largestcftid = 0;
	          rescaledcftpct = 0.0;
	          for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {
                      cftpct = outFINALPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] * 100.0 / totalcftpct;
	              cftpct = trunc(cftpct * CLMPCTTruncate) / CLMPCTTruncate;
		      outFINALPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = cftpct;
		      rescaledcftpct += cftpct;
		      if (cftpct > largestcftpct) {
		          largestcftpct = cftpct;
		          largestcftid = cftid;
                      }
	          }
 	          outFINALPCTCFTGrid[largestcftid][clm5lin * MAXCLMPIX + clm5pix] += (100.0 -  rescaledcftpct);
	      }
	  }
      }
  }

  return 0;
  
}


int writedescgrids(char *databasestr, char *seriesname, char *yearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  long outlin, outpix;
  long clm5lin, clm5pix;
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];

  for (cftid = 0; cftid < MAXCFTLUH3; cftid++) {

      for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
          for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
	      tempGrid[clm5lin * MAXCLMPIX + clm5pix] = outFINALPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix];
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
  char luh3descrawcftoutputdirname[1024];
  char clm5modisrawinputreferenceyearstr[256];
  char workrawclm5inputreferenceyearstr[256];
  char luh3descrawcftoutputreferenceyearstr[256];
  int currentyear;
  
  if(narg != 9){
        printf("Usage createalldescc4annCTSM53Deg025bin luh3descrawfile luh3descrawdir timeseriesrawfile timeseriesrawdir workrawfile workrawdir luh3desccropparams luh3descc4annpftnamelistfile\n");
        return 0;
  }
  
  readnamelistfile(argv[8]);

  readluh3descrawcftoutputlutfile(argv[1],argv[2],namelistluh3descrawcftoutput);
  readclm5modisrawinputlutfile(argv[3],argv[4],namelistclm5modisrawinput);
  readworkrawclm5inputlutfile(argv[5],argv[6],namelistworkrawclm5input);
  readcftparamfile(argv[7]);
  workrawclm5inputreferenceyear = atoi(namelistworkrawclm5inputreferenceyear);
  sprintf(workrawclm5inputreferenceyearstr,"%04d",workrawclm5inputreferenceyear);
  clm5modisrawinputreferenceyear = atoi(namelistclm5modisrawinputreferenceyear);
  sprintf(clm5modisrawinputreferenceyearstr,"%04d",clm5modisrawinputreferenceyear);
  clm5modisrawinputstartyear = atoi(namelistclm5modisrawinputstartyear);
  clm5modisrawinputendyear = atoi(namelistclm5modisrawinputendyear);

  initializeGrids();

  sprintf(workrawclm5inputdirname,"%s/%s",workrawclm5inputdb,workrawclm5inputname);
  readlandGrids(workrawclm5inputdirname,namelistworkrawclm5inputname,workrawclm5inputreferenceyearstr);
  readclm5referenceGrids(workrawclm5inputdirname,namelistworkrawclm5inputname,workrawclm5inputreferenceyearstr);    
  readclimvarseries(clm5modisrawinputdb,clm5modisrawinputname,"CLIM");
  readlaivarseries(clm5modisrawinputdb,clm5modisrawinputname,"FILLCLIM");
  
  genLUHLandGrids();

  sprintf(luh3descrawcftoutputdirname,"%s/%s",luh3descrawcftoutputdb,luh3descrawcftoutputname);
  for (currentyear = clm5modisrawinputstartyear; currentyear <= clm5modisrawinputendyear; currentyear++) {
      sprintf(clm5modisrawinputreferenceyearstr,"%04d",currentyear);
      sprintf(luh3descrawcftoutputreferenceyearstr,"%04d",currentyear);  
      clearGrids();
      readclm5cropGrids(clm5modisrawinputdb,clm5modisrawinputname,clm5modisrawinputreferenceyearstr,"C4ANN");
      genLUHDescriptorCropGrids("C4ANN");
      regroupDescriptorCropGrids("C4ANN");
      writedescgrids(luh3descrawcftoutputdirname,namelistluh3descrawcftoutputc4annname,luh3descrawcftoutputreferenceyearstr);
  }

  return 1;
  
}
