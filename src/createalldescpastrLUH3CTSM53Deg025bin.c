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
char namelistluh3descallpastrallpftname[1024];
char namelistluh3descrawpftoutput[1024];
char namelistluh3descrawpftoutputpastrname[1024];
char namelistworkrawclm5input[1024];
char namelistworkrawclm5inputname[1024];
char namelistworkrawclm5inputreferenceyear[1024];
char namelistworkrawextrappftinput[1024];
char namelistworkrawextrappftinputgrassname[1024];
char namelistworkrawextrappftinputreferenceyear[1024];
char namelistclm5modisrawinput[1024];
char namelistclm5modisrawinputreferenceyear[1024];

int clm5modisrawinputreferenceyear;
int workrawclm5inputreferenceyear;
int workrawextrappftinputreferenceyear;

char clm5modisrawinputname[1024];
char clm5modisrawinputdb[1024];
int clm5modisrawinputstartyear;
int clm5modisrawinputendyear;

char workrawclm5inputname[1024];
char workrawclm5inputdb[1024];
int workrawclm5inputstartyear;
int workrawclm5inputendyear;

char workrawextrappftinputname[1024];
char workrawextrappftinputdb[1024];
int workrawextrappftinputstartyear;
int workrawextrappftinputendyear;

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
float inCTSMPCTPFTGrid[MAXPFT][MAXCLMPIX * MAXCLMLIN];

float inCTSMGRASSEXTRAPPCTGRSC3Grid[MAXCLMPIX * MAXCLMLIN];
float inCTSMGRASSEXTRAPPCTGRSC4Grid[MAXCLMPIX * MAXCLMLIN];

float inGRASSEXTRAPPCTPFTGrid[MAXPFT][MAXCLMPIX * MAXCLMLIN];

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
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3descallpastrallpftname);
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3descrawpftoutput);
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3descrawpftoutputpastrname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5input);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5inputname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawclm5inputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawextrappftinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawextrappftinputgrassname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawextrappftinputreferenceyear);
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


int initializeGrids() {

  long clm5lin, clm5pix;
  long pftid, cftid;
  
  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
	  
          inCTSMGRASSEXTRAPPCTGRSC3Grid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;
          inCTSMGRASSEXTRAPPCTGRSC4Grid[clm5lin * MAXCLMPIX + clm5pix] = 0.0;

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


int readcurrworkgrassextrapGrids(char *databasestr,char *seriesname,char *referenceyearname) {

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
	      inGRASSEXTRAPPCTPFTGrid[pftid][clm5lin * MAXCLMPIX + clm5pix] = tempGrid[clm5lin * MAXCLMPIX + clm5pix];
	      if (pftid == 12 || pftid == 13) {
	          inCTSMGRASSEXTRAPPCTGRSC3Grid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
              }
	      if (pftid == 14) {
	          inCTSMGRASSEXTRAPPCTGRSC4Grid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
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
  float pastrtreepct, pastrgrasspct, pastrbarepct, pastrallpct;
  float addpct, removepct;
  float pastrndlevgtreepct, pastrndldectreepct, pastrbrdevgtreepct, pastrbrddectreepct, pastrshrpct, pastrgrsc3pct, pastrgrsc4pct;
  float pastralltreepct, pastrallgrasspct;
  float pastrndlevgtemptreepct, pastrndlevgborltreepct, pastrndldecborltreepct;
  float pastrbrdevgtroptreepct, pastrbrddectroptreepct, pastrbrdevgtemptreepct, pastrbrddectemptreepct, pastrbrddecborltreepct;
  float pastrshrevgtemppct, pastrshrdectemppct, pastrshrdecborlpct, pastrgrsc3arcpct, pastrgrsc3nonpct, pastrgrsc4nonpct;

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

	      pastrtreepct = 0.0;
	      pastrbarepct = inCTSMPCTPFTGrid[0][clm5lin * MAXCLMPIX + clm5pix];
              pastrbarepct = trunc(pastrbarepct * CLMPCTTruncate) / CLMPCTTruncate;
	      pastrgrasspct = 100.0 - pastrbarepct;
	      
	      pastrgrsc3pct = inCTSMGRASSEXTRAPPCTGRSC3Grid[clm5lin * MAXCLMPIX + clm5pix];
	      pastrgrsc4pct = inCTSMGRASSEXTRAPPCTGRSC4Grid[clm5lin * MAXCLMPIX + clm5pix];
              pastrallgrasspct = pastrgrsc3pct + pastrgrsc4pct;
	      
	      if (pastrallgrasspct > 0.0) {
	          pastrgrsc3pct = pastrgrsc3pct / pastrallgrasspct * pastrgrasspct;
                  pastrgrsc3pct = trunc(pastrgrsc3pct * CLMPCTTruncate) / CLMPCTTruncate;
	          pastrgrsc4pct = pastrgrsc4pct / pastrallgrasspct * pastrgrasspct;
                  pastrgrsc4pct = trunc(pastrgrsc4pct * CLMPCTTruncate) / CLMPCTTruncate;
              }
	      
              if (growdegdaysvalue <= 1000.0) {
                  pastrgrsc3arcpct = pastrgrsc3pct + pastrgrsc4pct;
                  pastrgrsc3nonpct = 0.0;
                  pastrgrsc4nonpct = 0.0;
              }
              else {
                  if (tempcoldestvalue > 22.0 && precipminvalue > 25.0) {
                      pastrgrsc3arcpct = 0.0;
                      pastrgrsc3nonpct = 0.0;
                      pastrgrsc4nonpct = pastrgrsc3pct + pastrgrsc4pct;
		  }
		  else {
                      pastrgrsc3arcpct = 0.0;
                      pastrgrsc3nonpct = pastrgrsc3pct;
                      pastrgrsc4nonpct = pastrgrsc4pct;
		  }
              }

	      pastrndlevgtemptreepct = 0.0;		  
	      pastrndlevgborltreepct = 0.0;		  
	      pastrndldecborltreepct = 0.0;		  
	      pastrbrdevgtroptreepct = 0.0;		  
	      pastrbrdevgtemptreepct = 0.0;		  
	      pastrbrddectroptreepct = 0.0;		  
	      pastrbrddectemptreepct = 0.0;		  
	      pastrbrddecborltreepct = 0.0;		  
	      pastrshrevgtemppct = 0.0;
	      pastrshrdectemppct = 0.0;
	      pastrshrdecborlpct = 0.0;
	      
	      outDESCMERGEDPCTPFTGrid[0][clm5lin * MAXCLMPIX + clm5pix] = pastrbarepct;
	      outDESCMERGEDPCTPFTGrid[1][clm5lin * MAXCLMPIX + clm5pix] = pastrndlevgtemptreepct;		  
	      outDESCMERGEDPCTPFTGrid[2][clm5lin * MAXCLMPIX + clm5pix] = pastrndlevgborltreepct;		  
	      outDESCMERGEDPCTPFTGrid[3][clm5lin * MAXCLMPIX + clm5pix] = pastrndldecborltreepct;		  
	      outDESCMERGEDPCTPFTGrid[4][clm5lin * MAXCLMPIX + clm5pix] = pastrbrdevgtroptreepct;		  
	      outDESCMERGEDPCTPFTGrid[5][clm5lin * MAXCLMPIX + clm5pix] = pastrbrdevgtemptreepct;		  
	      outDESCMERGEDPCTPFTGrid[6][clm5lin * MAXCLMPIX + clm5pix] = pastrbrddectroptreepct;		  
	      outDESCMERGEDPCTPFTGrid[7][clm5lin * MAXCLMPIX + clm5pix] = pastrbrddectemptreepct;		  
	      outDESCMERGEDPCTPFTGrid[8][clm5lin * MAXCLMPIX + clm5pix] = pastrbrddecborltreepct;		  
	      outDESCMERGEDPCTPFTGrid[9][clm5lin * MAXCLMPIX + clm5pix] = pastrshrevgtemppct;
	      outDESCMERGEDPCTPFTGrid[10][clm5lin * MAXCLMPIX + clm5pix] = pastrshrdectemppct;
	      outDESCMERGEDPCTPFTGrid[11][clm5lin * MAXCLMPIX + clm5pix] = pastrshrdecborlpct;
	      outDESCMERGEDPCTPFTGrid[12][clm5lin * MAXCLMPIX + clm5pix] = pastrgrsc3arcpct;
	      outDESCMERGEDPCTPFTGrid[13][clm5lin * MAXCLMPIX + clm5pix] = pastrgrsc3nonpct;
	      outDESCMERGEDPCTPFTGrid[14][clm5lin * MAXCLMPIX + clm5pix] = pastrgrsc4nonpct; 
	  }
      }
  }

  return 0;
  	      
}


int balanceMergedDescriptorVegGrids() {

  long clm5lin, clm5pix;
  long pftid, pastrmaxpftid;
  float landmask, landfrac, area, precipann, tempaverage;
  float pastrpftpct, pastrmaxpftpct, pastrallpftpct;
  float pastrremovepct, pastraddpct;

  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
          landmask = inLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix];
          landfrac = inLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix];
	  if (landfrac > 0.0) {
	  
	      pastrmaxpftid = -1;
	      pastrmaxpftpct = 0.0;
	      pastrallpftpct = 0.0;
	      
	      for (pftid = 0; pftid < MAXPFT; pftid++) {
	          pastrpftpct = outDESCMERGEDPCTPFTGrid[pftid][clm5lin * MAXCLMPIX + clm5pix];
		  pastrallpftpct += pastrpftpct;
		  if (pastrpftpct > pastrmaxpftpct) {
		      pastrmaxpftpct = pastrpftpct;
		      pastrmaxpftid = pftid;
		  }
              }
	      
	      if (pastrallpftpct > 100.0) {
	          pastrremovepct = pastrallpftpct - 100.0;
		  pastrmaxpftpct = outDESCMERGEDPCTPFTGrid[pastrmaxpftid][clm5lin * MAXCLMPIX + clm5pix];
		  if (pastrremovepct > pastrmaxpftpct) {
		      pastrremovepct = pastrmaxpftpct;
		  }
	          outDESCMERGEDPCTPFTGrid[pastrmaxpftid][clm5lin * MAXCLMPIX + clm5pix] -= pastrremovepct;
              }
	      
	      if (pastrallpftpct < 100.0) {
	          pastraddpct = 100.0 - pastrallpftpct;
	          outDESCMERGEDPCTPFTGrid[pastrmaxpftid][clm5lin * MAXCLMPIX + clm5pix] += pastraddpct;
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
  char workrawextrappftinputdirname[1024];
  char luh3descrawpftoutputdirname[1024];
  char clm5modisrawinputreferenceyearstr[256];
  char workrawclm5inputreferenceyearstr[256];
  char workrawextrappftinputreferenceyearstr[256];
  char luh3descrawpftoutputreferenceyearstr[256];
  int currentyear;
  
  if(narg != 8){
        printf("Usage createalldescpastrCTSM53Deg025bin luh3descrawfile luh3descrawdir timeseriesrawfile timeseriesrawdir workrawfile workrawdir luh3descpastrpftnamelistfile\n");
        return 0;
  }
  
  readnamelistfile(argv[7]);

  readluh3descrawpftoutputlutfile(argv[1],argv[2],namelistluh3descrawpftoutput);
  readclm5modisrawinputlutfile(argv[3],argv[4],namelistclm5modisrawinput);
  readworkrawclm5inputlutfile(argv[5],argv[6],namelistworkrawclm5input);
  readworkrawextrappftinputlutfile(argv[5],argv[6],namelistworkrawextrappftinput);
  clm5modisrawinputreferenceyear = atoi(namelistclm5modisrawinputreferenceyear);
  sprintf(clm5modisrawinputreferenceyearstr,"%04d",clm5modisrawinputreferenceyear);
  workrawclm5inputreferenceyear = atoi(namelistworkrawclm5inputreferenceyear);
  sprintf(workrawclm5inputreferenceyearstr,"%04d",workrawclm5inputreferenceyear);
  workrawextrappftinputreferenceyear = atoi(namelistworkrawextrappftinputreferenceyear);
  sprintf(workrawextrappftinputreferenceyearstr,"%04d",workrawextrappftinputreferenceyear);  
  currentyear = workrawclm5inputreferenceyear;

  initializeGrids();
  
  sprintf(workrawclm5inputdirname,"%s/%s",workrawclm5inputdb,workrawclm5inputname);
  sprintf(workrawextrappftinputdirname,"%s/%s",workrawextrappftinputdb,workrawextrappftinputname);
  readlandGrids(workrawclm5inputdirname,namelistworkrawclm5inputname,workrawclm5inputreferenceyearstr);
  readclm5referenceGrids(workrawclm5inputdirname,namelistworkrawclm5inputname,workrawclm5inputreferenceyearstr);    
  readcurrworkgrassextrapGrids(workrawextrappftinputdirname,namelistworkrawextrappftinputgrassname,workrawextrappftinputreferenceyearstr);
  readclimvarseries(clm5modisrawinputdb,clm5modisrawinputname,"CLIM");
  readlaivarseries(clm5modisrawinputdb,clm5modisrawinputname,"FILLCLIM");
  
  genLUHLandGrids();
  
  clearGrids();
  genMergedDescriptorVegGrids();
  balanceMergedDescriptorVegGrids();

  sprintf(luh3descrawpftoutputdirname,"%s/%s",luh3descrawpftoutputdb,luh3descrawpftoutputname);
  sprintf(luh3descrawpftoutputreferenceyearstr,"%04d",currentyear);  
  writedescmergedgrids(luh3descrawpftoutputdirname,namelistluh3descrawpftoutputpastrname,luh3descrawpftoutputreferenceyearstr);

  return 1; 
  
}
