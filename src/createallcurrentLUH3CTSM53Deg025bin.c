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
char namelistcurrentname[1024];
char namelistworkrawoutput[1024];
char namelistworkrawoutputname[1024];
char namelistluh3rawinput[1024];
char namelistluh3rawinputreferenceyear[1024];
char namelistclm5currentrawinput[1024];
char namelistclm5currentrawinputreferenceyear[1024];

int luh3rawinputreferenceyear;
int clm5currentrawinputreferenceyear;

char luh3rawinputname[1024];
char luh3rawinputdb[1024];
int luh3rawinputstartyear;
int luh3rawinputendyear;

char clm5currentrawinputname[1024];
char clm5currentrawinputdb[1024];
int clm5currentrawinputstartyear;
int clm5currentrawinputendyear;

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
float inPCTGLACIERGrid[MAXCLMPIX * MAXCLMLIN];
float inPCTLAKEGrid[MAXCLMPIX * MAXCLMLIN];
float inPCTWETLANDGrid[MAXCLMPIX * MAXCLMLIN];
float inPCTURBANGrid[MAXCLMPIX * MAXCLMLIN];
float inPCTNATVEGGrid[MAXCLMPIX * MAXCLMLIN];
float inPCTCROPGrid[MAXCLMPIX * MAXCLMLIN];
float inCURRENTPCTPFTGrid[MAXPFT][MAXCLMPIX * MAXCLMLIN];
float inCURRENTPCTCFTGrid[MAXCFT][MAXCLMPIX * MAXCLMLIN];

float inLUHPCTCROPGrid[MAXCLMPIX * MAXCLMLIN];

float intempwarmest[MAXCLMPIX * MAXCLMLIN];
float intempcoldest[MAXCLMPIX * MAXCLMLIN];
float ingrowdegdays[MAXCLMPIX * MAXCLMLIN];
float inprecipann[MAXCLMPIX * MAXCLMLIN];
float inprecipmin[MAXCLMPIX * MAXCLMLIN];
float inprecipwin[MAXCLMPIX * MAXCLMLIN];
float inprecipmaxtg22[MAXCLMPIX * MAXCLMLIN];

float inc3LAIGrid[MAXCLMPIX * MAXCLMLIN];
float inc4LAIGrid[MAXCLMPIX * MAXCLMLIN];
float inminLAIGrid[MAXCLMPIX * MAXCLMLIN];
float inmaxLAIGrid[MAXCLMPIX * MAXCLMLIN];


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
  fscanf(namelistfile,"%s %s",fieldname,namelistcurrentname);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawoutput);
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawoutputname);
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3rawinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistluh3rawinputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5currentrawinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5currentrawinputreferenceyear);
  
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
  fread(inPCTGLACIERGrid,sizeof(inPCTGLACIERGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTLAKE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPCTLAKEGrid,sizeof(inPCTLAKEGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTWETLAND.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPCTWETLANDGrid,sizeof(inPCTWETLANDGrid),1,globalbinfile);  
  fclose(globalbinfile);

  return 0;
  
}


int readclm5referenceGrids(char *databasestr,char *seriesname,char *yearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int pftid, cftid;
  char pftidstr[256], cftidstr[256];
  long clmlin, clmpix;
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTURBAN.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPCTURBANGrid,sizeof(inPCTURBANGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTNATVEG.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPCTNATVEGGrid,sizeof(inPCTNATVEGGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTCROP.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inPCTCROPGrid,sizeof(inPCTCROPGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      sprintf(pftidstr,"%02d",pftid);
      sprintf(globalbinfilename,"%s/%s/%s.PCTNATPFT%s.%s.dat",databasestr,seriesname,seriesname,pftidstr,yearname);
      printf("Reading: %s\n",globalbinfilename);
      globalbinfile = fopen(globalbinfilename,"r");
      fread(tempGrid,sizeof(tempGrid),1,globalbinfile);  
      fclose(globalbinfile);
      for (clmlin = 0; clmlin < MAXCLMLIN; clmlin++) {
          for (clmpix = 0; clmpix < MAXCLMPIX; clmpix++) {
	      inCURRENTPCTPFTGrid[pftid][clmlin * MAXCLMPIX + clmpix] = tempGrid[clmlin * MAXCLMPIX + clmpix];
	  }
      }
  }
  
  for (cftid = 0; cftid < MAXCFT; cftid++) {
      sprintf(cftidstr,"%02d",cftid);
      sprintf(globalbinfilename,"%s/%s/%s.PCTCFT%s.%s.dat",databasestr,seriesname,seriesname,cftidstr,yearname);
      printf("Reading: %s\n",globalbinfilename);
      globalbinfile = fopen(globalbinfilename,"r");
      fread(tempGrid,sizeof(tempGrid),1,globalbinfile);  
      fclose(globalbinfile);
      for (clmlin = 0; clmlin < MAXCLMLIN; clmlin++) {
          for (clmpix = 0; clmpix < MAXCLMPIX; clmpix++) {
	      inCURRENTPCTCFTGrid[cftid][clmlin * MAXCLMPIX + clmpix] = tempGrid[clmlin * MAXCLMPIX + clmpix];
	  }
      }
  }

  return 0;
  
}


int readluhcropGrids(char *databasestr,char *seriesname,char *yearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  int cftid;
  char cftidstr[256];
  long clmlin, clmpix;
  float cropfrac;
  
  sprintf(globalbinfilename,"%s/%s/%s.c3ann.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(tempGrid,sizeof(tempGrid),1,globalbinfile);  
  fclose(globalbinfile);

  for (clmlin = 0; clmlin < MAXCLMLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXCLMPIX; clmpix++) {
          cropfrac = tempGrid[clmlin * MAXCLMPIX + clmpix];
          if (cropfrac < 0.0 || cropfrac > 1.0) {
	      cropfrac = 0.0;
	  }
          inLUHPCTCROPGrid[clmlin * MAXCLMPIX + clmpix] = cropfrac * 100.0;
      }
  }

  sprintf(globalbinfilename,"%s/%s/%s.c4ann.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(tempGrid,sizeof(tempGrid),1,globalbinfile);  
  fclose(globalbinfile);

  for (clmlin = 0; clmlin < MAXCLMLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXCLMPIX; clmpix++) {
          cropfrac = tempGrid[clmlin * MAXCLMPIX + clmpix];
          if (cropfrac < 0.0 || cropfrac > 1.0) {
	      cropfrac = 0.0;
	  }
          inLUHPCTCROPGrid[clmlin * MAXCLMPIX + clmpix] += cropfrac * 100.0;
      }
  }

  sprintf(globalbinfilename,"%s/%s/%s.c3per.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(tempGrid,sizeof(tempGrid),1,globalbinfile);  
  fclose(globalbinfile);

  for (clmlin = 0; clmlin < MAXCLMLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXCLMPIX; clmpix++) {
          cropfrac = tempGrid[clmlin * MAXCLMPIX + clmpix];
          if (cropfrac < 0.0 || cropfrac > 1.0) {
	      cropfrac = 0.0;
	  }
          inLUHPCTCROPGrid[clmlin * MAXCLMPIX + clmpix] += cropfrac * 100.0;
      }
  }

  sprintf(globalbinfilename,"%s/%s/%s.c4per.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(tempGrid,sizeof(tempGrid),1,globalbinfile);  
  fclose(globalbinfile);

  for (clmlin = 0; clmlin < MAXCLMLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXCLMPIX; clmpix++) {
          cropfrac = tempGrid[clmlin * MAXCLMPIX + clmpix];
          if (cropfrac < 0.0 || cropfrac > 1.0) {
	      cropfrac = 0.0;
	  }
          inLUHPCTCROPGrid[clmlin * MAXCLMPIX + clmpix] += cropfrac * 100.0;
      }
  }

  sprintf(globalbinfilename,"%s/%s/%s.c3nfx.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(tempGrid,sizeof(tempGrid),1,globalbinfile);  
  fclose(globalbinfile);

  for (clmlin = 0; clmlin < MAXCLMLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXCLMPIX; clmpix++) {
          cropfrac = tempGrid[clmlin * MAXCLMPIX + clmpix];
          if (cropfrac < 0.0 || cropfrac > 1.0) {
	      cropfrac = 0.0;
	  }
          inLUHPCTCROPGrid[clmlin * MAXCLMPIX + clmpix] += cropfrac * 100.0;
      }
  }

  return 0;
  
}


int readclimvarseries(char *databasestr, char *seriesname, char *fileyear) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  
  sprintf(globalbinfilename,"%s/%s/%s.tempwarmest-%s.dat",databasestr,seriesname,seriesname,fileyear);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(intempwarmest,sizeof(intempwarmest),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.tempcoldest-%s.dat",databasestr,seriesname,seriesname,fileyear);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(intempcoldest,sizeof(intempcoldest),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.growdegdays-%s.dat",databasestr,seriesname,seriesname,fileyear);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(ingrowdegdays,sizeof(ingrowdegdays),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.precipann-%s.dat",databasestr,seriesname,seriesname,fileyear);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inprecipann,sizeof(inprecipann),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.precipmin-%s.dat",databasestr,seriesname,seriesname,fileyear);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inprecipmin,sizeof(inprecipmin),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.precipwin-%s.dat",databasestr,seriesname,seriesname,fileyear);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inprecipwin,sizeof(inprecipwin),1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s.precipmaxtg22-%s.dat",databasestr,seriesname,seriesname,fileyear);
  printf("Reading: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"r");
  fread(inprecipmaxtg22,sizeof(inprecipmaxtg22),1,globalbinfile);  
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
  fread(inminLAIGrid,sizeof(inminLAIGrid),1,globalbinlaifile);  
  fclose(globalbinlaifile);

  sprintf(globalbinlaifilename,"%s/%s/%s.MaxLai_1km.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinlaifilename);
  globalbinlaifile = fopen(globalbinlaifilename,"r");
  fread(inmaxLAIGrid,sizeof(inmaxLAIGrid),1,globalbinlaifile);  
  fclose(globalbinlaifile);

  sprintf(globalbinlaifilename,"%s/%s/%s.C3Lai_1km.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinlaifilename);
  globalbinlaifile = fopen(globalbinlaifilename,"r");
  fread(inc3LAIGrid,sizeof(inc3LAIGrid),1,globalbinlaifile);  
  fclose(globalbinlaifile);

  sprintf(globalbinlaifilename,"%s/%s/%s.C4Lai_1km.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Reading: %s\n",globalbinlaifilename);
  globalbinlaifile = fopen(globalbinlaifilename,"r");
  fread(inc4LAIGrid,sizeof(inc4LAIGrid),1,globalbinlaifile);  
  fclose(globalbinlaifile);

  return 0;

}


int genLUHLandMask() {

  long clmlin, clmpix;
  long luhlin, luhpix;
  
  for (luhlin = 0; luhlin < MAXLUHLIN; luhlin++) {
      for (luhpix = 0; luhpix < MAXLUHPIX; luhpix++) {
	  inLUHLANDMASKGrid[luhlin * MAXLUHPIX + luhpix] = 1.0;
      }
  }

  return 0;

}


float truncCLMValues(float CLMValueIn, float CLMMaxValue) {

  float CLMValueOut = 0.0;
  
  if (CLMValueIn >= CLMPCTThreshold) {
      CLMValueOut = ((float) ((int) (CLMValueIn / CLMPCTThreshold))) * CLMPCTThreshold;
  }
  
  if (CLMValueOut > CLMMaxValue) {
      CLMValueOut = CLMMaxValue;
  }
  
  return CLMValueOut;
    
}


int regroupcurrentpftgrids() {

  long clmlin, clmpix;
  int pftid;
  int resetextrap;
  float landfrac, vegfrac, natvegfrac, barefrac, treefrac;
  float clmvegpct, clmcroppct, luhcroppct;
  float newclmvegpct;
  float currentgrasspct, newgrasspct, currentotherpct, newotherpct;
  float totalpftpct, pftpct, pftlevel, maxextraplevel;
  float firstpftpct, secondpftpct, thirdpftpct, fourthpftpct, largestpftpct, rescaledpftpct;
  int pftcount, firstpftid, secondpftid, thirdpftid, fourthpftid, largestpftid;
  float minlaivalue, maxlaivalue, c3laivalue, c4laivalue;
  float tempcoldestvalue, tempwarmestvalue, growdegdaysvalue, precipannvalue, precipminvalue, precipwinvalue, precipmaxtg22value;
  float c3c4laivalue, c3laifrac, c4laifrac;
  float ndlevgtreepct, ndldectreepct, brdevgtreepct, brddectreepct, shrubpct, grasspct, barepct, oldbarepct;
  
  for (clmlin = 0; clmlin < MAXCLMLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXCLMPIX; clmpix++) {
          landfrac = inLANDFRACGrid[clmlin * MAXCLMPIX + clmpix];
          clmvegpct = inPCTNATVEGGrid[clmlin * MAXCLMPIX + clmpix];
          clmcroppct = inPCTCROPGrid[clmlin * MAXCLMPIX + clmpix];
	  luhcroppct = truncCLMValues(inLUHPCTCROPGrid[clmlin * MAXCLMPIX + clmpix],100.0);
	  if (landfrac > 0.0) {
	      if (clmcroppct > luhcroppct) {

	          minlaivalue = inminLAIGrid[clmlin * MAXCLMPIX + clmpix];
	          maxlaivalue = inmaxLAIGrid[clmlin * MAXCLMPIX + clmpix];
	          c3laivalue = inc3LAIGrid[clmlin * MAXCLMPIX + clmpix];
	          c4laivalue = inc4LAIGrid[clmlin * MAXCLMPIX + clmpix];
		  tempcoldestvalue = intempcoldest[clmlin * MAXCLMPIX + clmpix];
		  tempwarmestvalue = intempwarmest[clmlin * MAXCLMPIX + clmpix];
		  growdegdaysvalue = ingrowdegdays[clmlin * MAXCLMPIX + clmpix];
		  precipannvalue = inprecipann[clmlin * MAXCLMPIX + clmpix];
                  precipminvalue = inprecipmin[clmlin * MAXCLMPIX + clmpix];
                  precipwinvalue = inprecipwin[clmlin * MAXCLMPIX + clmpix];
		  precipmaxtg22value = inprecipmaxtg22[clmlin * MAXCLMPIX + clmpix];

	          newclmvegpct = clmvegpct + (clmcroppct - luhcroppct);
                  inPCTNATVEGGrid[clmlin * MAXCLMPIX + clmpix] = newclmvegpct;
                  inPCTCROPGrid[clmlin * MAXCLMPIX + clmpix] = luhcroppct;
		  
	          currentgrasspct = 0.0;
                  newgrasspct = currentgrasspct;
	          for (pftid = 12; pftid <= 14; pftid++) {
                      currentgrasspct += inCURRENTPCTPFTGrid[pftid][clmlin * MAXCLMPIX + clmpix];
                  }
		  if (newclmvegpct > 0.0) {
		      newgrasspct = currentgrasspct + (clmcroppct - luhcroppct) / newclmvegpct * 100.0;
		      if (newgrasspct > 100.0) {
		          newgrasspct = 100.0;
                      }
		  }

	          if (growdegdaysvalue <= 1000.0) {
	              inCURRENTPCTPFTGrid[12][clmlin * MAXCLMPIX + clmpix] = newgrasspct;
		      inCURRENTPCTPFTGrid[13][clmlin * MAXCLMPIX + clmpix] = 0.0;
		      inCURRENTPCTPFTGrid[14][clmlin * MAXCLMPIX + clmpix] = 0.0;
		  }
		  else {
		      if (tempwarmestvalue <= 22.0 || precipmaxtg22value <= 25.0) {
			  inCURRENTPCTPFTGrid[12][clmlin * MAXCLMPIX + clmpix] = 0.0;
			  inCURRENTPCTPFTGrid[13][clmlin * MAXCLMPIX + clmpix] = newgrasspct;
			  inCURRENTPCTPFTGrid[14][clmlin * MAXCLMPIX + clmpix] = 0.0;
	              }
		      else {
			  if (tempcoldestvalue > 22.0 && precipminvalue > 25.0) {
			      inCURRENTPCTPFTGrid[12][clmlin * MAXCLMPIX + clmpix] = 0.0;
			      inCURRENTPCTPFTGrid[13][clmlin * MAXCLMPIX + clmpix] = 0.0;
			      inCURRENTPCTPFTGrid[14][clmlin * MAXCLMPIX + clmpix] = newgrasspct;
			  }
			  else {
                              c3c4laivalue = c3laivalue + c4laivalue;
		              if (c3c4laivalue > 0.0) {
		                  c3laifrac = c3laivalue / c3c4laivalue;
	                          c4laifrac = 1.0 - c3laifrac;
			          inCURRENTPCTPFTGrid[12][clmlin * MAXCLMPIX + clmpix] = 0.0;
			          inCURRENTPCTPFTGrid[13][clmlin * MAXCLMPIX + clmpix] = c3laifrac * newgrasspct;
			          inCURRENTPCTPFTGrid[14][clmlin * MAXCLMPIX + clmpix] = c4laifrac * newgrasspct;
			      }
			      else {
			          inCURRENTPCTPFTGrid[12][clmlin * MAXCLMPIX + clmpix] = 0.0;
			          inCURRENTPCTPFTGrid[13][clmlin * MAXCLMPIX + clmpix] = 0.5 * newgrasspct;
			          inCURRENTPCTPFTGrid[14][clmlin * MAXCLMPIX + clmpix] = 0.5 * newgrasspct;
			      }
			  }
  		      }
		  }
		  
                  currentotherpct = 0.0;
	          for (pftid = 1; pftid <= 8; pftid++) {
                      currentotherpct += inCURRENTPCTPFTGrid[pftid][clmlin * MAXCLMPIX + clmpix];
                  } 
                  newotherpct = currentotherpct - (newgrasspct - currentgrasspct);
		  if (newotherpct < 0.0) {
		      newotherpct = 0.0;
		  }
		  
		  if (currentotherpct > 0.0) {
                      for (pftid = 1; pftid <= 8; pftid++) {
                          inCURRENTPCTPFTGrid[pftid][clmlin * MAXCLMPIX + clmpix] = newotherpct / currentotherpct * inCURRENTPCTPFTGrid[pftid][clmlin * MAXCLMPIX + clmpix];
                      }
		  }
		     	  
              }
	  }

	  if (landfrac > 0.0) {
              totalpftpct = 0.0;
              for (pftid = 0; pftid < MAXPFT; pftid++) {
                  totalpftpct += inCURRENTPCTPFTGrid[pftid][clmlin * MAXCLMPIX + clmpix];
              }
	      if (totalpftpct > 0.0) {
	          largestpftpct = 0.0;
	          largestpftid = 0;
	          rescaledpftpct = 0.0;
	          for (pftid = 0; pftid < MAXPFT; pftid++) {
                      pftpct = inCURRENTPCTPFTGrid[pftid][clmlin * MAXCLMPIX + clmpix] * 100.0 / totalpftpct;
	              pftpct = truncCLMValues(pftpct,100.0);
		      inCURRENTPCTPFTGrid[pftid][clmlin * MAXCLMPIX + clmpix] = pftpct;
		      rescaledpftpct += pftpct;
		      if (pftpct > largestpftpct) {
		          largestpftpct = pftpct;
		          largestpftid = pftid;
                      }
	          }
 	          inCURRENTPCTPFTGrid[largestpftid][clmlin * MAXCLMPIX + clmpix] += (100.0 -  rescaledpftpct);
	      }
	  }
      }
  }

  return 0;
  
}
		  
	      
int writecurrentgrids(char *databasestr, char *outputdir, char *seriesname, char *yearname) {

  FILE *globalbinfile;
  char globalbinfilename[256];
  long outlin, outpix;
  long clmlin, clmpix;
  int pftid, cftid;
  char pftidstr[256];
  char cftidstr[256];
 
  tempoutGrid = (float *) malloc(OUTDATASIZE);

  sprintf(globalbinfilename,"%s/%s/%s/%s.LANDMASK.%s.dat",databasestr,outputdir,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  for (outlin = 0; outlin < MAXOUTLIN; outlin++) {
      clmlin = outlin + OUTLATOFFSET;
      for (outpix = 0; outpix < MAXOUTPIX; outpix++) {
          clmpix = outpix + OUTLONOFFSET;
          tempoutGrid[outlin * MAXOUTPIX + outpix] = inLANDMASKGrid[clmlin * MAXCLMPIX + clmpix];
      }
  }
  fwrite(tempoutGrid,OUTDATASIZE,1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s/%s.LANDFRAC.%s.dat",databasestr,outputdir,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  for (outlin = 0; outlin < MAXOUTLIN; outlin++) {
      clmlin = outlin + OUTLATOFFSET;
      for (outpix = 0; outpix < MAXOUTPIX; outpix++) {
          clmpix = outpix + OUTLONOFFSET;
          tempoutGrid[outlin * MAXOUTPIX + outpix] = inLANDFRACGrid[clmlin * MAXCLMPIX + clmpix];
      }
  }
  fwrite(tempoutGrid,OUTDATASIZE,1,globalbinfile);  
  fclose(globalbinfile);

  sprintf(globalbinfilename,"%s/%s/%s/%s.AREA.%s.dat",databasestr,outputdir,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  for (outlin = 0; outlin < MAXOUTLIN; outlin++) {
      clmlin = outlin + OUTLATOFFSET;
      for (outpix = 0; outpix < MAXOUTPIX; outpix++) {
          clmpix = outpix + OUTLONOFFSET;
          tempoutGrid[outlin * MAXOUTPIX + outpix] = inAREAGrid[clmlin * MAXCLMPIX + clmpix];
      }
  }
  fwrite(tempoutGrid,OUTDATASIZE,1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s/%s.PCTGLACIER.%s.dat",databasestr,outputdir,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  for (outlin = 0; outlin < MAXOUTLIN; outlin++) {
      clmlin = outlin + OUTLATOFFSET;
      for (outpix = 0; outpix < MAXOUTPIX; outpix++) {
          clmpix = outpix + OUTLONOFFSET;
          tempoutGrid[outlin * MAXOUTPIX + outpix] = inPCTGLACIERGrid[clmlin * MAXCLMPIX + clmpix];
      }
  }
  fwrite(tempoutGrid,OUTDATASIZE,1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s/%s.PCTLAKE.%s.dat",databasestr,outputdir,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  for (outlin = 0; outlin < MAXOUTLIN; outlin++) {
      clmlin = outlin + OUTLATOFFSET;
      for (outpix = 0; outpix < MAXOUTPIX; outpix++) {
          clmpix = outpix + OUTLONOFFSET;
          tempoutGrid[outlin * MAXOUTPIX + outpix] = inPCTLAKEGrid[clmlin * MAXCLMPIX + clmpix];
      }
  }
  fwrite(tempoutGrid,OUTDATASIZE,1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s/%s.PCTWETLAND.%s.dat",databasestr,outputdir,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  for (outlin = 0; outlin < MAXOUTLIN; outlin++) {
      clmlin = outlin + OUTLATOFFSET;
      for (outpix = 0; outpix < MAXOUTPIX; outpix++) {
          clmpix = outpix + OUTLONOFFSET;
          tempoutGrid[outlin * MAXOUTPIX + outpix] = inPCTWETLANDGrid[clmlin * MAXCLMPIX + clmpix];
      }
  }
  fwrite(tempoutGrid,OUTDATASIZE,1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s/%s.PCTURBAN.%s.dat",databasestr,outputdir,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  for (outlin = 0; outlin < MAXOUTLIN; outlin++) {
      clmlin = outlin + OUTLATOFFSET;
      for (outpix = 0; outpix < MAXOUTPIX; outpix++) {
          clmpix = outpix + OUTLONOFFSET;
          tempoutGrid[outlin * MAXOUTPIX + outpix] = inPCTURBANGrid[clmlin * MAXCLMPIX + clmpix];
      }
  }
  fwrite(tempoutGrid,OUTDATASIZE,1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s/%s.PCTNATVEG.%s.dat",databasestr,outputdir,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  for (outlin = 0; outlin < MAXOUTLIN; outlin++) {
      clmlin = outlin + OUTLATOFFSET;
      for (outpix = 0; outpix < MAXOUTPIX; outpix++) {
          clmpix = outpix + OUTLONOFFSET;
          tempoutGrid[outlin * MAXOUTPIX + outpix] = inPCTNATVEGGrid[clmlin * MAXCLMPIX + clmpix];
      }
  }
  fwrite(tempoutGrid,OUTDATASIZE,1,globalbinfile);  
  fclose(globalbinfile);
  
  sprintf(globalbinfilename,"%s/%s/%s/%s.PCTCROP.%s.dat",databasestr,outputdir,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  for (outlin = 0; outlin < MAXOUTLIN; outlin++) {
      clmlin = outlin + OUTLATOFFSET;
      for (outpix = 0; outpix < MAXOUTPIX; outpix++) {
          clmpix = outpix + OUTLONOFFSET;
          tempoutGrid[outlin * MAXOUTPIX + outpix] = inPCTCROPGrid[clmlin * MAXCLMPIX + clmpix];
      }
  }
  fwrite(tempoutGrid,OUTDATASIZE,1,globalbinfile);  
  fclose(globalbinfile);
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      sprintf(pftidstr,"%02d",pftid);
      sprintf(globalbinfilename,"%s/%s/%s/%s.PCTNATPFT%s.%s.dat",databasestr,outputdir,seriesname,seriesname,pftidstr,yearname);
      printf("Writing: %s\n",globalbinfilename);
      globalbinfile = fopen(globalbinfilename,"w+");
      for (outlin = 0; outlin < MAXOUTLIN; outlin++) {
          clmlin = outlin + OUTLATOFFSET;
          for (outpix = 0; outpix < MAXOUTPIX; outpix++) {
              clmpix = outpix + OUTLONOFFSET;
              tempoutGrid[outlin * MAXOUTPIX + outpix] = inCURRENTPCTPFTGrid[pftid][clmlin * MAXCLMPIX + clmpix];
          }
      }
      fwrite(tempoutGrid,OUTDATASIZE,1,globalbinfile);  
      fclose(globalbinfile);
  }
  
  for (cftid = 0; cftid < MAXCFT; cftid++) {
      sprintf(cftidstr,"%02d",cftid);
      sprintf(globalbinfilename,"%s/%s/%s/%s.PCTCFT%s.%s.dat",databasestr,outputdir,seriesname,seriesname,cftidstr,yearname);
      printf("Writing: %s\n",globalbinfilename);
      globalbinfile = fopen(globalbinfilename,"w+");
      for (outlin = 0; outlin < MAXOUTLIN; outlin++) {
          clmlin = outlin + OUTLATOFFSET;
          for (outpix = 0; outpix < MAXOUTPIX; outpix++) {
              clmpix = outpix + OUTLONOFFSET;
              tempoutGrid[outlin * MAXOUTPIX + outpix] = inCURRENTPCTCFTGrid[cftid][clmlin * MAXCLMPIX + clmpix];
          }
      }
      fwrite(tempoutGrid,OUTDATASIZE,1,globalbinfile);  
      fclose(globalbinfile);
  }
  
  return 0;

}


int main(long narg, char **argv) {

  char luh3rawinputreferenceyearstr[256];
  char clm5currentrawreferenceyearstr[256];
  char workrawoutputreferenceyearstr[256];
  char workrawoutputdirname[1024];
  
  if(narg != 8){
        printf("Usage createalldesccurrentLUH3clm5Deg025bin luh3rawfile luh3rawdir timeseriesrawfile timeseriesrawdir workrawfile workrawdir generateallnamelist\n");
        return 0;
  }
  
  readnamelistfile(argv[7]);

  readluh3rawinputlutfile(argv[1],argv[2],namelistluh3rawinput);
  readclm5currentrawinputlutfile(argv[3],argv[4],namelistclm5currentrawinput);
  readworkrawoutputlutfile(argv[5],argv[6],namelistworkrawoutput);
  luh3rawinputreferenceyear = atoi(namelistluh3rawinputreferenceyear);
  sprintf(luh3rawinputreferenceyearstr,"%04d",luh3rawinputreferenceyear);
  clm5currentrawinputreferenceyear = atoi(namelistclm5currentrawinputreferenceyear);
  sprintf(clm5currentrawreferenceyearstr,"%04d",clm5currentrawinputreferenceyear);

  readlandGrids(clm5currentrawinputdb,clm5currentrawinputname,clm5currentrawreferenceyearstr);
  readclm5referenceGrids(clm5currentrawinputdb,clm5currentrawinputname,clm5currentrawreferenceyearstr);    
  readluhcropGrids(luh3rawinputdb,luh3rawinputname,luh3rawinputreferenceyearstr);
  readclimvarseries(clm5currentrawinputdb,clm5currentrawinputname,"CLIM");
  readlaivarseries(clm5currentrawinputdb,clm5currentrawinputname,"FILLCLIM");

  genLUHLandMask();
  regroupcurrentpftgrids();
  
  sprintf(workrawoutputreferenceyearstr,"%s",clm5currentrawreferenceyearstr);
  sprintf(workrawoutputdirname,"%s/%s",workrawoutputdb,workrawoutputname);
  writecurrentgrids(workrawoutputdb,workrawoutputname,namelistworkrawoutputname,clm5currentrawreferenceyearstr);

  return 1;
  
}
