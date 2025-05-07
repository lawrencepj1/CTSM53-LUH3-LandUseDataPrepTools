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
char namelistworkrawoutputname[1024];
char namelistworkrawinput[1024];
char namelistworkrawinputname[1024];
char namelistworkrawinputreferenceyear[1024];
char namelistclm5modisrawinput[1024];
char namelistclm5referenceyear[1024];
char namelistnonforesttreeclimatename[1024];
char namelistnonforestherbclimatename[1024];
char namelistnonforestbareclimatename[1024];

int workrawinputreferenceyear;

char clm5modisrawinputname[1024];
char clm5modisrawinputdb[1024];
int clm5modisrawinputstartyear;
int clm5modisrawinputendyear;

char workrawinputname[1024];
char workrawinputdb[1024];
int workrawinputstartyear;
int workrawinputendyear;

char workrawoutputname[1024];
char workrawoutputdb[1024];
int workrawoutputstartyear;
int workrawoutputendyear;

char nonforesttreeclimatename[1024];
char nonforesttreeclimatedb[1024];
int nonforesttreeclimatestartyear;
int nonforesttreeclimateendyear;

char nonforestherbclimatename[1024];
char nonforestherbclimatedb[1024];
int nonforestherbclimatestartyear;
int nonforestherbclimateendyear;

char nonforestbareclimatename[1024];
char nonforestbareclimatedb[1024];
int nonforestbareclimatestartyear;
int nonforestbareclimateendyear;

float nonforestpotentialtree[PRECIPCLASSES][TEMPCLASSES];
float nonforestpotentialherb[PRECIPCLASSES][TEMPCLASSES];
float nonforestpotentialbare[PRECIPCLASSES][TEMPCLASSES];

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
float inclm5PCTBAREGrid[MAXCLMPIX * MAXCLMLIN];

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
  fscanf(namelistfile,"%s %s",fieldname,namelistworkrawinputreferenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5modisrawinput);
  fscanf(namelistfile,"%s %s",fieldname,namelistclm5referenceyear);
  fscanf(namelistfile,"%s %s",fieldname,namelistnonforesttreeclimatename);
  fscanf(namelistfile,"%s %s",fieldname,namelistnonforestherbclimatename);
  fscanf(namelistfile,"%s %s",fieldname,namelistnonforestbareclimatename);
  
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
      itemsfound = fscanf(workrawoutputvaluefile,"%s %s %s %s %s",inworkrawoutputname, inworkrawoutputdb, inworkrawoutputstartyear, inworkrawoutputendyear);
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


int readnonforesttreeclimatelutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *nonforesttreeclimatevaluefile;
  int foundnonforesttreeclimate, itemsfound;
  char fullfilenamestr[1024], innonforesttreeclimatename[1024], innonforesttreeclimatedb[1024], innonforesttreeclimatefile[1024];
  char innonforesttreeclimatestartyear[1024], innonforesttreeclimateendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  nonforesttreeclimatevaluefile = fopen(fullfilenamestr,"r");
  
  foundnonforesttreeclimate = 0;

  while (foundnonforesttreeclimate == 0) {
      itemsfound = fscanf(nonforesttreeclimatevaluefile,"%s %s %s %s %s",innonforesttreeclimatename, innonforesttreeclimatedb, innonforesttreeclimatefile, innonforesttreeclimatestartyear, innonforesttreeclimateendyear);
      if (itemsfound != 5) {
          printf("Error: nonforesttreeclimate %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,innonforesttreeclimatename) == 0) {
          printf("Processing nonforesttreeclimate: %s\n",innonforesttreeclimatename);
	  sprintf(nonforesttreeclimatename,"%s",innonforesttreeclimatename);
          if (strcmp(innonforesttreeclimatedb,"<luh3modiscrurawdir>") == 0) {
              sprintf(nonforesttreeclimatedb,"%s%s",dirnamestr,innonforesttreeclimatefile);
	  }
	  else {
              sprintf(nonforesttreeclimatedb,"%s%s",innonforesttreeclimatedb,innonforesttreeclimatefile);
	  }
          nonforesttreeclimatestartyear = atoi(innonforesttreeclimatestartyear);
          nonforesttreeclimateendyear = atoi(innonforesttreeclimateendyear);
	  foundnonforesttreeclimate = 1;
      }
  }  
  
  return 0;

}


int readnonforestherbclimatelutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *nonforestherbclimatevaluefile;
  int foundnonforestherbclimate, itemsfound;
  char fullfilenamestr[1024], innonforestherbclimatename[1024], innonforestherbclimatedb[1024], innonforestherbclimatefile[1024];
  char innonforestherbclimatestartyear[1024], innonforestherbclimateendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  nonforestherbclimatevaluefile = fopen(fullfilenamestr,"r");
  
  foundnonforestherbclimate = 0;

  while (foundnonforestherbclimate == 0) {
      itemsfound = fscanf(nonforestherbclimatevaluefile,"%s %s %s %s %s",innonforestherbclimatename, innonforestherbclimatedb, innonforestherbclimatefile, innonforestherbclimatestartyear, innonforestherbclimateendyear);
      if (itemsfound != 5) {
          printf("Error: nonforestherbclimate %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,innonforestherbclimatename) == 0) {
          printf("Processing nonforestherbclimate: %s\n",innonforestherbclimatename);
	  sprintf(nonforestherbclimatename,"%s",innonforestherbclimatename);
          if (strcmp(innonforestherbclimatedb,"<luh3modiscrurawdir>") == 0) {
              sprintf(nonforestherbclimatedb,"%s%s",dirnamestr,innonforestherbclimatefile);
	  }
	  else {
              sprintf(nonforestherbclimatedb,"%s%s",innonforestherbclimatedb,innonforestherbclimatefile);
	  }
          nonforestherbclimatestartyear = atoi(innonforestherbclimatestartyear);
          nonforestherbclimateendyear = atoi(innonforestherbclimateendyear);
	  foundnonforestherbclimate = 1;
      }
  }  
  
  return 0;

}


int readnonforestbareclimatelutfile(char *filenamestr, char *dirnamestr, char *namestr) {

  FILE *nonforestbareclimatevaluefile;
  int foundnonforestbareclimate, itemsfound;
  char fullfilenamestr[1024], innonforestbareclimatename[1024], innonforestbareclimatedb[1024], innonforestbareclimatefile[1024];
  char innonforestbareclimatestartyear[1024], innonforestbareclimateendyear[1024];
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  nonforestbareclimatevaluefile = fopen(fullfilenamestr,"r");
  
  foundnonforestbareclimate = 0;

  while (foundnonforestbareclimate == 0) {
      itemsfound = fscanf(nonforestbareclimatevaluefile,"%s %s %s %s %s",innonforestbareclimatename, innonforestbareclimatedb, innonforestbareclimatefile, innonforestbareclimatestartyear, innonforestbareclimateendyear);
      if (itemsfound != 5) {
          printf("Error: nonforestbareclimate %s Not Found\n",namestr);
	  exit(0);
      }
      if (strcmp(namestr,innonforestbareclimatename) == 0) {
          printf("Processing nonforestbareclimate: %s\n",innonforestbareclimatename);
	  sprintf(nonforestbareclimatename,"%s",innonforestbareclimatename);
          if (strcmp(innonforestbareclimatedb,"<luh3modiscrurawdir>") == 0) {
              sprintf(nonforestbareclimatedb,"%s%s",dirnamestr,innonforestbareclimatefile);
	  }
	  else {
              sprintf(nonforestbareclimatedb,"%s%s",innonforestbareclimatedb,innonforestbareclimatefile);
	  }
          nonforestbareclimatestartyear = atoi(innonforestbareclimatestartyear);
          nonforestbareclimateendyear = atoi(innonforestbareclimateendyear);
	  foundnonforestbareclimate = 1;
      }
  }  
  
  return 0;

}


int readpotentialveglut(char *filenamestr, int veglutid) {

  FILE *climatevaluefile;
  int precipannindex, tempaverageindex;
  char fullfilenamestr[1024];
  char tempstrvalue[1024];
  char *temptokenvalue;
  float tempfloatvalue;
  
  sprintf(fullfilenamestr,"%s",filenamestr);
  printf("Reading %s\n",fullfilenamestr);
  climatevaluefile = fopen(fullfilenamestr,"r");
  
  for (precipannindex = 0; precipannindex < PRECIPCLASSES; precipannindex++) {
      fgets(tempstrvalue, 1024, climatevaluefile);
      temptokenvalue = strtok(tempstrvalue," ");
      for (tempaverageindex = 0; tempaverageindex < TEMPCLASSES; tempaverageindex++) {
	  tempfloatvalue = atof(temptokenvalue);
          switch (veglutid) {
              case 1 : 
                  nonforestpotentialtree[precipannindex][tempaverageindex] = tempfloatvalue;
                  break;
              case 2 : 
                  nonforestpotentialherb[precipannindex][tempaverageindex] = tempfloatvalue;
                  break;
              case 3 : 
                  nonforestpotentialbare[precipannindex][tempaverageindex] = tempfloatvalue;
                  break;
          }
          temptokenvalue = strtok(NULL," ");
      }
  }

  return 0;

}


int printpotentialveglut(int veglutid) {

  int precipannindex, tempaverageindex;
  int precipannindexcalc, tempaverageindexcalc;
  float precipannvalue, tempaveragevalue;
  float tempfloatvalue;


  for (precipannindex = 0; precipannindex < PRECIPCLASSES; precipannindex++) {
      precipannvalue = PRECIPCLASSSTART + ((float) precipannindex) * PRECIPCLASSSIZE;
      for (tempaverageindex = 0; tempaverageindex < TEMPCLASSES; tempaverageindex++) {
          tempaveragevalue = TEMPCLASSSTART + ((float) tempaverageindex) * TEMPCLASSSIZE;
          precipannindexcalc = (int) ((precipannvalue - PRECIPCLASSSTART) / PRECIPCLASSSIZE);
          if (precipannindexcalc < 0) {
              precipannindexcalc = 0;
          }
          if (precipannindexcalc >= PRECIPCLASSES) {
              precipannindexcalc = PRECIPCLASSES - 1;
          }
          tempaverageindexcalc = (int) ((tempaveragevalue - TEMPCLASSSTART) / TEMPCLASSSIZE);
          if (tempaverageindexcalc < 0) {
	      tempaverageindexcalc = 0;
          }
          if (tempaverageindexcalc >= TEMPCLASSES) {
              tempaverageindexcalc = TEMPCLASSES - 1;
          }
          switch (veglutid) {
              case 1 : 
                  printf("%f %f - %f\n",precipannvalue,tempaveragevalue,nonforestpotentialtree[precipannindexcalc][tempaverageindexcalc]);
                  break;
              case 2 : 
                  printf("%f %f - %f\n",precipannvalue,tempaveragevalue,nonforestpotentialherb[precipannindexcalc][tempaverageindexcalc]);
                  break;
              case 3 : 
                  printf("%f %f - %f\n",precipannvalue,tempaveragevalue,nonforestpotentialbare[precipannindexcalc][tempaverageindexcalc]);
                  break;
          }
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
          for (pftid = 0; pftid < MAXPFT; pftid++) {
              outPOTVEGPCTPFTGrid[pftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
	  }
          for (cftid = 0; cftid < MAXCFT; cftid++) {
              outPOTVEGPCTCFTGrid[cftid][clm5lin * MAXCLMPIX + clm5pix] = 0.0;
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
	      if (pftid >= 1 && pftid <= 8) {
	          inclm5PCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
              }
	      if (pftid >= 9 && pftid <= 14) {
	          inclm5PCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] += tempGrid[clm5lin * MAXCLMPIX + clm5pix];
              }
	      if (pftid == 0) {
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


int genLUHPotentialVegGrids() {

  long clm5lin, clm5pix;
  long pftid, cftid;
  float landmask, landfrac, area, precipann, tempaverage;
  float potvegtreepct, potvegherbpct, potvegbarepct;
  float currenttreepct, currentherbpct, currentbarepct;
  float climatemask;
  float newpotveg, availpotveg, removepotveg, addpotveg;
  int precipannindex, tempaverageindex;

  for (clm5lin = 0; clm5lin < MAXCLMLIN; clm5lin++) {
      for (clm5pix = 0; clm5pix < MAXCLMPIX; clm5pix++) {
          landmask = inLANDMASKGrid[clm5lin * MAXCLMPIX + clm5pix];
          landfrac = inLANDFRACGrid[clm5lin * MAXCLMPIX + clm5pix];
	  if (landfrac > 0.0) {
              precipann = inPRECIPANNGrid[clm5lin * MAXCLMPIX + clm5pix];
              tempaverage = inTEMPAVGGrid[clm5lin * MAXCLMPIX + clm5pix];
              precipannindex = (int) ((precipann - PRECIPCLASSSTART) / PRECIPCLASSSIZE);
              if (precipannindex < 0) {
                  precipannindex = 0;
              }
              if (precipannindex >= PRECIPCLASSES) {
                  precipannindex = PRECIPCLASSES - 1;
              }
              tempaverageindex = (int) ((tempaverage - TEMPCLASSSTART) / TEMPCLASSSIZE);
	      if (tempaverageindex < 0) {
	          tempaverageindex = 0;
              }
              if (tempaverageindex >= TEMPCLASSES) {
                  tempaverageindex = TEMPCLASSES - 1;
              }
	      potvegtreepct = nonforestpotentialtree[precipannindex][tempaverageindex];
	      potvegherbpct = nonforestpotentialherb[precipannindex][tempaverageindex];
	      potvegbarepct = nonforestpotentialbare[precipannindex][tempaverageindex];
	      
	      if (potvegtreepct < 0.0) {
	          potvegtreepct = 0.0;
              }
	      
	      if (potvegtreepct > 100.0) {
	          potvegtreepct = 100.0;
              }
	      
	      if (potvegherbpct < 0.0) {
	          potvegherbpct = 0.0;
              }
	      
	      if (potvegherbpct > 100.0) {
	          potvegherbpct = 100.0;
              }
	      
	      if (potvegbarepct < 0.0) {
	          potvegbarepct = 0.0;
              }
	      
	      if (potvegbarepct > 100.0) {
	          potvegbarepct = 100.0;
              }
	      
	      currenttreepct = inclm5PCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix];
	      currentherbpct = inclm5PCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix];
	      currentbarepct = inclm5PCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix];
	      
	      if (currenttreepct < 0.0) {
	          currenttreepct = 0.0;
              }
	      
	      if (currenttreepct > 100.0) {
	          currenttreepct = 100.0;
              }
	      
	      if (currentherbpct < 0.0) {
	          currentherbpct = 0.0;
              }
	      
	      if (currentherbpct > 100.0) {
	          currentherbpct = 100.0;
              }
	      
	      if (currentbarepct < 0.0) {
	          currentbarepct = 0.0;
              }
	      
	      if (currentbarepct > 100.0) {
	          currentbarepct = 100.0;
              }
	      
	      climatemask = precipann - 8.0 * (12.0 + tempaverage);
	      if (climatemask < 0.0 || tempaverage < 2.5 || currentbarepct > 30.0) {
	          if (potvegtreepct > currenttreepct) {
                      potvegtreepct = currenttreepct;
                      potvegherbpct = currentherbpct;
                      potvegbarepct = currentbarepct;
		  }
              }
	      
	      if (currentbarepct > potvegbarepct) {
	          removepotveg = currentbarepct - potvegbarepct;
		  availpotveg = potvegtreepct + potvegherbpct;
		  newpotveg = availpotveg - removepotveg;
		  if (availpotveg > 0.0) {
		      potvegtreepct = newpotveg / availpotveg * potvegtreepct;
		      potvegherbpct = newpotveg / availpotveg * potvegherbpct;
		      potvegbarepct = currentbarepct;
		  }
		  else {
		      potvegtreepct = 0.0;
		      potvegherbpct = 0.0;
		      potvegbarepct = currentbarepct;
                  }		      
              }
	      
	      if (currentbarepct < potvegbarepct) {
	          addpotveg = potvegbarepct - currentbarepct;
		  availpotveg = potvegtreepct + potvegherbpct;
		  newpotveg = availpotveg + addpotveg;
		  if (availpotveg > 0.0) {
		      potvegtreepct = newpotveg / availpotveg * potvegtreepct;
		      potvegherbpct = newpotveg / availpotveg * potvegherbpct;
		      potvegbarepct = currentbarepct;
		  }
		  else {
		      potvegtreepct += 0.5 * newpotveg;
		      potvegherbpct += 0.5 * newpotveg;
		      potvegbarepct = currentbarepct;
                  }		      
              }
	      
	      if (potvegtreepct < 0.0) {
	          potvegtreepct = 0.0;
              }
	      
	      if (potvegtreepct > 100.0) {
	          potvegtreepct = 100.0;
              }
	      
	      if (potvegherbpct < 0.0) {
	          potvegherbpct = 0.0;
              }
	      
	      if (potvegherbpct > 100.0) {
	          potvegherbpct = 100.0;
              }
	      
	      if (potvegbarepct < 0.0) {
	          potvegbarepct = 0.0;
              }
	      
	      if (potvegbarepct > 100.0) {
	          potvegbarepct = 100.0;
              }
	      
              outPOTVEGPCTTREEGrid[clm5lin * MAXCLMPIX + clm5pix] = potvegtreepct;
              outPOTVEGPCTHERBGrid[clm5lin * MAXCLMPIX + clm5pix] = potvegherbpct;
	      outPOTVEGPCTBAREGrid[clm5lin * MAXCLMPIX + clm5pix] = potvegbarepct;
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
  
  sprintf(globalbinfilename,"%s/%s/%s.PCTBARE.%s.dat",databasestr,seriesname,seriesname,yearname);
  printf("Writing: %s\n",globalbinfilename);
  globalbinfile = fopen(globalbinfilename,"w+");
  fwrite(outPOTVEGPCTBAREGrid,sizeof(outPOTVEGPCTBAREGrid),1,globalbinfile);  
  fclose(globalbinfile);
  
  return 0;

}


int main(long narg, char **argv) {

  char workrawinputdirname[1024];
  char workrawoutputdirname[1024];
  char workrawinputreferenceyearstr[256];
  char workrawoutputreferenceyearstr[256];
  
  if(narg != 8){
        printf("Usage createpotvegpftsclm552Deg025bin timeseriesrawfile timeseriesrawdir workrawfile workrawdir luh3modiscruclimatefile luh3modiscrurawdir generatenonforestpotveginitialvcfnamelist\n");
        return 0;
  }
  
  readnamelistfile(argv[7]);

  readclm5modisrawinputlutfile(argv[1],argv[2],namelistclm5modisrawinput);
  readworkrawinputlutfile(argv[3],argv[4],namelistworkrawinput);
  readworkrawoutputlutfile(argv[3],argv[4],namelistworkrawoutput);
  readnonforesttreeclimatelutfile(argv[5],argv[6],namelistnonforesttreeclimatename);
  readnonforestherbclimatelutfile(argv[5],argv[6],namelistnonforestherbclimatename);
  readnonforestbareclimatelutfile(argv[5],argv[6],namelistnonforestbareclimatename);
  workrawinputreferenceyear = atoi(namelistworkrawinputreferenceyear);
  sprintf(workrawinputreferenceyearstr,"%04d",workrawinputreferenceyear);
  
  readpotentialveglut(nonforesttreeclimatedb,1);
  readpotentialveglut(nonforestherbclimatedb,2);
  readpotentialveglut(nonforestbareclimatedb,3);

  initializeGrids();
  
  sprintf(workrawinputdirname,"%s/%s",workrawinputdb,workrawinputname);
  readlandGrids(workrawinputdirname,namelistworkrawinputname,workrawinputreferenceyearstr);
  readclm5referenceGrids(workrawinputdirname,namelistworkrawinputname,workrawinputreferenceyearstr);    
  readclimvarseries(clm5modisrawinputdb,clm5modisrawinputname,"CLIM");
  readlaivarseries(clm5modisrawinputdb,clm5modisrawinputname,"FILLCLIM");
  
  genLUHLandGrids();
  genLUHPotentialVegGrids();
  
  sprintf(workrawoutputreferenceyearstr,"%s",workrawinputreferenceyearstr);
  sprintf(workrawoutputdirname,"%s/%s",workrawoutputdb,workrawoutputname);
  writepotveggrids(workrawoutputdirname,namelistworkrawoutputname,workrawoutputreferenceyearstr);

  return 1;
  
}
