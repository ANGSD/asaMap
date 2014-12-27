#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "readplink.h"


/*
  output will be genotype[ind][nsites]
  g \in {0,1,2,3},counts of allele2, 3=NA/missing
*/
//modified from snpMatrix by clayton and hintak leung 2007
unsigned char **readbed(const char* file, int nrow,int ncol) {
  int i;
  const unsigned char recode[4] = {'\x01', '\x00', '\x02', '\x03'};
  const unsigned char mask = '\x03';


  FILE *in = fopen(file, "r");
  if (!in){
    fprintf(stderr,"Couldn't open input file: %s\n", file);
    exit(0);
  }
  unsigned char start[3];
  if (fread(start, 1, 3, in)!=3){
    fprintf(stderr,"Failed to read first 3 bytes");
    exit(0);
  }
  if (start[0]!='\x6C' || start[1]!='\x1B'){
    fprintf(stderr,"Input file does not appear to be a .bed file (%X, %X)", 
	   start[0], start[1]);
    exit(0);
  }
  /* Create output object */

  
  unsigned char** results =(unsigned char**) calloc(nrow,sizeof(unsigned char*));
  for(i=0;i<nrow;i++)
    results[i] =(unsigned char*) calloc(ncol,sizeof(unsigned char));

  /* Read in data */

  int snp_major = start[2];
  int part=0, ij=0, j=0;i=0;
  while (1) {
    unsigned char byte;
    if (!part) {
      if (feof(in) || !fread(&byte, 1, 1, in)) {
	printf("Unexpected end of file reached");
	exit(0);
      }
      part = 4;
    }
    unsigned char code = byte & mask;
    byte = byte >> 2;
    part--;
    unsigned char tmp = recode[code];
    if(tmp==0)
      results[i][j] = 3;
    else
      results[i][j] =tmp-1;

    assert(results[i][j]>=0 && results[i][j]<=3);

    if (snp_major) {
      ij++;
      i++;
      if (i==nrow) {
	i = part = 0;
	j++;
	if (j==ncol)
	  break;
      }
    }	
    else {
      ij += nrow;
      j++;
      if (j==ncol){
	j = part = 0;
	i++;
	if (i==nrow)
	  break;
	ij = i;
      }
    }
  }
  fclose(in);
  return results;
}

int nlines(const char *fname){
  FILE *in =NULL;
  if(!((in=fopen(fname,"r")))){
      fprintf(stderr,"Problem opening file: %s\n",fname);
      return 0;
  }
  int c;
  int n=0;
  while ( (c=fgetc(in)) != EOF ) 
    if ( c == '\n' )  n++;
  fclose(in);
  return n;
}

void kill_plink(plink *p){
  int i;
  for(i=0;i<p->x;i++)
    free(p->d[i]);
  
  free(p->d);
  free(p);
  p=NULL;
}

plink *readplink(char *str){
  plink *p =(plink*) malloc(sizeof( plink));

  char *fname =(char*) malloc(strlen(str)+5);
  strcpy(fname,str);strcpy(fname+strlen(fname),".bim");
  int ncol = nlines(fname);
  
  strcpy(fname+strlen(str),".fam");
  int nrow = nlines(fname);
  //  nrow=3;
  if(ncol==0||nrow==0)
    return 0;
  fprintf(stderr,"Done reading file: \'%s\' with dim ncol:%d\tnrow:%d\n",str,ncol,nrow);
  strcpy(fname+strlen(str),".bed");
  unsigned char **dat = readbed(fname,nrow,ncol);
  p->x=nrow;
  p->y=ncol;
  p->d=dat;
  free(fname);
  return p;
}


#ifdef __WITH_MAIN__

int main(int argc,char **argv){
  int i,j;
  if(argc==1){
    fprintf(stderr,"Supply prefix for plink binary files\n");
    return 0;
  }
  char *fname =(char*) malloc(strlen(argv[1])+5);
  strcpy(fname,argv[1]);strcpy(fname+strlen(fname),".bim");
  int ncol = nlines(fname);
  
  strcpy(fname+strlen(argv[1]),".fam");
  int nrow = nlines(fname);

  if(ncol==0||nrow==0)
    return 0;
  fprintf(stderr,"ncol:%d\tnrow:%d\n",ncol,nrow);
  strcpy(fname+strlen(argv[1]),".bed");
  unsigned char **dat = readbed(fname,nrow,ncol);
#if 1
  for(i=0;i<ncol;i++){
    for(j=0;j<nrow;j++)
      fprintf(stdout,"%d ",dat[j][i]);
    fprintf(stdout,"\n");
  }
#endif

  for(i=0;i<nrow;i++)
    free(dat[i]);
  free(dat);
  free(fname);
  return 0;
}
#endif
