/*
Lines project:
run simple scenario like:
7 samples times 5 snps
 ./line -p /home/thorfinn/line/data/plink -o tmp -c /home/thorfinn/line/data/sub.eig -s /home/thorfinn/line/data/sub.pheno -a /home/thorfinn/line/data/ad.sub -6 /home/thorfinn/line/data/sub.6 -f /home/thorfinn/line/data/adfreq.sub
 
 run all samples all sites like
./line -p /home/thorfinn/line/data/datplus_QCed_newIDs_IHIT_allSNPsPlus7xtraPGLU120sub -o tmp -c /home/thorfinn/line/data/datplus_QCed_newIDs_IHIT_allSNPsPlus7xtraPGLU120sub.myeigenvec  -s /home/thorfinn/line/data/pheno -a /home/thorfinn/line/data/admixProper.txt  -6 /home/thorfinn/line/data/6.txt -f /home/thorfinn/line/data/admixFreqs.txt

 */



#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "readplink.h"
#include "analysisFunction.h"
#include "anal.h"

void print_info(FILE *fp){
  fprintf(fp, "\n");
  fprintf(fp, "Usage: line  [options] \n");
  fprintf(fp, "Options:\n");
  fprintf(fp, "   -p <filename>       plink prefix filename\n");
  fprintf(fp, "   -o <filename>       outputfilename\n");
  fprintf(fp, "   -c <filename>       covariance matrix filename\n");
  fprintf(fp, "   -s <filename>       phenotypes\n");
  fprintf(fp, "   -a <filename>       admixprop (admixprop for source pop1)\n");
  fprintf(fp, "   -6 <filename>       the sex of the samples\n");
  fprintf(fp, "   -f <filename>       freqs\n");
  fprintf(fp, "\n");
}


int main(int argc,char **argv){
  char *pname = NULL;
  char *outname = NULL;
  char *covname = NULL;
  char *phename = NULL;
  char *adname = NULL;
  char *sexname = NULL;
  char *freqname=NULL;
  int n;
  while ((n = getopt(argc, argv, "p:o:c:s:a:6:f:")) >= 0) {
    switch (n) {
    case 'p': pname = strdup(optarg); break; 
    case 'o': outname = strdup(optarg); break;
    case 'c': covname = strdup(optarg); break;
    case 's': phename = strdup(optarg); break;
    case 'a': adname = strdup(optarg); break;
    case '6': sexname = strdup(optarg); break;
    case 'f': freqname = strdup(optarg); break;
    }
  }
  plink *p=readplink(pname);

#if 1
  if(!outname||!pname||!covname||!phename||!adname||!sexname||!freqname){
    print_info(stderr);
    return 1;
  }  
#endif
  
  Matrix<double> eig =getMatrix(covname);
  std::vector<double> pheno=getArray(phename);
  std::vector<double> adprop =getArray(adname);
  std::vector<double> sex =getArray(sexname);
  std::vector<double> f =getArray(freqname);
  

  wrap(p,pheno,adprop,sex,f);

  //cleanup
  kill_plink(p); 
  free(pname);
  free(outname);
  free(covname);
  free(phename);
  free(adname);
  free(sexname);
  free(freqname);
  for(int i=0;i<eig.x;i++)
    delete [] eig.matrix[i];
  delete [] eig.matrix;
  return 0;
}
