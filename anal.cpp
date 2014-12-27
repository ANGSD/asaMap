#include <cmath>
#include "readplink.c"
#include "anal.h"
#include "anders.h"
void wrap(const plink *plnk,const std::vector<double> &phe,const std::vector<double> &ad, const std::vector<double> &sex,const std::vector<double> &freq){
  //  fprintf(stderr,"\t-> plinkdim: x->%lu y->%lu\n",plnk->x,plnk->y);
  //return ;
  char **d = new char*[plnk->y];//transposed of plink->d. Not sure what is best, if we filterout nonmissing anyway.
  for(int i=0;i<plnk->y;i++){
    d[i] = new char[plnk->x];
    for(int j=0;j<plnk->x;j++)
      d[i][j]=plnk->d[j][i];
  }

  double **res=new double *[plnk->y];
  
  pars *p=init_pars(plnk->x,10);//we prep for threading. By using encapsulating all data need for a site in  struct called pars

  for(int y=0;y<plnk->y;y++){//loop over sites
    int cats[4] = {0,0,0,0};
    int cats2[4] = {0,0,0,0};
    res[y]=new double[47];
    for(int i=0;i<47;i++)
      res[y][i] = NAN; //<- didn't know this...
    
    for(int x=0;x<plnk->x;x++)
      cats[plnk->d[x][y]]++;
    for(int x=0;x<plnk->x;x++)//similar to above but with transposed plink matrix
      cats2[d[y][x]]++;
#if 0
    //print table
    //    fprintf(stdout,"[%d] %d %d %d %d freq:%f ",y,cats[0],cats[1],cats[2],cats[3],freq[y]);
    // fprintf(stdout,"[%d] %d %d %d %d freq:%f\n",y,cats2[0],cats2[1],cats2[2],cats2[3],freq[y]);
#endif
    
    //discard sites if missingness > 0.1
    if(cats[3]/(double)(cats[0]+cats[1],cats[2])>0.1){
      //fprintf(stderr,"skipping site[%d] due to excess missingness\n",y);
      res[y][0] = 1;
      continue;
    }
    
    //check that we observe atleast 10 obs of 2 different genotypes
    int n=0;
    for(int i=0;i<4;i++)
      if(cats[i]>10)
	n++;
    if(n<2){
      //fprintf(stderr,"skipping site[%d] due to categories filter\n",y);
      res[y][0]=2;
      continue;
    }
    //plug in freq
    res[y][1]=freq[y];res[y][2]=1-freq[y];

    if(freq[y]>0.999||freq[y]<0.001){
      //fprintf(stderr,"skipping site[%d] due to maf filter\n",y);
      res[y][0]=3;
      continue;
      
    }
    //set the pars, including the result vector
    set_pars(p,d[y],phe,ad,sex,freq);
    p->res=res[y];
    //find fit.
    getfit(p);

    //print it
    for(int i=0;i<47;i++)
      fprintf(stderr,"%f\t",res[y][i]);
    fprintf(stderr,"\n");
    
  }
  fprintf(stderr,"\t-> done\n");
  kill_pars(p,plnk->x);
  for(int i=0;i<plnk->y;i++)
    delete [] res[i];

  delete [] res;
  for(int i=0;i<plnk->y;i++)
    delete [] d[i];
  delete [] d;
  

}

pars *init_pars(int l,int ncov){
  pars *p=new pars;
  p->len=l;
  p->ncov=ncov;
  p->gs=new char[l];
  p->ys=new double[l];
  p->qvec=new double[l];
  p->mafs=new double[l];
  p->covs=new double*[l];
  p->sex=new double[l];
  for(int i=0;i<l;i++)
    p->covs[i]=new double[p->ncov];

  return p;
}

void kill_pars(pars *p,int l){
  delete [] p->gs;
  delete [] p->ys;
  delete [] p->qvec;
  delete [] p->mafs;
  delete [] p->sex;
  for(int i=0;i<l;i++)
    delete [] p->covs[i];
  delete [] p->covs;
  delete p;
  p=NULL;
}


void set_pars(pars*p,char *g,const std::vector<double> &phe,const std::vector<double> &ad, const std::vector<double> &sex,const std::vector<double> &freq){
  p->len=0;
  for(int i=0;i<phe.size();i++){
    if(g[i]!=3){
      p->gs[p->len] = g[i];
      p->ys[p->len] = phe[i];
      p->qvec[p->len] = ad[i];
      p->sex[p->len] = sex[i];
      p->mafs[p->len] = freq[i];
      p->len++;
    }
  }
  
}
