
#include <cmath>
#include <vector>
#include <zlib.h>
template<typename T>
struct Matrix {
  size_t x;
  size_t y;
  T **matrix;
};

double getMax(double a,double b, double c);
double addProtect2(double a,double b);
double addProtect3(double a,double b, double c);
Matrix<double> getMatrix(const char *name);

std::vector<double> getArray(const char *name);
int fexists(const char* str);
double sigm(double x);

void swapDouble (double& first, double& second);
int matinv( double x[], int n, int m, double space[]);
void deleteMatrix(Matrix<double> mat);
void printMatrix(Matrix<double> mat,FILE *file);
void logrescale(double *ary,int len);
void svd_inverse(double mat[],int xLen, int yLen);
std::vector<char*> getFilenames(const char * name,int nInd);
char *strpop(char **str,char split);

template <typename T>
T * allocArray(size_t len,T defval){
  T *ret= new T[len];
  for(size_t i=0;i<len;i++)
    ret[i]=defval;
  return ret;
  
  
}
template <typename T>
  T * allocArray(size_t len){
  T *ret= new T[len];
  return ret;
  
}

template <typename T>
T sum(const T *ary,size_t len){
  T s =0;
  for(size_t i=0;i<len ; i++)
    s+=ary[i];
  //  printf("sum:%f\n",s);
  return s;
}

void print_array(FILE *fp,double *ary,int len);
void print_array(FILE *fp,int *ary,int len);

double *readDouble(const char*fname,int hint);
int whichMax(double *d,int len);


size_t fsize(const char* fname);
int fexists(const char* str);//{///@param str Filename given as a string.
FILE *openFile(const char* a,const char* b);
gzFile openFileGz(const char* a,const char* b,const char *mode);
FILE *getFILE(const char*fname,const char* mode);
gzFile getGz(const char*fname,const char* mode);

int isNewer(const char *newer,const char *older);



