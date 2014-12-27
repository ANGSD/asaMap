#pragma once
#include <stdlib.h>
typedef struct{
  size_t x;
  size_t y;
  unsigned char **d;
}plink;

plink *readplink(char *str);
void kill_plink(plink *p);
