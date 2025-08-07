#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
{

  fprintf(stderr,"1");
  fprintf(stderr,"2");
  fprintf(stderr,"3");
  exit(1);

}

int *ivector(long nl, long nh)
{
  int *v;
  
  v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if(!v) nrerror("");
  return v-nl+NR_END;

}

void free_ivector(int *v, long nl, long nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}
