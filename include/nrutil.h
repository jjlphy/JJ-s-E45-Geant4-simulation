#ifndef nrutil
#define nrutil


void nrerror(char error_text[]);
int *ivector(long nl, long nh);
void free_ivector(int *v, long nl, long nh);


#endif
