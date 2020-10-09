#include "bspedupack.h"

long nloc(long p, long s, long n){
    /* Compute number of local components of processor s for vector
       of length n distributed cyclically over p processors.
       Also useful for block distribution with
            ceil(n/p)  components if s < (n mod p),
            floor(n/p) components otherwise.
    */

    return  (n+p-s-1)/p;

} /* end nloc */


/* The following functions can be used to allocate and deallocate
   vectors and matrices.  If not enough memory is available,
   one processor halts them all.  */

long *vecalloci(size_t n){
    /* This function allocates a vector of long
       integers of length n */
    long *pi;

    pi= malloc(MAX(n,1)*sizeof(long));
    if (pi==NULL)
        bsp_abort("vecalloci: not enough memory");
    return pi;

} /* end vecalloci */


bool *vecallocb(size_t n){
    /* This function allocates a vector of booleans of length n */
    bool *pb;

    pb= malloc(MAX(n,1)*sizeof(bool));
    if (pb==NULL)
        bsp_abort("vecallocb: not enough memory");
    return pb;

} /* end vecallocb */

double *vecallocd(size_t n){
    /* This function allocates a vector of doubles of length n */
    double *pd;

    pd= malloc(MAX(n,1)*sizeof(double));
    if (pd==NULL)
        bsp_abort("vecallocd: not enough memory");
    return pd;

} /* end vecallocd */

double complex *vecallocc(size_t n){
    /* This function allocates a vector of complex numbers
       of length n */
    double complex *pc;

    pc= malloc(MAX(n,1)*sizeof(double complex));
    if (pc==NULL)
        bsp_abort("vecallocc: not enough memory");
    return pc;

} /* end vecallocc */

Item *vecallocitem(size_t n){
    /* This function allocates a vector of items of length n */
    Item *pitem;

    pitem= malloc(MAX(n,1)*sizeof(Item));
    if (pitem==NULL)
        bsp_abort("vecallocitem: not enough memory");
    return pitem;

} /* end vecallocitem */

double **matallocd(size_t m, size_t n){
    /* This function allocates an m x n matrix of doubles */
    size_t i, m1, n1;
    double *pd, **ppd;

    m1= MAX(m,1);
    n1= MAX(n,1);
    ppd= malloc(m1*sizeof(double *));
    if (ppd==NULL)
        bsp_abort("matallocd: not enough memory");
    pd= malloc(m1*n1*sizeof(double)); 
    if (pd==NULL)
        bsp_abort("matallocd: not enough memory");
    ppd[0]= pd;
    for (i=1; i<m1; i++)
        ppd[i]= ppd[i-1]+n1;
     
    return ppd;

} /* end matallocd */

void vecfreei(long *pi){
    /* This function frees a vector of long integers */

    if (pi!=NULL)
        free(pi);

} /* end vecfreei */

void vecfreeb(bool *pb){
    /* This function frees a vector of booleans */

    if (pb!=NULL)
        free(pb);

} /* end vecfreeb */

void vecfreed(double *pd){
    /* This function frees a vector of doubles */

    if (pd!=NULL)
        free(pd);

} /* end vecfreed */

void vecfreec(double complex *pc){
    /* This function frees a vector of complex numbers */
    
    if (pc!=NULL)
        free(pc);
        
} /* end vecfreec */

void vecfreeitem(Item *pitem){
    /* This function frees a vector of items */

    if (pitem!=NULL)
        free(pitem);

} /* end vecfreeitem */


void matfreed(double **ppd){
    /* This function frees a matrix of doubles */

    if (ppd!=NULL){
        if (ppd[0]!=NULL)
            free(ppd[0]);
        free(ppd);
    }

} /* end matfreed */
