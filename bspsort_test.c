#include "bspedupack.h"
#include <sys/time.h>

/*  This is a test program which uses bspsort to sort an array
    of randomly generated doubles.
*/

long P; // number of processors requested
long N; // input length of array


int compare_doubles (const void *a, const void *b){
    /* Compares two doubles, for use in qsort or mergesort */

    if ( *(double *)a < *(double *)b )
        return -1;
    if ( *(double *)a > *(double *)b )
        return 1;
    return 0; // equality

} /* end compare_doubles */


int compare_items (const void *a, const void *b){
    /* Compares two items, for use in qsort or mergesort */

    Item A= *(Item *)a;
    Item B= *(Item *)b;

    if (A.weight < B.weight)
        return -1;
    if (A.weight > B.weight)
        return 1;
    if (A.index < B.index)
        return -1;
    if (A.index > B.index)
        return 1;
    return 0; //equality for both criteria

} /* end compare_items*/


void merge(char *x, char *tmp, long a, long b, long c, size_t size,
           int (*compare)(const void *,const void *)){

   /* This function merges two sorted parts X[a,b-1] and X[b,c]
      into one sorted part X[a,b], where a <= b <= c.
      Here, the indexing is such that X[i] = x[i*size],
      and X[i] is a pointer to the i-th item.
      The sort is by increasing order detemined by the comparison function
      compare which gives -1, 0, or 1 as output (similar to the one of the
      C function qsort).

      The function uses a temporary array of length >= c-a+1 to avoid
      repeated memory allocations.

      size is the size of the items to be sorted, e.g. sizeof(double).
   */

      if (a >= b || b > c) // one of the parts is empty
          return;

      long i= a*size; // index for range [a,b-1]
      long j= b*size; // index for range [b,c]
      long k= 0; // index for tmp
     
      /* Compare values as long as no part is empty */
      while (i < b*size && j <= c*size){
          if (compare(x+i, x+j) < 0 ){ // compares x[i] and x[j]
              memcpy (tmp+k, x+i, size);
              i += size;
          } else {
              memcpy (tmp+k, x+j, size); 
              j += size;
          }
          k += size;
      }

      /* Copy the remaining values */ 
      if (i >= b*size)
          memcpy (tmp+k, x+j, (c+1)*size - j);
      if (j > c*size)
          memcpy (tmp+k, x+i, b*size - i );

      /* Copy the values back into x */
      memcpy (x+(a*size), tmp, (c-a+1)*size);

} /* end merge */


void mergeparts(char *x, long *start, long p, size_t size,
                int (*compare)(const void *,const void *)){

    /* This function sorts the array x for
       the index range start[0] <= i < start[p].

       On input, subranges start[k] <= i < start[k+1] are assumed
       to have been sorted already, for k = 0,1,...,p-1.

      The sort is by increasing order detemined by the comparison function
      compare which gives -1, 0, or 1 as output (similar to the one of the
      C function qsort).
    */

    /* Initialize current data */
    long nparts= p; // current number of parts
    long *cstart = vecalloci(p+1); // current start
    for (long k=0; k<=p; k++)
        cstart[k]= start[k];
    
    /* Allocate sufficient temporary workspace */
    char *tmp = malloc( (start[p] - start[0])*size );
    if (tmp==NULL)
        bsp_abort("malloc: not enough memory");

    while (nparts > 1){
        
        /* Merge pairs of parts */
        for (long i=0; i<nparts/2; i++){
            merge(x,tmp,cstart[2*i],cstart[2*i+1],cstart[2*i+2]-1,
                  size,compare);
        }

        /* Adjust the starts and the number of parts */
        for (long i=0; i<nparts/2; i++)
            cstart[i] = cstart[2*i];
        if (nparts%2==1) // remaining unmerged part
            cstart[nparts/2]= cstart[nparts-1];
        nparts= (nparts+1)/2;
        cstart[nparts]= start[p];
    }

    free(tmp);
    vecfreei(cstart);

} /* end mergeparts */



void bspsort_test(){
    
    void bspsort(double *x, long n, long *nlout);
    bsp_begin(P);
    long p= bsp_nprocs();
    long s= bsp_pid();

    /* Obtain input length n */
    long n;
    if (s==0)
        n = N;
    
    bsp_push_reg(&n,sizeof(long));
    bsp_sync();

    bsp_get(0,&n,0,&n,sizeof(long));
    bsp_sync();
    bsp_pop_reg(&n);

    /* Set seed for random number generator, different for every processor
       and for every run, since it is based on the current time.  */
    struct timeval t1;
    gettimeofday(&t1, NULL);
    srand(t1.tv_usec * t1.tv_sec*s);

    /* Initialize array of doubles */ 
    long nl = nloc(p,s,n); // number of local elements 
    long nl0 = nloc(p,0,n); // maximum number of local elements of a processor
    double *x= vecallocd(2*nl0+p);
    bsp_push_reg(x, (2*nl0+p)*sizeof(double));

    /* Generate random number in range [0,1] as input */
    for (long i=0; i<nl; i++)
        x[i]= (double)rand() / (double)RAND_MAX;

    bsp_sync(); 
    double time0= bsp_time();

    long nlout= 0;
    bspsort(x,n,&nlout);
    bsp_sync();  
    double time1= bsp_time();

    /* Check if the values of x are indeed in increasing order */
    for (long i=0; i<nlout-1; i++){
        if (x[i] > x[i+1])
            bsp_abort("Proc %ld: output is not sorted correctly \n",s);
    }

    printf("Proc %ld: number of elements = %ld \n",s,nlout);
    fflush(stdout);
    if (nlout>0){
        printf("Proc %ld: first element = %.6lf last element = %.6lf \n",
                s,x[0],x[nlout-1]); fflush(stdout);
    }
    if (s==0)
        printf("This took only %.6lf seconds.\n", time1-time0);

    bsp_pop_reg(x);
    vecfreed(x);
    bsp_end();

} /* end bspsort_test */


int main(int argc, char **argv){

    bsp_init(bspsort_test, argc, argv);

    /* Sequential part */
    printf("How many processors do you want to use?\n"); fflush(stdout);
    scanf("%ld",&P);
    printf("Please enter array length N:\n"); fflush(stdout);
    scanf("%ld",&N);

    if (N < 1 || P < 1){
        printf("Error: input values  N or P are invalid.\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    /* Reduce P if needed for sample sort */
    while (P*P > N)
        P--;
    if(P > bsp_nprocs())
        P=  bsp_nprocs();
    printf("Using %ld processors\n", P); fflush(stdout);
        
    /* SPMD part */
    bspsort_test();

    /* Sequential part */
    exit(EXIT_SUCCESS);

} /* end main */
