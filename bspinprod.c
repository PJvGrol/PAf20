#include "bspedupack.h"

/*  This program computes the sum of the first n squares,
        sum = 1*1 + 2*2 + ... + n*n
    by computing the inner product of x=(1,2,...,n)^T and
    itself, for n>=0.
    The output should equal n*(n+1)*(2n+1)/6.
    The distribution of x is cyclic.
*/

long P; // number of processors requested

double bspip(long n, double *x, double *y){
    /* Compute inner product of vectors x and y
       of length n>=0 */

    long p= bsp_nprocs(); // p = number of processors obtained
    long s= bsp_pid();    // s = processor number

    double *Inprod= vecallocd(p);
    bsp_push_reg(Inprod,p*sizeof(double));
    bsp_sync();

    double inprod= 0.0;
    for (long i=0; i<nloc(p,s,n); i++)
        inprod += x[i]*y[i];
    
    for (long t=0; t<p; t++)
        bsp_put(t,&inprod,Inprod,s*sizeof(double),sizeof(double));
    
    bsp_sync();

    double alpha= 0.0;
    for (long t=0; t<p; t++)
        alpha += Inprod[t];
    
    bsp_pop_reg(Inprod);
    vecfreed(Inprod);

    return alpha;

} /* end bspip */

void bspinprod(){
    
    bsp_begin(P);
    long p= bsp_nprocs();
    long s= bsp_pid();
    long n;
    if (s==0){
        printf("Please enter n:\n"); fflush(stdout);
        scanf("%ld",&n);
        if(n<0)
            bsp_abort("Error in input: n is negative");
    }
    bsp_push_reg(&n,sizeof(long));
    bsp_sync();

    bsp_get(0,&n,0,&n,sizeof(long));
    bsp_sync();
    bsp_pop_reg(&n);

    long nl= nloc(p,s,n);
    double *x= vecallocd(nl);
    for (long i=0; i<nl; i++){
        long iglob= i*p+s;
        x[i]= iglob+1;
    }
    bsp_sync(); 
    double time0= bsp_time();

    double alpha= bspip(n,x,x);
    bsp_sync();  
    double time1= bsp_time();

    printf("Proc %ld: sum of squares up to %ld*%ld is %.lf\n",
            s,n,n,alpha); fflush(stdout);
    if (s==0){
        printf("This took only %.6lf seconds.\n", time1-time0);
        
        /* Compute exact output (number can become large) */
        long long sum= (n*(n+1)*(2*n+1)) / (long long)6;
        printf("n(n+1)(2n+1)/6 = %lld\n", sum); fflush(stdout);
    }

    vecfreed(x);
    bsp_end();

} /* end bspinprod */

int main(int argc, char **argv){

    bsp_init(bspinprod, argc, argv);

    /* Sequential part */
    printf("How many processors do you want to use?\n");
    fflush(stdout);

    scanf("%ld",&P);
    if (P > bsp_nprocs()){
        printf("Sorry, only %u processors available.\n",
                bsp_nprocs());
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    /* SPMD part */
    bspinprod();

    /* Sequential part */
    exit(EXIT_SUCCESS);

} /* end main */
