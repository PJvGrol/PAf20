#include "bspedupack.h"

/*  This is a test program which uses bspfft to perform 
    a Fast Fourier Transform and its inverse.

    The input vector is defined by x[j]=j+i, for 0 <= j < n.
    Here i= sqrt(-1).
 
    The output vector should equal the input vector,
    up to roundoff errors. Output is by triples (j, Re x[j], Im x[j]).
    Warning: don't rely on this test alone to check correctness. 
    (After all, deleting the main loop will give similar results ;) 

*/

#define NITERS 100 // Perform NITERS forward and backward transforms.
                   // A large NITERS helps to obtain accurate timings.
#define NPRINT 10  // Print at most NPRINT values per processor
#define MEGA 1000000.0

long P; // number of processors requested

void bspfft_test(){
    void bspfft(double complex *x, long n, bool forward, double complex *w,
                long *rho_np, long *rho_p);
    void bspfft_init(long n, double complex *w, long *rho_np, long *rho_p);

    bsp_begin(P);
    long p= bsp_nprocs();
    long s= bsp_pid();
    long n;

    bsp_push_reg(&n,sizeof(long));
    double *Error= vecallocd(p);
    bsp_push_reg(Error,p*sizeof(double));
    bsp_sync();

    if (s==0){
        printf("Please enter length n: \n");
        scanf("%ld",&n);
        if(n<2*p)
            bsp_abort("Error in input: n < 2p");
        for (long q=0; q<p; q++)
            bsp_put(q,&n,&n,0,sizeof(long));

        printf("FFT of vector of length %ld using %ld processors\n",n,p);
        printf("performing %d forward and %d backward transforms\n",
                NITERS,NITERS);
    }
    bsp_sync();

    /* Determine the number of computation supersteps, by computing
       the smallest integer t such that (n/p)^t >= p */

    long np= n/p;
    long t= 0;
    for (long c=1; c<p; c *= np)
        t++;

    /* Allocate, register,  and initialize vectors */
    double complex *w=  vecallocc((t+1)*np);
    double complex *x= vecallocc(np);
    bsp_push_reg(x,np*sizeof(double complex));
    long *rho_np= vecalloci(np);
    long *rho_p=  vecalloci(p);

    for (long j=0; j<np; j++){
        double jglob= j*p + s;
        x[j]= jglob + 1.0*I;
    }
    bsp_sync();
    double time0= bsp_time();


    /* Initialize the weight and bit reversal tables */
    for (long it=0; it<NITERS; it++)
        bspfft_init(n,w,rho_np,rho_p);
    bsp_sync();
    double time1= bsp_time();
  
    /* Perform the FFTs */
    for (long it=0; it<NITERS; it++){
        bspfft(x,n,true,w,rho_np,rho_p);
        bspfft(x,n,false,w,rho_np,rho_p);
    }
    bsp_sync();
    double time2= bsp_time();

    /* Compute the accuracy */
    double max_error= 0.0;
    for (long j=0; j<np; j++){
        double jglob= j*p + s;
        double error_re= fabs(creal(x[j])-jglob);
        double error_im= fabs(cimag(x[j])-1.0);
        double error= sqrt(error_re*error_re + error_im*error_im);
        if (error>max_error)
            max_error= error;
    }
    bsp_put(0,&max_error, Error, s*sizeof(double), sizeof(double));
    bsp_sync();

    if (s==0){
        max_error= 0.0;
        for (long q=0; q<p; q++){
            if (Error[q]>max_error)
                max_error= Error[q];
        }
    }

    printf("The computed solution is (<= %ld values per processor):\n",
           (long)NPRINT);
  
    for (long j=0; j<NPRINT && j<np; j++){
        long jglob= j*p + s;
        printf("Proc %ld: j=%ld Re= %lf Im= %lf \n",s,jglob,creal(x[j]),cimag(x[j]));
    }
    fflush(stdout);
    bsp_sync();

    /* Determine time and computing rate */
    if (s==0){
        printf("Time per initialization = %lf sec \n",
                (time1-time0)/NITERS);
        double ffttime= (time2-time1)/(2.0*NITERS);
        printf("Time per FFT = %lf sec \n", ffttime);

        double nflops= 5*n*log((double)n)/log(2.0) + 2*n;
        printf("Computing rate of FFT = %lf Mflop/s \n",
                nflops/(MEGA*ffttime));
        printf("Absolute error= %e \n", max_error);
        printf("Relative error= %e \n\n", max_error/n);
    }

  
    bsp_pop_reg(x);
    bsp_pop_reg(Error);
    bsp_pop_reg(&n);
    bsp_sync();

    vecfreei(rho_p);
    vecfreei(rho_np);
    vecfreec(x);
    vecfreec(w);
    vecfreed(Error);
    bsp_end();

} /* end bspfft_test */


int main(int argc, char **argv){
 
    bsp_init(bspfft_test, argc, argv);
 
    /* Sequential part */
    printf("How many processors do you want to use?\n"); fflush(stdout);
    scanf("%ld",&P);
    if (P > bsp_nprocs()){
        printf("Sorry, only %u processors available.\n", bsp_nprocs());
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
 
    /* SPMD part */
    bspfft_test();
 
    /* Sequential part */
    exit(EXIT_SUCCESS);
 
} /* end main */
