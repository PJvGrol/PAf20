#include "bspedupack.h"
 
/*  This program measures p, r, g, and l of a BSP computer
    using bsp_put for communication.
*/

#define NITERS 1000    // number of iterations
#define MAXN 1024      // maximum length of DAXPY
#define MAXH 2048      // maximum h in h-relation, must be >> p
#define MEGA 1000000.0

long P; // number of processors requested

void leastsquares(long h0, long h1, double *t,
                  double *g, double *l){
    /* This function computes the parameters g and l of the 
       linear function T(h)= g*h+l that best fits
       the data points (h,t[h]) with h0 <= h <= h1. */

    /* Compute sums:
        sumt  =  sum of t[h] over h0 <= h <= h1
        sumth =         t[h]*h
        nh    =         1
        sumh  =         h
        sumhh =         h*h     */

    double sumt= 0.0;
    double sumth= 0.0;
    for (long h=h0; h<=h1; h++){
        sumt  += t[h];
        sumth += t[h]*h;
    }
    long nh= h1-h0+1;
    double h00 = (double)h0*(h0-1);
    double h11 = (double)h1*(h1+1);
    double sumh= (h11 - h00)/2;  
    double sumhh= (h11*(2*h1+1) - h00*(2*h0-1))/6;

    /* Solve    nh*l +  sumh*g =  sumt 
              sumh*l + sumhh*g =  sumth */

    double a= nh/sumh;    // nh <= sumh

    /* Subtract a times second eqn from first eqn
       to obtain g */
    *g= (sumt-a*sumth)/(sumh-a*sumhh);

    /* Use second eqn to obtain l */
    *l= (sumth-sumhh* *g)/sumh;

} /* end leastsquares */

void bspbench(){
  
    /**** Determine p ****/
    bsp_begin(P);
    long p= bsp_nprocs(); // p = number of processors
    long s= bsp_pid();    // s = processor number
  
    double *Time= vecallocd(p);
    bsp_push_reg(Time,p*sizeof(double));
    double *dest= vecallocd(2*MAXH+p);
    bsp_push_reg(dest,(2*MAXH+p)*sizeof(double));
    bsp_sync();

    /**** Determine r ****/
    double x[MAXN], y[MAXN], z[MAXN];
    double r= 0.0;
    for (long n=1; n <= MAXN; n *= 2){
        /* Initialize scalars and vectors */
        double alpha= 1.0/3.0;
        double beta= 4.0/9.0;
        for (long i=0; i<n; i++){
            z[i]= y[i]= x[i]= (double)i;
        }
        /* Measure time of 2*NITERS DAXPY operations of length n */
        double time0= bsp_time();
        for (long iter=0; iter<NITERS; iter++){
            for (long i=0; i<n; i++)
                y[i] += alpha*x[i];
            for (long i=0; i<n; i++)
                z[i] -= beta*x[i];
        }
        double time1= bsp_time(); 
        double time= time1-time0; 
        bsp_put(0,&time,Time,s*sizeof(double),sizeof(double));
        bsp_sync();
    
        /* Processor 0 determines minimum, maximum, average
           computing rate */
        if (s==0){
            double mintime= Time[0];
            double maxtime= Time[0];
            for (long s1=1; s1<p; s1++){
                if (Time[s1] < mintime)
                    mintime= Time[s1];
                if (Time[s1] > maxtime)
                    maxtime= Time[s1];
            }
            if (mintime>0.0){
                /* Compute r = average computing rate in flop/s */
                long nflops= 4*NITERS*n;
                r= 0.0;
                for (long s1=0; s1<p; s1++)
                    r += nflops/Time[s1];
                r /= p; 
                printf("n= %5ld min= %7.3lf ",
                        n, nflops/(maxtime*MEGA));
                printf("max= %7.3lf av= %7.3lf Mflop/s ",
                       nflops/(mintime*MEGA), r/MEGA);
                /* Output for fooling too clever compilers */
                printf(" fool=%7.1lf\n",y[n-1]+z[n-1]);
            } else 
                printf("minimum time is 0\n");
            fflush(stdout);
        }
    }
    /* r is taken as the value at length MAXN */

    /**** Determine g and l ****/
    long destproc[MAXH], destindex[MAXH];
    double src[MAXH], t[MAXH+1]; // t[h] = time of h-relation 

    for (long h=0; h<=MAXH; h++){
        /* Initialize communication pattern */
        for (long i=0; i<h; i++){
            src[i]= (double)i;
            if (p==1){
                destproc[i]= 0;
                destindex[i]= i;
            } else {
                /* destination proc is one of the p-1 others */
                destproc[i]= (s+1 + i%(p-1)) %p;
                /* destination index is in my own part of dest */
                destindex[i]= s + (i/(p-1))*p;
            }
        }

        /* Measure time of NITERS h-relations */
        bsp_sync(); 
        double time0= bsp_time(); 
        for (long iter=0; iter<NITERS; iter++){
            for (long i=0; i<h; i++)
                bsp_put(destproc[i],&src[i],dest,
                        destindex[i]*sizeof(double),
                        sizeof(double));
            bsp_sync(); 
        }
        double time1= bsp_time();
        double time= time1-time0;
 
        /* Compute time of one h-relation */
        if (s==0){
            t[h]= (time*r)/NITERS; // time in flop units
            printf("Time of %5ld-relation = %.2lf"
                   " microsec = %8.0lf flops\n",
                   h, (time*MEGA)/NITERS, t[h]);
            fflush(stdout);
        }
    }

    if (s==0){
        double g0, l0, g, l; 
        printf("size double = %ld bytes\n",(long)sizeof(double));
        leastsquares(0,p,t,&g0,&l0); 
        printf("Range h=0 to p   : g= %.1lf, l= %.1lf\n",g0,l0);
        leastsquares(p,MAXH,t,&g,&l);
        printf("Range h=p to %ld: g= %.1lf, l= %.1lf\n",
               (long)MAXH, g,l);

        printf("The bottom line for this BSP computer is:\n");
        printf("p= %ld, r= %.3lf Mflop/s, g= %.1lf, l= %.1lf\n",
               p,r/MEGA,g,l);
        fflush(stdout);
    }
    bsp_pop_reg(dest);
    vecfreed(dest);
    bsp_pop_reg(Time);
    vecfreed(Time);
 
    bsp_end();

} /* end bspbench */

int main(int argc, char **argv){

    bsp_init(bspbench, argc, argv);
    printf("How many processors do you want to use?\n");
    fflush(stdout);
    scanf("%ld",&P);
    if (P > bsp_nprocs()){
        printf("Sorry, only %u processors available.\n",
               bsp_nprocs());
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
    bspbench();
    exit(EXIT_SUCCESS);

} /* end main */
