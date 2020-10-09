#include "bspedupack.h"

void butterfly_stage(double complex *x, long n, long k, 
                     bool forward, double complex *w){

    /* This sequential function computes one butterfly
       stage k of a complex vector x of length n, where k
       is even and n mod k = 0. The output overwrites x.

       The function uses a table w of k/2 complex weights
       which must have been initialized beforehand.

       Each pair (x[j], x[j+k/2]) with 0 <= j < n and 
       j mod k < k/2 is transformed into a pair
          (x[j] + weight*x[j+k/2], x[j] - weigth*x[j+k/2]),
       where weight = w[j mod k] in the forward case
       and its complex conjugate otherwise.

    */

    for (long r=0; r<n/k; r++){
        for (long j=0; j<k/2; j++){
            double complex weight;
            if (forward) {
                weight= w[j];
            } else {
                weight= conj(w[j]);
            }

            double complex tau= weight * x[r*k+j+k/2];
            x[r*k+j+k/2]= x[r*k+j] - tau;
            x[r*k+j] += tau;
        }
    }

} /* end butterfly_stage */


void ufft(double complex *x, long n, bool forward,
          double complex *w){

    /* This sequential function computes the unordered
       discrete Fourier transform of a complex vector x
       of length n, where n=2^m, m >= 0.
       The output overwrites x.

       If forward, then the forward unordered DFT is computed,
       and otherwise the backward unordered DFT.

       w is a table of complex weights of length n-1, 
       which must have been suitably initialized before
       calling this function.
    */

    long start= 0;
    for (long k=2; k<=n; k *=2){
        butterfly_stage(x,n,k,forward,&w[start]);
        start += k/2;
    }

} /* end ufft */


void permute(double complex *x, long n, long *sigma){

    /* This in-place sequential function permutes a complex
       vector x of length n >= 1 by the permutation sigma,
           y[j] = x[sigma[j]], 0 <= j < n.
       The output overwrites the vector x.

       sigma is a permutation of length n that must be 
       decomposable into disjoint swaps; an example is
       the bit-reversal permutation.
    */

    for (long j=0; j<n; j++){
        long sigmaj = sigma[j];
        if (j<sigmaj){
            /* swap components j and sigma[j] */
            double complex tmp= x[j];
            x[j]= x[sigmaj];
            x[sigmaj]= tmp;
        }
    }
           
} /* end permute */


void bitrev_init(long n, long *rho){

    /* This function initializes the bit-reversal permutation
       rho of length n, where n=2^m, m >= 0. */

    rho[0]= 0;
    for (long k=2; k<=n; k *=2){
        for (long j=0; j<k/2; j++){
            rho[j] *= 2;
            rho[j+k/2]= rho[j] + 1;
        }
    }
    
} /* end bitrev_init */


void bspredistr(double complex *x, long n, long c0, long c,
                bool rev, long *rho_p){

    /* This function redistributes the complex vector x
       of length n, from group-cyclic distribution over
       p processors with cycle c0 to cycle c, where
       c0, c, p, n are powers of two with
           1 <= c0 <= c <= p <= n.
       If rev, the function assumes the processor numbering
       is bit-reversed on input.
       rho_p is the bit-reversal permutation of length p.
    */

    long p= bsp_nprocs();
    long s= bsp_pid();
    long np= n/p;

    long j0, j2;
    if (rev) {
        j0= rho_p[s]%c0; 
        j2= rho_p[s]/c0; 
    } else {
        j0= s%c0;
        j2= s/c0;
    }

    long ratio= c/c0;
    long size=  (np >= ratio ? np/ratio : 1 );
    long npackets= np/size;
    double complex *tmp= vecallocc(size);

    for (long j=0; j<npackets; j++){
        long jglob= j2*c0*np + j*c0 + j0;
        long destproc=  (jglob/(c*np))*c + jglob%c;
        long destindex= (jglob%(c*np))/c;
        for (long r=0; r<size; r++)
            tmp[r]=   x[j+r*ratio]; 
         
        bsp_put(destproc,tmp,x,destindex*sizeof(double complex),
                size*sizeof(double complex));
    }
    bsp_sync();
    vecfreec(tmp);

} /* end bspredistr */


void bspfft(double complex *x, long n, bool forward,
            double complex *w, long *rho_np, long *rho_p){

    /* This parallel function computes the discrete
       Fourier transform of a complex vector x of length
       n=2^m, m >= 1.
       x must have been registered before calling this function.
       The number of processors p must be a power of two.

       The function uses a weight table w which stores
       in sequence all the weights used locally in the
       butterfly stages of the FFT.

       The function uses two bit-reversal permutations:
           rho_np of length n/p,
           rho_p of length p. 
       The weight table and bit-reversal permutations 
       must have been initialized before calling this function.

       If forward, then the DFT is computed,
           y[k] = sum j=0 to n-1 exp(-2*pi*i*k*j/n)*x[j],
       for 0 <= k < n.
       Otherwise, the inverse DFT is computed,
           y[k] = (1/n) sum j=0 to n-1 exp(+2*pi*i*k*j/n)*x[j],
       for 0 <= k < n.
       Here, i=sqrt(-1). The output vector y overwrites x.
    */

    long p= bsp_nprocs();
    long np= n/p;
    long c= 1;
    bool rev= true;

    /* Perform a local ordered FFT of length n/p.
       This part can be replaced easily by your favourite
       sequential FFT */
    permute(x,np,rho_np);
    ufft(x,np,forward,w);

    long k= 2*np;
    long start= np-1; // start of current weights in w
    while (c < p){
        long c0= c;
        c= ( np*c <= p ? np*c : p);
        bspredistr(x,n,c0,c,rev,rho_p);
        rev= false;

        while (k <= np*c){
            butterfly_stage(x,np,k/c,forward,&w[start]);
            start += k/(2*c);
            k *= 2;
        }
    }

    /* Normalize the inverse FFT */
    if (!forward){
        double ninv= 1 / (double)n;
        for (long j=0; j<np; j++)
            x[j] *= ninv;
    }

} /* end bspfft */


void bspfft_init(long n, double complex *w, long *rho_np,
                 long *rho_p){

    /* This parallel function initializes all the tables
       used in the FFT. */

    long p= bsp_nprocs();
    long s= bsp_pid();
    long np= n/p;

    /* Initialize bit-reversal tables */
    bitrev_init(np,rho_np);
    bitrev_init(p,rho_p);

    /* Initialize weight tables */
    long c= 1;
    long k= 2;
    long start= 0;

    while (c <= p){
        /* Initialize weights for superstep with cycle c */
        long j0= s%c;
        while (k <= np*c){
            /* Initialize k/(2c) weights for stage k */
            double theta= -2.0 * M_PI / (double)k;

            for (long j=0; j<k/(2*c); j++, start++){
                double jtheta= (j0 + j*c) * theta; 
                w[start]= cexp(jtheta*I);
            }

            k *= 2;
        }
        if (c < p)
            c= ( np*c <= p ? np*c : p);
        else
            c= 2*p; //done
    }

} /* end bspfft_init */
