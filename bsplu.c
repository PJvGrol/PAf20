#include "bspedupack.h"

#define EPS 1.0e-15
#define BLOCK 16

void matmat_tall_skinny(double **A, double **B, double **C,
                        long m, long n, long b,
                        long i0, long j0){

    /* This function multiplies the m by b matrix A and the
       transpose of the n by b matrix B and subtracts the
       result A(B^T) from the submatrix C(i0:i0+m-1, j0:j0+n-1).

       The function is written in a cache-friendly way for
       tall-and-skinny matrices A and B (b << m,n), using
       a block size b.
    */

    if (m < 1 || n < 1 || b < 1)
        return;

    for (long ia=0; ia<m; ia+=b){
        long imax= MIN(ia+b,m); 
        for (long jb=0; jb<n; jb+=b){
            long jmax= MIN(jb+b,n); 

            /* Multiply a block from A with a block from B */
            for (long i=ia; i<imax; i++){
                for (long j=jb; j<jmax; j++){
                    double sum= 0.0;
                    for (long k=0; k<b; k++)
                        sum += A[i][k]*B[j][k];
                    C[i0+i][j0+j] -= sum;
                }
            }
        }
    }

} /* end matmat_tall_skinny */


void bsp_permute_rows(long M, long *Src, long *Dest, long nperm,
                      double *pa, long nc, long jc, long jc1){

    /* This function permutes the rows of matrix A such that
       local row Src[i] moves to row Dest[i], for 0 <= i < nperm.
       This is only done for local column indices jc <= j < jc1.
       M = number of processor rows,
       nc = length of local rows of matrix A,
       A is stored as a 1D array pa, which must have been
       registered previously.
    */

    if (jc1 <= jc)
        return;

    long pid= bsp_pid();
    long t= pid/M;

    for (long i=0; i<nperm; i++){
        /* Store row Src[i] of A in row r=Dest[i] on P(r%M,t) */
        long r= Dest[i];
        bsp_put(r%M+t*M,&pa[(Src[i]/M)*nc+jc],pa,
                ((r/M)*nc+jc)*sizeof(double),
                (jc1-jc)*sizeof(double));
    }

} /* end bsp_permute_rows */


void bsp_broadcast(double *x, long n, long src, long s0,
                   long stride, long p0, long phase){

    /* Broadcast the vector x of length n from processor src to
       processors s0+t*stride, 0 <= t < p0. Here n, p0 >= 1.
       The vector x must have been registered previously.
       Processors are numbered in 1D fashion.

       phase = phase of two-phase broadcast (0 or 1)
       Only one phase is performed, without synchronization.
    */

    if (n < 1 || p0 < 1)
        return;

    long s= bsp_pid(); // processor number
    long b= (n%p0==0 ?  n/p0 : n/p0+1); // broadcast block size

    if ((phase==0 && s==src) ||
        (phase==1 && s0 <= s && s < s0+p0*stride &&
                               (s-s0)%stride==0)){
        /* Participate */

        for (long t=0; t<p0; t++){
            long dest= s0+t*stride;
            long t1 = (phase==0 ? t : (s-s0)/stride);
                 // in phase 1: s= s0+t1*stride 

            long nbytes= MIN(b,n-t1*b)*sizeof(double);
            if (nbytes>0 && dest!=src)
                bsp_put(dest,&x[t1*b],x,
                        t1*b*sizeof(double),nbytes);
        }
    }

} /* end bsp_broadcast */


void bsplu(long M, long N, long n, long *pi, double **a){

    /* Compute LU decomposition of n by n matrix A with partial
       row pivoting. The processors are numbered in 2D fashion.
       Program text for P(s,t) = processor s+t*M,
       with 0 <= s < M and 0 <= t < N.
       A is distributed by the M by N cyclic distribution.
       pi stores the output permutation, with a cyclically
       distributed copy in every processor column.
    */

    void swaps2perm(long M, long k0, long b, long *R,
                    long *Src, long *Dest, long *nperm);

    long pid= bsp_pid();
    long s= pid%M;
    long t= pid/M;

    long nr=  nloc(M,s,n); // number of local rows
    long nc=  nloc(N,t,n); // number of local columns

    long r; /* index of pivot row */
    bsp_push_reg(&r,sizeof(long));

    /* Set pointer for 1D access to A */
    double *pa= NULL;
    if (nr>0)
        pa= a[0];
    bsp_push_reg(pa,nr*nc*sizeof(double));

    /* Initialize permutation vector pi */
    for (long i=0; i<nr; i++)
        pi[i]= i*M+s; /* global row index */
    bsp_push_reg(pi,nr*sizeof(long));

    /* Allocate memory for BLOCK columns of L and U^T */
    double **UT= matallocd(nc,BLOCK);
    double **L= matallocd(nr,BLOCK);

    long *R= vecalloci(BLOCK); // global indices of pivot rows
    long *Src= vecalloci(2*BLOCK); // source indices of rows
                                   // to be swapped
    long *Dest= vecalloci(2*BLOCK); // destination indices

    double *Uk=  vecallocd(nc); bsp_push_reg(Uk,nc*sizeof(double));
    double *Lk=  vecallocd(nr); bsp_push_reg(Lk,nr*sizeof(double));
    double *Max= vecallocd(M); bsp_push_reg(Max,M*sizeof(double));
    long *Imax=  vecalloci(M); bsp_push_reg(Imax,M*sizeof(long));

    bsp_sync();

    long k0= 0;
    long b= MIN(BLOCK,n); // algorithmic block size

    for (long k=0; k<n; k++){

        /****** Superstep (0)/(1) ******/
        long kr=   nloc(M,s,k); // first local row
                                // with global index >= k
        long kr1=  nloc(M,s,k+1);     // >= k+1
        long k0rb= nloc(M,s,k0+b);    // >= k0+b
        long kc=   nloc(N,t,k); // first local column
        long kc1=  nloc(N,t,k+1);
        long k0c=  nloc(N,t,k0); 
        long k0cb= nloc(N,t,k0+b);

        if (k%N==t){   /* k=kc*N+t */
            /* Search for local absolute maximum
               in column k of A */
            double absmax= 0.0;
            long imax= DUMMY;
            for (long i=kr; i<nr; i++){
                if (fabs(a[i][kc])>absmax){
                    absmax= fabs(a[i][kc]);
                    imax= i;
                }
            } 

            double max= 0.0;
            if (absmax>0.0)
                max= a[imax][kc];

            /* Broadcast value and local index
               of maximum to P(*,t) */
            for (long s1=0; s1<M; s1++){
                bsp_put(s1+t*M,&max,Max,
                        s*sizeof(double),sizeof(double));
                bsp_put(s1+t*M,&imax,Imax,
                        s*sizeof(long),sizeof(long));
            }
        }
        bsp_sync();
        
        /****** Superstep (2)/(3) ******/
        if (k%N==t){
            /* Determine global absolute maximum (redundantly) */
            double absmax= 0.0;
            long smax= DUMMY;
            for (long s1=0; s1<M; s1++){
                if (fabs(Max[s1])>absmax){
                    absmax= fabs(Max[s1]);
                    smax= s1;
                }
            }
            if (absmax > EPS){
                long imax= Imax[smax];
                r= imax*M+smax; // global index of pivot row
                double pivot= Max[smax];
                for (long i=kr; i<nr; i++)
                    a[i][kc] /= pivot;
                if (s==smax)
                    a[imax][kc]= pivot; // Restore value of pivot
            } else {
                bsp_abort("bsplu at stage %d:"
                          " matrix is singular\n",k);
            }
        }
        /* Obtain index of pivot row */
        bsp_get(s+(k%N)*M,&r,0,&r,sizeof(long));
        bsp_sync();


        /****** Superstep (4)/(5) ******/
        R[k-k0]= r; /* store index of pivot row */

        long nperm= 0;
        long Src2[2], Dest2[2];
        if (k%M==s && r!=k){
            /* Store pi(k) in pi(r) on P(r%M,t) */
            bsp_put(r%M+t*M,&pi[k/M],pi,(r/M)*sizeof(long),
                    sizeof(long));
            Src2[nperm]= k; Dest2[nperm]= r; nperm++;
        }
        if (r%M==s && r!=k){
            bsp_put(k%M+t*M,&pi[r/M],pi,(k/M)*sizeof(long),
                    sizeof(long));
            Src2[nperm]= r; Dest2[nperm]= k; nperm++;
        }
        /* Swap rows k and r for columns in range k0..k0+b-1 */
        bsp_permute_rows(M,Src2,Dest2,nperm,pa,nc,k0c,k0cb);
        bsp_sync();

        /****** Superstep (6) ******/
        /* Phase 0 of two-phase broadcasts */
        if (k%N==t){ 
            /* Store new column k in Lk */
            for (long i=kr1; i<nr; i++)     
                Lk[i-kr1]= a[i][kc];
        }
        if (k%M==s){ 
            /* Store new row k in Uk for columns
               in range k+1..k0+b-1 */
            for (long j=kc1; j<k0cb; j++)
                Uk[j-kc1]= a[kr][j];
        }
        bsp_broadcast(Lk,nr-kr1,s+(k%N)*M,s,M,N,0);
        bsp_broadcast(Uk,k0cb-kc1,(k%M)+t*M,t*M,1,M,0);
        bsp_sync();
       
        /****** Superstep (7) ******/
        /* Phase 1 of two-phase broadcasts */
        bsp_broadcast(Lk,nr-kr1,s+(k%N)*M,s,M,N,1); 
        bsp_broadcast(Uk,k0cb-kc1,(k%M)+t*M,t*M,1,M,1);
        bsp_sync();
      
        /****** Superstep (0') ******/
        /* Update of A for columns in range k+1..k0+b-1 */
        for (long i=kr1; i<nr; i++){
            for (long j=kc1; j<k0cb; j++)
                a[i][j] -= Lk[i-kr1]*Uk[j-kc1];
        }

        /* Store column k in L for rows in range k0+b..n-1 */
        for (long i=k0rb; i<nr; i++)
            L[i-k0rb][k-k0]= Lk[i-kr1];

        /* Perform postponed operations of stages k0,...,k0+b-1 */ 
        if (k==k0+b-1){
            /****** Superstep (8) ******/
            nperm=0;
            swaps2perm(M,k0,b,R,Src,Dest,&nperm);
            bsp_permute_rows(M,Src,Dest,nperm,pa,nc,0,k0c);
            bsp_permute_rows(M,Src,Dest,nperm,pa,nc,k0cb,nc);

            for (long q=0; q<nperm; q++){
                if (Src[q] >= k0+b){ 
                    /* Obtain row Src[q] of L for columns
                       in range k0..k0+b-1 again */ 
                    for (long j=k0; j<k0+b; j++){
                        long tj= j%N;
                        long ncj= nloc(N,tj,n);
                        long i= Src[q]/M;
                        bsp_get(s+tj*M,pa,
                                (i*ncj+j/N)*sizeof(double),
                                &(L[i-k0rb][j-k0]),sizeof(double));
                    }
                }
            }
            bsp_sync();

            if (k < n-1){
                for (long q=0; q<b; q++){
                    /****** Superstep (9) ******/
                    if ((k0+q)%M==s){
                        /* Store row k0+q in Uk for columns
                           in range k0+b..n-1 */
                        for (long j=k0cb; j<nc; j++)
                            Uk[j-k0cb]= a[(k0+q)/M][j];
                    }
                    bsp_broadcast(Uk,nc-k0cb,
                                  (k0+q)%M+t*M,t*M,1,M,0);
                    bsp_sync();
    
                    /****** Superstep (10) ******/
                    bsp_broadcast(Uk,nc-k0cb,(k0+q)%M+t*M,
                                  t*M,1,M,1);
    
                    /* Obtain column k0+q of L 
                       in range k0+q+1..k0+b-1 */
                    long k0rq1= nloc(M,s,k0+q+1);
                    for (long i= k0rq1; i<k0rb; i++){
                        long tj= (k0+q)%N;
                        long ncj= nloc(N,tj,n);
                        bsp_get(s+tj*M,pa,
                                (i*ncj+(k0+q)/N)*sizeof(double),
                                &Lk[i-k0rq1],sizeof(double));
                    }
                    bsp_sync();
    
                    /****** Superstep (9') ******/
                    /* Store row k0+q in UT for columns
                       in range k0+b..n-1 */
                    for (long j=0; j<nc-k0cb; j++)
                        UT[j][q]= Uk[j]; // Store as a column of UT
    
                    /* Update of A in rows k0+q+1..k0+b-1 */
                    for (long i= k0rq1; i<k0rb; i++){
                        for (long j=k0cb; j<nc; j++)
                            a[i][j] -= Lk[i-k0rq1]*Uk[j-k0cb];
                    }
                }
    
                /* Rank-b update of the matrix A */
                matmat_tall_skinny(L,UT,a, nr-k0rb, nc-k0cb,
                                   b, k0rb, k0cb);
    
                /* Prepare for next block */
                k0 += b;
                b= MIN(BLOCK,n-k0);
            }
        }
    }

    bsp_pop_reg(Imax); vecfreei(Imax);
    bsp_pop_reg(Max);  vecfreed(Max);
    bsp_pop_reg(Lk);   vecfreed(Lk);
    bsp_pop_reg(Uk);   vecfreed(Uk);
    bsp_pop_reg(pi); bsp_pop_reg(pa); bsp_pop_reg(&r);
    vecfreei(Dest); vecfreei(Src); vecfreei(R);
    matfreed(L); matfreed(UT);

} /* end bsplu */
