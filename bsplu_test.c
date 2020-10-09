#include "bspedupack.h"

/*  This is a test program which uses bsplu to decompose an n by n
    matrix A into triangular factors L and U, with partial row pivoting.
    The decomposition is A(pi(i),j)=(LU)(i,j), for 0 <= i,j < n,
    where pi is a permutation.
 
    The input matrix A is a row-rotated version of a matrix B:
        the matrix B is defined by: B(i,j)= 0.5*i+1     if i<=j 
                                            0.5*j+0.5      i>j,
        the matrix A is defined by: A(i,j)= B((i-1) mod n, j).

    This should give as output: 
        the matrix L given by: L(i,j)= 0   if i<j,
                                     = 1      i=j,
                                     = 0.5    i>j.
        the matrix U given by: U(i,j)= 1   if i<=j,
                                     = 0      i>j.
        the permutation pi given by:  pi(i)= (i+1) mod n.

    Output of L and U is in triples (i,j,L\U(i,j)):
        (i,j,0.5) for i>j
        (i,j,1)   for i<=j
    Output of pi is in pairs (i,pi(i))
        (i,(i+1) mod n) for all i.

    The matrix A is constructed such that the pivot choice is unique.
    In stage k of the LU decomposition, row k is swapped with row r=k+1.
    For the M by N cyclic distribution this forces a row swap
    between processor rows. 
*/

long M, N;

void swaps2perm(long M, long k0, long b, long *R,
                long *Src, long *Dest, long *nperm){

    /* This function determines the row permutation to be done
       based on the pivoting swaps from stages k0,k0+1,...,k0+b-1.
       M = number of processor rows,
       R[i]= pivot row of stage k0+i, for 0 <= i < b.

       The output states that row Src[i] has to be permuted into
       row Dest[i], for 0 <= i < nperm. The source rows are local.
       All row indices are given as global indices.
    */

    long pid= bsp_pid();
    long s= pid%M;

    /* Determine local rows i >= k0+b to be permuted */
    long count= 0;
    for (long i=0; i<b; i++){
        if (R[i]%M==s && R[i] >= k0+b ){
            /* Check if R[i] already occurs */
            bool duplicate= false;
            for (long j=0; j<count; j++){
                if (R[i]==Src[j]){
                    duplicate= true;
                    break;
                }
            }
            if (!duplicate){
                Src[count]= R[i]; 
                count++;
            }
        }
    }

    /* Determine local rows i with k0 <= i < k0+b to be permuted */
    for (long i=0; i<b; i++){
        if ((k0+i)%M==s && R[i]!=k0+i){
            Src[count]= k0+i;
            count++;
        }
    }
    *nperm= count;

    /* Compute source and destination for rows to be permuted */
    for (long j=0; j<count; j++){
        long dest= Src[j];
        for (long i=0; i<b; i++){
            /* Compute effect of swap(k0+i,R[i]) on dest */
            if (dest==k0+i){
                dest= R[i];
            } else if (dest==R[i]){
                dest= k0+i;
            }
        }
        Dest[j]= dest;
    }

} /* end swaps2perm */


void bsplu_test(){

    void bsplu(long M, long N, long n, long *pi, double **a);
  
    bsp_begin(M*N);
    long p= bsp_nprocs(); /* p=M*N */
    long pid= bsp_pid();

    bsp_push_reg(&M,sizeof(long));
    bsp_push_reg(&N,sizeof(long));
    long n; /* matrix size */
    bsp_push_reg(&n,sizeof(long));
    bsp_sync();

    if (pid==0){
        printf("Please enter matrix size n:\n");
        scanf("%ld",&n);
        for (long q=0; q<p; q++){
            bsp_put(q,&M,&M,0,sizeof(long));
            bsp_put(q,&N,&N,0,sizeof(long));
            bsp_put(q,&n,&n,0,sizeof(long));
        }
    }
    bsp_sync();
    bsp_pop_reg(&n); /* not needed anymore */
    bsp_pop_reg(&N);
    bsp_pop_reg(&M);

    /* Compute 2D processor numbering from 1D numbering */
    long s= pid%M;  /* 0 <= s < M */
    long t= pid/M;  /* 0 <= t < N */

    /* Allocate and initialize matrix */
    long nlr=  nloc(M,s,n); /* number of local rows */
    long nlc=  nloc(N,t,n); /* number of local columns */
    double **a= matallocd(nlr,nlc);
    long *pi= vecalloci(nlr);
  
    if (s==0 && t==0){
        printf("LU decomposition of %ld by %ld matrix\n",n,n);
        printf("using the %ld by %ld cyclic distribution\n",M,N);
    }
    for (long i=0; i<nlr; i++){
        long iglob= i*M+s;    /* Global row index in A */
        iglob= (iglob-1+n)%n; /* Global row index in B */
        for (long j=0; j<nlc; j++){
            long jglob= j*N+t;     /* Global column index in A and B */
            a[i][j]= (iglob<=jglob ? 0.5*iglob+1 : 0.5*(jglob+1) );
        }
    }
  
    if (s==0 && t==0)
        printf("Start of LU decomposition\n");
    bsp_sync();
    double time0= bsp_time();
 
    bsplu(M,N,n,pi,a);
    bsp_sync();
    double time1= bsp_time();
 
    if (s==0 && t==0){
        printf("End of LU decomposition\n");
        printf("This took only %.6lf seconds.\n", time1-time0);
        printf("\nThe output permutation is:\n"); fflush(stdout);
    }

    if (t==0){
        bool pi_error= false;
        for (long i=0; i<nlr; i++){
            long iglob= i*M+s;
            printf("i=%ld, pi=%ld, proc=(%ld,%ld)\n",iglob,pi[i],s,t);
            if (pi[i] != (iglob+1)%n)
                pi_error = true;
        }
        if (pi_error)
            printf("Error in pi found by proc=(%ld,%ld)\n",s,t);
        fflush(stdout);
    }
    bsp_sync();

    if (s==0 && t==0){  
        printf("\nThe output matrix is:\n"); fflush(stdout);
    }   
    double maxerror= 0.0;
    for (long i=0; i<nlr; i++){
        long iglob= i*M+s;
        for (long j=0; j<nlc; j++){
            long jglob= j*N+t;
            printf("i=%ld, j=%ld, a=%f, proc=(%ld,%ld)\n",
                   iglob,jglob,a[i][j],s,t);
            double error= (iglob<=jglob ? a[i][j] - 1.0 : a[i][j] - 0.5);
            if (fabs(error) > maxerror)
                maxerror= fabs(error);
        }
    }
    printf("Maximum error of proc=(%ld,%ld) = %lf\n",s,t,maxerror);

    vecfreei(pi);
    matfreed(a);

    bsp_end();

} /* end bsplu_test */


int main(int argc, char **argv){
 
    bsp_init(bsplu_test, argc, argv);

    printf("Please enter number of processor rows M:\n");
    scanf("%ld",&M);
    printf("Please enter number of processor columns N:\n");
    scanf("%ld",&N);
    if (M*N > bsp_nprocs()){
        printf("Sorry, only %u processors available.\n", bsp_nprocs()); 
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    bsplu_test();
    exit(EXIT_SUCCESS);

} /* end main */
