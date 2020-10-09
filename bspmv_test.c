#include "bspedupack.h"
#include "bspsparse_input.h"

/* This is a test program which uses bspmv to multiply a 
   sparse matrix A and a dense vector u to obtain a dense vector v.
   The sparse matrix and its distribution are read from an input file.
   The dense vector v is initialized by the test program.
   The distribution of v is read from an input file.
   The distribution of u is read from another input file.

   The output vector is defined by
       u[i]= (sum: 0 <= j < n: a[i][j]*v[j]).
*/

#define NITERS 1000
#define NPRINT 10  // Print at most NPRINT values per processor

long P;

void bspmv_init(long n, long nrows, long ncols,
                long nv, long nu, long *rowindex, long *colindex,
                long *vindex, long *uindex,
                long *srcprocv, long *srcindv,
                long *destprocu, long *destindu);

void bspmv(long n, long nz, long nrows, long ncols,
           double *a, long *inc,
           long *srcprocv, long *srcindv,
           long *destprocu, long *destindu,
           long nv, long nu, double *v, double *u);

void bspinputvec(const char *filename,
                 long *pn, long *pnv, long **pvindex){
  
    /* This function reads the distribution of a dense
       vector v from the input file and initializes the
       corresponding local index array.
       The input consists of one line
           n p    (number of components, processors)
       followed by n lines in the format
           i proc (index, processor number),
       where i=1,2,...,n.
       
       Input:
       filename is the name of the input file.

       Output:
       n is the global length of the vector.
       nv is the local length.
       vindex[i] is the global index corresponding to
                 the local index i, 0 <= i < nv.
    */

    long p= bsp_nprocs(); // p = number of processors obtained
    long s= bsp_pid(); // processor number

    /* Initialize fp and Nv */
    FILE *fp = NULL;
    long *Nv = NULL; // Nv[q] = number of components of P(q)

    bsp_push_reg(pn,sizeof(long));
    bsp_push_reg(pnv,sizeof(long));
    bsp_sync();

    if (s==0){
        /* Open the file and read the header */
        long procsv;
        fp=fopen(filename,"r");
        fscanf(fp,"%ld %ld\n", pn, &procsv);
        if(procsv!=p)
            bsp_abort("Error: p not equal to p(vec)\n"); 
        for (long q=0; q<p; q++)
            bsp_put(q,pn,pn,0,sizeof(long));
    }
    bsp_sync();

    /* The owner of the global index i and its
       local index are stored in temporary arrays 
       which are distributed cyclically. */
    long n = *pn;
    long np= nloc(p,s,n);
    long *tmpproc= vecalloci(np);
    long *tmpind= vecalloci(np);
    bsp_push_reg(tmpproc,np*sizeof(long));
    bsp_push_reg(tmpind,np*sizeof(long));
    bsp_sync();

    if (s==0){
        /* Allocate component counters */
        Nv= vecalloci(p);
        for (long q=0; q<p; q++)
            Nv[q]= 0;
    }

    long b= (n%p==0 ? n/p : n/p+1); // block size for vector read
    for (long q=0; q<p; q++){
        if(s==0){
            /* Read the owners of the vector components from
               file and put owner and local index into a
               temporary location. This is done n/p
               components at a time to save memory  */
            for (long k=q*b; k<(q+1)*b && k<n; k++){
                long i, proc, ind;
                fscanf(fp,"%ld %ld\n", &i, &proc);
                /* Convert index and processor number
                   to ranges 0..n-1 and 0..p-1,
                   assuming they were in 1..n and 1..p */
                i--;  
                proc--;
                ind= Nv[proc];
                if(i!=k)
                    bsp_abort("Error: i not equal to index \n");
                bsp_put(i%p,&proc,tmpproc,(i/p)*sizeof(long),sizeof(long));
                bsp_put(i%p,&ind,tmpind,(i/p)*sizeof(long),sizeof(long));
                Nv[proc]++;
            }
        }
        bsp_sync();
    }

    if(s==0){
        for (long q=0; q<p; q++)
            bsp_put(q,&Nv[q],pnv,0,sizeof(long));
    }
    bsp_sync();

    /* Store the components at their final destination */
    long *vindex= vecalloci(*pnv);  
    bsp_push_reg(vindex,(*pnv)*sizeof(long));
    bsp_sync();

    for (long k=0; k<np; k++){
        long globk= k*p+s;
        bsp_put(tmpproc[k],&globk,vindex,
                tmpind[k]*sizeof(long),sizeof(long));
    }
    bsp_sync();

    if (s==0)
        fclose(fp);
    bsp_pop_reg(vindex);
    bsp_pop_reg(tmpind);
    bsp_pop_reg(tmpproc);
    bsp_pop_reg(pnv);
    bsp_pop_reg(pn);
    bsp_sync();

    vecfreei(Nv);
    vecfreei(tmpind);
    vecfreei(tmpproc);

    *pvindex= vindex;

} /* end bspinputvec */


void bspmv_test(){
    
    bsp_begin(P);
    long p= bsp_nprocs(); /* p=P */
    long s= bsp_pid();
    
    /* Input of sparse matrix into triple storage */
    long n, nz, *ia, *ja;
    double *a;
    double suma = bspinput2triple(&n,&nz,&ia,&ja,&a);

    /* Convert data structure to incremental compressed row storage */
    long nrows, ncols, *rowindex, *colindex, *start;
    triple2icrs(n,nz,ia,ja,a,&nrows,&ncols,&rowindex,&colindex,&start);
    
    /* Read vector distributions */
    char vfilename[STRLEN], ufilename[STRLEN]; // only used by P(0)
    if (s==0){
        printf("Please enter the filename of the v-vector distribution\n");
        scanf("%s",vfilename);
    }
    long nv, *vindex;
    bspinputvec(vfilename,&n,&nv,&vindex);

    if (s==0){ 
        printf("Please enter the filename of the u-vector distribution\n");
        scanf("%s",ufilename);
    }
    long nu, *uindex;
    bspinputvec(ufilename,&n,&nu,&uindex);

    if (s==0){
        printf("Sparse matrix-vector multiplication");
        printf(" using %ld processors\n",p);
    }

    /* Initialize input vector v */
    double *v= vecallocd(nv);
    for (long i=0; i<nv; i++){
        v[i]= 1.0;
    }
    double *u= vecallocd(nu);
    
    if (s==0){
        printf("Initialization for matrix-vector multiplications\n");
        fflush(stdout);
    }
    bsp_sync(); 
    double time0= bsp_time();
    
    long *srcprocv= vecalloci(ncols);
    long *srcindv= vecalloci(ncols);
    long *destprocu= vecalloci(nrows);
    long *destindu= vecalloci(nrows);
    bspmv_init(n,nrows,ncols,nv,nu,
               rowindex,colindex,vindex,uindex,
               srcprocv,srcindv,destprocu,destindu);

    if (s==0){
        printf("Start of %ld matrix-vector multiplications.\n",
               (long)NITERS);
        fflush(stdout);
    }
    bsp_sync(); 
    double time1= bsp_time();
    
    for (long iter=0; iter<NITERS; iter++)
        bspmv(n,nz,nrows,ncols,a,ia,srcprocv,srcindv,
              destprocu,destindu,nv,nu,v,u);
    bsp_sync();
    double time2= bsp_time();
    
    if (s==0){
        printf("End of matrix-vector multiplications.\n");
        printf("Initialization took only %.6lf seconds.\n",time1-time0);
        printf("Each matvec took only %.6lf seconds.\n", 
                      (time2-time1)/(double)NITERS);
        printf("The computed solution is (<= %ld values per processor):\n",
                (long)NPRINT);
        fflush(stdout);
    }

    for (long i=0; i < NPRINT && i<nu; i++){
        long iglob= uindex[i];
        printf("Proc %ld: i=%ld, u=%lf \n",s,iglob,u[i]);
    }

    /* Check error by computing
           sum (u[i] : 0 <= i < n) = sum (a[i][j] : 0 <= i, j < n) */
  
    double *SumU= vecallocd(p);
    bsp_push_reg(SumU,p*sizeof(double));
    bsp_sync();

    double sumu= 0.0;
    for (long i=0; i<nu; i++)
        sumu += u[i];
    bsp_put(0,&sumu, SumU, s*sizeof(double), sizeof(double));
    bsp_sync();

    if (s==0){
        double totalU= 0.0;
        for (long t=0; t<p; t++)
            totalU += SumU[t];
        printf("Sum(u)=%lf sum(A)=%lf checksum error = %lf\n",
               totalU,suma,fabs(totalU-suma));
    }
    bsp_pop_reg(SumU);
    bsp_sync();
    vecfreed(SumU); 

    vecfreei(destindu); vecfreei(destprocu); 
    vecfreei(srcindv);  vecfreei(srcprocv); 
    vecfreed(u);        vecfreed(v);
    vecfreei(uindex);   vecfreei(vindex);
    vecfreei(rowindex); vecfreei(colindex);
    vecfreei(start);    vecfreei(ja);
    vecfreei(ia);       vecfreed(a);

    bsp_end();
    
} /* end bspmv_test */


int main(int argc, char **argv){

    bsp_init(bspmv_test, argc, argv);

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
    bspmv_test();

    /* Sequential part */
    exit(EXIT_SUCCESS);

} /* end main */
