#include "bspedupack.h"

void bspmv(long n, long nz, long nrows, long ncols,
           double *a, long *inc,
           long *srcprocv, long *srcindv,
           long *destprocu, long *destindu,
           long nv, long nu, double *v, double *u){

    /* This function multiplies a sparse matrix A with a
       dense vector v, giving a dense vector u=Av.
       A is n by n, and u,v are vectors of length n. 
       A, u, and v are distributed arbitrarily on input.
       They are all accessed using local indices, but the local
       matrix indices may differ from the local vector indices.
       The local matrix nonzeros are stored in an incremental
       compressed row storage (ICRS) data structure defined by
       nz, nrows, ncols, a, inc.

       All rows and columns in the local data structure
       are nonempty.
      
       n is the global size of the matrix A.
       nz is the number of local nonzeros.
       nrows is the number of local rows.
       ncols is the number of local columns.

       a[k] is the numerical value of the k'th local nonzero
           of the sparse matrix A, 0 <= k < nz.
       inc[k] is the increment in the local column index 
           of the k'th local nonzero, compared to the column
           index of the (k-1)th nonzero if this nonzero is
           in the same row; otherwise, ncols is added
           to the difference.
           Here, the column index of the -1'th nonzero is 0.

       srcprocv[j] is the source processor of the component in v
           corresponding to local column j, 0 <= j < ncols.
       srcindv[j] is the local index on the source processor
           of the component in v corresponding to local column j.
       destprocu[i] is the destination processor of the partial sum
           corresponding to local row i, 0 <= i < nrows.
       destindu[i] is the local index in the vector u on
           the destination processor corresponding to local row i.
    
       nv is the number of local components of the input vector v.
       nu is the number of local components of the output vector u.
       v[k] is the k'th local component of v, 0 <= k < nv.
       u[k] is the k'th local component of u, 0 <= k < nu.
    */

    /* Allocate, register, and set tag size */
    double *vloc= vecallocd(ncols);
    bsp_push_reg(v,nv*sizeof(double));
    bsp_size_t tagsize= sizeof(long);
    bsp_set_tagsize(&tagsize);
    bsp_sync();

    /****** Superstep (0). Fanout ******/
    long s= bsp_pid();
    for (long j=0; j<ncols; j++){
        if (srcprocv[j] == s)
            vloc[j]= v[srcindv[j]];
        else 
            bsp_get(srcprocv[j],v,srcindv[j]*sizeof(double),
                    &vloc[j],sizeof(double)); 
    }
    bsp_sync();

    /****** Superstep (1)/(2) ******/
    /* Local matrix-vector multiplication and fanin */
    double *pa= a;
    long *pinc= inc;
    double *pvloc= vloc;
    double *pvloc_end= pvloc + ncols;

    for (long i=0; i<nu; i++)
        u[i]= 0.0;

    pvloc += *pinc;
    for (long i=0; i<nrows; i++){
        double sum= 0.0;
        while (pvloc<pvloc_end){
            sum += (*pa) * (*pvloc);
            pa++; 
            pinc++;
            pvloc += *pinc;
        }
        if (destprocu[i] == s)
            u[destindu[i]]= sum;
        else
            bsp_send(destprocu[i],&destindu[i],
                     &sum,sizeof(double)); 
        pvloc -= ncols;
    }
    bsp_sync();

    /****** Superstep (3) ******/
    /* Summation of nonzero partial sums */

    bsp_nprocs_t nsums;
    bsp_size_t nbytes, status;
    long  i;
    bsp_qsize(&nsums,&nbytes);
    bsp_get_tag(&status,&i);
    for (long k=0; k<nsums; k++){
        /* status != -1, but its value is not used */
        double sum;
        bsp_move(&sum,sizeof(double)); 
        u[i] += sum;
        bsp_get_tag(&status,&i);
    }

    bsp_pop_reg(v);
    vecfreed(vloc);

} /* end bspmv */


void bspmv_init(long n, long nrows, long ncols,
                long nv, long nu, long *rowindex, long *colindex,
                long *vindex, long *uindex,
                long *srcprocv, long *srcindv,
                long *destprocu, long *destindu){

    /* This function initializes the communication data structure
       needed for multiplying a sparse matrix A with a dense
       vector v, giving a dense vector u=Av.

       Input: the arrays rowindex, colindex, vindex, uindex,
       containing the global indices corresponding to
       the local indices of the matrix and the vectors.

       Output: initialized arrays srcprocv, srcindv,
       destprocu, destindu containing the processor number
       and the local index on the remote processor of vector
       components corresponding to local matrix columns and rows. 
      
       n, nrows, ncols, nv, nu are the same as in bspmv.

       rowindex[i] is the global index of local row i,
                   0 <= i < nrows.
       colindex[j] is the global index of local column j,
                   0 <= j < ncols.
       vindex[j] is the global index of local component v[j],
                   0 <= j < nv.
       uindex[i] is the global index of local component u[i],
                   0 <= i < nu.

       srcprocv, srcindv, destprocu, destindu
                   are the same as in bspmv.
    */

    long p= bsp_nprocs(); // p = number of processors obtained
    long s= bsp_pid();    // s = processor number
  
    /* Allocate and register temporary arrays */
    long np= nloc(p,s,n);
    long *tmpprocv=vecalloci(np);
    long *tmpindv=vecalloci(np);
    long *tmpprocu=vecalloci(np);
    long *tmpindu=vecalloci(np);

    bsp_push_reg(tmpprocv,np*sizeof(long));
    bsp_push_reg(tmpindv,np*sizeof(long));
    bsp_push_reg(tmpprocu,np*sizeof(long));
    bsp_push_reg(tmpindu,np*sizeof(long));
    bsp_sync();

    /* Write my announcement into temporary arrays */
    for (long j=0; j<nv; j++){
        long jglob= vindex[j];
        /* Use the cyclic distribution */
        bsp_put(jglob%p,&s,tmpprocv,
               (jglob/p)*sizeof(long),sizeof(long));
        bsp_put(jglob%p,&j,tmpindv,
               (jglob/p)*sizeof(long),sizeof(long));
    }
    for (long i=0; i<nu; i++){
        long iglob= uindex[i];
        bsp_put(iglob%p,&s,tmpprocu,
               (iglob/p)*sizeof(long),sizeof(long));
        bsp_put(iglob%p,&i,tmpindu,
               (iglob/p)*sizeof(long),sizeof(long));
    }
    bsp_sync();

    /* Read announcements from temporary arrays */
    for (long j=0; j<ncols; j++){
        long jglob= colindex[j];
        bsp_get(jglob%p,tmpprocv,(jglob/p)*sizeof(long),
                &srcprocv[j],sizeof(long));
        bsp_get(jglob%p,tmpindv, (jglob/p)*sizeof(long),
                &srcindv[j], sizeof(long));
    }
    for (long i=0; i<nrows; i++){
        long iglob= rowindex[i];
        bsp_get(iglob%p,tmpprocu,(iglob/p)*sizeof(long),
                &destprocu[i],sizeof(long));
        bsp_get(iglob%p,tmpindu, (iglob/p)*sizeof(long),
                &destindu[i], sizeof(long));
    }
    bsp_sync();

    /* Deregister temporary arrays */
    bsp_pop_reg(tmpindu); bsp_pop_reg(tmpprocu);
    bsp_pop_reg(tmpindv); bsp_pop_reg(tmpprocv);

    /* Free temporary arrays */
    vecfreei(tmpindu); vecfreei(tmpprocu);          
    vecfreei(tmpindv); vecfreei(tmpprocv);   

    bsp_sync();
} /* end bspmv_init */
