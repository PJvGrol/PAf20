#include "bspedupack.h"
#include "bspsparse_input.h"

/* This is a test program which uses bspmatch to match vertices in 
   a weighted undirected graph with edge weights > 0. 
   The graph G and its distribution are read from an input file
   in the format of a sparse adjacency matrix A,
   where nonzero a[i][j] represents the weight of the edge (i,j).

   The matrix must have an empty diagonal (a[i][i]=0 for all i)
   and it must be symmetric (a[i][j]=a[j][i] for all i,j),
   with both entries present.

   The matrix must have been distributed by a 1D row distribution,
   which also determines the vertex distribution.
   
   The output is the number of matches, the total weight of the matches,
   and a number of matched vertex pairs.
*/

#define NPRINT 10  // Print at most NPRINT matches per processor

long P;

void bspmatch(long nvertices, long nedges, long nhalo, long *v0, long *v1,
              long *destproc, double *weight, long *weight1,
              long *Adj, long *Start, long *degree, long maxops,
              long *nmatch, long *match, long *nsteps, long *nops);


void bspmatch_test(){
    
    bsp_begin(P);

    /***** Part 0: prepare input *****/

    long p= bsp_nprocs(); // p=P
    long s= bsp_pid();
    
    /* Input of sparse matrix into triple storage */
    long n, nz, *ia, *ja;
    double *weight;
    double suma= bspinput2triple(&n,&nz,&ia,&ja,&weight);

    /* Sequential part */
    long maxops=0;
    bsp_push_reg(&maxops,sizeof(long));
    bsp_sync();

    if (s==0){
        printf("Please enter the maximum number of operations per superstep\n");
        printf("    (0 if no maximum)\n"); fflush(stdout);
        scanf("%ld",&maxops);
        if (maxops>0)
           printf("Maximum number of operations per superstep = %ld\n\n", maxops);
        else
           printf("No maximum number of operations per superstep\n\n");
        for (long t=0; t<p; t++)
            bsp_put(t,&maxops,&maxops,0,sizeof(long));
    }

    /* Convert data structure to incremental compressed row storage (ICRS),
       where ia contains the local column index increments,
       ja the local column indices, and Start the starting points of rows. */

    long nrows, ncols, *rowindex, *colindex, *Start;
    triple2icrs(n,nz,ia,ja,weight,&nrows,&ncols,&rowindex,&colindex,&Start);

    vecfreei(ia); // increments are not needed

    /* Translate to graph language. Here, nz is the number of edges
       including symmetric duplicates. */
    long nvertices= nrows; // number of local vertices (with degree > 0)
    long *v0= vecalloci(nz);
    long *v1= vecalloci(nz);
    long *weight1= vecalloci(nz);
    long *degree= vecalloci(nvertices);
    
    /* Search for the local row (vertex) corresponding to each global column.
      We use that both rowindex and colindex are ordered by increasing
      global index. */ 

    long *rowvertex= vecalloci(ncols);
    long r=0; //local row index
    for (long j=0; j<ncols; j++){
        long jglob= colindex[j];
        while (r<nrows && rowindex[r]<jglob)
           r++;
        if (r<nrows && rowindex[r]==jglob)
            rowvertex[j]= r;     // local vertex
        else
            rowvertex[j]= DUMMY; // nonlocal vertex
    }

    /* Initialize v0, v1, weight1, degree */
    long nhalo= 0; // number of halo edges
    double sum_max= 0.0; // sum of the maximum weights of the local rows
    for (long i=0; i<nrows; i++){
        degree[i]= Start[i+1] - Start[i];
        double maxw= 0.0; // maximum weight of row i
        for (long k=Start[i]; k<Start[i+1]; k++){
            v0[k]= i; // local index of row vertex
            v1[k]= rowvertex[ja[k]]; // local index of column vertex
            if (v1[k] == DUMMY){ // halo edge
                nhalo++;
                weight1[k]= rowindex[i] + colindex[ja[k]];
            } else
                weight1[k]= 2*n + rowindex[i] + colindex[ja[k]];
            if (weight[k] > maxw)
                maxw= weight[k];
        }
        sum_max += maxw;
    }
    if ((nz - nhalo)%2==1)
        bsp_abort("Error on input: nz-nhalo is odd\n");
    long nedges= (nz - nhalo)/2;

    /* Determine new edge numbers */
    long *new=vecalloci(nz);
    long count= 0;
    for (long k=0; k<nz; k++){
        if (v1[k]!=DUMMY && v0[k]<v1[k]){ // local edges are registered first,
           new[k]= count;                 // and only once
           count++;
        }
    }
    for (long k=0; k<nz; k++){
        if (v1[k]==DUMMY){ // halo edges second
           new[k]= count;
           count++;
        }
    }

    /* Insert edges into adjacency lists */
    long *Adj= vecalloci(nz);
    long *free= vecalloci(nrows);
 
    for (long i=0; i<nrows; i++)
        free[i]= Start[i]; // first free position for row i

    for (long k=0; k<nz; k++){
        if (v1[k]==DUMMY){
            Adj[free[v0[k]]]= new[k];
            free[v0[k]]++;
        } else if (v0[k] < v1[k]){ // insert edge in two lists
            Adj[free[v0[k]]]= new[k]; free[v0[k]]++; 
            Adj[free[v1[k]]]= new[k]; free[v1[k]]++; 
        }
    }
    vecfreei(free);
   
    /*  Copy the edges according to the new numbering */
    long nedges_tot= nedges + nhalo;
    long *v0new= vecalloci(nedges_tot);
    long *v1new= vecalloci(nedges_tot);
    long *janew= vecalloci(nedges_tot);
    double *weightnew= vecallocd(nedges_tot);
    long *weight1new= vecalloci(nedges_tot);

    for (long k=0; k<nz; k++){
        if (v1[k]==DUMMY || v0[k] < v1[k]){
            v0new[new[k]]= v0[k];
            janew[new[k]]= ja[k]; // keep a copy for output printing
            if (v1[k]==DUMMY){
                v1new[new[k]]= ja[k]; // register the local column index
            } else
                v1new[new[k]]= v1[k]; 
            weightnew[new[k]]= weight[k];
            weight1new[new[k]]= weight1[k];
        }
    }
    vecfreei(new);
    vecfreei(weight1);
    vecfreed(weight);
    vecfreei(v1);
    vecfreei(v0);
    vecfreei(ja);
    
    /* Initialize communication data structure */
    long np= nloc(p,s,n);
    long *tmpproc=vecalloci(np); // temporary array for storing the owners
                                 // of the vertices
    bsp_push_reg(tmpproc,np*sizeof(long));

    /* Initialize owner to 0 as a default for vertices with degree 0 */
    for (long i=0; i<np; i++)
        tmpproc[i]= 0;

    /* Set tagsize for communication of halo edge numbers */
    bsp_size_t tagsize= sizeof(indexpair);
    bsp_set_tagsize(&tagsize);
    bsp_sync();

    /* Announce my vertices (rows) */
    for (long i=0; i<nrows; i++){
        long iglob= rowindex[i];
        bsp_put(iglob%p,&s,tmpproc,(iglob/p)*sizeof(long),sizeof(long));
    }
    bsp_sync();

    /* Determine the owner proc[j] of vertex colindex[j],
       for each local column j, by reading the announcements */
    long *proc= vecalloci(ncols);
    for (long j=0; j<ncols; j++){
        if (rowvertex[j]==DUMMY){
            long jglob= colindex[j];
            bsp_get(jglob%p,tmpproc,(jglob/p)*sizeof(long),
                    &(proc[j]),sizeof(long));
        } else
            proc[j]= s;
    }
    vecfreei(rowvertex);
    bsp_sync();
    bsp_pop_reg(tmpproc);

    /* Initialize destproc[e-nedges]= owner of halo edge e 
       and send the local edge number e to this owner */
    long *destproc= vecalloci(nhalo);
    indexpair tag;
    for (long e=nedges; e<nedges_tot; e++){
        destproc[e-nedges]= proc[v1new[e]]; // local column index
        tag.i= rowindex[v0new[e]]; // global row index
        tag.j= colindex[v1new[e]]; // global column index
        bsp_send(destproc[e-nedges], &tag, &e, sizeof(long));
    }
    bsp_sync();

    vecfreei(proc); 
    vecfreei(tmpproc); 

    /* Receive triples (i,j,e), where i,j are global indices and
       e is the local edge number on the sending processor */

    bsp_nprocs_t nmsg; // total number of messages received
    bsp_size_t nbytes;    // total size in bytes received
    bsp_qsize(&nmsg,&nbytes);
    if (nmsg != nhalo)
        bsp_abort("Error: number of messages <> nhalo\n");

    long *imsg= vecalloci(nmsg);
    long *jmsg= vecalloci(nmsg);
    long *emsg= vecalloci(nmsg);
    for (long k=0; k<nmsg; k++){
        bsp_size_t status; // not used
        bsp_get_tag(&status, &tag);
        imsg[k]= tag.i;
        jmsg[k]= tag.j;
        bsp_move(&(emsg[k]), sizeof(long));
    }

    /* Sort nmsg edges by primary key j and secondary key i,
       using as radix the smallest power of two >= sqrt(n).
       The div and mod operations are cheap for powers of two.
       A radix of about sqrt(n) minimizes memory and time. */

    long radix;
    for (radix=1; radix*radix<n; radix *= 2)
       ;

    sort(n,nmsg,imsg,jmsg,emsg,radix,MOD); // imsg goes first
    sort(n,nmsg,imsg,jmsg,emsg,radix,DIV);
    sort(n,nmsg,jmsg,imsg,emsg,radix,MOD); // jmsg goes first
    sort(n,nmsg,jmsg,imsg,emsg,radix,DIV); 

    /* Couple the local halo edges with the remote halo edges.
       The local halo edges e = (i,j) have been sorted
       in the CRS data structure by primary key i and secondary key j,
       with i and j global indices.

       A remote edge (j,i) corresponds to a local edge (i,j).
       For this reason, the received edges have been sorted
       by primary key j and secondary key i. */

    for (long e=nedges; e<nedges_tot; e++)
        v1new[e]= emsg[e-nedges];

    vecfreei(emsg);
    vecfreei(jmsg);
    vecfreei(imsg);

    long nmatch= 0;          // number of matches found
    long *match= vecalloci(nvertices); // matches found
    long nsteps= 0;  // number of (mixed) supersteps taken
    long nops= 0;    // number of elemental operations carried out

    /***** Part 1: run matching *****/
    if (s==0){
        printf("Start of graph matching using %ld processors\n",p);
        fflush(stdout);
    }

    bsp_sync(); // an extra sync to make sure the receive buffer
                // is empty before bspmatch starts

    bsp_sync(); 
    double time0= bsp_time();

    bspmatch(nvertices, nedges, nhalo, v0new, v1new,
             destproc, weightnew, weight1new,
             Adj, Start, degree, maxops, &nmatch, match, &nsteps, &nops);

    bsp_sync(); 
    double time1= bsp_time();

    /***** Part 2: print output *****/
    
    if (s==0){
        printf("End of the graph matching.\n");
        printf("Matching took only %.6lf seconds.\n\n",time1-time0);
        printf("The computed solution is (<= %ld values per processor, 1-based):\n",
                (long)NPRINT);
        fflush(stdout);
    }

    /* Print up to NPRINT matches, including external edges found twice.
       The 0-based vertex numbers are converted to 1-based for output. */
    for (long k=0; k < NPRINT && k<nmatch; k++){
        long e= match[k];
        long iglob= rowindex[v0new[e]];
        long jglob= colindex[janew[e]];
        printf("Proc %ld: found match (i,j)= (%ld, %ld) with weight=%lf \n",
                   s,iglob+1,jglob+1,weightnew[e]);
        fflush(stdout);
    }

    /* Compute the total weight, number of matches, and number of
       operations done */
    double *SumW= vecallocd(p);
    long *SumM= vecalloci(p);
    long *SumO= vecalloci(p);
    double *Sum_Max= vecallocd(p);
    bsp_push_reg(SumW,p*sizeof(double));
    bsp_push_reg(SumM,p*sizeof(long));
    bsp_push_reg(SumO,p*sizeof(long));
    bsp_push_reg(Sum_Max,p*sizeof(double));
    bsp_sync();

    double sumw= 0.0;
    long nmatch_int= 0;
    long nmatch_ext= 0;
    for (long k=0; k<nmatch; k++){
        long e= match[k];
        if (e < nedges){
            sumw += 2.0*weightnew[e]; // weight of internal match
            nmatch_int++;             // is counted twice
        } else {
            sumw += weightnew[e]; // weight of external match
            nmatch_ext++;         // is counted on two processors
        }
    }
    long summ= 2*nmatch_int + nmatch_ext;
    printf("Proc %ld: found %ld internal matches and %ld external matches\n",
           s, nmatch_int, nmatch_ext);
    fflush(stdout);

    bsp_put(0,&sumw, SumW, s*sizeof(double), sizeof(double));
    bsp_put(0,&summ, SumM, s*sizeof(long), sizeof(long));
    bsp_put(0,&nops, SumO, s*sizeof(long), sizeof(long));
    bsp_put(0,&sum_max, Sum_Max, s*sizeof(double), sizeof(double));
    bsp_sync();

    if (s==0){
        double totalW= 0.0;
        long totalM= 0;
        long totalO= 0;
        double totalMax= 0.0;
        for (long t=0; t<p; t++){
            totalW += SumW[t];
            totalM += SumM[t];
            totalO += SumO[t];
            totalMax += Sum_Max[t];
            printf("Proc %ld: number of operations performed = %ld\n", t, SumO[t]);
        }
        totalW /= 2.0;    // every match was counted twice
        if (totalM%2==1)
            printf("Error: total external match count is odd\n");
        totalM /= 2;  
        totalMax /= 2.0;  // every vertex contributes the maximum possible weight
                          // of a half-edge
        suma /= 2.0;      // every edge occurs twice in adjacency matrix

        /* Print the total amount of work, measured as a number of operations,
           and other statistics */
        printf("\nTotal number of operations = %ld\n", totalO);
        printf("Total number of supersteps = %ld\n", nsteps);
        printf("Total number of matches = %ld\n", totalM);
        printf("Total matching weight= %lf\n", totalW);

        /* Compute two upper bounds on the matching weight as a sanity check:
           half the sum of the maximum weights of the vertices
           and the total edge weight.*/
        printf("Upper bound on matching weight= %lf\n", totalMax);
        printf("Total edge weight= %lf\n", suma);
        if (totalW > totalMax)
            printf("Error: matching weight > 0.5 * sum of maximum weights\n");
        if (totalW > suma)
            printf("Error: matching weight > total edge weight\n");

    }
    bsp_pop_reg(Sum_Max);
    bsp_pop_reg(SumO);
    bsp_pop_reg(SumM);
    bsp_pop_reg(SumW);
    bsp_sync();

    vecfreed(Sum_Max);    vecfreei(SumO);
    vecfreei(SumM);       vecfreed(SumW);       vecfreei(match);
    vecfreei(destproc);   vecfreei(weight1new); vecfreed(weightnew);
    vecfreei(janew);      vecfreei(v1new);      vecfreei(v0new);
    vecfreei(Adj);        vecfreei(degree);     vecfreei(Start);
    vecfreei(rowindex);   vecfreei(colindex);

    bsp_end();
    
} /* end bspmv_test */


int main(int argc, char **argv){

    bsp_init(bspmatch_test, argc, argv);

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
    bspmatch_test();

    /* Sequential part */
    exit(EXIT_SUCCESS);

} /* end main */
