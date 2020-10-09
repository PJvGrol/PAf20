#include "bspedupack.h"

#define PROPOSE 0 // possible tags
#define ACCEPT 1
#define REJECT 2

#define SPLITMIN 5 // minimum array size
                   // for splitting to be worthwhile

bool heavier(long e0, long e1, double *weight, long *weight1){

    /* This function checks whether edge e0 is heavier than edge e1
       using weight as a primary criterion and weight1 for breaking
       ties. */

    if (e0 == DUMMY)
        return false;

    if (e1 == DUMMY || weight[e0] > weight[e1] || 
       (weight[e0] == weight[e1] && weight1[e0] > weight1[e1]))
        return true;

    return false;

} /* end heavier */


void find_alive (long v, long *Adj, long nedges, double *weight, 
                 long *weight1, long *v0, long *v1, bool *Alive, 
                 long *Suitor, long lo, long *d){

    /* This function finds the highest index i of a living edge
       e=Adj[i] in the adjacency list of vertex v, with range
       [lo,lo+dv-1], where dv is the degree of v. It removes dead
       edges along the way and reduces the degree accordingly.

       Input:
           v is the local vertex number,
           Adj = adjacency list of vertex v, of length d[v],
           nedges is the number of internal (nonhalo) edges,
           weight[e] is the weight of edge e,
               where for internal edges, e  = (v0[e],v1[e]),
           weight1[e] is a secondary weight used for breaking ties,
           Alive[e] is a boolean stating whether edge e is alive,
           Suitor[u] = edge number of suitor of local vertex u.

       Output:
           d[v] is the new vertex degree, d[v]>=0,
           Alive[Adj[lo+d[v]-1]] holds if d[v]>0.
    */

    for (long i= lo+d[v]-1; i>=lo; i--){
        long e = Adj[i];
            
        /* Kill internal edge if it cannot be matched any more */
        if (e < nedges){
            /* Determine other endpoint of edge e */
            long u = (v0[e]==v ? v1[e] : v0[e]);
            if ((d[u]==0 && Suitor[v]!=e) ||
                heavier(Suitor[u], e, weight, weight1))
                Alive[e] = false;
        }

        if (Alive[e]){
           return;
        } else {
            d[v]--;
        }
    }

} /* end find_alive */ 


long find_splitter (long v, long *Adj, long nedges,
                    double *weight, long *weight1,
                    long *v0, long *v1, bool *Alive,
                    bool *Splitter, long *Suitor,
                    long lo, long *d){

    /* This function finds the highest splitter r of the adjacency
       list of vertex v, with range [lo,lo+d[v]-1], where d[v] is
       the degree of v. If no splitter exists, the function
       returns r=lo.

       The function removes dead edges along the way and reduces
       the degree accordingly. On output, it holds that
           lo <= r < lo+d[v].

       Input and output parameters are the same as for find_alive,
       except for an extra boolean array Splitter (of the same
       length as Adj) which denotes whether an index is a splitter.
       
    */

    for (long r= lo+d[v]-1; r>=lo; r--){
        long e = Adj[r];

        /* Kill internal edge if it cannot be matched any more */
        if (e < nedges){
            /* Determine other endpoint of edge e */
            long u = (v0[e]==v ? v1[e] : v0[e]);
            if ((d[u]==0 && Suitor[v]!=e) ||
                heavier(Suitor[u], e, weight, weight1))
                Alive[e] = false;
        }

        if (!Alive[e]){
            Adj[r] = Adj[lo+d[v]-1]; // copy highest alive
            d[v]--;
        }

        if (Splitter[r])
            return r;
    }

    return lo; // in case no splitter was found

} /* end find_splitter */


void swap (long *x, long i, long j){

    /* This function swaps x[i] and x[j] */

    long tmp = x[i];
    x[i] = x[j];
    x[j] = tmp;

} /* end swap */


void find_pref(long *Adj, double *weight, long *weight1,
               long lo, long hi){

    /* This function finds the maximum weight[Adj[i]] with
       lo <= i <= hi and swaps Adj[i] and Adj[hi]. Here, weight1
       is a secondary weight criterion used for breaking ties. */

    if (hi <= lo )
        return;

    long imax= hi;
    long emax= Adj[hi];
    for (long i=lo; i<hi; i++){
        if (heavier(Adj[i], emax, weight, weight1)){
            imax= i;
            emax= Adj[i];
        }
    }
    swap(Adj, imax, hi);

    return;

} /* end find_pref */


long split_adj(long *Adj, double *weight, long *weight1,
               long lo, long hi){

    /* This function splits the range [lo,hi] with lo <= hi
       of the adjacency list Adj and returns a splitter r,
       such that
           weight[Adj[i]] < weight[Adj[j]], for all i < r <= j.
    */

    if (hi-lo >= 2){
        long piv = (lo+hi)/2; // a simple random-like pivot

        /* Swap Adj[piv] and Adj[hi] */
        swap(Adj, piv, hi);
        long epiv= Adj[hi];
        long i = lo; // first free position
                     // for values < pivot weight
        for (long j=lo; j<hi; j++){
            /* Loop invariant: lo <= i <= j < hi,
                               weight[Adj[lo:i-1]] < wpiv,
                               weight[Adj[i:j-1]] >= wpiv */

           if (heavier(epiv, Adj[j], weight, weight1)){
                swap(Adj, i, j);
                i++;
           }
        }
        swap(Adj, i, hi);
        return i;

    } else if (hi-lo == 1){
        if (heavier(Adj[lo], Adj[hi], weight, weight1))
            swap(Adj, lo, hi);
        return hi;
    } else if (hi==lo){
        return hi;
    }

    return DUMMY; // in case hi < lo

} /* end split_adj */


void push(long v, long nvertices, long q_lo, long *nq, long *Q,
          long *Pref){

    /* This function pushes vertex v onto the queue
       and sets its preference to DUMMY.

       nvertices= number of vertices,
       q_lo= start of the queue,
       nq= number of queue entries,
       Q=  array of size nvertices storing the queue,
       Pref= array of size nvertices storing preferences.
    */

    long q_hi= q_lo + (*nq); // first free position
    if (q_hi >= nvertices)
        q_hi -= nvertices;
    Q[q_hi]= v;
    (*nq)++;
    Pref[v]= DUMMY;

} /* end push */


long pop(long nvertices, long *q_lo, long *nq, long *Q){

    /* This function pops a vertex v from the queue */

    long i= *q_lo; 
    (*q_lo)++;
    if (*q_lo >= nvertices)
        *q_lo -= nvertices;
    (*nq)--;

    return Q[i];

} /* end pop */


void reject_suitor(long v, long e, long q_lo, long *nq, long *Q,
                   long nvertices, long nedges,
                   long *v0, long *v1, long *destproc,
                   bool *Alive, long *Pref){

    /* This function rejects suitor e of vertex v.

       q_lo, nq, Q, nvertices, Pref are the same as
       in the push function.

       nedges= number of internal edges,
       For local edge e < nedges:
           v0[e], v1[e] are the vertices of edge e.
       For halo edge e >= nedges:
           v0[e]= the local vertex of edge e,
           v1[e]= the local edge number on the remote processor
                  for edge e.
       destproc[e-nedges]= the remote processor that
           owns halo edge e. The shift by nedges is to
           save memory and store destproc only for halo edges.
       Alive[e]= boolean that says whether edge e is still alive.
    */

    if (e == DUMMY)
        return;

    if (e < nedges){
        /* Determine other endpoint of edge e */
        long x = (v0[e]==v ? v1[e] : v0[e]);
        push(x, nvertices, q_lo, nq, Q, Pref);
    } else {
        long tag= REJECT;
        bsp_send(destproc[e-nedges], &tag, &(v1[e]), sizeof(long));
    }
    Alive[e]= false;

} /* end reject_suitor */


void bsp_process_recvd_msgs(long q_lo, long *nq, long *Q,
                            long *nmatch, long *match,
                            long nvertices, long nedges,
                            long *v0, long *v1, long *destproc,
                            double *weight, long *weight1, 
                            bool *Alive, long *Suitor, long *Pref,
                            long *degree){

    /* This function processes the messages received at the start
       of a superstep. The messages can be of three types:
       a proposal for a match (with a tag 0), acceptance (tag 1),
       or rejection (tag 2).

       In case of a proposal, the proposer either becomes
       the new suitor, a match, or the proposal is rejected. 
       An acceptance message leads to the registration of a match. 
       A rejection leads to the proposer being pushed
       back onto the queue. 

       q_lo, nq, Q, nvertices, Pref are the same
           as in the push function.
       v0, v1, destproc, Alive are the same as in the
           reject_suitor function.
       nmatch= number of matches registered so far,
       match[i]= edge number of match i,
       weight[e]= weight of edge e,
       weight1[e]= secondary weight of edge e,
                   used for breaking ties,
       Suitor[v]= edge corresponding to the suitor of vertex v,
       degree[v]= the degree of vertex v.
    */

    bsp_nprocs_t nmessages; // total number of messages received
    bsp_size_t nbytes;      // total size in bytes received
    bsp_qsize(&nmessages,&nbytes);

    for (long i=0; i<nmessages; i++){
        bsp_size_t status; // not used
        long tag;      // tag of received message
        bsp_get_tag(&status, &tag);

        long e; // local number of halo edge e received 
        bsp_move(&e, sizeof(long));
        long v= v0[e]; // local vertex of halo edge e

        if (tag==PROPOSE){

            /* Register a match if the preference is mutual */
            if (e==Pref[v]){ 
                /* v has proposed to e in the previous superstep.
                   No need to send an accept. */
                match[*nmatch]= e;
                (*nmatch)++;
                degree[v]= 0;
            }

            /* Assign new suitor */
            long e0= Suitor[v]; // previous suitor
            if (heavier(e, e0, weight, weight1)){
                Suitor[v]= e;
                reject_suitor(v,e0,q_lo,nq,Q,nvertices,nedges,
                              v0,v1,destproc,Alive,Pref);
            } else { // reject the proposal
                long tag_new= REJECT;
                bsp_send(destproc[e-nedges], &tag_new,
                         &(v1[e]), sizeof(long));
                Alive[e]= false;
            }

        } else if (tag==ACCEPT){
            match[*nmatch]= e;
            (*nmatch)++;
            degree[v]= 0;

            /* Reject previous suitor */
            long e0= Suitor[v];
            Suitor[v]= e; // so future proposers know
            reject_suitor(v,e0,q_lo,nq,Q,nvertices,nedges,
                          v0,v1,destproc,Alive,Pref);

        } else if (tag==REJECT){
            push(v, nvertices, q_lo, nq, Q, Pref);
            Alive[e]= false;
        }
    }

} /* end bsp_process_recvd_msgs */


void bspmatch(long nvertices, long nedges, long nhalo,
              long *v0, long *v1, long *destproc,
              double *weight, long *weight1, long *Adj,
              long *Start, long *degree, long maxops,
              long *nmatch, long *match, long *nsteps,
              long *nops){

    /* This function matches the vertices of a graph
       by a local domination algorithm. This guarantees
       obtaining at least half the maximum possible weight.

       Input: 
           nvertices is the number of local vertices,
           nedges is the number of internal (nonhalo) edges,
           nhalo is the number of external (halo) edges,
           degree[v] = degree of local vertex v,
                       0 <= v < nvertices,
           weight[e] = weight of edge e,
           weight1[e] = secondary weight of edge e,
                        used for breaking ties.

           For internal edges, 0 <= e < nedges:
           v0[e] = local vertex number of minimum
                   endpoint of edge e,
           v1[e] = local vertex number of maximum
                   endpoint of edge e.
           
           For external edges, nedges <= e < nedges+nhalo:
           v0[e] = vertex number of the local endpoint of edge e,
           v1[e] = local edge number on remote processor of edge e,
           destproc[e-nedges] = destination processor that owns
                   the nonlocal vertex of edge e.  
    
           Adj is an array of length 2*nedges+nhalo, which contains
           the adjacency lists of all the local vertices. 
           The list of local vertex v contains the local
           indices of the edges incident to v.
           This list is initially stored in positions
           Start[v]..Start[v+1]-1.
           Adj may change during the matching.

           maxops is the maximum number of elemental
           operations carried out in a superstep, which
           is a measure of the amount of work.
           Once maxops is reached, a synchronization is called.
           maxops=0 means no maximum number is imposed.

       Output:
           nmatch is the number of matches found locally,
           match[i] = local edge number of match i,
                      0 <= i < nmatch,
           nsteps is the number of mixed supersteps taken,
                   each consisting of processing the received
                   messages and then setting preferences,
           nops is the total number of elemental operations carried
                   out locally.
    */

    long p= bsp_nprocs(); // p = number of processors obtained
    long s= bsp_pid();    // s = processor number

    /* Allocate, register, and set tag size */
    long *Done= vecalloci(p);
    bsp_push_reg(Done,p*sizeof(long));

    bsp_size_t tagsize= sizeof(long);
    bsp_set_tagsize(&tagsize);
    bsp_sync();

    /* Initialize vertices */
    long *Suitor= vecalloci(nvertices);
    long *Pref= vecalloci(nvertices);
    long *Q= vecalloci(nvertices);
    for (long v=0; v<nvertices; v++){
        Suitor[v]= DUMMY;
        Pref[v]= DUMMY;
        Q[v]= v;
    }

    /* Initialize edges */
    bool *Alive= vecallocb(nedges+nhalo);
    bool *Splitter= vecallocb(2*nedges+nhalo); // same size as Adj
    for (long e=0; e<nedges+nhalo; e++)
        Alive[e]= true;
    for (long i=0; i<2*nedges+nhalo; i++)
        Splitter[i]= false;

    long q_lo=0;        // start of the queue
    long nq= nvertices; // number of vertices in the queue
    *nmatch= 0;         // number of matches registered
    *nsteps= 0;         // number of supersteps taken
    *nops= 0;           // number of operations

    bool alldone= false;

    while (!alldone){
        long nops_step= 0; // number of operations
                           // of this superstep

        /* Initialize all processors to not done yet */
        for (long t=0; t<p; t++)
            Done[t]= false;

        /* Determine if the local processor is done
           for this superstep */
        bsp_nprocs_t nmessages; // total number of messages 
        bsp_size_t nbytes;      // total size in bytes
        bsp_qsize(&nmessages,&nbytes);

        if (nmessages==0 && nq==0){
            long done= true;
            for (long t=0; t<p; t++)
                bsp_put(t,&done,Done,s*sizeof(long),sizeof(long));
        } else {
            bsp_process_recvd_msgs(q_lo,&nq,Q,nmatch,
               match,nvertices,nedges,v0,v1,destproc,
               weight,weight1,Alive,Suitor,Pref,degree);

            while (nq > 0 && (maxops==0 || nops_step < maxops)){
                long v= pop(nvertices, &q_lo, &nq, Q);

                /* Find highest living edge */
                if (degree[v] > 0){
                    long degree_old= degree[v]; // it may decrease
                    find_alive (v,Adj,nedges,weight,weight1,v0,v1,
                                Alive,Suitor,Start[v],degree);
                    nops_step += degree_old-degree[v]+1;
                }

                /* Find highest splitter r */
                long r= Start[v];
                if (degree[v] > 0){
                    long hi= Start[v]+degree[v]-1;
                    r= find_splitter (v,Adj,nedges,weight,weight1,
                                       v0,v1,Alive,Splitter,Suitor,
                                       Start[v],degree);
                    nops_step += hi-r+1; 
                }

                /* Find preference */
                if (degree[v] > 0){
                    // Start[v] <= r <= hi
                    long hi= Start[v]+degree[v]-1;
                    find_pref(Adj, weight, weight1, r, hi); 
                    nops_step += hi-r+1;
                    long e= Adj[hi];
                    Pref[v]= e;

                    /* Remove the top entry from
                       the adjacency list of v */
                    degree[v]--;

                    /* Register a match or propose */
                    if (e==Suitor[v]){
                        match[*nmatch]= e;
                        (*nmatch)++;
                        if (e < nedges){ // internal edge
                            degree[v0[e]]= 0;
                            degree[v1[e]]= 0;
                        } else {
                            degree[v]= 0;
                            long tag= ACCEPT;
                            bsp_send(destproc[e-nedges], &tag,
                                     &(v1[e]), sizeof(long)); 
                        }
                    } else if (e >= nedges){
                        long tag= PROPOSE;
                        bsp_send(destproc[e-nedges], &tag,
                                 &(v1[e]), sizeof(long));
                    }

                    /* Replace the previous suitor
                       of the preference */
                    if (e < nedges){
                        /* Determine other endpoint of edge e */
                        long u = (v0[e]==v ? v1[e] : v0[e]);
                        long e0= Suitor[u]; 
                        Suitor[u]= e;
                        reject_suitor(u,e0,q_lo,&nq,
                                      Q,nvertices,nedges,
                                      v0,v1,destproc,Alive,Pref);
                    }

                    /* Split the top part of the adjacency list
                       of v if the part is large enough */ 
                    hi= Start[v]+degree[v]-1; // new hi
                    if (hi-r >= SPLITMIN){
                        long r_new= split_adj(Adj, weight,
                                              weight1, r, hi);
                        nops_step += hi-r+1;
                        Splitter[r_new]= true;
                    }
                }
            }
        }
        *nops += nops_step;
        (*nsteps)++; 
        bsp_sync();

        /* Determine if the algorithm has terminated */
        alldone= true;
        for (long t=0; t<p; t++){
            if(Done[t] == false){
                alldone= false;
                break;
            }
        }
    }
    bsp_pop_reg(Done);

    vecfreeb(Splitter); vecfreeb(Alive);
    vecfreei(Q);        vecfreei(Pref);
    vecfreei(Suitor);   vecfreei(Done);

} /* end bspmatch */
