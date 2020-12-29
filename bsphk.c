#include "stdlib.h"
#include "time.h"
#include "bspedupack.h"

/*  This program computes the primes up to a given number
    using a prime number sieve.
    The distribution of the sieve is by block
*/

long M; // number of vertices
long N; // approximation of matrix denseness (1=~ 20%, 2=~50%, 3=~ 80%)
long P; // number of processors requested

void bsphk(){
    
    bsp_begin(P);

    double startTime = bsp_time();

    long p= bsp_nprocs(); // p = number of processors
    long s= bsp_pid();    // s = processor number
    long m= M;
    long n= N;
    long counter= 0;

    // storing matrix as vector
    // with 0..(m-1) the connections from 1 in u to 1..m in v
    long *edges= vecalloci(m * m);
    bsp_push_reg(edges, m * m * sizeof(long));

    // u stores for the u nodes with v nodes they are connected to.
    long *u = vecalloci(m);
    long *v = vecalloci(m);
    bsp_push_reg(u, m * sizeof(long));
    bsp_push_reg(v, m * sizeof(long));

    bsp_sync();

    for (long i = 0; i < m * m; i++){
        edges[i] = 0;
    }

    if (s == 0){
        srand(time(NULL));
        
        for (long r = 0; r < m; r++){
            for (long c = 0; c < m; c++){
                long temp = rand();

                if (n == 1){
                    if (temp % 5 == 4){
                        edges[r * m + c] = 1;
                        counter++;
                    }
                }
                else if (n == 2){
                    if (temp % 2 == 1){
                        edges[r * m + c] = 1;
                        counter++;
                    }
                }
                else{
                    if (temp % 5 > 0){
                        edges[r * m + c] = 1;
                        counter++;
                    }
                }
            }
        }

        for (int i = 0; i < p; i++){
            bsp_put(i, edges, edges, 0, m * m * sizeof(long));
        }
    }

    bsp_sync();
    bsp_pop_reg(edges);

    for (long i = 0; i < m; i++){
        u[i] = -1;
        v[i] = -1;
    }

    long *ul = vecalloci(m);
    long *vl = vecalloci(m);

    bool *visitedV = vecallocb(m);

    bool *currentLayerU = vecallocb(m);
    bool *currentLayerV = vecallocb(m);

    bool *currentLayersU = vecallocb(p * m);
    bool *currentLayersV = vecallocb(p * m);
    bool *bfsDones = vecallocb(p);
    bsp_push_reg(currentLayersU, p * m * sizeof(bool));
    bsp_push_reg(currentLayersV, p * m * sizeof(bool));
    bsp_push_reg(bfsDones, p * sizeof(bool));
    bsp_sync();

    bool *layerUHasEdges = vecallocb(m);
    bool *layerVHasEdges = vecallocb(m);

    for (long i = 0; i < m; i++){
        layerUHasEdges[i] = false;
        layerVHasEdges[i] = false;
    }

    bool *layerVFree = vecallocb(m);

    long *bfsLayers = vecalloci(2 * m);
    long *bfsResult = vecalloci(m * m);

    long *candidatePath = vecalloci(2 * m);

    long previousMatchings = 0;

    for (long i = 0; i < m; i++){
        for (long j = 0; j < m; j++){
            if (edges[m * i + j] == 1){
                layerUHasEdges[i] = true;
                layerVHasEdges[j] = true;
            }
        }
    }

    long maximumMatchings = 0;
    long maximumUMatchings = 0;
    long maximumVMatchings = 0;

    for (long i = 0; i < m; i++){
        if (layerUHasEdges[i]){
            maximumUMatchings++;
        }
        if (layerVHasEdges[i]){
            maximumVMatchings++;
        }
    }

    maximumMatchings = maximumUMatchings;
    if (maximumVMatchings < maximumMatchings){
        maximumMatchings = maximumVMatchings;
    }
    
    if (true){
        bool done = false;

        while (!done){
            // Start with BFS
            // We can do an inital round to get a starting position
            bool bfsDone = false;
            long k = 0;
            long currentIndex = 0;
            long newMatchings = 0;

            bfsLayers[k] = currentIndex;

            for (long i = 0; i < m; i++){
                ul[i] = u[i];
                vl[i] = v[i];
                layerVFree[i] = false;
                visitedV[i] = false;
                currentLayerU[i] = false;
                currentLayerV[i] = false;
            }
            
            for (long i = s; i < m; i+=p){
                // check locally, then communitate
                if (ul[i] == -1 && layerUHasEdges[i]){
                    // bfsResult[currentIndex] = i;
                    currentLayerU[i] = true;
                    // currentIndex++;
                }
            }
            
            // communication step
            // use which are visited to recreate bfsResult
            for (long i = 0; i < p; i++){
                bsp_put(i, currentLayerU, currentLayersU, s * m * sizeof(bool), m * sizeof(bool));
            }
            bsp_sync();

            for (long i = 0; i < m; i++){
                bool currentVertexVisited = currentLayerU[i];
                for (long j = 0; j < p; j++){
                    currentVertexVisited = currentVertexVisited || currentLayersU[j * m + i];
                }
                currentLayerU[i] = currentVertexVisited;
                if (currentVertexVisited){
                    bfsResult[currentIndex] = i;
                    currentIndex++;
                }
            }

            k++;
            bfsLayers[k] = currentIndex;

            for (long i = 0; i < m; i++){
                if (currentLayerU[i]){
                    for (long j = s; j < m; j+=p){
                        if (edges[i * m + j] == 1 && !currentLayerV[j]){
                            // bfsResult[currentIndex] = j;
                            currentLayerV[j] = true;
                            // currentIndex++;

                            if (v[j] == -1){
                                // As soon as we find any free vertex in V we are done
                                // but we still need to complete all other searches in this level of BFS
                                // layerVFree[j] = true;
                                bfsDone = true;
                            }
                        }
                    }
                }
                currentLayerU[i] = false;
            }
            
            // communication step
            // use which are visited to recreate bfsResult
            for (long i = 0; i < p; i++){
                bsp_put(i, currentLayerV, currentLayersV, s * m * sizeof(bool), m * sizeof(bool));
                bsp_put(i, &bfsDone, bfsDones, s * sizeof(bool), sizeof(bool));
            }
            bsp_sync();

            for (long i = 0; i < p; i++){
                bfsDone = bfsDone || bfsDones[i];
            }

            printf("Proc %d has bsfDone %d\n", s, bfsDone);

            for (long i = 0; i < m; i++){
                bool currentVertexVisited = currentLayerV[i];
                for (long j = 0; j < p; j++){
                    currentVertexVisited = currentVertexVisited || currentLayersV[j * m + i];
                }
                currentLayerV[i] = currentVertexVisited;
                if (currentVertexVisited){
                    if (bfsDone){
                        if (v[i] == -1){
                            layerVFree[i] = true;
                            bfsResult[currentIndex] = i;
                            currentIndex++;
                        }
                    }
                    else{
                        bfsResult[currentIndex] = i;
                        currentIndex++;
                    }
                }
            }

            k++;
            bfsLayers[k] = currentIndex;
            
            for (long i = 0; i < m; i++){
                visitedV[i] = currentLayerV[i];
            }

            bool bfsFailed = false;
            
            while (!bfsDone && !bfsFailed){
                // We were previously on the V side of things, so now in U
                // This means we can simply go from V to U via matches
                long previousLayerIndex = bfsLayers[k - 1];
                
                if (k % 2 == 0){
                    for (long i = s; i < m; i+=p){
                        if (currentLayerV[i]){
                            // Since v[i] contains the index of the vertex in U that vertex i in V is connected to, we can use that
                            // Additionally since we always have matched edges in this direction and unmatched in reverse
                            // We need not verify that we haven't visited this vertex in U yet, since we certainly won't have
                            // (As long as we properly verify we don't revisit vertices in V)
                            // bfsResult[currentIndex] = v[i];
                            currentLayerU[v[i]] = true;
                            // currentIndex++;
                        }
                    }
            
                    // communication step
                    // use which are visited to recreate bfsResult
                    for (long i = 0; i < p; i++){
                        bsp_put(i, currentLayerU, currentLayersU, s * m * sizeof(bool), m * sizeof(bool));
                    }
                    bsp_sync();

                    for (long i = 0; i < m; i++){
                        currentLayerV[i] = false;
                        bool currentVertexVisited = currentLayerU[i];
                        for (long j = 0; j < p; j++){
                            currentVertexVisited = currentVertexVisited || currentLayersU[j * m + i];
                        }
                        currentLayerU[i] = currentVertexVisited;
                        if (currentVertexVisited){
                            bfsResult[currentIndex] = i;
                            currentIndex++;
                        }
                    }
                }
                // We were previously on the U side of things, so now in V
                // This means we go from U to V via unmatched edges
                else{
                    for (long i = 0; i < m; i++){
                        if (currentLayerU[i]){
                            // We need to find edges, verify these are unmatched and verify we don't visit already seen vertices in V.
                            // edges[i...i + m - 1] contains all edges from vertex i in U to vertices in V
                            for (long j = s; j < m; j+=p){
                                if (!visitedV[j] && u[i] != j && edges[i * m + j] == 1){
                                    // bfsResult[currentIndex] = j;
                                    currentLayerV[j] = true;
                                    // currentIndex++;

                                    if (v[j] == -1){
                                        // As soon as we find any free vertex in V we are done
                                        // but we still need to complete all other searches in this level of BFS
                                        // layerVFree[j] = true;
                                        bfsDone = true;
                                    }
                                }
                            }
                        }
                    }
            
                    // communication step
                    // use which are visited to recreate bfsResult
                    for (long i = 0; i < p; i++){
                        bsp_put(i, currentLayerV, currentLayersV, s * m * sizeof(bool), m * sizeof(bool));
                        bsp_put(i, &bfsDone, bfsDones, s * sizeof(bool), sizeof(bool));
                    }
                    bsp_sync();

                    for (long i = 0; i < p; i++){
                        bfsDone = bfsDone || bfsDones[i];
                    }

                    for (long i = 0; i < m; i++){
                        currentLayerU[i] = false;
                        bool currentVertexVisited = currentLayerV[i];
                        for (long j = 0; j < p; j++){
                            currentVertexVisited = currentVertexVisited || currentLayersV[j * m + i];
                        }
                        currentLayerV[i] = currentVertexVisited;
                        if (currentVertexVisited){
                            if (bfsDone){
                                if (v[i] == -1){
                                    layerVFree[i] = true;
                                    bfsResult[currentIndex] = i;
                                    currentIndex++;
                                }
                            }
                            else{
                                bfsResult[currentIndex] = i;
                                currentIndex++;
                            }
                        }
                    }
                }

                k++;
                bfsLayers[k] = currentIndex;
            
                for (long i = 0; i < m; i++){
                    visitedV[i] = visitedV[i] || currentLayerV[i];
                }

                bool canGoOn = false;

                for (long i = 0; i < m; i++){
                    canGoOn = canGoOn || currentLayerU[i] || currentLayerV[i];
                }

                printf("Proc %d has done a round of BFS\n", s);

                if (!canGoOn){
                    bfsFailed = true;
                }
            }

            printf("Proc %d has index %d\n", s, currentIndex);

            // We've now completed the BFS, we use DFS to find paths.

            // for (long i = 0; i < k; i++){
            //     printf("Layer %d:\n", i);
            //     long currentLayerIndex = bfsLayers[i];
            //     long nextLayerIndex = bfsLayers[i + 1];
            //     for (long j = currentLayerIndex; j < nextLayerIndex; j++){
            //         printf("Vertex %d\n", bfsResult[j]);
            //     }
            // }
            
            // Note that we've stored the vertices in U that lead to free vertices in V
            // As well as the free vertices in V that we've hit
            if (!bfsFailed){
                for (long i = 0; i < m; i++){
                    // As a first step we look at the step from the free vertex in V
                    // to a vertex in U that may be free
                    // This vertex needs to not have changed status yet
                    // We can make use of the layers that we've found, but still need to check for
                    // any changes that may have happened with a previous dfs
                    long candidatePathIndex = 0;
                    if (layerVFree[i]){
                        // k-2 for we already have layerVFree to check layer k-1
                        long uIndex = 0;
                        long vIndex = i;
                        candidatePath[0] = vIndex;
                        candidatePathIndex = 1;
                        bool dfsDone = false;
                        for (long j = k - 2; j >= 0; j--){
                            long layerStartIndex = bfsLayers[j];
                            long nextLayerIndex = bfsLayers[j + 1];
                            for (long l = layerStartIndex; l < nextLayerIndex; l++){
                                if (j % 2 == 0){
                                    uIndex = bfsResult[l];
                                    vIndex = candidatePath[candidatePathIndex - 1];
                                    if (ul[uIndex] != u[uIndex] || edges[uIndex * m + vIndex] != 1){
                                        continue;
                                    }
                                    candidatePath[candidatePathIndex] = uIndex;
                                    candidatePathIndex++;
                                    if (ul[uIndex] == -1){
                                        dfsDone = true;
                                    }
                                    break;
                                }
                                else{
                                    uIndex = candidatePath[candidatePathIndex - 1];
                                    vIndex = bfsResult[l];
                                    if (vl[vIndex] != v[vIndex] || vl[vIndex] != uIndex || ul[uIndex] != vIndex){
                                        continue;
                                    }
                                    candidatePath[candidatePathIndex] = vIndex;
                                    candidatePathIndex++;
                                    break;
                                }
                            }
                            if (dfsDone){
                                break;
                            }
                        }

                        // make required changes
                        if (dfsDone){
                            for (long j = 0; j < candidatePathIndex; j++){
                                if (j % 2 == 0){
                                    vl[candidatePath[j]] = candidatePath[j + 1];
                                }
                                else{
                                    ul[candidatePath[j]] = candidatePath[j - 1];
                                }
                            }
                        }

                        for (long j = 0; j <= candidatePathIndex; j++){
                            candidatePath[j] = -1;
                        }
                    }
                }
            }

            bool error = false;

            // printf("MAXIMUMMATCHINGS: %d\n", maximumMatchings);

            // printf("INITIAL MATCHING\n");

            // for (long i = 0; i < m; i++){
            //     printf("Vertex %d in U connected to vertex %d in V\n", i, u[i]);
            //     if (edges[m * i + ul[i]] != 1){
            //         error = true;
            //     }
            // }


            for (long i = 0; i < m; i++){
                u[i] = ul[i];
                v[i] = vl[i];
            }
            
            for (long i = 0; i < m; i++){
                if (u[i] != -1){
                    newMatchings++;
                }
            }

            if (newMatchings == previousMatchings || newMatchings == maximumMatchings){
                printf("DONE!\n");
                done = true;
                if (newMatchings != maximumMatchings){
                    printf("MATRIX\n");
                    for (long i = 0; i < m; i++){
                        for (long j = 0; j < m; j++){
                            if (edges[i * m + j] == 1){
                                printf("%d, %d connected\n", i, j);
                            }
                        }
                    }
                    printf("NEW MATCHING\n");

                    for (long i = 0; i < m; i++){
                        printf("Vertex %d in U connected to vertex %d in V\n", i, ul[i]);
                        if (edges[m * i + ul[i]] != 1 && ul[i] != -1){
                            error = true;
                        }
                    }
                }
            }
            
            previousMatchings = newMatchings;

            long maxBfsIndex = bfsLayers[k];

            for (long i = 0; i <= k; i++){
                bfsLayers[i] = -1;
            }
            for (long i = 0; i <= maxBfsIndex; i++){
                bfsResult[i] = -1;
            }

            printf("FINISH ROUND\n");
            if (error){
                printf("ERROR!\n");
            }
        }
    }
    else{

    }

    bsp_pop_reg(u);
    bsp_pop_reg(v);
    vecfreei(edges);
    vecfreei(u);
    vecfreei(v);
    bsp_end();

} /* end bsphk */

int main(int argc, char **argv){

    bsp_init(bsphk, argc, argv);

    /* Sequential part for M*/
    printf("How many vertices on either side of the bipartite graph?\n");
    fflush(stdout);

    scanf("%ld",&M);

    /* Sequential part for N*/
    printf("How many edges per vertex as percentage of possible connections?\n1~20%%, 2~50%%, 3~80%%.\n");
    fflush(stdout);

    scanf("%ld",&N);
    
    /* Sequential part for P*/
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
    bsphk();

    /* Sequential part */
    exit(EXIT_SUCCESS);

} /* end main */
