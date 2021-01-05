#include "stdlib.h"
#include "time.h"
#include "bspedupack.h"

/*  This program uses the Hopcroft-Karp algorithm
    to produce a maximum matching for a bipartite graph
*/

long M; // number of U nodes
long N; // number of V nodes
long D; // approximation of matrix density (1=~ 20%, 2=~50%, 3=~ 80%)
long P; // number of processors requested

void bsphk(){
    
    bsp_begin(P);

    long p= bsp_nprocs(); // p = number of processors
    long s= bsp_pid();    // s = processor number
    long m= M;
    long n= N;
    long d= D;
    long counter= 0;

    bsp_push_reg(&m, sizeof(long));
    bsp_push_reg(&n, sizeof(long));
    bsp_sync();

    if (s == 0){
        for (int i = 0; i < p; i++){
            bsp_put(i, &m, &m, 0, sizeof(long));
            bsp_put(i, &n, &n, 0, sizeof(long));
        }
    }

    bsp_sync();
    bsp_pop_reg(&m);
    bsp_pop_reg(&n);

    // n = m;

    // storing matrix as vector
    // with 0..(n-1) the connections from 1 in u to 1..n in v
    long *edges= vecalloci(m * n);
    bsp_push_reg(edges, m * n * sizeof(long));

    // u stores for the u nodes with v nodes they are connected to.
    long *u = vecalloci(m);
    long *v = vecalloci(n);

    bsp_sync();

    for (long i = 0; i < m * n; i++){
        edges[i] = 0;
    }

    if (s == 0){
        srand(time(NULL));
        
        for (long r = 0; r < m; r++){
            for (long c = 0; c < n; c++){
                long temp = rand();

                if (d == 1){
                    if (temp % 5 == 4){
                        edges[r * n + c] = 1;
                        counter++;
                    }
                }
                else if (d == 2){
                    if (temp % 2 == 1){
                        edges[r * n + c] = 1;
                        counter++;
                    }
                }
                else{
                    if (temp % 5 > 0){
                        edges[r * n + c] = 1;
                        counter++;
                    }
                }
            }
        }

        for (int i = 0; i < p; i++){
            bsp_put(i, edges, edges, 0, m * n * sizeof(long));
        }
    }

    bsp_sync();
    bsp_pop_reg(edges);

    for (long i = 0; i < m; i++){
        u[i] = -1;
    }

    for (long i = 0; i < n; i++){
        v[i] = -1;
    }

    // ul and vl are the local/in progress versions of u and v
    long *ul = vecalloci(m);
    long *vl = vecalloci(n);

    bool *visitedVerticesV = vecallocb(n);

    bool *currentVerticesU = vecallocb(m);
    bool *currentVerticesV = vecallocb(n);

    bool *currentVerticesUs = vecallocb(p * m);
    bool *currentVerticesVs = vecallocb(p * n);
    bool *bfsDones = vecallocb(p);
    long *paths = vecalloci(p * (m + n));
    long *pathIndices = vecalloci(p);
    bool *pathsFound = vecallocb(p);
    bsp_push_reg(currentVerticesUs, p * m * sizeof(bool));
    bsp_push_reg(currentVerticesVs, p * n * sizeof(bool));
    bsp_push_reg(bfsDones, p * sizeof(bool));
    bsp_push_reg(paths, p * (m + n) * sizeof(long));
    bsp_push_reg(pathIndices, p * sizeof(long));
    bsp_push_reg(pathsFound, p * sizeof(bool));
    bsp_sync();

    bool *layerUHasEdges = vecallocb(m);
    bool *layerVHasEdges = vecallocb(n);

    for (long i = 0; i < m; i++){
        layerUHasEdges[i] = false;
    }

    for (long i = 0; i < n; i++){
        layerVHasEdges[i] = false;
    }

    bool *finalVerticesV = vecallocb(n);

    long *bfsLayers = vecalloci(m + n);
    long *bfsResult = vecalloci(m * n);

    long *path = vecalloci(m + n);

    long oldMatchingCount = 0;

    for (long i = 0; i < m; i++){
        for (long j = 0; j < n; j++){
            if (edges[n * i + j] == 1){
                layerUHasEdges[i] = true;
                layerVHasEdges[j] = true;
            }
        }
    }

    long maxMatchingCount = 0;
    long maximumUMatchings = 0;
    long maximumVMatchings = 0;

    for (long i = 0; i < m; i++){
        if (layerUHasEdges[i]){
            maximumUMatchings++;
        }
    }

    for (long i = 0; i < n; i++){
        if (layerVHasEdges[i]){
            maximumVMatchings++;
        }
    }

    maxMatchingCount = maximumUMatchings;
    if (maximumVMatchings < maxMatchingCount){
        maxMatchingCount = maximumVMatchings;
    }

    bool done = false;

    long roundCounter = 0;

    long bfsSuperSteps = 0;
    long dfsSuperSteps = 0;

    long *nrPathsFound = vecalloci(p);

    for (long i = 0; i < p; i++){
        nrPathsFound[i] = 0;
    }

    double outerLoopStartTime = bsp_time();

    while (!done){
        double startTime = bsp_time();

        bool bfsDone = false;

        long layer = 0;
        long index = 0;

        for (long i = 0; i < m; i++){
            ul[i] = u[i];
            currentVerticesU[i] = false;
        }

        for (long i = 0; i < n; i++){
            vl[i] = v[i];
            currentVerticesV[i] = false;
            visitedVerticesV[i] = false;
            finalVerticesV[i] = false;
        }

        // Layer 0 and 1 of BFS, cyclic distr.

        for (long i = s; i < m; i += p){
            if (layerUHasEdges[i] && ul[i] == -1){
                currentVerticesU[i] = true;

                for (long j = 0; j < n; j++){
                    if (edges[i * n + j] == 1){
                        currentVerticesV[j] = true;

                        if (vl[j] == -1){
                            bfsDone = true;
                        }
                    }
                }
            }
        }

        for (long i = 0; i < p; i++){
            bsp_put(i, currentVerticesU, currentVerticesUs, s * m * sizeof(bool), m * sizeof(bool));
            bsp_put(i, currentVerticesV, currentVerticesVs, s * n * sizeof(bool), n * sizeof(bool));
            bsp_put(i, &bfsDone, bfsDones, s * sizeof(bool), sizeof(bool));
        }

        bsp_sync();
        bfsSuperSteps++;

        bfsDone = false;

        for (long i = 0; i < p; i++){
            if (bfsDones[i]){
                bfsDone = true;
            }
        }

        bfsLayers[layer] = index;

        for (long i = 0; i < m; i++){
            bool vertexVisited = false;

            for (long j = 0; j < p; j++){
                if (currentVerticesUs[j * m + i]){
                    vertexVisited = true;
                }
            }

            if (vertexVisited){
                currentVerticesU[i] = true;
                bfsResult[index] = i;
                index++;
            }
        }

        layer++;
        bfsLayers[layer] = index;

        for (long i = 0; i < m; i++){
            currentVerticesU[i] = false;
        }

        for (long i = 0; i < n; i++){
            bool vertexVisited = false;

            for (long j = 0; j < p; j++){
                if (currentVerticesVs[j * n + i]){
                    vertexVisited = true;
                }
            }

            if (vertexVisited){
                visitedVerticesV[i] = true;

                if (bfsDone && vl[i] == -1){
                    finalVerticesV[i] = true;
                    currentVerticesV[i] = true;
                    bfsResult[index] = i;
                    index++;
                }
                else if (!bfsDone){
                    currentVerticesV[i] = true;
                    bfsResult[index] = i;
                    index++;
                }
            }
        }

        layer++;
        bfsLayers[layer] = index;

        // Value of layer is now 2
        
        bool bfsCanContinue = true;

        double initialBfsTime = bsp_time();

        while (!bfsDone && bfsCanContinue){
            // layer % 2 == 0 means we're in V and need to get to U over a matched edge
            if (layer % 2 == 0){
                for (long i = s; i < n; i += p){
                    if (currentVerticesV[i]){
                        currentVerticesU[vl[i]] = true;
                    }
                }

                for (long i = 0; i < p; i++){
                    bsp_put(i, currentVerticesU, currentVerticesUs, s * m * sizeof(bool), m * sizeof(bool));
                }

                bsp_sync();
                bfsSuperSteps++;

                for (long i = 0; i < n; i++){
                    currentVerticesV[i] = false;
                }

                for (long i = 0; i < m; i++){
                    bool vertexVisited = false;

                    for (long j = 0; j < p; j++){
                        if (currentVerticesUs[j * m + i]){
                            vertexVisited = true;
                        }
                    }

                    if (vertexVisited){
                        currentVerticesU[i] = true;
                        bfsResult[index] = i;
                        index++;
                    }
                }
            }
            else{
                for (long i = s; i < m; i += p){
                    if (currentVerticesU[i]){
                        for (long j = 0; j < n; j++){
                            if (!visitedVerticesV[j] && edges[i * n + j] == 1 && ul[i] != j){
                                currentVerticesV[j] = true;

                                if (vl[j] == -1){
                                    bfsDone = true;
                                }
                            }
                        }
                    }
                }

                for (long i = 0; i < p; i++){
                    bsp_put(i, currentVerticesV, currentVerticesVs, s * n * sizeof(bool), n * sizeof(bool));
                    bsp_put(i, &bfsDone, bfsDones, s * sizeof(bool), sizeof(bool));
                }

                bsp_sync();
                bfsSuperSteps++;

                bfsDone = false;

                for (long i = 0; i < p; i++){
                    if (bfsDones[i]){
                        bfsDone = true;
                    }
                }

                for (long i = 0; i < m; i++){
                    currentVerticesU[i] = false;
                }

                for (long i = 0; i < n; i++){
                    bool vertexVisited = false;

                    for (long j = 0; j < p; j++){
                        if (currentVerticesVs[j * n + i]){
                            vertexVisited = true;
                        }
                    }

                    if (vertexVisited){
                        visitedVerticesV[i] = true;

                        if (bfsDone && vl[i] == -1){
                            finalVerticesV[i] = true;
                            currentVerticesV[i] = true;
                            bfsResult[index] = i;
                            index++;
                        }
                        else if (!bfsDone){
                            currentVerticesV[i] = true;
                            bfsResult[index] = i;
                            index++;
                        }
                    }
                }
            }

            layer++;
            bfsLayers[layer] = index;

            bfsCanContinue = false;

            for (long i = 0; i < m; i++){
                if (currentVerticesU[i]){
                    bfsCanContinue = true;
                    break;
                }
            }

            for (long i = 0; i < n; i++){
                if (currentVerticesV[i]){
                    bfsCanContinue = true;
                    break;
                }
            }
        }

        double finalBfsTime = bsp_time();

        bool allDfsDone = false;

        while (!allDfsDone){
            bool pathFound = false;
            long pathIndex = 0;
            for (long i = s; i < n; i += p){
                if (!pathFound && finalVerticesV[i]){
                    long uIndex = 0;
                    long vIndex = i;

                    path[0] = vIndex;
                    pathIndex = 1;
                    
                    for (long j = layer - 2; j > -1; j--){
                        long startIndex = bfsLayers[j];
                        long endIndex = bfsLayers[j + 1];

                        for (long k = startIndex; k < endIndex; k++){
                            if (j % 2 == 0){
                                uIndex = bfsResult[k];
                                vIndex = path[pathIndex - 1];

                                if (ul[uIndex] == u[uIndex] && edges[uIndex * n + vIndex] == 1){
                                    path[pathIndex] = uIndex;
                                    pathIndex++;

                                    if (ul[uIndex] == -1){
                                        pathFound = true;
                                    }

                                    break;
                                }
                            }
                            else{
                                uIndex = path[pathIndex - 1];
                                vIndex = bfsResult[k];

                                if (vl[vIndex] == v[vIndex] && vl[vIndex] == uIndex){
                                    path[pathIndex] = vIndex;
                                    pathIndex++;

                                    break;
                                }
                            }
                        }

                        if (pathFound){
                            break;
                        }
                    }

                    if (pathFound){
                        break;
                    }
                    else{
                        finalVerticesV[i] = false;
                    }
                }
            }

            for (long i = 0; i < p; i++){
                bsp_put(i, path, paths, s * (m + n) * sizeof(long), (m + n) * sizeof(long));
                bsp_put(i, &pathIndex, pathIndices, s * sizeof(long), sizeof(long));
                bsp_put(i, &pathFound, pathsFound, s * sizeof(bool), sizeof(bool));
            }

            bsp_sync();
            dfsSuperSteps++;

            allDfsDone = true;

            for (long i = 0; i < p; i++){
                if (pathsFound[i]){
                    allDfsDone = false;
                    nrPathsFound[i]++;

                    long pathStartIndex = i * (m + n);
                    bool pathAccepted = true;
                    
                    for (long j = 0; j < pathIndices[i] && pathAccepted; j++){
                        if (j % 2 == 0){
                            if (vl[paths[pathStartIndex + j]] != v[paths[pathStartIndex + j]]){
                                pathAccepted = false;
                            }
                        }
                        else{
                            if (ul[paths[pathStartIndex + j]] != u[paths[pathStartIndex + j]]){
                                pathAccepted = false;
                            }
                        }
                    }

                    if (pathAccepted){
                        finalVerticesV[paths[pathStartIndex]] = false;

                        for (long j = 0; j < pathIndices[i]; j++){
                            if (j % 2 == 0){
                                vl[paths[pathStartIndex + j]] = paths[pathStartIndex + j + 1];
                            }
                            else{
                                ul[paths[pathStartIndex + j]] = paths[pathStartIndex + j - 1];
                            }
                        }
                    }
                }
            }

            for (long i = 0 ; i < pathIndex + 1; i++){
                path[i] = -1;
            }
        }

        bool finalDfsTime = bsp_time();

        roundCounter++;

        long newMatchingCount = 0;

        for (long i = 0; i < m; i++){
            u[i] = ul[i];

            if (u[i] != -1){
                newMatchingCount++;
            }
        }

        for (long i = 0; i < n; i++){
            v[i] = vl[i];
        }

        if (newMatchingCount == oldMatchingCount || newMatchingCount == maxMatchingCount){
            done = true;

            if (s == 0){
                bool error = false;

                printf("DONE!\n");

                if (s == 0 && (true || newMatchingCount != maxMatchingCount)){
                    printf("A total of %ld rounds of HK were required for the solution\n", roundCounter);

                    printf("A total of %ld BFS communication related supersteps occurred\n", bfsSuperSteps);
                    printf("A total of %ld DFS communication related supersteps occurred\n", dfsSuperSteps);

                    for (long i = 0; i < p; i++){
                        printf("Proc %ld found a total of %ld augmenting paths\n", i, nrPathsFound[i]);
                    }

                    printf("MATRIX\n");

                    printf("%ld connections out of %ld possible\n", counter, m * n);

                    // for (long i = 0; i < m; i++){
                    //     for (long j = 0; j < n; j++){
                    //         if (edges[i * n + j] == 1){
                    //             printf("%ld, %ld connected\n", i, j);
                    //         }
                    //     }
                    // }

                    printf("MATCHING\n");

                    printf("Found a maximum matching of size %ld\n", newMatchingCount);

                    // printf("NEW MATCHING\n");

                    // for (long i = 0; i < m; i++){

                    //     printf("Vertex %ld in U connected to vertex %ld in V\n", i, ul[i]);

                    //     if (edges[n * i + ul[i]] != 1 && ul[i] != -1){
                    //         error = true;
                    //     }
                    // }
                }
                if (error){
                    printf("ERROR!\n");
                }
            }
        }

        oldMatchingCount = newMatchingCount;

        long maxBfsIndex = bfsLayers[layer];

        for (long i = 0; i < layer + 1; i++){
            bfsLayers[i] = -1;
        }

        for (long i = 0; i < maxBfsIndex + 1; i++){
            bfsResult[i] = -1;
        }

        for (long i = 0; i < p; i++){
            bfsDones[i] = false;
            pathsFound[i] = false;
        }
    }

    bsp_pop_reg(currentVerticesUs);
    bsp_pop_reg(currentVerticesVs);
    bsp_pop_reg(bfsDones);
    bsp_pop_reg(paths);
    bsp_pop_reg(pathIndices);
    bsp_pop_reg(pathsFound);   

    vecfreei(edges);
    vecfreei(u);
    vecfreei(v);
    vecfreei(ul);
    vecfreei(vl);
    vecfreeb(visitedVerticesV);
    vecfreeb(currentVerticesU);
    vecfreeb(currentVerticesV);
    vecfreeb(currentVerticesUs);
    vecfreeb(currentVerticesVs);
    vecfreeb(bfsDones);
    vecfreei(paths);
    vecfreei(pathIndices);
    vecfreeb(pathsFound);
    vecfreeb(layerUHasEdges);
    vecfreeb(layerVHasEdges);
    vecfreeb(finalVerticesV);
    vecfreei(bfsLayers);
    vecfreei(bfsResult);
    vecfreei(path);

    bsp_end();

} /* end bsphk */

int main(int argc, char **argv){

    bsp_init(bsphk, argc, argv);

    /* Sequential part for M*/
    printf("How many vertices on the left side of the bipartite graph?\n");
    fflush(stdout);

    scanf("%ld",&M);

    /* Sequential part for N*/
    printf("How many vertices on the right side of the bipartite graph?\n");
    fflush(stdout);

    scanf("%ld",&N);

    /* Sequential part for N*/
    printf("How many edges per vertex as percentage of possible connections?\n1~20%%, 2~50%%, 3~80%%.\n");
    fflush(stdout);

    scanf("%ld",&D);
    
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
