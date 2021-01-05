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

    // ul and vl are the local/in progress versions of u and v
    long *ul = vecalloci(m);
    long *vl = vecalloci(m);

    bool *visitedVerticesV = vecallocb(m);

    bool *currentVerticesU = vecallocb(m);
    bool *currentVerticesV = vecallocb(m);

    bool *currentVerticesUs = vecallocb(p * m);
    bool *currentVerticesVs = vecallocb(p * m);
    bool *bfsDones = vecallocb(p);
    long *paths = vecalloci(p * m * 2);
    long *pathIndices = vecalloci(p);
    bool *pathsFound = vecallocb(p);
    bsp_push_reg(currentVerticesUs, p * m * sizeof(bool));
    bsp_push_reg(currentVerticesVs, p * m * sizeof(bool));
    bsp_push_reg(bfsDones, p * sizeof(bool));
    bsp_push_reg(paths, p * m * 2 * sizeof(long));
    bsp_push_reg(pathIndices, p * sizeof(long));
    bsp_push_reg(pathsFound, p * sizeof(bool));
    bsp_sync();

    bool *layerUHasEdges = vecallocb(m);
    bool *layerVHasEdges = vecallocb(m);

    for (long i = 0; i < m; i++){
        layerUHasEdges[i] = false;
        layerVHasEdges[i] = false;
    }

    bool *finalVerticesV = vecallocb(m);

    long *bfsLayers = vecalloci(2 * m);
    long *bfsResult = vecalloci(m * m);

    long *path = vecalloci(2 * m);

    long oldMatchingCount = 0;

    for (long i = 0; i < m; i++){
        for (long j = 0; j < m; j++){
            if (edges[m * i + j] == 1){
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
        if (layerVHasEdges[i]){
            maximumVMatchings++;
        }
    }

    maxMatchingCount = maximumUMatchings;
    if (maximumVMatchings < maxMatchingCount){
        maxMatchingCount = maximumVMatchings;
    }

    printf("Starting HK\n");

    bool done = false;

    while (!done){
        bool bfsDone = false;

        long layer = 0;
        long index = 0;

        for (long i = 0; i < m; i++){
            ul[i] = u[i];
            vl[i] = v[i];
            currentVerticesU[i] = false;
            currentVerticesV[i] = false;
            visitedVerticesV[i] = false;
            finalVerticesV[i] = false;
        }

        // Layer 0 and 1 of BFS, cyclic distr.

        for (long i = s; i < m; i += p){
            if (layerUHasEdges[i] && ul[i] == -1){
                currentVerticesU[i] = true;

                for (long j = 0; j < m; j++){
                    if (edges[i * m + j] == 1){
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
            bsp_put(i, currentVerticesV, currentVerticesVs, s * m * sizeof(bool), m * sizeof(bool));
            bsp_put(i, &bfsDone, bfsDones, s * sizeof(bool), sizeof(bool));
        }

        bsp_sync();

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

            bool vertexVisited = false;

            for (long j = 0; j < p; j++){
                if (currentVerticesVs[j * m + i]){
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

        printf("Round Zero BFS Done\n");

        while (!bfsDone && bfsCanContinue){
            // layer % 2 == 0 means we're in V and need to get to U over a matched edge
            if (layer % 2 == 0){
                for (long i = s; i < m; i += p){
                    if (currentVerticesV[i]){
                        currentVerticesU[vl[i]] = true;
                    }
                }

                for (long i = 0; i < p; i++){
                    bsp_put(i, currentVerticesU, currentVerticesUs, s * m * sizeof(bool), m * sizeof(bool));
                }

                bsp_sync();

                for (long i = 0; i < m; i++){
                    currentVerticesV[i] = false;

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
                        for (long j = 0; j < m; j++){
                            if (!visitedVerticesV[j] && edges[i * m + j] == 1 && ul[i] != j){
                                currentVerticesV[j] = true;

                                if (vl[j] == -1){
                                    bfsDone = true;
                                }
                            }
                        }
                    }
                }

                for (long i = 0; i < p; i++){
                    bsp_put(i, currentVerticesV, currentVerticesVs, s * m * sizeof(bool), m * sizeof(bool));
                    bsp_put(i, &bfsDone, bfsDones, s * sizeof(bool), sizeof(bool));
                }

                bsp_sync();

                for (long i = 0; i < p; i++){
                    if (bfsDones[i]){
                        bfsDone = true;
                    }
                }

                for (long i = 0; i < m; i++){
                    currentVerticesU[i] = false;

                    bool vertexVisited = false;

                    for (long j = 0; j < p; j++){
                        if (currentVerticesVs[j * m + i]){
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
                if (currentVerticesU[i] || currentVerticesV[i]){
                    bfsCanContinue = true;
                    break;
                }
            }
        }

        printf("BFS Done\n");

        bool allDfsDone = false;

        while (!allDfsDone){
            bool pathFound = false;
            long pathIndex = 0;
            for (long i = s; i < m; i += p){
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

                                if (ul[uIndex] == u[uIndex] && edges[uIndex * m + vIndex] == 1){
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
                bsp_put(i, path, paths, s * 2 * m * sizeof(long), 2 * m * sizeof(long));
                bsp_put(i, &pathIndex, pathIndices, s * sizeof(long), sizeof(long));
                bsp_put(i, &pathFound, pathsFound, s * sizeof(bool), sizeof(bool));
            }

            bsp_sync();

            allDfsDone = true;

            for (long i = 0; i < p; i++){
                if (pathsFound[i]){
                    allDfsDone = false;

                    long pathStartIndex = i * 2 * m;
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

        printf("DFS Done\n");

        long newMatchingCount = 0;

        for (long i = 0; i < m; i++){
            u[i] = ul[i];
            v[i] = vl[i];

            if (u[i] != -1){
                newMatchingCount++;
            }
        }

        if (newMatchingCount == oldMatchingCount || newMatchingCount == maxMatchingCount){
            done = true;

            if (s == 0){
                bool error = false;

                printf("DONE!\n");

                if (s == 0 && (true || newMatchingCount != maxMatchingCount)){

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
