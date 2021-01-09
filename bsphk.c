#include "stdlib.h"
#include "time.h"
#include "bspedupack.h"

/*  This program uses the Hopcroft-Karp algorithm
    to produce a maximum matching for a bipartite graph
*/

long F; // 0 = no improv, 1 = in-process, 2 = between-rounds, 3 = both
long M; // number of U nodes
long N; // number of V nodes
long D; // approximation of matrix density (1=~ 20%, 2=~50%, 3=~ 80%)
long P; // number of processors requested
long R; // Frequency of postprocessing

void bsphk(){
    
    bsp_begin(P);

    long p= bsp_nprocs(); // p = number of processors
    long s= bsp_pid();    // s = processor number
    long f= F;
    long r= R;
    long m= M;
    long n= N;
    long d= D;
    long counter= 0;

    bsp_push_reg(&f, sizeof(long));
    bsp_push_reg(&m, sizeof(long));
    bsp_push_reg(&n, sizeof(long));
    bsp_push_reg(&r, sizeof(long));
    bsp_sync();

    if (s == 0){
        for (int i = 0; i < p; i++){
            bsp_put(i, &f, &f, 0, sizeof(long));
            bsp_put(i, &m, &m, 0, sizeof(long));
            bsp_put(i, &n, &n, 0, sizeof(long));
            bsp_put(i, &r, &r, 0, sizeof(long));
        }
    }

    bsp_sync();
    bsp_pop_reg(&f);
    bsp_pop_reg(&m);
    bsp_pop_reg(&n);
    bsp_pop_reg(&r);

    bool preProcessing = false;
    bool postProcessing = false;

    if (f == 1){
        preProcessing = true;
    }
    if (f == 2){
        postProcessing = true;
    }
    if (f == 3){
        preProcessing = true;
        postProcessing = true;
    }

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

    bool *visitedVerticesU = vecallocb(m);
    bool *visitedVerticesV = vecallocb(n);

    bool *disallowedVerticesU = vecallocb(m);
    bool *disallowedVerticesV = vecallocb(n);

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


    long newMatchingCount = 0;
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
        disallowedVerticesU[i] = false;
        if (layerUHasEdges[i]){
            maximumUMatchings++;
        }
    }

    for (long i = 0; i < n; i++){
        disallowedVerticesV[i] = false;
        if (layerVHasEdges[i]){
            maximumVMatchings++;
        }
    }

    maxMatchingCount = maximumUMatchings;
    if (maximumVMatchings < maxMatchingCount){
        maxMatchingCount = maximumVMatchings;
    }

    bool done = false;

    long roundCounter = 1;

    long bfsSuperSteps = 0;
    long bfsPostSuperSteps = 0;
    long dfsSuperSteps = 0;

    double totalPostProcessingTime = 0;

    long *nrPathsFound = vecalloci(p);

    long preMatchingCount = 0;

    long totalPaths = 0;
    long acceptedPaths = 0;
    long rejectedPaths = 0;

    for (long i = 0; i < p; i++){
        nrPathsFound[i] = 0;
    }

    double preProcessingStart = bsp_time();

    if (preProcessing){
        // Perform HK within the boundaries of the thread
        bool preProcessingDone = false;

        while (!preProcessingDone){
            bool preBfsDone = false;

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

                    for (long j = s; j < n; j += p){
                        if (edges[i * n + j] == 1){
                            currentVerticesV[j] = true;

                            if (vl[j] == -1){
                                preBfsDone = true;
                            }
                        }
                    }
                }
            }

            bfsLayers[layer] = index;

            for (long i = 0; i < m; i++){
                bool vertexVisited = currentVerticesU[i];

                if (vertexVisited){
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
                bool vertexVisited = currentVerticesV[i];

                if (vertexVisited){
                    if (preBfsDone && vl[i] == -1){
                        finalVerticesV[i] = true;
                        bfsResult[index] = i;
                        index++;
                    }
                    else if (!preBfsDone){
                        bfsResult[index] = i;
                        index++;
                    }
                }
            }

            layer++;
            bfsLayers[layer] = index;

            bool preBfsCanContinue = true;

            while (!preBfsDone && preBfsCanContinue){
                // layer % 2 == 0 means we're in V and need to get to U over a matched edge
                if (layer % 2 == 0){
                    for (long i = s; i < n; i += p){
                        if (currentVerticesV[i] && vl[i] % p == s){
                            currentVerticesU[vl[i]] = true;
                            bfsResult[index] = i;
                            index++;
                        }
                    }

                    for (long i = 0; i < n; i++){
                        currentVerticesV[i] = false;
                    }
                }
                else{
                    for (long i = s; i < m; i += p){
                        if (currentVerticesU[i]){
                            for (long j = s; j < n; j += p){
                                if (!visitedVerticesV[j] && edges[i * n + j] == 1 && ul[i] != j){
                                    currentVerticesV[j] = true;
                                    visitedVerticesV[j] = true;
                                    bfsResult[index] = j;
                                    index++;

                                    if (vl[j] == -1){
                                        preBfsDone = true;
                                        finalVerticesV[j] = true;
                                    }
                                }
                            }
                        }
                    }

                    for (long i = 0; i < m; i++){
                        currentVerticesU[i] = false;
                    }
                }

                layer++;
                bfsLayers[layer] = index;

                preBfsCanContinue = false;

                for (long i = s; i < m; i += p){
                    if (currentVerticesU[i]){
                        preBfsCanContinue = true;
                        break;
                    }
                }

                for (long i = s; i < n; i += p){
                    if (currentVerticesV[i]){
                        preBfsCanContinue = true;
                        break;
                    }
                }
            }

            bool preDfsDone = false;

            long preAugmentingPaths = 0;
            long preAugmentingPathsAccepted = 0;
            long preAugmentingPathsRejected = 0;

            while (!preDfsDone){
                bool prePathFound = false;
                long prePathIndex = 0;
                for (long i = s; i < n; i += p){
                    if (!prePathFound && finalVerticesV[i]){
                        long uIndex = 0;
                        long vIndex = i;

                        path[0] = vIndex;
                        prePathIndex = 1;
                        
                        for (long j = layer - 2; j > -1; j--){
                            long startIndex = bfsLayers[j];
                            long endIndex = bfsLayers[j + 1];

                            for (long k = startIndex; k < endIndex; k++){
                                if (j % 2 == 0){
                                    uIndex = bfsResult[k];
                                    vIndex = path[prePathIndex - 1];

                                    if (ul[uIndex] == u[uIndex] && edges[uIndex * n + vIndex] == 1){
                                        path[prePathIndex] = uIndex;
                                        prePathIndex++;

                                        if (ul[uIndex] == -1){
                                            prePathFound = true;
                                        }

                                        break;
                                    }
                                }
                                else{
                                    uIndex = path[prePathIndex - 1];
                                    vIndex = bfsResult[k];

                                    if (vl[vIndex] == v[vIndex] && vl[vIndex] == uIndex){
                                        path[prePathIndex] = vIndex;
                                        prePathIndex++;

                                        break;
                                    }
                                }
                            }

                            if (prePathFound){
                                break;
                            }
                        }

                        if (prePathFound){
                            break;
                        }
                        else{
                            finalVerticesV[i] = false;
                        }
                    }
                }

                preDfsDone = true;

                if (prePathFound){
                    preAugmentingPaths++;
                    totalPaths++;
                    preDfsDone = false;
                    nrPathsFound[s]++;

                    long pathStartIndex = 0;
                    bool pathAccepted = true;
                    
                    for (long j = 0; j < prePathIndex && pathAccepted; j++){
                        if (j % 2 == 0){
                            if (vl[path[pathStartIndex + j]] != v[path[pathStartIndex + j]]){
                                pathAccepted = false;
                            }
                        }
                        else{
                            if (ul[path[pathStartIndex + j]] != u[path[pathStartIndex + j]]){
                                pathAccepted = false;
                            }
                        }
                    }

                    if (pathAccepted){
                        acceptedPaths++;
                        preAugmentingPathsAccepted++;
                        finalVerticesV[path[pathStartIndex]] = false;

                        for (long j = 0; j < prePathIndex; j++){
                            if (j % 2 == 0){
                                vl[path[pathStartIndex + j]] = path[pathStartIndex + j + 1];
                            }
                            else{
                                ul[path[pathStartIndex + j]] = path[pathStartIndex + j - 1];
                            }
                        }
                    }
                    else{
                        rejectedPaths++;
                        preAugmentingPathsRejected++;
                    }
                }

                for (long i = 0 ; i < prePathIndex + 1; i++){
                    path[i] = -1;
                }
            }

            newMatchingCount = 0;

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
                preProcessingDone = true;
            }

            oldMatchingCount = newMatchingCount;

            long maxBfsIndex = bfsLayers[layer];

            for (long i = 0; i < layer + 1; i++){
                bfsLayers[i] = -1;
            }

            for (long i = 0; i < maxBfsIndex + 1; i++){
                bfsResult[i] = -1;
            }
        }

        long *us = vecalloci(p * m);
        long *vs = vecalloci(p * n);

        bsp_push_reg(us, p * m * sizeof(long));
        bsp_push_reg(vs, p * n * sizeof(long));
        bsp_push_reg(nrPathsFound, p * sizeof(long));
        bsp_sync();

        long tempNrPathsFound = nrPathsFound[s];

        for (long i = 0; i < p; i++){
            bsp_put(i, u, us, s * m * sizeof(long), m * sizeof(long));
            bsp_put(i, v, vs, s * n * sizeof(long), n * sizeof(long));
            bsp_put(i, &tempNrPathsFound, nrPathsFound, s * sizeof(long), sizeof(long));
        }

        bsp_sync();

        for (long i = 0; i < p; i++){
            for (long j = i; j < m; j += p){
                u[j] = us[i * m + j];
            }
            for (long j = i; j < n; j += p){
                v[j] = vs[i * n + j];
            }
        }

        bsp_pop_reg(us);
        bsp_pop_reg(vs);

        vecfreei(us);
        vecfreei(vs);

        bsp_sync();

        if (s == 0){
            long tempMatchingCount = 0;
            for (long i = 0; i < m; i++){
                if (u[i] > -1){
                    tempMatchingCount++;
                }
            }
            preMatchingCount = tempMatchingCount;
            oldMatchingCount = tempMatchingCount;
        }
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
            visitedVerticesU[i] = disallowedVerticesU[i];
        }

        for (long i = 0; i < n; i++){
            vl[i] = v[i];
            currentVerticesV[i] = false;
            visitedVerticesV[i] = disallowedVerticesV[i];
            finalVerticesV[i] = false;
        }

        // Layer 0 and 1 of BFS, cyclic distr.

        for (long i = s; i < m; i += p){
            if (layerUHasEdges[i] && ul[i] == -1 && !visitedVerticesU[i]){
                currentVerticesU[i] = true;

                for (long j = 0; j < n; j++){
                    if (edges[i * n + j] == 1 && !visitedVerticesV[j]){
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
                    if (currentVerticesV[i] && !visitedVerticesU[vl[i]]){
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

        long augmentingPaths = 0;
        long augmentingPathsAccepted = 0;
        long augmentingPathsRejected = 0;
        long dfsRoundCounter = 0;

        while (!allDfsDone && bfsDone){
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

            if (dfsRoundCounter % 2 == 0){
                for (long i = 0; i < p; i++){
                    if (pathsFound[i]){
                        augmentingPaths++;
                        totalPaths++;
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
                            acceptedPaths++;
                            augmentingPathsAccepted++;
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
                        else{
                            rejectedPaths++;
                            augmentingPathsRejected++;
                        }
                    }
                }
            }
            else{
                for (long i = p - 1; i > -1; i--){
                    if (pathsFound[i]){
                        augmentingPaths++;
                        totalPaths++;
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
                            acceptedPaths++;
                            augmentingPathsAccepted++;
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
                        else{
                            rejectedPaths++;
                            augmentingPathsRejected++;
                        }
                    }
                }
            }   

            for (long i = 0 ; i < pathIndex + 1; i++){
                path[i] = -1;
            }

            dfsRoundCounter++;
        }

        if (s == 0){
            printf("In round %ld a total of %ld augmenting paths were found of length %ld\n", roundCounter, augmentingPaths, layer);
            printf("Of these paths %ld were accepted and therefore %ld rejected\n", augmentingPathsAccepted, augmentingPathsRejected);
        }

        double finalDfsTime = bsp_time();

        newMatchingCount = 0;

        for (long i = 0; i < m; i++){
            u[i] = ul[i];

            if (u[i] != -1){
                newMatchingCount++;
            }
        }

        for (long i = 0; i < n; i++){
            v[i] = vl[i];
        }

        double postProcessingStart = bsp_time();

        if (postProcessing && roundCounter % r == 0){
            // This is where the fun begins. Unrestricted BFS
            bool searchFinished = false;

            bool newVerticesDiscovered = false;

            long initialIndex = s;
            long finalIndex = initialIndex;

            while (!searchFinished){
                long layer = 0;
                newVerticesDiscovered = false;

                for (long i = 0; i < m; i++){
                    currentVerticesU[i] = false;
                    visitedVerticesU[i] = disallowedVerticesU[i];
                }

                for (long i = 0; i < n; i++){
                    currentVerticesV[i] = false;
                    visitedVerticesV[i] = disallowedVerticesV[i];
                }

                for (long i = initialIndex; i < m; i += p){
                    finalIndex = i + p;
                    if (layerUHasEdges[i] && ul[i] == -1 && !visitedVerticesU[i]){
                        visitedVerticesU[i] = true;

                        for (long j = 0; j < n; j++){
                            if (!visitedVerticesV[j] && edges[i * n + j] == 1){
                                newVerticesDiscovered = true;
                                currentVerticesV[j] = true;

                                if (vl[j] == -1){
                                    bfsDone = true;
                                }
                            }
                        }
                        if (newVerticesDiscovered){
                            break;
                        }
                    }
                }

                initialIndex = finalIndex;

                bool internalBfsDone = false;

                while (!internalBfsDone && newVerticesDiscovered){
                    newVerticesDiscovered = false;
                    if (layer % 2 == 0){
                        for (long i = 0; i < n; i++){
                            if (currentVerticesV[i]){
                                for (long j = 0; j < m; j++){
                                    if (!visitedVerticesU[j] && edges[j * n + i] == 1){
                                        newVerticesDiscovered = true;
                                        currentVerticesU[j] = true;
                                        visitedVerticesU[i] = true;
                                    }
                                }
                            }
                            currentVerticesV[i] = false;
                        }
                    }
                    else{
                        for (long i = 0; i < m; i++){
                            if (currentVerticesU[i]){
                                for (long j = 0; j < n; j++){
                                    if (!visitedVerticesV[j] && edges[i * n + j] == 1){
                                        newVerticesDiscovered = true;
                                        currentVerticesV[j] = true;
                                        visitedVerticesV[j] = true;

                                        if (vl[j] == -1){
                                            internalBfsDone = true;
                                        }
                                    }
                                }
                            }
                            currentVerticesU[i] = false;
                        }
                    }

                    layer++;
                }

                if (internalBfsDone){
                    for (long i = 0; i < m; i++){
                        visitedVerticesU[i] = false;
                    }

                    for (long i = 0; i < n; i++){
                        visitedVerticesV[i] = false;
                    }
                }

                for (long i = 0; i < p; i++){
                    bsp_put(i, visitedVerticesU, currentVerticesUs, s * m * sizeof(bool), m * sizeof(bool));
                    bsp_put(i, visitedVerticesV, currentVerticesVs, s * n * sizeof(bool), n * sizeof(bool));
                    bsp_put(i, &internalBfsDone, bfsDones, s * sizeof(bool), sizeof(bool));
                }

                bsp_sync();

                searchFinished = true;

                for (long i = 0; i < p; i++){
                    if (bfsDones[i]){
                        searchFinished = false;
                    }
                }

                for (long i = 0; i < m; i++){
                    bool vertexDisallowed = false;
                    for (long j = 0; j < p; j++){
                        if (visitedVerticesU[j * m + i]){
                            vertexDisallowed = true;
                        }
                    }
                    disallowedVerticesU[i] = vertexDisallowed;
                }

                for (long i = 0; i < n; i++){
                    bool vertexDisallowed = false;
                    for (long j = 0; j < p; j++){
                        if (visitedVerticesV[j * n + i]){
                            vertexDisallowed = true;
                        }
                    }
                    disallowedVerticesV[i] = vertexDisallowed;
                }
            }
            
            bool anyVerticesAllowed = false;

            for (long i = 0; i < m; i++){
                if (!disallowedVerticesU[i]){
                    anyVerticesAllowed = true;
                    break;
                }
            }

            if (!anyVerticesAllowed){
                for (long i = 0; i < n; i++){
                    if (!disallowedVerticesV[i]){
                        anyVerticesAllowed = true;
                        break;
                    }
                }
            }

            if (!anyVerticesAllowed){
                done = true;
            }
        }

        double endTime = bsp_time();

        if (postProcessing){
            totalPostProcessingTime += (endTime - postProcessingStart);
        }

        if (s == 0){
            printf("In round %ld, %ld augmenting paths were found, resulting in %ld new matchings\n", roundCounter, augmentingPaths, newMatchingCount - oldMatchingCount);
        }

        if (newMatchingCount == oldMatchingCount || newMatchingCount == maxMatchingCount || done){
            done = true;

            if (s == 0){
                bool error = false;

                printf("DONE!\n");

                if (s == 0 && (true || newMatchingCount != maxMatchingCount)){
                    printf("A total of %ld rounds of HK were required for the solution\n", roundCounter);

                    printf("A total of %ld BFS communication related supersteps occurred\n", bfsSuperSteps);
                    printf("A total of %ld DFS communication related supersteps occurred\n", dfsSuperSteps);

                    printf("In total %ld augmenting paths were found, %ld accepted and %ld rejected\n", totalPaths, acceptedPaths, rejectedPaths);

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
                    if (preProcessing){
                        printf("In preprocessing a matching of size %ld was found\n", preMatchingCount);
                    }

                    printf("Found a maximum matching of size %ld\n", newMatchingCount);

                    printf("Final timing\n");
                    if (f == 0){
                        printf("In total the algorith took %.6lf seconds\n", endTime - outerLoopStartTime);
                    }
                    if (f == 1){
                        printf("In total the preprocessing step took %.6lf seconds\n", outerLoopStartTime - preProcessingStart);
                        printf("In total the algorith took %.6lf seconds\n", endTime - preProcessingStart);
                    }
                    if (f == 2){
                        printf("In total postprocessing took %.6lf seconds\n", totalPostProcessingTime);
                        printf("In total the algorith took %.6lf seconds\n", endTime - outerLoopStartTime);
                    }
                    if (f == 3){
                        printf("In total the preprocessing step took %.6lf seconds\n", outerLoopStartTime - preProcessingStart);
                        printf("In total postprocessing took %.6lf seconds\n", totalPostProcessingTime);
                        printf("In total the algorith took %.6lf seconds\n", endTime - preProcessingStart);
                    }

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

            break;
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

        roundCounter++;
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
    vecfreei(nrPathsFound);

    bsp_end();

} /* end bsphk */

int main(int argc, char **argv){

    bsp_init(bsphk, argc, argv);

    /* Sequential part for F*/
    printf("Which performance improvements? 0 for none, 1 for preprocessing, 2 for postprocessing, 3 for both\n");
    fflush(stdout);

    scanf("%ld",&F);

    if (F > 1){
        /* Sequential part for F*/
        printf("Select the frequency of postprocessing, where 1 means every round, 2 every 2 rounds etc.\n");
        fflush(stdout);

        scanf("%ld",&R);
    }

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
