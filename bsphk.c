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

        // printf("MATRIX\n");
        // for (long i = 0; i < m; i++){
        //     for (long j = 0; j < m; j++){
        //         if (edges[i * m + j] == 1){
        //             printf("%d, %d connected\n", i, j);
        //         }
        //     }
        // }
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

    bool *layerUFreeV = vecallocb(m);

    long *bfstrack = vecalloci(m * m);
    
    if (p == 1){
        bool done = false;

        while (!done){
            // start with bfs
            bool bfsDone = false;
            long k = 0;
            bfstrack[0] = -1;
            long previousBlockIndex = 0;
            long nextBlockIndex = 1;
            long blockCountIndex = 2;
            long index = 3; // vertex index

            for (long i = 0; i < m; i++){
                ul[i] = u[i];
                vl[i] = v[i];
                layerUFreeV[i] = false;
                visitedV[i] = false;
                currentLayerU[i] = false;
                currentLayerV[i] = false;
            }

            long blockCount = 0;
            
            for (long i = 0; i < m; i++){
                if (ul[i] == -1){
                    long edgeCountIndex = index + 1;
                    long edgeCount = 0;
                    long tempIndex = index + 2;
                    bool connected = false;
                    for (long j = 0; j < m; j++){
                        if (edges[i * m + j] == 1){
                            connected = true;
                            bfstrack[edgeCountIndex - 1] = i;
                            if (vl[j] == -1){
                                layerUFreeV[i] = true;
                                bfsDone = true;
                            }
                            currentLayerV[j] = true;
                            bfstrack[tempIndex] = j;
                            edgeCount++;
                            tempIndex++;
                        }
                    }
                    if (connected){
                        blockCount++;
                        index = tempIndex;
                    }
                    bfstrack[edgeCountIndex] = edgeCount;
                }
            }

            for (long i = 0; i < m; i++){
                visitedV[i] = currentLayerV[i];
            }

            bfstrack[index] = previousBlockIndex;
            bfstrack[nextBlockIndex] = index;
            bfstrack[blockCountIndex] = blockCount;

            previousBlockIndex = index;
            nextBlockIndex = index + 1;
            blockCountIndex = index + 2;
            index+=3;
            k++;

            blockCount = 0;
            
            while (!bfsDone){
                for (long i = 0; i < m; i++){
                    if (k % 2 == 1){
                        if (currentLayerV[i] && v[i] > -1){
                            blockCount++;
                            bfstrack[index] = i;
                            bfstrack[index + 1] = 1;
                            bfstrack[index + 2] = v[i];
                            currentLayerU[v[i]] = true;
                            index += 3;
                        }
                    }
                    else{
                        if (currentLayerU[i]){
                            long edgeCountIndex = index + 1;
                            long edgeCount = 0;
                            long tempIndex = index + 2;
                            bool connected = false;
                            for (long j = 0; j < m; j++){
                                if (edges[i * m + j] == 1 && u[i] != j && !visitedV[j]){
                                    connected = true;
                                    bfstrack[edgeCountIndex - 1] = i;
                                    if (vl[j] == -1){
                                        layerUFreeV[i] = true;
                                        bfsDone = true;
                                    }
                                    currentLayerV[j] = true;
                                    bfstrack[tempIndex] = j;
                                    edgeCount++;
                                    tempIndex++;
                                }
                            }
                            if (connected){
                                blockCount++;
                                index = tempIndex;
                            }
                            bfstrack[edgeCountIndex] = edgeCount;
                        }
                    }
                }

                for (long i = 0; i < m; i++){
                    if (k % 2 == 1){
                        currentLayerV[i] = false;
                    }
                    else{
                        currentLayerU[i] = false;
                        visitedV[i] = visitedV[i] || currentLayerV[i];
                    }
                }

                bfstrack[index] = previousBlockIndex;
                bfstrack[nextBlockIndex] = index;
                bfstrack[blockCountIndex] = blockCount;

                bfsDone = currentLayerU[0] == currentLayerV[0];

                for (long i = 1; i < m; i++){
                    bfsDone = bfsDone && (currentLayerU[i] == currentLayerV[i]);
                }

                if (!bfsDone){
                    previousBlockIndex = index;
                    nextBlockIndex = index + 1;
                    blockCountIndex = index + 2;
                    index+=3;
                    k++;
                }

                blockCount = 0;
            }

            // printf("DATA/\n");
            // for (long i = 0; i < 100 && i < previousBlockIndex; i++){
            //     printf("%d\n", bfstrack[i]);
            // }

            done = true;

            // bfs is done now
            // previousBlockIndex now contains where the last block starts
            long blockIndex = previousBlockIndex + 2;

            // note: below doesn't work as we can only write a path when we know that it ends

            // there'll never be more than m free vertices in V anyway
            // for (long i = 0; i < m && blockIndex < index; i++){
            //     // uIndex is the index of the vertex in U
            //     long uIndex = bfstrack[blockIndex];
            //     // how many connections to V
            //     long nrConn = bfstrack[blockIndex + 1];
            //     blockIndex += 2;
            //     if (!layerUFreeV[uIndex]){
            //         blockIndex += nrConn;
            //         continue;
            //     }
            //     for (long j = 0; j < nrConn; j++){
            //         if (vl[bfstrack[blockIndex + j]] == -1){
            //             long vIndex = bfstrack[blockIndex + j];
            //             // we found the non-connected vertex that we need to connect
            //             // now we need dfs for the path back.
            //             bool dfsDone = false;
            //             long oldVIndex = -1;

            //             if (ul[uIndex] == -1){
            //                 ul[uIndex] = vIndex;
            //                 vl[vIndex] = uIndex;
            //                 dfsDone = true;
            //             }
            //             else{
            //                 if (ul[uIndex] != u[uIndex]){
            //                     dfsDone = true;
            //                 }
            //                 else{
            //                     oldVIndex = ul[uIndex];
            //                     ul[uIndex] = vIndex;
            //                     vl[vIndex] = uIndex;
            //                 }
            //             }

            //             long tempPreviousBlockIndex = previousBlockIndex;

            //             while (!dfsDone){
            //                 vIndex = oldVIndex;
            //                 tempPreviousBlockIndex = bfstrack[bfstrack[tempPreviousBlockIndex]];
            //                 long tempNextBlockIndex = bfstrack[tempPreviousBlockIndex + 1];
            //                 long tempUIndex = bfstrack[tempPreviousBlockIndex + 2];
            //                 long tempNrConn = bfstrack[tempPreviousBlockIndex + 3];

            //                 long tempIndex = tempUIndex + 2;

            //                 while (tempIndex < tempNextBlockIndex){
            //                     if (tempIndex == tempUIndex + 2 + tempNrConn){
            //                         tempUIndex = bfstrack[tempIndex];
            //                         tempNrConn = bfstrack[tempIndex +1];
            //                         tempIndex += 2;
            //                     }

            //                     if (bfstrack[tempIndex] == vIndex){
            //                         if (ul[tempUIndex] == -1){
            //                             ul[tempUIndex] = vIndex;
            //                             vl[vIndex] = tempUIndex;
            //                             dfsDone = true;
            //                             break;
            //                         }
            //                         else{
            //                             if (ul[uIndex] != u[uIndex]){
            //                                 dfsDone = true;
            //                                 break;
            //                             }
            //                             else{
            //                                 oldVIndex = ul[uIndex];
            //                                 ul[uIndex] = vIndex;
            //                                 vl[vIndex] = uIndex;
            //                             }
            //                         }
            //                     }
            //                     tempIndex++;
            //                 }
            //             }
            //         }
            //     }
            //     blockIndex += nrConn;
            // }


        }
    }
    else{

    }


    // /*  A vector in which to store prime candidates (pc)  */
    // long *pcs = vecalloci(P);
    // bsp_push_reg(pcs, P * sizeof(long));
    // bsp_sync();
    
    // for (long i = 0; i < P; i++){
    //     pcs[i] = 2;
    // }

    // /*  This rounds down by default  */
    // long size = N / P;

    // /*  So we get the ceiling by adding 1 if the remainder is not zero  */
    // if (N % P != 0){
    //     size++;
    // }

    // /*  Local storage of all found prime numbers  */
    // bool *primes = vecallocb(size);

    // /*  Set all numbers except 1 as a prime number  */
    // for (long i = 0; i < size; i++){
    //     primes[i] = true;
    // }

    // if (s == 0){
    //     primes[0] = false;
    // }

    // /*  Our prime candidate, we know that the first one will be 2
    //     All processors will start with this candidate
    // */
    // long globalPc = 2;

    // /*  The actual values of the first and last number for each processor  */
    // long startValue = s * size + 1;
    // long endValue = (s + 1) * size;

    // /*  Correct endvalue for the last processor  */
    // if (endValue > N){
    //     endValue = N;
    // }

    // /*  For any pc it suffices to start checking from pc^2, since the
    //     previous rounds of the sieve will already have multiplied our pc
    //     with all values lower than pc. This also means that we can stop checking
    //     for multiples once we know that globalPc^2 > endValue
    // */

    // double endSetup = bsp_time(); 

    // while (globalPc * globalPc <= N){

    //     /*  The index from which we can start striking duplicates as non-prime
    //         note that we need to remove one from the index, since C is zero-based
    //     */

    //     long value = 0;

    //     /*  We need to find the smallest multiple that will be in range.
    //         This cuts down on the inital duplications until we hit the startvalue.
    //         If the found smallest multiple is smaller than globalPc^2, we use that.
    //     */
    //     long temp = startValue / globalPc;
    //     if (startValue % globalPc != 0){
    //         temp++;
    //     }
    //     value = temp * globalPc;

    //     if (value < globalPc * globalPc){
    //         value = globalPc * globalPc;
    //     }

    //     while (value <= endValue){
    //         /*  To ensure we don't check values below the range  */
    //         if (value < startValue){
    //             value += globalPc;
    //         }
    //         /*  And also not outside of the range  */
    //         else if (value <= endValue){
    //             /*  Get the actual index of the value we're checking  */
    //             long actualIndex = value - startValue;
    //             primes[actualIndex] = false;
    //             value += globalPc;
    //         }
    //     }
            
    //     long newPc;

    //     if (globalPc == pcs[s]){
    //         /*  For each processor find the smallest prime candidate  */
    //         for (long i = 0; i < size; i++){
    //             if (primes[i] && startValue + i > globalPc){
    //                 newPc = startValue + i;
    //                 break;
    //             }
    //         }

    //         /*  Each processor puts their pc in pcs on their index  */
    //         for (long i = 0; i < P; i++) {
    //             bsp_put(i, &newPc, pcs, s * sizeof(long), sizeof(long));
    //         }
    //     }
        
    //     bsp_sync();
        
    //     newPc = -1;
        
    //     /*  Check for the smallest pc that increases globalPc  */
    //     for (long i = 0; i < P; i++){
    //         if (pcs[i] > globalPc && (pcs[i] < newPc || newPc == -1)){
    //             newPc = pcs[i];
    //         }
    //     }

    //     globalPc = newPc;
    // }

    // double endTime = bsp_time();

    // printf("Setup took %.6lf seconds\n", endSetup - startTime);
    // printf("Calculating primes up to %ld on %ld processors took only %.6lf seconds.\n", N, P, endTime-startTime);
    
    // // for (long j = 0; j < size; j++){
    // //     /*  Check on endvalue for cases where P does not divide N  */
    // //     if (primes[j] && startValue + j <= endValue){
    // //         printf("Processor %ld found prime number %ld\n", s, startValue + j);
    // //     }
    // // }

    // bsp_pop_reg(pcs);
    // vecfreei(pcs);
    // vecfreeb(primes);

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
