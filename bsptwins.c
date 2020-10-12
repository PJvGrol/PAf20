#include "bspedupack.h"

/*  This program computes the primes up to a given number
    using a prime number sieve.
    The distribution of the sieve is by block
*/

long N; // number up until which to calculate prime numbers
long P; // number of processors requested

void bsptwins(){
    
    bsp_begin(P);

    double startTime = bsp_time();

    long p= bsp_nprocs(); // p = number of processors
    long s= bsp_pid();    // s = processor number

    /*  A vector in which to store prime candidates (pc)  */
    long *pcs = vecalloci(P);
    bsp_push_reg(pcs, P * sizeof(long));
    bsp_sync();
    
    for (long i = 0; i < P; i++){
        pcs[i] = 2;
    }

    /*  This rounds down by default  */
    long size = N / P;

    /*  So we get the ceiling by adding 1 if the remainder is not zero  */
    if (N % P != 0){
        size++;
    }

    /*  Local storage of all found prime numbers  */
    bool *primes = vecallocb(size);

    /*  Set all numbers except 1 as a prime number  */
    for (long i = 0; i < size; i++){
        primes[i] = true;
    }

    if (s == 0){
        primes[0] = false;
    }

    /*  Our prime candidate, we know that the first one will be 2
        All processors will start with tbspsievehis candidate
    */
    long globalPc = 2;

    /*  The actual values of the first and last number for each processor  */
    long startValue = s * size + 1;
    long endValue = (s + 1) * size;

    /*  Correct endvalue for the last processor  */
    if (endValue > N){
        endValue = N;
    }

    /*  For any pc it suffices to start checking from pc^2, since the
        previous rounds of the sieve will already have multiplied our pc
        with all values lower than pc. This also means that we can stop checking
        for multiples once we know that pc^2 > endValue
    */

    double endSetup = bsp_time(); 

    while (globalPc * globalPc <= N){

        /*  The index from which we can start striking duplicates as non-prime
            note that we need to remove one from the index, since C is zero-based
        */

        long value = 0;

        /*  We need to find the smallest multiple that will be in range.
            This cuts down on the inital duplications until we hit the startvalue.
            If the found smallest multiple is smaller than globalPc^2, we use that.
        */
        long temp = startValue / globalPc;
        if (startValue % globalPc != 0){
            temp++;
        }
        value = temp * globalPc;

        if (value < globalPc * globalPc){
            value = globalPc * globalPc;
        }

        while (value <= endValue){
            /*  To ensure we don't check values below the range  */
            if (value < startValue){
                value += globalPc;
            }
            /*  And also not outside of the range  */
            else if (value <= endValue){
                /*  Get the actual index of the value we're checking  */
                long actualIndex = value - startValue;
                primes[actualIndex] = false;
                value += globalPc;
            }
        }

        long newPc;

        /*  For each processor find the smallest prime candidate  */
        for (long i = 0; i < size; i++){
            if (primes[i] && startValue + i > globalPc){
                newPc = startValue + i;
                break;
            }
        }

        /*  Each processor puts their pc in pcs on their index  */
        for (long i = 0; i < P; i++) {
            bsp_put(i, &newPc, pcs, s * sizeof(long), sizeof(long));
        }
        
        bsp_sync();
        
        newPc = -1;
        
        /*  Check for the smallest pc that increases globalPc  */
        for (long i = 0; i < P; i++){
            if (pcs[i] > globalPc && (pcs[i] < newPc || newPc == -1)){
                newPc = pcs[i];
            }
        }

        globalPc = newPc;
    }

    double endPrimes = bsp_time();

    // for (long j = 0; j < size; j++){
    //     /*  Check on endvalue for cases where P does not divide N  */
    //     if (primes[j] && startValue + j <= endValue){
    //         printf("Processor %ld found prime number %ld\n", s, startValue + j);
    //     }
    // }

    bsp_pop_reg(pcs);
    vecfreei(pcs);

    /*  The lowest prime on the current processor  */
    long lowestPrime;

    for (long i = 0; i < size; i++){
        if (primes[i]){
            lowestPrime = startValue + i;
            break;
        }
    }

    /*  The lowest prime on the next processor
        This uses that we have a block distribution
        where the next processor will have the next values
    */
    long nextPrime;

    bsp_push_reg(&nextPrime, sizeof(long));
    bsp_sync();

    if (s > 0){
        bsp_put(s - 1, &lowestPrime, &nextPrime, 0, sizeof(long));
    }

    bsp_sync();
    bool *twins = vecallocb(size);

    for (long i = 0; i < size; i++){
        twins[i] = false;
    }

    for (long i = 0; i < size; i++){
        /*  We can check for twin primes using that they are (k * 6) +-1
            The only pair that is not found with it is (3,5), so an exception for that
        */
        if ((startValue + i < 5 || (startValue + i + 1) % 6 == 0) && primes[i]){
            if ((i + 2 <= endValue && primes[i + 2]) || startValue + i + 2 == nextPrime){
                /*  We track only the lower of the two primes that form a twin
                    This allows us to use nextPrime from the nextblock to check for
                    twins in cases where the lower twin is endValue (-1)
                */
                twins[i] = true;
            }
        }
    }

    double endTime = bsp_time();

    printf("Setup took %.6lf seconds\n", endSetup - startTime);
    printf("From setup to calculating all primes took %.6lf seconds\n", endPrimes - endSetup);
    printf("Calculating twin primes up to %ld on %ld processors took only %.6lf seconds.\n", N, P, endTime-startTime);
    
    
    for (long i = 0; i < size; i++){
        if (twins[i]){
            printf("Processor %ld found twin primes (%ld, %ld)\n", s, startValue + i, startValue + i + 2);
        }
    }

    bsp_pop_reg(nextPrime);
    vecfreeb(primes);
    vecfreeb(twins);

    bsp_end();

} /* end bsptwins */

int main(int argc, char **argv){

    /* Sequential part for N*/
    printf("Up until which number do you want to find twin primes?\n");
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
    bsp_init(bsptwins, argc, argv);
    bsptwins();

    /* Sequential part */
    exit(EXIT_SUCCESS);

} /* end main */
