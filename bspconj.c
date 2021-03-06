#include "bspedupack.h"

/*  This program computes the primes up to a given number
    using a prime number sieve.
    The distribution of the sieve is by block
*/

long N; // number up until which to calculate prime numbers
long P; // number of processors requested

void bspconj(){
    
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
        All processors will start with this candidate
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

    bsp_pop_reg(pcs);
    vecfreei(pcs);

    bool *allPrimes = vecallocb(P * size);
    bsp_push_reg(allPrimes, P * size * sizeof(bool));
    bsp_sync();

    for (long i = 0; i < P; i++) {
        bsp_put(i, primes, allPrimes, s * size * sizeof(bool), size * sizeof(bool));
    }

    bsp_sync();
    vecfreeb(primes);

    long isSumSize = size / 2;

    if (size % 2 != 0){
        isSumSize++;
    }
    
    long *isSumValues = vecalloci(isSumSize);

    for (long i = 0; i < isSumSize; i++){
        isSumValues[i] = -1;
    }

    long isSumStartValue = startValue;

    if (startValue % 2 != 0){
        isSumStartValue++;
    }

    long isSumEndValue = endValue;

    if (endValue % 2 != 0){
        isSumEndValue--;
    }

    long isSumValue = isSumStartValue;

    for (isSumValue; isSumValue <= isSumEndValue; isSumValue += 2){
        /*  What we make use of here is that if an even number can be written
            as n = p + q with p,q prime, then given that n = 2m we know that m
            must be equidistant of p and q
        */

        long halfSumValue = (isSumValue / 2);
        long halfSumIndex = halfSumValue - 1;

        /*  To stay equidistant we simply subtract and add the same value to the halfvalue  */
        for (long i = 0; i < halfSumValue; i++){
            if (allPrimes[halfSumIndex - i] && allPrimes[halfSumIndex + i]){
                isSumValues[(isSumValue - isSumStartValue) / 2] = halfSumValue - i;
                break;
            }
        }
    }

    double endTime = bsp_time();

    printf("Setup took %.6lf seconds\n", endSetup - startTime);
    printf("From setup to calculating all primes took %.6lf seconds\n", endPrimes - endSetup);
    printf("Verifying conjecture up to %ld on %ld processors took only %.6lf seconds.\n", N, P, endTime-startTime);
    

    // for (long i = 0; i < isSumSize; i++){
    //     if (isSumValues[i] == -1 && isSumStartValue <= isSumEndValue){
    //         if (isSumStartValue + (i * 2) > 2){
    //             printf("Processor %ld did not find a twinpair for value %ld\n", s, isSumStartValue + (i * 2));
    //         }
    //     }
    //     else if(isSumStartValue + (i * 2) > 2){
    //         printf("Processor %ld found primesum %ld + %ld for even number %ld\n", s, isSumValues[i], isSumStartValue + (i * 2) - isSumValues[i], isSumStartValue + (i * 2));
    //     }
    // }

    bsp_pop_reg(allPrimes);
    vecfreeb(allPrimes);
    vecfreei(isSumValues);

    bsp_end();

} /* end bspconj */

int main(int argc, char **argv){

    bsp_init(bspconj, argc, argv);

    /* Sequential part for N*/
    printf("Up until which number do you want to verify Goldbach's conjecture?\n");
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
    bspconj();

    /* Sequential part */
    exit(EXIT_SUCCESS);

} /* end main */
