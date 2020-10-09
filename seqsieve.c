#include "bspedupack.h"

/*  This program computes the primes up to a given number
    using a prime number sieve.
    The distribution of the sieve is by block
*/

long N; // number up until which to calculate prime numbers
long P; // number of processors requested

void seqsieve(){
    /*  Sequentially sieve for prime numbers  */
    bsp_begin(1);

    /*  Boolean array of length N  */
    bool *primes = vecallocb(N);

    /*  Set all numbers except 1 as a prime number  */
    for (long i = 0; i < N; i++){
        primes[i] = true;
    }

    primes[0] = false;
    
    /*  Our prime candidate (pc), we know that the first one will be 2  */
    long pc = 2;

    /*  For any pc it suffices to start checking from pc^2, since the
        previous rounds of the sieve will already have multiplied our pc
        with all values lower than pc. This also means that we can stop checking
        for multiples once we know that pc^2 > N
    */
    while (pc * pc <= N){
        /*  The index from which we can start striking duplicates as non-prime
            note that we need to remove one from the index, since C is zero-based
        */
        long index = (pc * pc) - 1;

        /*  Here we also use N as the length of our array  */
        while (index < N){
            primes[index] = false;
            index += pc;
        }

        /*  Starting from our current pc, the next number marked
            as prime will be our next pc
        */
        for (long i = pc; i < N; i++){
            if (primes[i]){
                pc = i + 1;
                break;
            }
        }
    }

    for (long i = 0; i < N; i++){
        if (primes[i]){
            printf("%ld\n", i + 1);
        }
    }

    vecfreeb(primes);
    bsp_end();
} /* end seqsieve */

int main(int argc, char **argv){

    /* Sequential part for N*/
    printf("Up until which number do you want to find prime numbers?\n");
    fflush(stdout);

    scanf("%ld",&N);
    bsp_init(seqsieve, argc, argv);
    seqsieve();

    /* Sequential part */
    exit(EXIT_SUCCESS);

} /* end main */
