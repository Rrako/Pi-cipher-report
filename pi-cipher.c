#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

int w= 32;

// ARX Mu transformation constants
const unsigned long muConst[4] = {0xF0E8E4E2, 0xE1D8D4D2, 0xD1CCCAC9, 0xC6C5C3B8};

// ARX Nu transformation constants
const unsigned long nuConst[4] = {0xB4B2B1AC, 0xAAA9A6A5, 0xA39C9A99, 0x9695938E};

// Rotation vectors used in µ and ν
const unsigned int rou[4] = {5, 11, 17, 23};
const unsigned int rov[4] = {3, 10, 19, 29};

// Round constants for π32-Cipher
unsigned long int C[24] = {0x8D8B8778, 0x7472716C, 0x6A696665, 0x635C5A59,
                           0x5655534E, 0x4D4B473C, 0x3A393635, 0x332E2D2B,
                           0x271E1D1B, 0x170FF0E8, 0xE4E2E1D8, 0xD4D2D1CC,
                           0xCAC9C6C5, 0xC3B8B4B2, 0xB1ACAAA9, 0xA6A5A39C,
                           0x9A999695, 0x938E8D8B, 0x87787472, 0x716C6A69,
                           0x6665635C, 0x5A595655, 0x534E4D4B, 0x473C3A39};

unsigned long int ROTL32(unsigned long int value, unsigned int shift) {
    //obtaining the left shifted value
    unsigned long ls = value<<shift;
    //obtaining right shift value 
    unsigned long rs = value>>(32-shift);
    //obtaining XOR of right shift and left shift to get circular left shift
    return(ls^rs);
}

unsigned long longadd(unsigned long a, unsigned long b) {return (a+b)%((long long int)pow(2,w));}

void tuplecopy(int Xi, unsigned long int *X, int Yi, unsigned long int *Y) {
    for (int i=0; i<4; i++) {
        X[Xi + i] = Y[Yi + i];
    }
}

void ARX32(unsigned long int X[4], unsigned long int Y[4], unsigned long int Z[4]) {
    //temporary variables for transformation
    unsigned long T[12];
    
    //µ–transformation for X
        //addition and rotation
        T[0] = ROTL32(longadd(longadd(longadd(muConst[0], X[0]), X[1]), X[2]), rou[0]);
        T[1] = ROTL32(longadd(longadd(longadd(muConst[1], X[0]), X[1]), X[3]), rou[1]);
        T[2] = ROTL32(longadd(longadd(longadd(muConst[2], X[0]), X[2]), X[3]), rou[2]);
        T[3] = ROTL32(longadd(longadd(longadd(muConst[3], X[1]), X[2]), X[3]), rou[3]);

        //XOR
        T[4] = T[0]^T[1]^T[3];
        T[5] = T[0]^T[1]^T[2];
        T[6] = T[1]^T[2]^T[3];
        T[7] = T[0]^T[2]^T[3];
    
    //ν–transformation for Y
        //addition and rotation
        T[0] = ROTL32(longadd(longadd(longadd(nuConst[0], Y[0]), Y[2]), Y[3]), rov[0]);
        T[1] = ROTL32(longadd(longadd(longadd(nuConst[1], Y[1]), Y[2]), Y[3]), rov[1]);
        T[2] = ROTL32(longadd(longadd(longadd(nuConst[2], Y[0]), Y[1]), Y[2]), rov[2]);
        T[3] = ROTL32(longadd(longadd(longadd(nuConst[3], Y[0]), Y[1]), Y[3]), rov[3]);

        //XOR
        T[8] = T[1]^T[2]^T[3];
        T[9] = T[0]^T[2]^T[3];
        T[10] = T[0]^T[1]^T[3];
        T[11] = T[0]^T[1]^T[2];
    
    //σ–transformation for both µ(X) and ν(Y)
        Z[3] = longadd(T[4], T[8]);
        Z[0] = longadd(T[5], T[9]);
        Z[1] = longadd(T[6], T[10]);
        Z[2] = longadd(T[7], T[11]);
}

// Inputs: C(Tuple), J(Array), N(No. of I), I1,...,IN(Tuples)
void E1(unsigned long int C[4], unsigned long int *I, int N, unsigned long int *J) {
    //counter variable(s)
    int i;

    // output array of ARX32 sub-function
    unsigned long int Z[4];

    // the 4-tuples that are iterated through values of i
    unsigned long int Ii[4];
    unsigned long int Ji[4];
    
    // --- initial condition ---
    //fetching I1 from I
    tuplecopy(0, Ii, 0, I);

    //calculating C * I1 (store in Z)
    ARX32(C, Ii, Z);

    //store Z to J1
    tuplecopy(0, J, 0, Z);

    // --- looping for each value i = 2,...,N ---
    for (i=2; i<=N; i++) {

        //fetching Ji-1
        tuplecopy(0, Ji, (i-2)*4, J);

        //fetching Ii
        tuplecopy(0, Ii, (i-1)*4, I);

        //calculating Ji-1 * Ii
        ARX32(Ji, Ii, Z);

        //store Z to Ji
        tuplecopy((i-1)*4, J, 0, Z);
    }
}

// Inputs: C(Tuple), J(Array), N(No. of I), I1,...,IN(Tuples)
void E2(unsigned long int C[4], unsigned long int *I, int N, unsigned long int *J) {
    //counter variable(s)
    int i;

    // output array of ARX32 sub-function
    unsigned long int Z[4];

    // the 4-tuples that are iterated through values of i
    unsigned long int Ii[4];
    unsigned long int Ji[4];
    
    // --- initial condition ---
    //fetching IN from I
    tuplecopy(0, Ii, (N-1)*4, I);

    //calculating IN * C (store in Z)
    ARX32(Ii, C, Z);

    //store Z to JN
    tuplecopy((N-1)*4, J, 0, Z);

    // --- looping for each value i = 1,...,N-1 ---
    for (i=1; i<=N-1; i++) {

        //fetching IN-i
        tuplecopy(0, Ii, (N-1-i)*4, I);

        //fetching JN-i+1
        tuplecopy(0, Ji, (N-i)*4, J);

        //calculating IN-i * JN-i+1
        ARX32(Ii, Ji, Z);

        //store Z to JN-i
        tuplecopy((N-1-i)*4, J, 0, Z);
    }
}

// R is No. of rounds, C is constants, N is No. of I blocks, I is internal states, J is modified internal states
void pi(int R, unsigned long int *C, int N, unsigned long int *I, unsigned long int *J) {
    //R defines the number of rounds of pi, and therefore the number of C tuples (C=2*R)

    //counter variable(s)
    int i;

    //constant tuples
    unsigned long int Cx[4], Cy[4];

    // exception case
    if (R<1 || R>3) {
        printf("error: R value exception");
        return;
    } else {
        // for each round of pi
        for (i=0; i<R; i++) {
            
            //fetching constant tuples from C list
            tuplecopy(0, Cx, (i*8),  C);
            tuplecopy(0, Cy, (i*8)+4, C);
            
            //execute one round of pi
            E1(Cx, I, N, J);
            E2(Cy, J, N, I);
        }
    }
}
