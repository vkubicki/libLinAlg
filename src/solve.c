/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
*/

/*
 *  Project:    libLinAlg
 *  Created on: April 20, 2016
 *  Authors:    Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#include <math.h>
#include <stdio.h>

#include "solve.h"

void matrixCholesky(const matrix* A,
                    matrix* B) {
    indexCoeff nI, nJ;
    matrixGetDim(A, &nI, &nJ);
    for (int i = 0; i < nI; ++i) {
        for (int j = 0; j <= i; ++j) {
            real s = 0.;
            for (int k = 0; k < j; k++) {
                s += matrixGetCoeff(B, i, k) * matrixGetCoeff(B, j, k);
            }
            matrixSetCoeff(B, i, j, (i == j) ?
                                    sqrt(matrixGetCoeff(A, i, i) - s) :
                                    (1.0 / matrixGetCoeff(B, j, j) * (matrixGetCoeff(A, i, j) - s)));
        }
    }
}

void forwardSubstitutionSolve(const matrix* L,
                              const vector* y,
                              vector* x) {
    indexCoeff nI = vectorGetDim(y);
    real sumPrevious;
    for (int i = 0; i < nI; ++i) {
        sumPrevious = 0.;
        for (int c = 0; c < i; ++c) { // part based on previously computed coefficients of x
            sumPrevious += matrixGetCoeff(L, i, c) * vectorGetCoeff(x, c);
        }
        vectorSetCoeff(x, i, (vectorGetCoeff(y, i) - sumPrevious) / matrixGetCoeff(L, i, i));
    }
}

void backwardSubstitutionSolve(const matrix* L,
                               const vector* y,
                               vector* x) {
    indexCoeff nI = vectorGetDim(y);
    real sumPrevious;
    for (int i = nI - 1; i > -1; --i) {
        sumPrevious = 0.;
        for (int c = i + 1; c < nI; ++c) { // part based on previously computed coefficients of x
            sumPrevious += matrixGetCoeff(L, c, i) * vectorGetCoeff(x, c);
        }
        vectorSetCoeff(x, i, (vectorGetCoeff(y, i) - sumPrevious) / matrixGetCoeff(L, i, i));
    }
}

void solveCholesky(const matrix* A,
                   const vector* y,
                   vector* x) {
    indexCoeff nI, nJ;
    matrixGetDim(A, &nI, &nJ);

    matrix L; // Cholesky decomposition of A, lower triangular matrix
    real LRaw[nI * nI];
    matrixCreateRaw(nI, nI, LRaw, &L);
    matrixCholesky(A, &L);

    vector z; // stores z so that L * z = rhs, by forward substitution
    real zRaw[nI];
    vectorCreateRaw(nJ, zRaw, &z);

    forwardSubstitutionSolve(&L, y, &z);
    backwardSubstitutionSolve(&L, &z, x); // backward substitution to compute x, which solves the least square problem
}

void leastSquare(const matrix* A,
                 const vector* y,
                 vector* x) {
    indexCoeff nI, nJ;
    matrixGetDim(A, &nI, &nJ);

    matrix AGram; // lower triangular component of the Gram matrix of A
    real AGramRaw[nJ * nJ];
    matrixCreateRaw(nJ, nJ, AGramRaw, &AGram);
    matrixGramLower(A, &AGram);

    matrix L; // Cholesky decomposition of AGram, lower triangular matrix
    real LRaw[nJ * nJ];
    matrixCreateRaw(nJ, nJ, LRaw, &L);
    matrixCholesky(&AGram, &L);

    vector rhs; // stores right hand side vector, A^t * y
    real rhsRaw[nJ];
    vectorCreateRaw(nJ, rhsRaw, &rhs);
    matrixTransposeVectorProduct(A, y, &rhs);

    vector z;
    real zRaw[nJ];
    vectorCreateRaw(nJ, zRaw, &z);

    forwardSubstitutionSolve(&L, &rhs, &z); // stores z so that L * z = rhs, by forward substitution
    backwardSubstitutionSolve(&L, &z, x); // backward substitution to compute x, which solves the least square problem
#ifdef LINDEBUG_OFF
    printf("AGram:\n");
    matrixPrint(&AGram);
    printf("L:\n");
    matrixPrint(&L);
    printf("rhs:\n");
    vectorPrint(&rhs);
    printf("z:\n");
    vectorPrint(&z);
    printf("x:\n");
    vectorPrint(x);
#endif
}
