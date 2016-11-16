/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
*/

/*
 *  Project:    libLinAlg
 *  Created on: April 19, 2016
 *  Authors:    Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

matrix* matrixCreateMalloc(indexCoeff nI,
                           indexCoeff nJ) {
    matrix* matOut;

    real* data = malloc(nI * nJ * sizeof(real));
    matOut = malloc(sizeof(matrix));
    matOut->nI = nI;
    matOut->nJ = nJ;
    matOut->data = data;

    return matOut;
}

void matrixFree(matrix* mat) {
    free(mat->data);
    free(mat);
}

void matrixCreateRaw(indexCoeff nI,
                     indexCoeff nJ,
                     real* raw,
                     matrix* mat) {
    mat->nI = nI;
    mat->nJ = nJ;
    mat->data = raw;
}

void matrixSetCoeff(matrix* mat,
                    indexCoeff i,
                    indexCoeff j,
                    real val) {
    mat->data[i * mat->nJ + j] = val;
}

void matrixSetConstant(real val,
                       matrix* mat) {
    indexCoeff nI, nJ;
    matrixGetDim(mat, &nI, &nJ);
    for (indexCoeff i = 0; i < nI; ++i) {
        for (indexCoeff j = 0; j < nJ; ++j) {
            matrixSetCoeff(mat, i, j, val);
        }
    }
}

void matrixColumnSetConstant(real val,
                             indexCoeff col,
                             matrix* mat) {
    indexCoeff nI, nJ;
    matrixGetDim(mat, &nI, &nJ);
    for (indexCoeff i = 0; i < nI; ++i) {
        matrixSetCoeff(mat, i, col, val);
    }
}

void matrixSetReal(const real* data,
                   matrix* mat) {
    indexCoeff nI, nJ;
    matrixGetDim(mat,
                 &nI,
                 &nJ);
    for (indexCoeff i = 0; i < nI; ++i) {
        for (indexCoeff j = 0; j < nJ; ++j) {
            matrixSetCoeff(mat,
                           i,
                           j,
                           data[i * mat->nJ + j]);
        }
    }
}

real matrixGetCoeff(const matrix* mat,
                    indexCoeff i,
                    indexCoeff j) {
    return mat->data[i * mat->nJ + j];
}

void matrixGetDim(const matrix* mat,
                  indexCoeff* nI,
                  indexCoeff* nJ) {
    *nI = mat->nI;
    *nJ = mat->nJ;
}

void matrixProduct(const matrix* A,
                   const matrix* B,
                   matrix* C) {
    indexCoeff Ai, Aj, Bi, Bj;
    matrixGetDim(A, &Ai, &Aj);
    matrixGetDim(B, &Bi, &Bj);

    real sum;
    for (indexCoeff i = 0; i < Ai; ++i) {
        for (indexCoeff j = 0; j < Bj; ++j) {
            sum = 0.;
            for (indexCoeff c = 0; c < Aj; ++c) {
                sum += matrixGetCoeff(A,
                                      i,
                                      c) * matrixGetCoeff(B,
                                                          c,
                                                          j);
            }
            matrixSetCoeff(C,
                           i,
                           j,
                           sum);
        }
    }
}

void matrixVectorProduct(const matrix* A,
                         const vector* x,
                         vector* y) {
    indexCoeff nI, nJ;
    matrixGetDim(A, &nI, &nJ);
    real sum;

    for (indexCoeff i = 0; i < nI; ++i) {
        sum = 0.;
        for (indexCoeff c = 0; c < nJ; ++c) {
            sum += matrixGetCoeff(A, i, c) * vectorGetCoeff(x, c);
        }
        vectorSetCoeff(y, i, sum);
    }
}

void matrixTransposeVectorProduct(const matrix* A,
                                  const vector* x,
                                  vector* y) {
    indexCoeff nI, nJ;
    matrixGetDim(A, &nI, &nJ);
    real sum;

    for (indexCoeff j = 0; j < nJ; ++j) {
        sum = 0.;
        for (indexCoeff c = 0; c < nI; ++c) {
            sum += matrixGetCoeff(A, c, j) * vectorGetCoeff(x, c);
        }
        vectorSetCoeff(y, j, sum);
    }
}

void matrixPrint(const matrix* mat) {
    indexCoeff nI, nJ;
    matrixGetDim(mat,
                 &nI,
                 &nJ);

    for (indexCoeff i = 0; i < nI; ++i) {
        printf("%f", matrixGetCoeff(mat, i, 0));
        for (indexCoeff j = 1; j < nJ; ++j) {
            printf(" %f", matrixGetCoeff(mat, i, j));
        }
        printf("\n");
    }
}

int matrixApproxEqual(const matrix* A,
                      const matrix* B,
                      real epsilon) {
    int approxEqual = 1;
    indexCoeff Ai, Aj, Bi, Bj;
    matrixGetDim(A, &Ai, &Aj);
    matrixGetDim(B, &Bi, &Bj);

    if (Ai != Bi || Aj != Bj) {
        return 0;
    }

    for (indexCoeff i = 0; i < Ai; ++i) {
        for (indexCoeff j = 0; j < Aj; ++j) {
            if (fabs(matrixGetCoeff(A, i, j) -
                     matrixGetCoeff(B, i, j))
                > epsilon) {
                return 0;
            }
        }
    }
    return approxEqual;
}

int matrixLowerTrigApproxEqual(const matrix* A,
                               const matrix* B,
                               real epsilon) {
    int approxEqual = 1;
    indexCoeff Ai, Aj, Bi, Bj;
    matrixGetDim(A, &Ai, &Aj);
    matrixGetDim(B, &Bi, &Bj);

    if (Ai != Bi || Aj != Bj) {
        return 0;
    }

    for (indexCoeff i = 0; i < Ai; ++i) {
        for (indexCoeff j = 0; j <= i; ++j) {
            if (fabs(matrixGetCoeff(A, i, j) -
                     matrixGetCoeff(B, i, j))
                > epsilon) {
                return 0;
            }
        }
    }
    return approxEqual;
}

void matrixTranspose(const matrix* A,
                     matrix* B) {
    indexCoeff Ai, Aj;
    matrixGetDim(A, &Ai, &Aj);

    for (indexCoeff i = 0; i < Ai; ++i) {
        for (indexCoeff j = 0; j < Aj; ++j) {
            matrixSetCoeff(B, j, i, matrixGetCoeff(A, i, j));
        }
    }
}

void matrixGramLower(const matrix* A,
                     matrix* B) {
    indexCoeff Ai, Aj;
    matrixGetDim(A, &Ai, &Aj);
    real sum;

    for (indexCoeff i = 0; i < Aj; ++i) { // Gram matrix is an Aj x Aj matrix
        for (indexCoeff j = 0; j <= i; ++j) {
            sum = 0.;
            for (int c = 0; c < Ai; ++c) {
                sum += matrixGetCoeff(A, c, i) * matrixGetCoeff(A, c, j);
            }
            matrixSetCoeff(B, i, j, sum);
        }
    }
}

void matrixCopyVectorToColumn(const vector* vec,
                              indexCoeff firstInd,
                              indexCoeff lastInd,
                              indexCoeff col,
                              matrix* mat) {
    indexCoeff nInd = lastInd - firstInd + 1;
    for (indexCoeff i = 0; i < nInd; ++i) {
        matrixSetCoeff(mat,
                       i,
                       col,
                       vectorGetCoeff(vec, firstInd + i)); // column data is contiguous but not row data, thus preventing a convenient use of memcpy
    }
}

void matrixCopyVectorToColumnCyclic(const vector* vec,
                                    indexCoeff firstInd,
                                    indexCoeff lastInd,
                                    indexCoeff col,
                                    matrix* mat) {
    indexCoeff nInd = lastInd - firstInd + 1;
    for (indexCoeff i = 0; i < nInd; ++i) {
        matrixSetCoeff(mat,
                       i,
                       col,
                       vectorGetCoeffCycle(vec, firstInd + i)); // The only difference with the non cycle implementation is the vectorGetCoeffCycle value access
    }
}
