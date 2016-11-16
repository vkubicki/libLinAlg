/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
*/

/*
 *  Project:    libLinAlg
 *  Created on: April 19, 2016
 *  Authors:    Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#ifndef MATRIX_H
#define MATRIX_H

#include "typedef.h"
#include "vector.h"

/** Storage in matrix is row-major ordered */
typedef struct {
    indexCoeff nI;
    indexCoeff nJ;
    real* data;
} matrix;

matrix* matrixCreateMalloc(indexCoeff nI,
                           indexCoeff nJ);

void matrixFree(matrix* mat);

void matrixCreateRaw(indexCoeff nI,
                     indexCoeff nJ,
                     real* raw,
                     matrix* mat);

void matrixSetCoeff(matrix* mat,
                    indexCoeff i,
                    indexCoeff j,
                    real val);

void matrixSetConstant(real val,
                       matrix* mat);

void matrixColumnSetConstant(real val,
                             indexCoeff col,
                             matrix* mat);

void matrixSetReal(const real* data,
                   matrix* mat);

real matrixGetCoeff(const matrix* mat,
                    indexCoeff i,
                    indexCoeff j);

void matrixGetDim(const matrix* mat,
                  indexCoeff* nI,
                  indexCoeff* nJ);

/** C = A * B */
void matrixProduct(const matrix* A,
                   const matrix* B,
                   matrix* C);

void matrixVectorProduct(const matrix* A,
                         const vector* x,
                         vector* y);

void matrixTransposeVectorProduct(const matrix* A,
                                  const vector* x,
                                  vector* y);

void matrixPrint(const matrix* mat);

int matrixApproxEqual(const matrix* A,
                      const matrix* B,
                      real epsilon);

int matrixLowerTrigApproxEqual(const matrix* A,
                               const matrix* B,
                               real epsilon);

void matrixTranspose(const matrix* A,
                     matrix* B);

/** Compute the lower part of the Gram matrix.*/
void matrixGramLower(const matrix* A,
                     matrix* B);

void matrixCopyVectorToColumn(const vector* vec,
                              indexCoeff firstInd,
                              indexCoeff lastInd,
                              indexCoeff col,
                              matrix* mat);

void matrixCopyVectorToColumnCyclic(const vector* vec,
                                    indexCoeff firstInd,
                                    indexCoeff lastInd,
                                    indexCoeff col,
                                    matrix* mat);

#endif
