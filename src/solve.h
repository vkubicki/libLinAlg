/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
*/

/*
 *  Project:    libLinAlg
 *  Created on: April 20, 2016
 *  Authors:    Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#ifndef SOLVE_H
#define SOLVE_H

#include "matrix.h"
#include "vector.h"

/** Note that only the lower triangular values are set to speed up computation. The behaviour on the remainder is undefined. */
void matrixCholesky(const matrix* A,
                    matrix* B);

/** L is a lower triangular matrix, as output by matrixCholesky for example. The function solves L * x = y
 * by forward substitution. */
void forwardSubstitutionSolve(const matrix* L,
                              const vector* y,
                              vector* x);

/** L is a lower triangular matrix, as output by matrixCholesky for example. The function solves L^t * x = y
 * by forward substitution. Note that the transpose of L is multiplied by x, so that the result of matrixCholesky
 * could be used directly. */
void backwardSubstitutionSolve(const matrix* L,
                               const vector* y,
                               vector* x);

void solveCholesky(const matrix* A,
                   const vector* y,
                   vector* x);

void leastSquare(const matrix* A,
                 const vector* y,
                 vector* x);

#endif
