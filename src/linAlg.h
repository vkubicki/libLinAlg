/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
*/

/*
 *  Project:    libLinAlg
 *  Created on: March 9, 2016
 *  Authors:    Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#ifndef LINALG_H
#define LINALG_H

#include "vector.h"

/** Add real b to vector a by modifying vector a in place */
void vecAddReal(real realB,
                vector* vecA);

void vecAddVec(const vector* vecIn,
               vector* vecOut);

void vecMulReal(real realB,
                vector* vecA);

real vecDotVec(const vector* vecA,
               const vector* vecB);

real vecNorm(const vector* vecA);

real vecSquaredNorm(const vector* vecA);

#endif
