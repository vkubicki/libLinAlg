/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
*/

/*
 *  Project:    libLinAlg
 *  Created on: March 9, 2016
 *  Authors:    Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#include <math.h>

#include "linAlg.h"

void vecAddReal(real realB,
                vector* vecA) {
    indexCoeff nI = vectorGetDim(vecA);
    for (indexCoeff i = 0; i < nI; ++i) {
        vectorSetCoeff(vecA,
                       i,
                       vectorGetCoeff(vecA, i) + realB);
    }
}

void vecMulReal(real realB,
                vector* vecA) {
    indexCoeff nI = vectorGetDim(vecA);
    for (indexCoeff i = 0; i < nI; ++i) {
        vectorSetCoeff(vecA,
                       i,
                       vectorGetCoeff(vecA, i) * realB);
    }
}

void vecAddVec(const vector* vecIn,
               vector* vecOut) {
    indexCoeff nI = vectorGetDim(vecIn);
    for (indexCoeff i = 0; i < nI; ++i) {
        vectorSetCoeff(vecOut,
                       i,
                       vectorGetCoeff(vecOut, i) + vectorGetCoeff(vecIn, i));
    }
}

real vecDotVec(const vector* vecA,
               const vector* vecB) {
    indexCoeff nI = vectorGetDim(vecA);
    real dot = 0.;
    for (indexCoeff i = 0; i < nI; ++i) {
        dot += vectorGetCoeff(vecA, i) * vectorGetCoeff(vecB, i);
    }
    return dot;
}

real vecNorm(const vector* vecA) {
    return sqrt(vecDotVec(vecA, vecA));
}

real vecSquaredNorm(const vector* vecA) {
    return vecDotVec(vecA, vecA);
}
