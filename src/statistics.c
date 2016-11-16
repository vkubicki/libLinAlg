/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
*/

/*
 *  Project:    libLinAlg
 *  Created on: March 9, 2016
 *  Authors:    Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#include "linAlg.h"
#include "statistics.h"

real vecMean(const vector* vecIn) {
    real mean = 0.;
    indexCoeff nI = vectorGetDim(vecIn);

    for (indexCoeff i = 0; i < nI; ++i) {
        mean += vectorGetCoeff(vecIn, i);
    }

    mean /= nI;
    return mean;
}

real corrCoef(const vector* vecInA,
              const vector* vecInB) {
    real coeff;
    indexCoeff nI = vectorGetDim(vecInA);

    vector centeredA;
    real rawCenteredA[nI];
    vectorCreateRaw(nI, rawCenteredA, &centeredA);
    vectorCopy(vecInA, &centeredA);
    vecAddReal(-vecMean(vecInA), &centeredA);

    vector centeredB;
    real rawCenteredB[nI];
    vectorCreateRaw(nI, rawCenteredB, &centeredB);
    vectorCopy(vecInB, &centeredB);
    vecAddReal(-vecMean(vecInB), &centeredB);

    coeff = vecDotVec(&centeredA, &centeredB);
    coeff /= vecNorm(&centeredA) * vecNorm(&centeredB);

    return coeff;
}

void decorr(const vector* noise,
            vector* signal) {
    indexCoeff nI = vectorGetDim(signal);

    vector centeredNoise;
    real rawCenteredNoise[nI];
    vectorCreateRaw(nI, rawCenteredNoise, &centeredNoise);
    vectorCopy(noise, &centeredNoise);

    vecAddReal(-vecMean(noise), &centeredNoise);
    real scaling = vecDotVec(&centeredNoise, signal) / vecSquaredNorm(&centeredNoise);
    vecMulReal(-scaling, &centeredNoise);
    vecAddVec(&centeredNoise, signal);
}

void decorr3(const vector* noise0,
             const vector* noise1,
             const vector* noise2,
             vector* signal) {
    decorr(noise0, signal);
    decorr(noise1, signal);
    decorr(noise2, signal);
}

void smooth(const vector* signalRaw,
            int windowSize,
            vector* signalSmooth) {
    indexCoeff sizeIn = vectorGetDim(signalRaw);
    indexCoeff sizeOut = vectorGetDim(signalSmooth);
    real wReal = (real) windowSize;
    for (int i = 0; i < sizeOut; ++i) {
        real average = 0.;
        for (int iP = i; iP < i + windowSize; ++iP) {
            average += vectorGetCoeff(signalRaw, iP);
        }
        average /= wReal;
        vectorSetCoeff(signalSmooth,
                       i,
                       average);
    }
}

void smoothFromRaw(const data* dataIn,
                   indexCoeff inSize,
                   indexCoeff smoothSize,
                   vector* smoothSignal) {
    vector vDataIn;
    real rDataIn[inSize];
    vectorCreateRaw(inSize, rDataIn, &vDataIn);
    vectorSetDataOffset(dataIn, 0, &vDataIn);

    smooth(&vDataIn,
           inSize,
           smoothSignal);
}
