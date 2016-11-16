/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
*/

/*
 *  Project:    libLinAlg
 *  Created on: March 9, 2016
 *  Authors:    Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#ifndef STATISTICS_H
#define STATISTICS_H

#include "vector.h"

real vecMean(const vector* vecIn);

real corrCoef(const vector* vecInA,
              const vector* vecInB);

void decorr(const vector* noise,
            vector* signal);

void decorr3(const vector* noise0,
             const vector* noise1,
             const vector* noise2,
             vector* signal);

/** Forward smoothing by averaging. */
void smooth(const vector* signalRaw,
            int windowSize,
            vector* signalSmooth);

void smoothFromRaw(const data* dataIn,
                   indexCoeff inSize,
                   indexCoeff smoothSize,
                   vector* smoothSignal);

#endif
