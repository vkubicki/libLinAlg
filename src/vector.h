/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
*/

/*
 *  Project:    libLinAlg
 *  Created on: March 8, 2016
 *  Authors:    Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#ifndef VECTOR_H
#define VECTOR_H

#include "typedef.h"

typedef struct {
    indexCoeff nI_;
    real* data_;
} vector;

vector* vectorCreateMalloc(indexCoeff size);

void vectorCreateRaw(indexCoeff size,
                     real* raw,
                     vector* vec);

void vectorSetCoeff(vector* vec,
                    indexCoeff n,
                    real val);

void vectorSetConstant(real val,
                       vector* vec);

void vectorSetReal(const real* dataIn,
                   vector* vec);

void vectorSetData(const data* dataIn,
                   vector* vec);

void vectorSetDataOffset(const data* dataIn,
                         integer offset,
                         vector* vec);

real vectorGetCoeff(const vector* vec,
                    indexCoeff n);

void vectorSetCoeffCycle(vector* vec,
                         indexCoeff n,
                         real val);

real vectorGetCoeffCycle(const vector* vec,
                         indexCoeff n);

void vectorCopy(const vector* vecB,
                vector* vecA);

void vectorNegShift(integer shift,
                    vector* vec);

void vectorTailCopy(const vector* vecB,
                    indexCoeff sizeOut,
                    vector* vecA);

void vectorTailCopyCycle(const vector* vecB,
                         indexCoeff lastIndex,
                         indexCoeff sizeOut,
                         vector* vecA);

void vectorPartialCopy(const vector* vecIn,
                       indexCoeff firstIndIn,
                       indexCoeff lastIndIn,
                       indexCoeff firstIndOut,
                       vector* vecOut);

void vectorPartialCopyToArray(const vector* vecIn,
                              indexCoeff firstIndIn,
                              indexCoeff lastIndIn,
                              indexCoeff firstIndOut,
                              data* dataOut);

void vectorFree(vector* vecOut);

indexCoeff vectorGetDim(const vector* vec);

void vectorPrint(const vector* vec);

int vectorApproxEqual(const vector* vec1,
                      const vector* vec2,
                      real epsilon);

/** Add n Elements to the left of an malloc'ed vector */
void vectorAddLeft (indexCoeff nElement,
                    vector* vec);

/** Resize an malloc'ed vector, and preserve the most right coefficients possible */
void vectorResize (indexCoeff nElement,
                   vector* vec);

void vectorResizeNoCopy(indexCoeff newSize,
                        vector* vec);

real vectorMin(const vector* vec);

real vectorMax(const vector* vec);

#endif
