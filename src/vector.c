/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
*/

/*
 *  Project:    libLinAlg
 *  Created on: March 8, 2016
 *  Authors:    Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "base.h"
#include "vector.h"

vector* vectorCreateMalloc(indexCoeff size) {
    vector* vecOut;

    real* data = malloc(size * sizeof(real));
    vecOut = malloc(sizeof(vector));
    vecOut->nI_ = size;
    vecOut->data_ = data;

    return vecOut;
}

void vectorCreateRaw(indexCoeff size,
                     real* raw,
                     vector* vec) {
    vec->nI_ = size;
    vec->data_ = raw;
}

void vectorSetCoeff(vector* vec,
                    indexCoeff n,
                    real val) {
    vec->data_[n] = val;
}

void vectorSetConstant(real val,
                       vector* vec) {
    indexCoeff nI = vectorGetDim(vec);
    for (indexCoeff i = 0; i < nI; ++i) {
        vectorSetCoeff(vec,
                       i,
                       val);
    }
}

void vectorSetReal(const real* dataIn,
                   vector* vec) {
    indexCoeff nI = vectorGetDim(vec);
    for (indexCoeff i = 0; i < nI; ++i) {
        vectorSetCoeff(vec, i, dataIn[i]);
    }
}

void vectorSetData(const data* dataIn,
                   vector* vec) {
    indexCoeff nI = vectorGetDim(vec);
    for (indexCoeff i = 0; i < nI; ++i) {
        vectorSetCoeff(vec, i, dataIn[i]);
    }
}

void vectorSetDataOffset(const data* dataIn,
                         integer offset,
                         vector* vec) {
    indexCoeff nI = vectorGetDim(vec);
    for (indexCoeff i = 0; i < nI; ++i) {
        vectorSetCoeff(vec,
                       i,
                       dataIn[i + offset]);
    }
}

real vectorGetCoeff(const vector* vec,
                    indexCoeff n) {
    return vec->data_[n];
}

void vectorSetCoeffCycle(vector* vec,
                         indexCoeff n,
                         real val) {
    indexCoeff nI = vectorGetDim(vec);
    indexCoeff localN = n % nI;
#ifdef LINDEBUG_OFF
    printf("vectorSetCoeffCycle localN: %zu\n", localN);
#endif
    vec->data_[localN] = val;
}

real vectorGetCoeffCycle(const vector* vec,
                         indexCoeff n) {
    indexCoeff nI = vectorGetDim(vec);
    indexCoeff localN = n % nI;
#ifdef LINDEBUG_OFF
    printf("vectorGetCoeffCycle, nI: %zu, n: %d, localN: %zu\n", nI, n, localN);
#endif
    return vec->data_[localN];
}

void vectorCopy(const vector* vecB,
                vector* vecA) {
    indexCoeff nI = vectorGetDim(vecA);
    for (indexCoeff i = 0; i < nI; ++i) {
        vectorSetCoeff(vecA,
                       i,
                       vectorGetCoeff(vecB,
                                      i));
    }
}

void vectorNegShift(integer shift,
                    vector* vec) { // shift is provided as a negative value
    indexCoeff dim = vectorGetDim(vec);
    if (-shift < dim) {
        memmove(vec->data_,
                vec->data_ - shift,
                (vectorGetDim(vec) + shift) * sizeof(real));
    }
}

void vectorTailCopy(const vector* vecB,
                    indexCoeff sizeOut,
                    vector* vecA) {
    indexCoeff nA = vectorGetDim(vecA);
    indexCoeff nB = vectorGetDim(vecB);
    indexCoeff minN = min(nA, nB);
    minN = min(sizeOut, minN);

    indexCoeff fromB = nB - minN;
    indexCoeff toA = nA - minN;

    memcpy(vecA->data_ + toA,
           vecB->data_ + fromB,
           minN * sizeof(real));
}

void vectorTailCopyCycle(const vector* vecB,
                         indexCoeff lastIndex,
                         indexCoeff sizeOut,
                         vector* vecA) {
    indexCoeff nA = vectorGetDim(vecA);
    indexCoeff nB = vectorGetDim(vecB);
    indexCoeff minN = min(nA, nB);
    minN = min(sizeOut, minN);

    indexCoeff fromB = lastIndex - minN + 1;
    indexCoeff toA = nA - minN;

    indexCoeff localB;

    for (indexCoeff i = 0; i < minN; ++i) {
        localB = mod(fromB + i, nB); // mod is necessary instead of % to handle negative indices correctly
        vectorSetCoeff(vecA,
                       toA + i,
                       vectorGetCoeff(vecB,
                                      localB));
    }
}

void vectorPartialCopy(const vector* vecIn,
                       indexCoeff firstIndIn,
                       indexCoeff lastIndIn,
                       indexCoeff firstIndOut,
                       vector* vecOut) {
    indexCoeff nCopy = lastIndIn - firstIndIn + 1;
    for (indexCoeff i = 0;
         i < nCopy;
         ++i) {
        vectorSetCoeff(vecOut,
                       firstIndOut + i,
                       vectorGetCoeff(vecIn,
                                      firstIndIn + i));
    }
}

void vectorPartialCopyToArray(const vector* vecIn,
                              indexCoeff firstIndIn,
                              indexCoeff lastIndIn,
                              indexCoeff firstIndOut,
                              data* dataOut) {
    indexCoeff nCopy = lastIndIn - firstIndIn + 1;
    for (indexCoeff i = 0;
         i < nCopy;
         ++i) {
        dataOut[firstIndOut + i] = vectorGetCoeff(vecIn,
                                                  firstIndIn + i);
    }
}

void vectorFree(vector* vecOut) {
    free(vecOut->data_);
    free(vecOut);
}

indexCoeff vectorGetDim(const vector* vec) {
    return vec->nI_;
}

void vectorPrint(const vector* vec) {
    printf("%f", vec->data_[0]);
    for (indexCoeff i = 1; i < vec->nI_; ++i) {
        printf(" %f", vectorGetCoeff(vec,
                                     i));
    }
    printf("\n");
}

int vectorApproxEqual(const vector* vec1,
                      const vector* vec2,
                      real epsilon) {
    indexCoeff approxEqual = 1;
    indexCoeff nI = vectorGetDim(vec1);
    for (indexCoeff i = 0; i < nI; ++i) {
        if (fabs(vectorGetCoeff(vec1,
                                i) -
                 vectorGetCoeff(vec2,
                                i)) > epsilon) {
            return 0;
        }
    }
    return approxEqual;
}

void vectorAddLeft(indexCoeff nElement,
                   vector* vec) {
    indexCoeff currSize = vectorGetDim(vec);
    indexCoeff newSize = currSize + nElement;
    real* newData = malloc(newSize * sizeof(real));
    memmove(newData + nElement, vec->data_, currSize * sizeof(real)); // no memory overlap, memmove will be faster than memcopy
    free(vec->data_);
    vec->nI_ = newSize;
    vec->data_ = newData;
}

void vectorResize(indexCoeff newSize,
                  vector* vec) {
    indexCoeff currSize = vectorGetDim(vec);
    indexCoeff delta = newSize - currSize;
    real* newData = malloc(newSize * sizeof(real));

    if (newSize > currSize) {
        memmove(newData + delta, vec->data_        , currSize * sizeof(real)); // no memory overlap, memmove will be faster than memcopy
        free(vec->data_);
        vec->nI_ = newSize;
        vec->data_ = newData;
    }
    else if (newSize < currSize) {
        memmove(newData        , vec->data_ - delta, (currSize + delta) * sizeof(real)); // no memory overlap, memmove will be faster than memcopy
        free(vec->data_);
        vec->nI_ = newSize;
        vec->data_ = newData;
    }
}

void vectorResizeNoCopy(indexCoeff newSize,
                        vector* vec) {
    free(vec->data_);
    vec->nI_ = newSize;
    vec->data_ = malloc(newSize * sizeof(real));
}

real vectorMin(const vector* vec) {
    indexCoeff dim = vectorGetDim(vec);
    real min = vectorGetCoeff(vec, 0);

    for (indexCoeff c = 1; c < dim; ++c) {
        real val = vectorGetCoeff(vec, c);
        if (val < min) {
            min = val;
        }
    }

    return min;
}

real vectorMax(const vector* vec) {
    indexCoeff dim = vectorGetDim(vec);
    real max = vectorGetCoeff(vec, 0);

    for (indexCoeff c = 1; c < dim; ++c) {
        real val = vectorGetCoeff(vec, c);
        if (max < val) {
            max = val;
        }
    }

    return max;
}
