/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
 */

/*
 *  Project:    libLinAlg
 *  Created on: April 20, 2016
 *  Author:     Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#include "../UTest.h"
#include "math.h"

void testVectorCreateMalloc() {
    indexCoeff n = 5;
    vector* vec = vectorCreateMalloc(n);

    for (indexCoeff i = 0; i < n; ++i) {
        vectorSetCoeff(vec, i, i);
    }

    int correctVector = 1;
    for (int i = 0; i < n; ++i) {
        if (vectorGetCoeff(vec, i) != i) {
            correctVector = 0;
        }
    }

    TEST_ASSERT_EQUAL_INT(correctVector, 1);

    vectorFree(vec);
}

void testVectorCreateRaw() {
    indexCoeff n = 5;

    vector vec;
    real raw[n];
    vectorCreateRaw(n, raw, &vec);

    for (indexCoeff i = 0; i < n; ++i) {
        vectorSetCoeff(&vec, i, i);
    }

    int correctVector = 1;
    for (int i = 0; i < n; ++i) {
        if (vectorGetCoeff(&vec, i) != i) {
            correctVector = 0;
        }
    }

    TEST_ASSERT_EQUAL_INT(correctVector, 1);
}

void testVecMean() {
    indexCoeff n = 5;

    vector vec;
    real raw[] = {3, -1, 4, 12, -8};
    vectorCreateRaw(n, raw, &vec);

    TEST_ASSERT_EQUAL_FLOAT(vecMean(&vec), 2);
}

void testVecAddReal() {
    indexCoeff n = 5;
    real offset = 12;

    vector vec;
    real raw[] = {3, -1, 4, 12, -8};
    vectorCreateRaw(n, raw, &vec);

    real meanBefore = vecMean(&vec);
    vecAddReal(offset, &vec);
    real meanAfter = vecMean(&vec);

    TEST_ASSERT_EQUAL_FLOAT(meanBefore + offset, meanAfter);
}

void testVecMulReal() {
    indexCoeff n = 5;
    real offset = 12;

    vector vec;
    real raw[] = {3, -1, 4, 12, -8};
    vectorCreateRaw(n, raw, &vec);

    real meanBefore = vecMean(&vec);
    vecMulReal(offset, &vec);
    real meanAfter = vecMean(&vec);

    TEST_ASSERT_EQUAL_FLOAT(meanBefore * offset, meanAfter);
}

void testVecDotVec() {
    indexCoeff n = 2;

    vector vec0;
    real raw0[] = {1, 1};
    vectorCreateRaw(n, raw0, &vec0);

    vector vec1;
    real raw1[] = {2, 2};
    vectorCreateRaw(n, raw1, &vec1);

    vector vec2;
    real raw2[] = {-1, 1};
    vectorCreateRaw(n, raw2, &vec2);

    TEST_ASSERT_EQUAL_FLOAT(4, vecDotVec(&vec0, &vec1));
    TEST_ASSERT_EQUAL_FLOAT(0, vecDotVec(&vec0, &vec2));
}

void testVectorCopy() {
    indexCoeff n = 5;

    vector vecA;
    real rawA[] = {3, -1, 4, 12, -8};
    vectorCreateRaw(n, rawA, &vecA);

    vector vecB;
    real rawB[n];
    vectorCreateRaw(n, rawB, &vecB);

    vectorCopy(&vecA, &vecB);

    TEST_ASSERT_EQUAL_FLOAT(vecMean(&vecA), vecMean(&vecB));
}

void testVectorNegShiftA() {
    indexCoeff n = 5;
    indexCoeff shift = -2;

    vector aComputed;
    real aComputedRaw[] = {3, -1, 4, 12, -8};
    vectorCreateRaw(n, aComputedRaw, &aComputed);

    vector aExpected;
    real aExpectedRaw[] = {4, 12, -8, 12, -8};
    vectorCreateRaw(n, aExpectedRaw, &aExpected);

    vectorNegShift(shift, &aComputed);

    TEST_ASSERT_TRUE(vectorApproxEqual(&aExpected, &aComputed, epsilon));
}

void testVectorNegShiftB() {
    indexCoeff n = 5;
    indexCoeff shift = -5;

    vector aComputed;
    real aComputedRaw[] = {3, -1, 4, 12, -8};
    vectorCreateRaw(n, aComputedRaw, &aComputed);

    vector aExpected;
    real aExpectedRaw[] = {3, -1, 4, 12, -8};
    vectorCreateRaw(n, aExpectedRaw, &aExpected);

    vectorNegShift(shift, &aComputed);

    TEST_ASSERT_TRUE(vectorApproxEqual(&aExpected, &aComputed, epsilon));
}

void testVectorTailCopyA() {
    indexCoeff n = 5;
    indexCoeff sizeOut = 3;

    vector vecA;
    real rawA[] = {3, -1, 4, 12, -8};
    vectorCreateRaw(n, rawA, &vecA);

    vector vecB;
    real rawB[sizeOut];
    vectorCreateRaw(sizeOut, rawB, &vecB);

    vectorTailCopy(&vecA,
                   sizeOut,
                   &vecB);

    TEST_ASSERT_EQUAL_FLOAT((4. + 12. - 8.) / 3., vecMean(&vecB));
}

void testVectorTailCopyB() {
    indexCoeff nSrc = 5;
    indexCoeff nDst = 3;

    vector src;
    real srcData[] = {3, -1, 4, 12, -8};
    vectorCreateRaw(nSrc, srcData, &src);

    vector destComputed;
    real destComputedData[nDst];
    vectorCreateRaw(nDst, destComputedData, &destComputed);

    vector destExpected;
    real destExpectedData[] = {4, 12, -8};
    vectorCreateRaw(nDst, destExpectedData, &destExpected);

    vectorTailCopy(&src, nSrc, &destComputed);

    TEST_ASSERT_EQUAL_INT(true, vectorApproxEqual(&destExpected, &destComputed, epsilon));
}

void testVectorTailCopyCycleB() {
    indexCoeff nSrc = 5;
    indexCoeff nDst = 3;

    vector src;
    real srcData[] = {3., -1., 4., 12., -8.};
    vectorCreateRaw(nSrc, srcData, &src);

    vector destComputed;
    real destComputedData[nDst];
    vectorCreateRaw(nDst, destComputedData, &destComputed);

    vector destExpected;
    real destExpectedData[] = {-8., 3., -1.};
    vectorCreateRaw(nDst, destExpectedData, &destExpected);

    vectorTailCopyCycle(&src, 1, nSrc, &destComputed);

    TEST_ASSERT_EQUAL_INT(true, vectorApproxEqual(&destExpected, &destComputed, epsilon));
}

void testVectorPartialCopy() {
    indexCoeff n = 5;
    indexCoeff sizeOut = 3;

    vector vecA;
    real rawA[] = {3, -1, 4, 12, -8};
    vectorCreateRaw(n, rawA, &vecA);

    vector vecB;
    real rawB[sizeOut];
    vectorCreateRaw(sizeOut, rawB, &vecB);

    vectorPartialCopy(&vecA,
                      2,
                      4,
                      0,
                      &vecB);

    TEST_ASSERT_EQUAL_FLOAT((4. + 12. - 8.) / 3., vecMean(&vecB));
}

void testVectorSetData() {
    indexCoeff sizeVec = 8;
    data val = 12;

    vector vComputed;
    real rComputed[sizeVec];
    vectorCreateRaw(sizeVec, rComputed, &vComputed);

    vector vExpected;
    real rExpected[sizeVec];
    vectorCreateRaw(sizeVec, rExpected, &vExpected);
    vectorSetConstant(12, &vExpected);

    data raw[sizeVec];
    for (int i = 0; i < sizeVec; ++i) {
        raw[i] = val;
    }

    vectorSetDataOffset(raw,
                  0,
                  &vComputed);

    TEST_ASSERT_TRUE(vectorApproxEqual(&vExpected, &vComputed, 1e-8));
}

void testVectorApproxEqual() {
    indexCoeff n = 4;

    vector vAccX;
    double rAccX[] = {0., 0.03, 0.01, 0.05};
    vectorCreateRaw(n, rAccX, &vAccX);

    vector vAccY;
    double rAccY [n];
    vectorCreateRaw(n, rAccY, &vAccY);
    vectorCopy(&vAccX, &vAccY);
    vecAddReal(1e-8, &vAccY);

    vector vAccZ;
    double rAccZ [n];
    vectorCreateRaw(n, rAccZ, &vAccZ);
    vectorCopy(&vAccX, &vAccZ);
    vecAddReal(1e-3, &vAccZ);

    TEST_ASSERT_TRUE(vectorApproxEqual(&vAccX, &vAccY, 1e-4));
    TEST_ASSERT_FALSE(vectorApproxEqual(&vAccX, &vAccZ, 1e-4));
}

void testVectorAddLeft() {
    indexCoeff indexInitialSize = 3;
    indexCoeff indexFinalSize = 5;

    vector* VComputed = vectorCreateMalloc(indexInitialSize);
    double VComputedData[] = {4., 87., -12.};
    vectorSetReal(VComputedData, VComputed);

    double VExpectedData[] = {-1, -1, 4, 87, -12};

    vectorAddLeft(indexFinalSize - indexInitialSize, VComputed);

    bool isCorrect = true;
    for (int i = indexFinalSize - indexInitialSize; i < indexFinalSize; ++i) {
        if (fabs(vectorGetCoeff(VComputed, i) - VExpectedData[i]) > epsilon) {
            isCorrect = false;
        }
    }

    TEST_ASSERT_EQUAL_INT(true, isCorrect);

    vectorFree(VComputed);
}

void testVectorResizeLarger() {
    indexCoeff initialSize = 3;
    indexCoeff finalSize = 5;

    vector* VComputed = vectorCreateMalloc(initialSize);
    double VComputedData[] = {4., 87., -12.};
    vectorSetReal(VComputedData, VComputed);

    double VExpectedData[] = {-1, -1, 4, 87, -12};

    vectorResize(finalSize, VComputed);

    bool isCorrect = true;
    for (int i = finalSize - initialSize; i < finalSize; ++i) {
        if (fabs(vectorGetCoeff(VComputed, i) - VExpectedData[i]) > epsilon) {
            isCorrect = false;
        }
    }

    TEST_ASSERT_EQUAL_INT(true, isCorrect);

    vectorFree(VComputed);
}

void testVectorResizeSameSize() {
    indexCoeff initialSize = 3;
    indexCoeff finalSize = 3;

    vector* VComputed = vectorCreateMalloc(initialSize);
    double VComputedData[] = {4., 87., -12.};
    vectorSetReal(VComputedData, VComputed);

    double VExpectedData[] = {4, 87, -12};

    vectorResize(finalSize, VComputed);

    bool isCorrect = true;
    for (int i = finalSize - initialSize; i < finalSize; ++i) {
        if (fabs(vectorGetCoeff(VComputed, i) - VExpectedData[i]) > epsilon) {
            isCorrect = false;
        }
    }

    TEST_ASSERT_EQUAL_INT(true, isCorrect);

    vectorFree(VComputed);
}

void testVectorResizeSmaller() {
    indexCoeff initialSize = 5;
    indexCoeff finalSize = 3;

    vector* VComputed = vectorCreateMalloc(initialSize);
    double VComputedData[] = {4., 87., -12., 45., 16.};
    vectorSetReal(VComputedData, VComputed);

    double VExpectedData[] = {-12., 45., 16.};

    vectorResize(finalSize, VComputed);

    bool isCorrect = true;
    for (int i = finalSize - initialSize; i < finalSize; ++i) {
        if (fabs(vectorGetCoeff(VComputed, i) - VExpectedData[i]) > epsilon) {
            isCorrect = false;
        }
    }

    TEST_ASSERT_EQUAL_INT(true, isCorrect);

    vectorFree(VComputed);
}

void testVectorMinMax() {
    indexCoeff dim = 7;

    vector vec;
    real vecData[] = {-12., -5., 99., 134., 50., -75., 5.};
    vectorCreateRaw(7, vecData, &vec);

    TEST_ASSERT_EQUAL_FLOAT(-75., vectorMin(&vec));
    TEST_ASSERT_EQUAL_FLOAT(134., vectorMax(&vec));
}
