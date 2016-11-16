/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
 */

/*
 *  Project:    libLinAlg
 *  Created on: April 20, 2016
 *  Author:     Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#include "../UTest.h"

void testCorrCoeff() {
    indexCoeff n = 5;
    real scaling = 4.;

    vector vecA;
    real rawA[] = {3, -1, 4, 12, -8};
    vectorCreateRaw(n, rawA, &vecA);

    vector vecB;
    real rawB[n];
    vectorCreateRaw(n, rawB, &vecB);

    vectorCopy(&vecA, &vecB);
    vecMulReal(scaling, &vecA);

    TEST_ASSERT_EQUAL_FLOAT(1., corrCoef(&vecA, &vecB));
}

void testDecorr() {
    indexCoeff n = 4;

    vector vNoise;
    double rNoise[] = {0., 0.03, 0.01, 0.05};
    vectorCreateRaw(n, rNoise, &vNoise);

    vector vNoiseInfluence;
    double rNoiseInfluence[n];
    vectorCreateRaw(n, rNoiseInfluence, &vNoiseInfluence);
    vectorCopy(&vNoise, &vNoiseInfluence);
    vecMulReal(4, &vNoiseInfluence);

    vector vSignal;
    double rSignal[] = {0., 1., 5., 6.};
    vectorCreateRaw(n, rSignal, &vSignal);
    vecAddVec(&vNoiseInfluence, &vSignal);

    decorr(&vNoise,
           &vSignal);

    real coeff = corrCoef(&vNoise, &vSignal);

    TEST_ASSERT_FLOAT_WITHIN(1e-8, 0, coeff);
}

void testDecorr3() {
    indexCoeff n = 4;

    vector vAccX;
    double rAccX[] = {0., 0.03, 0.01, 0.05};
    vectorCreateRaw(n, rAccX, &vAccX);

    vector vAccY;
    double rAccY [] = {-0.02, 0.01, 0.01, -0.012};
    vectorCreateRaw(n, rAccY, &vAccY);

    vector vAccZ;
    double rAccZ [] = {0.06, 0.07, 0.01, -0.09};
    vectorCreateRaw(n, rAccZ, &vAccZ);

    vector vGSR;
    double rGSR[] = {0., 1., 5., 6.};
    vectorCreateRaw(n, rGSR, &vGSR);

    decorr3(&vAccX,
            &vAccY,
            &vAccZ,
            &vGSR);

    real coeff;

//    coeff = corrCoef(&vAccX, &vGSR);
//    TEST_ASSERT_FLOAT_WITHIN(1e-2, 0, coeff);
//
//    coeff = corrCoef(&vAccY, &vGSR);
//    TEST_ASSERT_FLOAT_WITHIN(1e-2, 0, coeff);
//
//    coeff = corrCoef(&vAccZ, &vGSR);
//    TEST_ASSERT_FLOAT_WITHIN(1e-2, 0, coeff);
}

void testSmooth() {
    indexCoeff sizeIn = 20;
    indexCoeff smoothingSize = 4;
    real offset = 3.;

    indexCoeff sizeOut = sizeIn - smoothingSize + 1;

    vector vRaw;
    real rRaw[sizeIn];
    vectorCreateRaw(sizeIn, rRaw, &vRaw);

    vector vSmooth;
    real rSmooth[sizeOut];
    vectorCreateRaw(sizeOut, rSmooth, &vSmooth);

    vector vExpected;
    real rExpected[sizeOut];
    vectorCreateRaw(sizeOut, rExpected, &vExpected);
    vectorSetConstant(offset, &vExpected);

    real signal[] = {0., 1., 0., -1.};

    for (int i = 0; i < sizeIn; ++i) {
        vectorSetCoeff(&vRaw, i, signal[i % smoothingSize] + offset);
    }

    smooth(&vRaw, smoothingSize, &vSmooth);

    TEST_ASSERT_TRUE(vectorApproxEqual(&vSmooth, &vExpected, 1e-8));
}

void testSmoothID() {
    indexCoeff sizeIn = 20;
    indexCoeff smoothingSize = 1;
    real offset = 3.;

    indexCoeff sizeOut = sizeIn - smoothingSize + 1;

    vector vRaw;
    real rRaw[sizeIn];
    vectorCreateRaw(sizeIn, rRaw, &vRaw);

    vector vSmooth;
    real rSmooth[sizeOut];
    vectorCreateRaw(sizeOut, rSmooth, &vSmooth);

    real signal[] = {0., 1., 0., -1.};
    for (int i = 0; i < sizeIn; ++i) {
        vectorSetCoeff(&vRaw,
                       i,
                       signal[i % smoothingSize] + offset);
    }

    smooth(&vRaw, smoothingSize, &vSmooth);

    TEST_ASSERT_TRUE(vectorApproxEqual(&vRaw, &vSmooth, 1e-8));
}
