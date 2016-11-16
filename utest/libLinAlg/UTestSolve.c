/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
 */

/*
 *  Project:    libLinAlg
 *  Created on: April 20, 2016
 *  Author:     Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#include "../UTest.h"

void testMatrixCholesky() {
    matrix* A = matrixCreateMalloc(3, 3);
    matrix* computedRes = matrixCreateMalloc(3, 3);
    matrix* expectedRes = matrixCreateMalloc(3, 3);

    real dataA[] = {4., 12., -16.,
                    12., 37., -43.,
                    -16., -43., 98.};

    real dataExpectedRes[] = {2., 0., 0.,
                              6., 1., 0.,
                              -8., 5., 3.};

    matrixSetReal(dataA, A);
    matrixSetReal(dataExpectedRes, expectedRes);

    matrixCholesky(A, computedRes);

    TEST_ASSERT_EQUAL_INT(1, matrixLowerTrigApproxEqual(expectedRes, computedRes, epsilon));

    matrixFree(A);
    matrixFree(computedRes);
    matrixFree(expectedRes);
}

void testMatrixCholesky2() {
    indexCoeff n = 3;

    matrix AGram;
    real AGramData[] = {21.   , 0.  , 0.  ,
                        13030., 8085328., 0.,
                        13016., 8076651., 8067996.};
    matrixCreateRaw(n, n, AGramData, &AGram);

    matrix LComputed;
    real LComputedData[n * n];
    matrixCreateRaw(n, n, LComputedData, &LComputed);

    matrix LComputedTranspose;
    real LComputedTransposeData[n * n];
    matrixCreateRaw(n, n, LComputedTransposeData, &LComputedTranspose);

    matrix LLt;
    real LLtData[n * n];
    matrixCreateRaw(n, n, LLtData, &LLt);

    matrixCholesky(&AGram, &LComputed);
    matrixTranspose(&LComputed, &LComputedTranspose);
    matrixProduct(&LComputed, &LComputedTranspose, &LLt);

    TEST_ASSERT_EQUAL_INT(1, matrixLowerTrigApproxEqual(&AGram, &LLt, epsilon));
}

void testForwardSubstitutionSolve() {
    matrix* L = matrixCreateMalloc(3, 3);
    vector* x = vectorCreateMalloc(3);
    vector* computedY = vectorCreateMalloc(3);
    vector* expectedY = vectorCreateMalloc(3);

    real dataL[] = {4., 0., 0.,
                    12., 37., 0.,
                    -16., -43., 98.};

    real dataY[] = {2., -4., 27.};

    matrixSetReal(dataL, L);
    vectorSetReal(dataY, expectedY);

    forwardSubstitutionSolve(L, expectedY, x);
    matrixVectorProduct(L, x, computedY);

    TEST_ASSERT_EQUAL_INT(1,vectorApproxEqual(expectedY, computedY, epsilon));

    matrixFree(L);
    vectorFree(x);
    vectorFree(computedY);
    vectorFree(expectedY);
}

void testBackwardSubstitutionSolve() {
    matrix* L = matrixCreateMalloc(3, 3);
    vector* x = vectorCreateMalloc(3);
    vector* computedY = vectorCreateMalloc(3);
    vector* expectedY = vectorCreateMalloc(3);

    real dataL[] = {4., 0., 0.,
                    12., 37., 0.,
                    -16., -43., 98.};

    real dataY[] = {2., -4., 27.};

    matrixSetReal(dataL, L);
    vectorSetReal(dataY, expectedY);

    backwardSubstitutionSolve(L, expectedY, x);
    matrixTransposeVectorProduct(L, x, computedY);

    TEST_ASSERT_EQUAL_INT(1,vectorApproxEqual(expectedY, computedY, epsilon));

    matrixFree(L);
    vectorFree(x);
    vectorFree(computedY);
    vectorFree(expectedY);
}

void testSolveCholesky() {
    matrix* A = matrixCreateMalloc(3, 3);
    matrix* LA = matrixCreateMalloc(3, 3);
    vector* xExpected = vectorCreateMalloc(3);
    vector* xComputed = vectorCreateMalloc(3);
    vector* y = vectorCreateMalloc(3);

    real dataA[] = {4., 12., -16.,
                    12., 37., -43.,
                    -16., -43., 98.};
    matrixSetReal(dataA, A);

    real dataLA[] = {4., 0., 0.,
                    12., 37., 0.,
                    -16., -43., 98.};
    matrixSetReal(dataLA, LA);

    real dataXExpected[] = {1., 12., -45.};
    vectorSetReal(dataXExpected, xExpected);

    matrixVectorProduct(A, xExpected, y);

    solveCholesky(LA, y, xComputed); // lower triangular component of A is used in the solver, to prove that only this part is needed, and not the whole matrix A

    TEST_ASSERT_EQUAL_INT(1,vectorApproxEqual(xExpected, xComputed, epsilon));

    matrixFree(A);
    matrixFree(LA);
    vectorFree(xExpected);
    vectorFree(xComputed);
    vectorFree(y);
}

void testMatrixGramLower() {
    int nObservation = 3;
    int nCoeff = 2;

    matrix* A = matrixCreateMalloc(nObservation, nCoeff);
    matrix* lowerGramExpected = matrixCreateMalloc(nCoeff, nCoeff);
    matrix* lowerGramComputed = matrixCreateMalloc(nCoeff, nCoeff);

    real dataA[] = {1., 2.,
                    3., 4.,
                    5., 6.};
    matrixSetReal(dataA, A);

    real dataLowerGramExpected[] = {35., 44.,
                                    44., 56};
    matrixSetReal(dataLowerGramExpected, lowerGramExpected);

    matrixGramLower(A, lowerGramComputed);

    TEST_ASSERT_EQUAL_INT(1, matrixLowerTrigApproxEqual(lowerGramExpected, lowerGramComputed, epsilon));

    matrixFree(A);
    matrixFree(lowerGramExpected);
    matrixFree(lowerGramComputed);
}

void testLeastSquare() {
    int nObservation = 40;
    int nCoeff = 2;

    matrix* A = matrixCreateMalloc(nObservation, nCoeff);
    vector* xComputed = vectorCreateMalloc(nCoeff);
    vector* xExpected = vectorCreateMalloc(nCoeff);
    vector* y = vectorCreateMalloc(nObservation);

    for (int c = 0; c < nCoeff; ++c) {
        vectorSetCoeff(xExpected, c, pow(c + 2, 2));
    }

    for (int i = 0; i < nObservation; ++i) {
        matrixSetCoeff(A, i, 0, i     );
        matrixSetCoeff(A, i, 1, i + 10);
        vectorSetCoeff(y, i,    i       * vectorGetCoeff(xExpected, 0)
                             + (i + 10) * vectorGetCoeff(xExpected, 1));
    }

    leastSquare(A, y, xComputed);

    matrixFree(A);
    vectorFree(xComputed);
    vectorFree(xExpected);
    vectorFree(y);
}

void testLeastSquare2() {
    indexCoeff nObs = 21;
    indexCoeff modelOrder = 3;

    matrix PhiLearn;
    real PhiLearnData[] = {1., 613., 612.,
                           1., 613., 613.,
                           1., 613., 613.,
                           1., 614., 613.,
                           1., 614., 614.,
                           1., 615., 614.,
                           1., 616., 615.,
                           1., 618., 616.,
                           1., 620., 618.,
                           1., 622., 620.,
                           1., 622., 622.,
                           1., 624., 622.,
                           1., 624., 624.,
                           1., 624., 624.,
                           1., 625., 624.,
                           1., 625., 625.,
                           1., 625., 625.,
                           1., 625., 625.,
                           1., 626., 625.,
                           1., 626., 626.,
                           1., 626., 626.};
    matrixCreateRaw(nObs, modelOrder, PhiLearnData, &PhiLearn);

    vector rhs;
    real rhsData[] = {613.,
                      613.,
                      614.,
                      614.,
                      615.,
                      616.,
                      618.,
                      620.,
                      622.,
                      622.,
                      624.,
                      624.,
                      624.,
                      625.,
                      625.,
                      625.,
                      625.,
                      626.,
                      626.,
                      626.,
                      625.};
    vectorCreateRaw(nObs, rhsData, &rhs);

    vector thetaComputed;
    real thetaComputedData[modelOrder];
    vectorCreateRaw(modelOrder, thetaComputedData, &thetaComputed);

    vector thetaExpected;
    real thetaExpectedData[] = {32.22329,
                                 1.14652,
                                -0.19775}; // computed using Matlab
    vectorCreateRaw(modelOrder, thetaExpectedData, &thetaExpected);

    leastSquare(&PhiLearn, &rhs, &thetaComputed);

    TEST_ASSERT_EQUAL_INT(1, vectorApproxEqual(&thetaExpected, &thetaComputed, 1e-3));
}
