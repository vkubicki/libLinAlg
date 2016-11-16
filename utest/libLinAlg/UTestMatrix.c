/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
 */

/*
 *  Project:    libLinAlg
 *  Created on: April 20, 2016
 *  Author:     Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#include "../UTest.h"

void testMatrixProduct() {
    matrix* A = matrixCreateMalloc(3, 2);
    matrix* B = matrixCreateMalloc(2, 1);
    matrix* computedC = matrixCreateMalloc(3, 1);
    matrix* expectedC = matrixCreateMalloc(3, 1);

    real dataA[] = {1., 2.,
                    3., 4.,
                    5., 6.};

    real dataB[] = {7.,
                    8.};

    real dataExpectedC[] = {23.,
                            53.,
                            83.};

    matrixSetReal(dataA, A);
    matrixSetReal(dataB, B);
    matrixSetReal(dataExpectedC, expectedC);

    matrixProduct(A, B, computedC);

    TEST_ASSERT_EQUAL_INT(1, matrixApproxEqual(expectedC, computedC, epsilon));

    matrixFree(A);
    matrixFree(B);
    matrixFree(computedC);
    matrixFree(expectedC);
}

void testMatrixVectorProduct() {
    matrix* A = matrixCreateMalloc(2, 3);
    vector* x = vectorCreateMalloc(3);
    vector* yExpected = vectorCreateMalloc(2);
    vector* yComputed = vectorCreateMalloc(2);

    real AData[] = {1., 2., 3.,
                    4., 5., 6.};
    matrixSetReal(AData, A);

    real xData[] = {7.,
                    8.,
                    9.};
    vectorSetReal(xData, x);

    real yExpectedData[] = {50., 122.};
    vectorSetReal(yExpectedData, yExpected);

    matrixVectorProduct(A, x, yComputed);

    TEST_ASSERT_EQUAL_INT(1, vectorApproxEqual(yExpected, yComputed, epsilon));

    matrixFree(A);
    vectorFree(x);
    vectorFree(yExpected);
    vectorFree(yComputed);
}

void testMatrixTransposeVectorProduct() {
    matrix* A = matrixCreateMalloc(3, 2);
    vector* x = vectorCreateMalloc(3);
    vector* yExpected = vectorCreateMalloc(2);
    vector* yComputed = vectorCreateMalloc(2);

    real AData[] = {1., 4.,
                    2., 5.,
                    3., 6.};
    matrixSetReal(AData, A);

    real xData[] = {7.,
                    8.,
                    9.};
    vectorSetReal(xData, x);

    real yExpectedData[] = {50., 122.};
    vectorSetReal(yExpectedData, yExpected);

    matrixTransposeVectorProduct(A, x, yComputed);

    TEST_ASSERT_EQUAL_INT(1, vectorApproxEqual(yExpected, yComputed, epsilon));

    matrixFree(A);
    vectorFree(x);
    vectorFree(yExpected);
    vectorFree(yComputed);
}

void testCopyVectorToColumn() {
    indexCoeff rowSize = 3;
    indexCoeff colSize = 2;
    indexCoeff vecSize = 9;
    indexCoeff firstInd = 3;
    indexCoeff lastInd = 5;

    matrix* AComputed = matrixCreateMalloc(rowSize, colSize);
    real AComputedData[] = {1., 4.,
                            2., 5.,
                            3., 6.};
    matrixSetReal(AComputedData, AComputed);

    vector* x = vectorCreateMalloc(vecSize);
    real xData[] = {24.,
                    9.,
                    13.,
                    -7.,
                    22.,
                    -9.,
                    -87.,
                    7.,
                    9.};
    vectorSetReal(xData, x);

    matrix* AExpected = matrixCreateMalloc(rowSize, colSize);
    real AExpectedData[] = {-7., 4.,
                            22., 5.,
                            -9., 6.};
    matrixSetReal(AExpectedData, AExpected);

    matrixCopyVectorToColumn(x, firstInd, lastInd, 0, AComputed);

    TEST_ASSERT_EQUAL_INT(1, matrixApproxEqual(AExpected, AComputed, epsilon));

    matrixFree(AComputed);
    matrixFree(AExpected);
    vectorFree(x);
}

void testCopyVectorToColumnCyclic() {
    indexCoeff rowSize = 3;
    indexCoeff colSize = 2;
    indexCoeff vecSize = 9;

    matrix* AComputed = matrixCreateMalloc(rowSize, colSize);
    real AComputedData[] = {1., 4.,
                            2., 5.,
                            3., 6.};
    matrixSetReal(AComputedData, AComputed);

    vector* x = vectorCreateMalloc(vecSize);
    real xData[] = {24.,
                    9.,
                    13.,
                    -7.,
                    22.,
                    -9.,
                    -87.,
                    7.,
                    9.};
    vectorSetReal(xData, x);

    matrix* AExpected = matrixCreateMalloc(rowSize, colSize);
    real AExpectedData[] = {1., 7.,
                            2., 9.,
                            3., 24.};
    matrixSetReal(AExpectedData, AExpected);

    matrixCopyVectorToColumnCyclic(x, 16, 18, 1, AComputed);

    TEST_ASSERT_EQUAL_INT(1, matrixApproxEqual(AExpected, AComputed, epsilon));

    matrixFree(AComputed);
    matrixFree(AExpected);
    vectorFree(x);
}
