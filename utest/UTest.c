/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
 */

/*
 *  Project:    libLinAlg
 *  Created on: March 9, 2016
 *  Author:     Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#include "UTest.h"

extern void testMinMax(void);

extern void testMatrixProduct(void);
extern void testMatrixVectorProduct(void);
extern void testMatrixTransposeVectorProduct(void);
extern void testCopyVectorToColumn(void);
extern void testCopyVectorToColumnCyclic(void);

extern void testMatrixCholesky(void);
extern void testMatrixCholesky2(void);
extern void testForwardSubstitutionSolve(void);
extern void testBackwardSubstitutionSolve(void);
extern void testSolveCholesky(void);
extern void testMatrixGramLower(void);
extern void testLeastSquare(void);
extern void testLeastSquare2(void);

extern void testVectorCreateMalloc(void);
extern void testVectorCreateRaw(void);
extern void testVecMean(void);
extern void testVecAddReal(void);
extern void testVecMulReal(void);
extern void testVecDotVec(void);
extern void testVectorCopy(void);
extern void testVectorNegShiftA(void);
extern void testVectorNegShiftB(void);
extern void testVectorTailCopyA(void);
extern void testVectorTailCopyB(void);
extern void testVectorTailCopyCycleB(void);
extern void testVectorPartialCopy(void);
extern void testVectorSetData(void);
extern void testVectorAddLeft(void);
extern void testVectorResizeSmaller(void);
extern void testVectorResizeSameSize(void);
extern void testVectorResizeLarger(void);
extern void testVectorMinMax(void);

extern void testCorrCoeff(void);
extern void testDecorr(void);
extern void testDecorr3(void);
extern void testVectorApproxEqual(void);
extern void testSmooth(void);
extern void testSmoothID(void);

extern void testCompleteProcessing(void);
extern void testComputeStabA(void);
extern void testComputeStabB(void);

extern void testReadCsv(void);
extern void testWriteCsv(void);

int main(void) {
  UnityBegin("");

  RUN_TEST(testMinMax);

  RUN_TEST(testMatrixProduct);
  RUN_TEST(testMatrixVectorProduct);
  RUN_TEST(testMatrixTransposeVectorProduct);
  RUN_TEST(testCopyVectorToColumn);
  RUN_TEST(testCopyVectorToColumnCyclic);

  RUN_TEST(testMatrixCholesky);
  RUN_TEST(testMatrixCholesky2);
  RUN_TEST(testForwardSubstitutionSolve);
  RUN_TEST(testBackwardSubstitutionSolve);
  RUN_TEST(testSolveCholesky);
  RUN_TEST(testMatrixGramLower);
  RUN_TEST(testLeastSquare);
  RUN_TEST(testLeastSquare2);

  RUN_TEST(testVectorCreateMalloc);
  RUN_TEST(testVectorCreateRaw);
  RUN_TEST(testVecMean);
  RUN_TEST(testVecAddReal);
  RUN_TEST(testVecMulReal);
  RUN_TEST(testVecDotVec);
  RUN_TEST(testVectorCopy);
  RUN_TEST(testVectorNegShiftA);
  RUN_TEST(testVectorNegShiftB);
  RUN_TEST(testVectorTailCopyA);
  RUN_TEST(testVectorTailCopyB);
  RUN_TEST(testVectorTailCopyCycleB);
  RUN_TEST(testVectorPartialCopy);
  RUN_TEST(testVectorSetData);
  RUN_TEST(testVectorAddLeft);
  RUN_TEST(testVectorResizeSmaller);
  RUN_TEST(testVectorResizeSameSize);
  RUN_TEST(testVectorResizeLarger);
  RUN_TEST(testVectorMinMax);

  RUN_TEST(testCorrCoeff);
  RUN_TEST(testDecorr);
  RUN_TEST(testDecorr3);
  RUN_TEST(testVectorApproxEqual);
  RUN_TEST(testSmooth);
  RUN_TEST(testSmoothID);

  RUN_TEST(testReadCsv);
  RUN_TEST(testWriteCsv);

  return (UnityEnd());
}
