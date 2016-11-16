/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
 */

/*
 *  Project:    libLinAlg
 *  Created on: March 9, 2016
 *  Author:     Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#include "../UTest.h"

void testReadCsv() {
    int nCoeff = 4;
    char name[] = "data/small.csv";

    vector* vecRead     = vectorCreateMalloc(nCoeff);
    readVectorCsv(name, nCoeff, vecRead);

    vector* vecExpected = vectorCreateMalloc(nCoeff);
    for (int i = 0; i < nCoeff; ++i) {
        vectorSetCoeff(vecExpected, i, i);
    }

    TEST_ASSERT_TRUE(vectorApproxEqual(vecRead, vecExpected, 1e-8));

    vectorFree(vecRead);
    vectorFree(vecExpected);
}

void testWriteCsv() {
    int nCoeff = 4;
    char name[] = "data/smallOutput.csv";

    vector* vecExpected = vectorCreateMalloc(nCoeff);
    for (int i = 0; i < nCoeff; ++i) {
        vectorSetCoeff(vecExpected, i, i);
    }
    writeVectorCsv(name,
                   vecExpected);

    vector* vecRead     = vectorCreateMalloc(nCoeff);
    readVectorCsv(name, nCoeff, vecRead);

    TEST_ASSERT_TRUE(vectorApproxEqual(vecRead, vecExpected, 1e-8));

    vectorFree(vecRead);
    vectorFree(vecExpected);
}
