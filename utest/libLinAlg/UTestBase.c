/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
 */

/*
 *  Project:    libLinAlg
 *  Created on: April 26, 2016
 *  Author:     Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#include "../UTest.h"

void testMinMax() {
    TEST_ASSERT_EQUAL_INT(4, max(3, 4));
    TEST_ASSERT_EQUAL_INT(4, max(4, 3));
    TEST_ASSERT_EQUAL_INT(4, max(4, 4));

    TEST_ASSERT_EQUAL_INT(3, min(3, 4));
    TEST_ASSERT_EQUAL_INT(3, min(4, 3));
    TEST_ASSERT_EQUAL_INT(4, min(4, 4));
}
