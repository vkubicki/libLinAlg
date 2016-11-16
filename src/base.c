/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
*/

/*
 *  Project:    libLinAlg
 *  Created on: April 26, 2016
 *  Authors:    Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#ifndef BASE_H
#define BASE_H

#include "typedef.h"

integer max(integer a,
            integer b) {
    return (a < b) ? b: a;
}

integer min(integer a,
            integer b) {
    return (a < b) ? a: b;
}

indexCoeff mod(integer a,
               integer b) {
    integer r = a % b;
    return r < 0 ? r + b : r;
}

#endif
