/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
*/

/*
 *  Project:    libLinAlg
 *  Created on: April 1, 2016
 *  Authors:    Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#ifndef IO_H
#define IO_H

#include "vector.h"

void readVectorCsv(const char* fileName,
                   int sizeVec,
                   vector* vec);

void writeVectorCsv(const char* fileName,
                    const vector* vec);
#endif
