/*--------------------------------------------------------------------*/
/*  Copyright (C) Inria 2016
*/

/*
 *  Project:    libLinAlg
 *  Created on: April 1, 2016
 *  Authors:    Vincent KUBICKI <vincent.kubicki@inria.fr>
 **/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "IO.h"

void readVectorCsv(const char* fileName,
                   int sizeVec,
                   vector* vec) {
    vec->nI_ = sizeVec;
    int bufferSize = 1024;
    FILE *inputFile = fopen(fileName, "rb");

    char line[bufferSize];
    for (int i = 0; i < sizeVec; ++i) {
        fgets(line, bufferSize, inputFile);
        real val = strtod(line, NULL);
        vectorSetCoeff(vec, i, val);
    }

    fclose(inputFile);
}

void writeVectorCsv(const char* fileName,
                    const vector* vec) {
    int nI = vectorGetDim(vec);
    FILE *outFile = fopen(fileName, "wb");
    real val;
    for(int i = 0; i < nI; ++i) {
        val = vectorGetCoeff(vec, i);
        fprintf(outFile, "%f\n", val);
    }


    fclose(outFile);
}
