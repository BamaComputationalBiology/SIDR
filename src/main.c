//
//  main.c
//  SIDR 2.0
//
//  Created by Adam Case on 01/31/2021.
//  Copyright © 2021 Adam Case.
//  All rights reserved.
//

#include "param.h"
#include "pipeline.h"
#include "pylink.h"


int main(int argc, char *argv[]) {
//    printf("starting\n");
    Py_Initialize();
//    printf("after initializing\n"); 

// Calls param.c reads in arguments, gets file paths, etc
    PARAM param;
    initPARAM(&param, argc, argv);
//    printf("after initPARAM main.ci\n");

    int ret = 0;

// Calls pipeline.c
    ret |= run_pipeline(&param);
//    printf("after run_pipeline, main.c \n");
    ret |= run_analysis(&param);
//    printf("after run_analysis, main.c \n");

    freePARAM(&param);
//    printf("after freePARAM main.c\n");

    return ret;

}
