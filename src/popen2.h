#ifndef _POPEN_2_H
#define _POPEN_2_H

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#define READ 0
#define WRITE 1

pid_t
popen2(const char *command, const char* pars, int *infp, int *outfp);

#endif
