
#ifndef DEFAULT_HEADER_INCLUDED
#define DEFAULT_HEADER_INCLUDED 1

#include "Python.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#if defined(__STDC__) || defined(__ANSI_CPP__)  /*###################*/
#define STRINGIFY(NAME) # NAME
#else
#define STRINGIFY(NAME) "NAME"
#endif

/* IF YOU USE CALLOC, YOUR FUNCTION MUST DEFINE A HANDLER WITH LABEL
   handle_malloc_failure:
   THIS HANDLER SHOULD RELEASE ANY TEMPORARILY ALLOCATED MEMORY AND
   return AN ERROR CODE SIGNALING ITS CALLER THAT IT WAS UNABLE TO 
   ALLOCATE THE DESIRED RESOURCE...  */
#define MALLOC_FAILURE_ACTION goto handle_malloc_failure

/* USE OF THIS MACRO WILL HANDLE MALLOC FAILURES BY
   RAISING THE APPROPRIATE PYTHON EXCEPTIONS
   INSTEAD OF CRASHING PYTHON OR FORCING THE PROCESS TO ABORT... */
#define CALLOC(memptr,N,ATYPE) \
  if ((N)<=0) {\
    char errstr[1024]; \
    sprintf(errstr,"%s, line %d: *** invalid memory request: %s[%d].\n",\
              __FILE__,__LINE__,STRINGIFY(memptr),(N));   \
    PyErr_SetString(PyExc_ValueError,errstr); \
    MALLOC_FAILURE_ACTION;\
  }\
  else if (NULL == ((memptr)=(ATYPE *)calloc((size_t)(N),sizeof(ATYPE))))  { \
    char errstr[1024]; \
    sprintf(errstr,"%s, line %d: memory request failed: %s[%d].\n",\
              __FILE__,__LINE__,STRINGIFY(memptr),(N));   \
    PyErr_SetString(PyExc_MemoryError,errstr); \
    MALLOC_FAILURE_ACTION;\
  }

#define FREE(P) if (P) {free(P);(P)=NULL;}

#endif
