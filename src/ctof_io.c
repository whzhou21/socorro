// Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
// of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

/* Reflect C low-level binary I/O functions to Fortran */


/*#include <sys/types.h>*/
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>

#include "ctof_io.h"

#ifdef NO_UNDERSCORE
# define FTN(x) x
#elifdef DOUBLE_UNDERSCORE
# define FTN(x) x##__
#else
# define FTN(x) x##_
#endif

FILE *gfd;


void FTN(swapbytes)(char *baseptr, longint *dsize, longint *nelem);

/*// ***************** open/close routines *************/

FILE *FTN(openf)(char *name, char *mode)
{
  FILE *fd;
  
    gfd = fopen(name,mode);
fd = gfd;
/*    printf("%d =fopen('%s','%s')\n",fd,name,mode); */
    return(fd);
  /*return(fopen(name, mode));*/
}

int FTN(closef)(FILE **fd)
{
  return(fclose(*fd));
}


/**************** move/pos routines ************/


void FTN(getposf)(FILE **fd, longlongint *pos)
{
  *pos = ftell(*fd);

}

void FTN(setposf)(FILE **fd, longlongint *pos, longint *err)
{
  longlongint fp;

  fp = *pos;
  *err = fseek(*fd, fp, SEEK_SET);
}


void FTN(rewindf)(FILE **fd)
{
  rewind(*fd);
}


longlongint FTN(cfilesize)(const char *fname) 
{
  struct stat buf;
  longlongint fs;

  stat(fname, &buf);

  fs = buf.st_size;
  return(fs);
}



/**************** read/write routines ************/

void FTN(readf)(void *ptr, longint *fsize, longint *fnmemb, FILE **fd, int *swap, longint *err)
{
  size_t size, nmemb;

  size = *fsize;  nmemb = *fnmemb;
  *err = fread(ptr, size, nmemb, *fd);

  if ((*swap > 0) && (*fsize > 1)) FTN(swapbytes)(ptr, fsize, fnmemb);
}

void FTN(writef)(const void *ptr, longint *fsize, longint *fnmemb, FILE **fd, longint *err)
{
  size_t size, nmemb;

  size = *fsize;  nmemb = *fnmemb;
/*  printf("writef('%1s',%d, %d, %d)\n",ptr, size, nmemb, *fd);*/
  *err = fwrite(ptr, size, nmemb, *fd);
}

/************** Determine byte ordering *******/

int FTN(machinebyteorder)() {
  int UNKNOWN_ENDIAN = 0;
  int BIG_ENDIAN     = 1;
  int LITTLE_ENDIAN  = 2;
  int MIDDLE_ENDIAN  = 3;

  longint num=0x12345678;
  unsigned char *cptr = (unsigned char *)&num;
  
  int byteorder;

  if ( *cptr==0x12     && *(cptr+1)==0x34 &&
       *(cptr+2)==0x56 && *(cptr+3)==0x78 )
    byteorder = BIG_ENDIAN;
  else if ( *cptr==0x78     && *(cptr+1)==0x56 &&
            *(cptr+2)==0x34 && *(cptr+3)==0x12 )
    byteorder = LITTLE_ENDIAN;
  else if ( *cptr==0x34     && *(cptr+1)==0x12 &&
            *(cptr+2)==0x78 && *(cptr+3)==0x56 )
    byteorder = MIDDLE_ENDIAN;
  else
    byteorder = UNKNOWN_ENDIAN;

  return(byteorder);
}
     
/************** Swap bytes *******/

void FTN(swapbytes)(char *baseptr, longint *dsize, longint *nelem) {
  char tmp;
  char *ptr0, *ptr1, *ptr;
  int  i, j, half;
  
  ptr = baseptr;
  half = *dsize/2;
  for (i=0; i<*nelem; i++) {
    
    ptr0 = ptr;
    ptr1 = ptr + *dsize - 1;

    for (j=0; j<half; j++) {
      tmp = *ptr0;
      *ptr0 = *ptr1;
      *ptr1 = tmp;

      ptr0++;
      ptr1--;
    }

    ptr = ptr + *dsize;
  }

}

/* --- */
