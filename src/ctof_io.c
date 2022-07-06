/* ---------------------------------------------------------------------------------------------------------------------------------
   Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.
   See the README file in the top-level directory.

   Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).
   This software is distributed uner the modified Berkeley Software Distribution (BSD) License.
   Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.
--------------------------------------------------------------------------------------------------------------------------------- */

/* Reflect C low-level binary I/O functions to Fortran */

#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>

/* ------------------------------------------------------------------ */

typedef int longint;                                /* 4 byte integer */
typedef long longlongint;                           /* 8 byte integer */

/* ------------------------------------------------------------------ */

#ifdef   NO_UNDERSCORE
# define FTN(x) x
#elifdef DOUBLE_UNDERSCORE
# define FTN(x) x##__
#else
# define FTN(x) x##_
#endif

FILE *gfd;

void FTN(swapbytes)(char *baseptr, longint *dsize, longint *nelem);

/* ------------------------------------------------------------------ */

// Routine to open a file

FILE *FTN(openf)(char *name, char *mode)
{
  FILE *fd;

  gfd = fopen(name,mode);
  fd = gfd;

  return(fd);
}

/* ------------------------------------------------------------------ */

// Routine to close a file

int FTN(closef)(FILE **fd)
{
  return(fclose(*fd));
}

/* ------------------------------------------------------------------ */

// Routine to get the current file position

void FTN(getposf)(FILE **fd, longlongint *pos)
{
  *pos = ftell(*fd);
}

/* ------------------------------------------------------------------ */

// Routine to set the current file position

void FTN(setposf)(FILE **fd, longlongint *pos, longint *err)
{
  longlongint fp;

  fp = *pos;
  *err = fseek(*fd, fp, SEEK_SET);
}

/* ------------------------------------------------------------------ */

// Routine to set the current file position to the beginning of the file

void FTN(rewindf)(FILE **fd)
{
  rewind(*fd);
}

/* ------------------------------------------------------------------ */

// Routine get the size of the current file

longlongint FTN(cfilesize)(const char *fname) 
{
  struct stat buf;
  longlongint fs;

  stat(fname, &buf);

  fs = buf.st_size;
  return(fs);
}

/* ------------------------------------------------------------------ */

// Routine to read the current file

void FTN(readf)(void *ptr, longint *fsize, longint *fnmemb, FILE **fd, int *swap, longint *err)
{
  size_t size, nmemb;

  size = *fsize;
  nmemb = *fnmemb;
  *err = fread(ptr, size, nmemb, *fd);

  if ((*swap > 0) && (*fsize > 1)) FTN(swapbytes)(ptr, fsize, fnmemb);
}

/* ------------------------------------------------------------------ */

// Routine to write to the current file

void FTN(writef)(const void *ptr, longint *fsize, longint *fnmemb, FILE **fd, longint *err)
{
  size_t size, nmemb;

  size = *fsize;  nmemb = *fnmemb;
  *err = fwrite(ptr, size, nmemb, *fd);
}

/* ------------------------------------------------------------------ */

// Routine to determine byte ordering

int FTN(machinebyteorder)()
{
  int UNKNOWN_ENDIAN = 0;
  int BIG_ENDIAN     = 1;
  int LITTLE_ENDIAN  = 2;
  int MIDDLE_ENDIAN  = 3;
  int byteorder;
  longint num=0x12345678;
  unsigned char *cptr = (unsigned char *)&num;

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

/* ------------------------------------------------------------------ */

// Routine to swap bytes     

void FTN(swapbytes)(char *baseptr, longint *dsize, longint *nelem)
{
  char tmp;
  char *ptr0, *ptr1, *ptr;
  int  i, j, half;

  ptr = baseptr;
  half = *dsize/2;
  for (i=0; i<*nelem; i++)
  {
    ptr0 = ptr;
    ptr1 = ptr + *dsize - 1;
    for (j=0; j<half; j++)
    {
      tmp = *ptr0;
      *ptr0 = *ptr1;
      *ptr1 = tmp;
      ptr0++;
      ptr1--;
    }
    ptr = ptr + *dsize;
  }
}

/* ------------------------------------------------------------------ */
