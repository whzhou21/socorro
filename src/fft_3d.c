/* ---------------------------------------------------------------------------------------------------------------------------------
   Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.
   See the README file in the top-level directory.

   Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).
   This software is distributed uner the modified Berkeley Software Distribution (BSD) License.
   Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.
--------------------------------------------------------------------------------------------------------------------------------- */

/* 3d FFT and data remapping functions

   Original author: Steve Plimpton, Sandia National Laboratories
   Original code: Parallel FFT Package - 1998, 1999
*/

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/* ------------------------------------------------------------------ */

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#ifdef  NO_UNDERSCORE
#define FTN(x) x
#elif   DOUBLE_UNDERSCORE
#define FTN(x) x##__
#else
#define FTN(x) x##_
#endif

#if (MPI_VERSION == 1)
#define FOR_MPI_COMM MPI_Comm
#else
#define FOR_MPI_COMM MPI_Fint
#endif

#if !defined(PACK_POINTER) && !defined(PACK_MEMCPY)
#define PACK_ARRAY
#endif

#define FFT_FFTW
#define FFT_PRECISION 2

/* ------------------------------------------------------------------ */

/* loop counters for doing a pack/unpack */

struct pack_plan_3d {
  int nfast;                 /* # of elements in fast index */
  int nmid;                  /* # of elements in mid index */
  int nslow;                 /* # of elements in slow index */
  int nstride_line;          /* stride between successive mid indices */
  int nstride_plane;         /* stride between successive slow indices */
  int nqty;                  /* # of values/element */
};

/* function prototypes */

void pack_3d(double *, double *, struct pack_plan_3d *);
void unpack_3d(double *, double *, struct pack_plan_3d *);
void unpack_3d_permute1_1(double *, double *, struct pack_plan_3d *);
void unpack_3d_permute1_2(double *, double *, struct pack_plan_3d *);
void unpack_3d_permute1_n(double *, double *, struct pack_plan_3d *);
void unpack_3d_permute2_1(double *, double *, struct pack_plan_3d *);
void unpack_3d_permute2_2(double *, double *, struct pack_plan_3d *);
void unpack_3d_permute2_n(double *, double *, struct pack_plan_3d *);

/* details of how to do a 3d remap */

struct remap_plan_3d {
  double *sendbuf;                  /* buffer for MPI sends */
  double *scratch;                  /* scratch buffer for MPI recvs */
  void (*pack)();                   /* which pack function to use */
  void (*unpack)();                 /* which unpack function to use */
  int *send_offset;                 /* extraction loc for each send */
  int *send_size;                   /* size of each send message */
  int *send_stride;
  int *send_proc;                   /* proc to send each message to */
int send_mallocd;                  /* boolean that tracks whether send the send arrays were malloc'd */
  struct pack_plan_3d *packplan;    /* pack plan for each send message */
  int *recv_offset;                 /* insertion loc for each recv */
  int *recv_size;                   /* size of each recv message */
  int *recv_stride;
  int *recv_proc;                   /* proc to recv each message from */
  int *recv_bufloc;                 /* offset in scratch buf for each recv */
int recv_mallocd;                  /* boolean that tracks whether  the receive arrays were malloc'd */
  MPI_Request *request;             /* MPI request for each posted recv */
  struct pack_plan_3d *unpackplan;  /* unpack plan for each recv message */
  int nrecv;                        /* # of recvs from other procs */
  int nsend;                        /* # of sends to other procs */
  int self;                         /* whether I send/recv with myself */
  int howmany;
  int in_stride;
  int out_stride;
  int memory;                       /* user provides scratch space or not */
  MPI_Comm comm;                    /* group of procs performing remap */
};

/* collision between 2 regions */

struct extent_3d {
  int ilo,ihi,isize;
  int jlo,jhi,jsize;
  int klo,khi,ksize;
};

/* function prototypes */

void remap_3d_cfunc(double *, double *, double *, struct remap_plan_3d *);
struct remap_plan_3d *remap_3d_create_plan_cfunc(MPI_Comm, 
  int, int, int, int, int, int,	int, int, int, int, int, int,
  int, int, int, int, int);
void remap_3d_destroy_plan_cfunc(struct remap_plan_3d *);
int remap_3d_collide(struct extent_3d *, 
		     struct extent_3d *, struct extent_3d *);

/* machine specifics */

#ifdef T3E_KLUDGE

#define remap_3d_ REMAP_3D
#define remap_3d_create_plan_ REMAP_3D_CREATE_PLAN
#define remap_3d_destroy_plan_ REMAP_3D_DESTROY_PLAN
#define MPI_DOUBLE_PRECISION MPI_DOUBLE

#endif

/* ------------------------------------------------------------------------- */

/* Data types for single-precision complex */

#if FFT_PRECISION == 1

#ifdef FFT_SGI
#include <fft.h>
typedef complex FFT_DATA;
#define FFT_1D cfft1d
#define FFT_1D_INIT cfft1di
#endif

#ifdef FFT_INTEL
typedef struct {
  float re;
  float im;
} FFT_DATA;
#define FFT_1D cfft1d_
#define FFT_1D_INIT cfft1d_
#endif

#ifdef FFT_DEC
typedef struct {
  float re;
  float im;
} FFT_DATA;
#define FFT_1D cfft_
#endif

#ifdef FFT_T3E
#include <complex.h>
typedef complex single FFT_DATA;
#define FFT_1D GGFFT
#define FFT_1D_INIT GGFFT
#endif

#ifdef FFT_FFTW
#include <fftw.h>
typedef FFTW_COMPLEX FFT_DATA;
#endif

#endif

/* ------------------------------------------------------------------------- */

/* Data types for double-precision complex */

#if FFT_PRECISION == 2

#ifdef FFT_SGI
#include <fft.h>
typedef zomplex FFT_DATA;
#define FFT_1D zfft1d
#define FFT_1D_INIT zfft1di
#endif

#ifdef FFT_INTEL
typedef struct {
  double re;
  double im;
} FFT_DATA;
#define FFT_1D zfft1d_
#define FFT_1D_INIT zfft1d_
#endif

#ifdef FFT_DEC
typedef struct {
  double re;
  double im;
} FFT_DATA;
#define FFT_1D zfft_
#endif

#ifdef FFT_T3E
#include <complex.h>
typedef complex double FFT_DATA;
#define FFT_1D CCFFT
#define FFT_1D_INIT CCFFT
#endif

#ifdef FFT_FFTW
#include <fftw.h>
typedef FFTW_COMPLEX FFT_DATA;
#endif

#endif

/* ------------------------------------------------------------------------- */

/* details of how to do a 3d FFT */

struct fft_plan_3d {
  struct remap_plan_3d *pre_plan;       /* remap from input -> 1st FFTs */
  struct remap_plan_3d *mid1_plan;      /* remap from 1st -> 2nd FFTs */
  struct remap_plan_3d *mid2_plan;      /* remap from 2nd -> 3rd FFTs */
  struct remap_plan_3d *post_plan;      /* remap from 3rd FFTs -> output */
  FFT_DATA *copy;                   /* memory for remap results (if needed) */
  FFT_DATA *scratch;                /* scratch space for remaps */
  int total1,total2,total3;         /* # of 1st,2nd,3rd FFTs (times length) */
  int length1,length2,length3;      /* length of 1st,2nd,3rd FFTs */
  int pre_target;                   /* where to put remap results */
  int mid1_target,mid2_target;
  int scaled;                       /* whether to scale FFT results */
  int fftsize;                      /* # of values to rescale */
  int howmany;                      /* number of contiguous ffts to perform */
  double norm;                      /* normalization factor for rescaling */
                                    /* system specific 1d FFT info */
#ifdef FFT_SGI
  FFT_DATA *coeff1;
  FFT_DATA *coeff2;
  FFT_DATA *coeff3;
#endif
#ifdef FFT_INTEL
  FFT_DATA *coeff1;
  FFT_DATA *coeff2;
  FFT_DATA *coeff3;
#endif
#ifdef FFT_T3E
  double *coeff1;
  double *coeff2;
  double *coeff3;
  double *work1;
  double *work2;
  double *work3;
#endif
#ifdef FFT_FFTW
  fftw_plan plan_fast_forward;
  fftw_plan plan_fast_backward;
  fftw_plan plan_mid_forward;
  fftw_plan plan_mid_backward;
  fftw_plan plan_slow_forward;
  fftw_plan plan_slow_backward;
#endif
};

/* function prototypes */

void fft_3d_cfunc(FFT_DATA *, FFT_DATA *, int, struct fft_plan_3d *);
struct fft_plan_3d *fft_3d_create_plan_cfunc(MPI_Comm, int, int, int,
  int, int, int, int, int, int, int, int, int, int, int, int,
  int, int, int, int *);
void fft_3d_destroy_plan_cfunc(struct fft_plan_3d *);
void factor(int, int *, int *);
void bifactor(int, int *, int *);

/* ------------------------------------------------------------------- */
/* F77 wrappers on fft */

void FTN(fft_3d)(FFT_DATA *in, FFT_DATA *out,
	     int *flag, struct fft_plan_3d **plan)

{
  fft_3d_cfunc(in,out,*flag,*plan);
}

void FTN(fft_3d_ip)(FFT_DATA *inout,
	     int *flag, struct fft_plan_3d **plan)

{
  fft_3d_cfunc(inout,inout,*flag,*plan);
}

/* ------------------------------------------------------------------- */
/* F77 wrapper on fft_create_plan */

void FTN(fft_3d_create_plan)(
       FOR_MPI_COMM *for_comm, int *nfast, int *nmid, int *nslow,
       int *in_ilo, int *in_ihi, int *in_jlo, int *in_jhi,
       int *in_klo, int *in_khi,
       int *out_ilo, int *out_ihi, int *out_jlo, int *out_jhi,
       int *out_klo, int *out_khi,
       int *howmany,
       int *scaled, int *permute, int *nbuf,
       struct fft_plan_3d **plan)

{
  int me;
  MPI_Comm comm;

#if (MPI_VERSION == 1)      //Pass the comm correctly
  comm = *for_comm;
#else
  comm = MPI_Comm_f2c(*for_comm);
#endif


/* convert F77 indices into C */

  *plan = fft_3d_create_plan_cfunc(comm,*nfast,*nmid,*nslow,
			     *in_ilo-1,*in_ihi-1,*in_jlo-1,*in_jhi-1,
			     *in_klo-1,*in_khi-1,
			     *out_ilo-1,*out_ihi-1,*out_jlo-1,*out_jhi-1,
			     *out_klo-1,*out_khi-1,
			     *howmany,
			     *scaled, *permute, nbuf);
  if (*plan == NULL) {
    MPI_Comm_rank(comm,&me);
    printf("ERROR: FFT 3d plan is NULL on proc %d\n",me);
  }
}

/* ------------------------------------------------------------------- */
/* F77 wrapper on fft_destroy_plan */

void FTN(fft_3d_destroy_plan)(struct fft_plan_3d **plan)

{
  fft_3d_destroy_plan_cfunc(*plan);
}

/* parallel FFT functions - 1998, 1999

   Steve Plimpton, MS 1111, Dept 9221, Sandia National Labs
   (505) 845-7873
   sjplimp@cs.sandia.gov
*/

/* ------------------------------------------------------------------- */
/* Data layout for 3d FFTs:

   data set of Nfast x Nmid x Nslow elements is owned by P procs
   on input, each proc owns a subsection of the elements
   on output, each proc will own a (possibly different) subsection
   my subsection must not overlap with any other proc's subsection,
     i.e. the union of all proc's input (or output) subsections must
     exactly tile the global Nfast x Nmid x Nslow data set
   when called from C, all subsection indices are 
     C-style from 0 to N-1 where N = Nfast or Nmid or Nslow
   when called from F77, all subsection indices are 
     F77-style from 1 to N where N = Nfast or Nmid or Nslow
   a proc can own 0 elements on input or output
     by specifying hi index < lo index
   on both input and output, data is stored contiguously on a processor
     with a fast-varying, mid-varying, and slow-varying index
*/
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
/* Perform 3d FFT */

/* Arguments:

   in           starting address of input data on this proc
   out          starting address of where output data for this proc
                  will be placed (can be same as in)
   flag         1 for forward FFT, -1 for inverse FFT
   plan         plan returned by previous call to fft_3d_create_plan
*/

void fft_3d_cfunc(FFT_DATA *in, FFT_DATA *out, int flag, struct fft_plan_3d *plan)
{
  int i,total,length,offset,num;
  double norm;
  int hackloop,trashme;
  FFT_DATA *data,*copy;

/* system specific constants */

#ifdef FFT_DEC
  char c = 'C';
  char f = 'F';
  char b = 'B';
  int one = 1;
#endif
#ifdef FFT_T3E
  int isys = 0;
  double scalef = 1.0;
#endif

/* pre-remap to prepare for 1st FFTs if needed
   copy = loc for remap result */

  if (plan->pre_plan) {
    if (plan->pre_target == 0)
      copy = out;
    else
      copy = plan->copy;
    remap_3d_cfunc((double *) in,
	     (double *) copy,
	     (double *) plan->scratch,
	     plan->pre_plan);
    data = copy;
  }
  else 
    data = in;

/* 1d FFTs along fast axis */
  total = plan->total1*plan->howmany;
  length = plan->length1;

#ifdef FFT_SGI
  for (offset = 0; offset < total; offset += length)
    FFT_1D(flag,length,&data[offset],1,plan->coeff1);
#endif
#ifdef FFT_INTEL
  for (offset = 0; offset < total; offset += length)
    FFT_1D(&data[offset],&length,&flag,plan->coeff1);
#endif
#ifdef FFT_DEC
  if (flag == -1)
    for (offset = 0; offset < total; offset += length)
      FFT_1D(&c,&c,&f,&data[offset],&data[offset],&length,&one);
  else
    for (offset = 0; offset < total; offset += length)
      FFT_1D(&c,&c,&b,&data[offset],&data[offset],&length,&one);
#endif
#ifdef FFT_T3E
  for (offset = 0; offset < total; offset += length)
    FFT_1D(&flag,&length,&scalef,&data[offset],&data[offset],plan->coeff1,
	   plan->work1,&isys);
#endif
#ifdef FFT_FFTW
  if (flag == -1)
    fftw(plan->plan_fast_forward,total/length,data,1,length,NULL,0,0);
  else
    fftw(plan->plan_fast_backward,total/length,data,1,length,NULL,0,0);
#endif

/* 1st mid-remap to prepare for 2nd FFTs
   copy = loc for remap result */

  if (plan->mid1_target == 0)
    copy = out;
  else
    copy = plan->copy;
  remap_3d_cfunc((double *) data,
	   (double *) copy,
	   (double *) plan->scratch,
	   plan->mid1_plan);
  data = copy;

/* 1d FFTs along mid axis */
  total = plan->total2*plan->howmany;
  length = plan->length2;

#ifdef FFT_SGI
  for (offset = 0; offset < total; offset += length)
    FFT_1D(flag,length,&data[offset],1,plan->coeff2);
#endif
#ifdef FFT_INTEL
  for (offset = 0; offset < total; offset += length)
    FFT_1D(&data[offset],&length,&flag,plan->coeff2);
#endif
#ifdef FFT_DEC
  if (flag == -1)
    for (offset = 0; offset < total; offset += length)
      FFT_1D(&c,&c,&f,&data[offset],&data[offset],&length,&one);
  else
    for (offset = 0; offset < total; offset += length)
      FFT_1D(&c,&c,&b,&data[offset],&data[offset],&length,&one);
#endif
#ifdef FFT_T3E
  for (offset = 0; offset < total; offset += length)
    FFT_1D(&flag,&length,&scalef,&data[offset],&data[offset],plan->coeff2,
	   plan->work2,&isys);
#endif
#ifdef FFT_FFTW
  if (flag == -1)
    fftw(plan->plan_mid_forward,total/length,data,1,length,NULL,0,0);
  else
    fftw(plan->plan_mid_backward,total/length,data,1,length,NULL,0,0);
#endif

/* 2nd mid-remap to prepare for 3rd FFTs
   copy = loc for remap result */

  if (plan->mid2_target == 0) 
    copy = out;
  else 
    copy = plan->copy;
  remap_3d_cfunc((double *) data,
	   (double *) copy,
	   (double *) plan->scratch,
	   plan->mid2_plan);
  data = copy;

/* 1d FFTs along slow axis */
  total = plan->total3*plan->howmany;
  length = plan->length3;

#ifdef FFT_SGI
  for (offset = 0; offset < total; offset += length)
    FFT_1D(flag,length,&data[offset],1,plan->coeff3);
#endif
#ifdef FFT_INTEL
  for (offset = 0; offset < total; offset += length)
    FFT_1D(&data[offset],&length,&flag,plan->coeff3);
#endif
#ifdef FFT_DEC
  if (flag == -1)
    for (offset = 0; offset < total; offset += length)
      FFT_1D(&c,&c,&f,&data[offset],&data[offset],&length,&one);
  else
    for (offset = 0; offset < total; offset += length)
      FFT_1D(&c,&c,&b,&data[offset],&data[offset],&length,&one);
#endif
#ifdef FFT_T3E
  for (offset = 0; offset < total; offset += length)
    FFT_1D(&flag,&length,&scalef,&data[offset],&data[offset],plan->coeff3,
	   plan->work3,&isys);
#endif
#ifdef FFT_FFTW
  if (flag == -1)
    fftw(plan->plan_slow_forward,total/length,data,1,length,NULL,0,0);
  else
    fftw(plan->plan_slow_backward,total/length,data,1,length,NULL,0,0);
#endif

/* post-remap to put data in output format if needed
   destination is always out */

  if (plan->post_plan)
    remap_3d_cfunc((double *) data, (double *) out, (double *) plan->scratch,
	     plan->post_plan);


/* scaling if required */

#ifndef FFT_T3E
  if (flag == 1 && plan->scaled) {
    norm = plan->norm;
    num = plan->fftsize*plan->howmany;
    for (i = 0; i < num; i++) {
      out[i].re *= norm;
      out[i].im *= norm;
    }
  }
#endif

#ifdef FFT_T3E
  if (flag == 1 && plan->scaled) {
    norm = plan->norm;
    num = plan->fftsize*plan->howmany;
    for (i = 0; i < num; i++)
      out[i] *= (norm,norm);
  }
#endif

}
/****
void fft_3d_cfunc(FFT_DATA *in, FFT_DATA *out, int flag, struct fft_plan_3d *plan) {
  int i;
  for (i=0;i<plan->howmany;i++) {
    fft_3d_looper(in+i*plan->fftsize,out+i*plan->fftsize,flag,plan);
  }
}
***/

/* ------------------------------------------------------------------- */
/* Create plan for performing a 3d FFT */

/* Arguments:

   comm                 MPI communicator for the P procs which own the data
   nfast,nmid,nslow     size of global 3d matrix
   in_ilo,in_ihi        input bounds of data I own in fast index
   in_jlo,in_jhi        input bounds of data I own in mid index
   in_klo,in_khi        input bounds of data I own in slow index
   out_ilo,out_ihi      output bounds of data I own in fast index
   out_jlo,out_jhi      output bounds of data I own in mid index
   out_klo,out_khi      output bounds of data I own in slow index
   howmany              number of contiguous ffts to perform
   scaled               0 = no scaling of result, 1 = scaling
   permute              permutation in storage order of indices on output
                          0 = no permutation
			  1 = permute once = mid->fast, slow->mid, fast->slow
			  2 = permute twice = slow->fast, fast->mid, mid->slow
   nbuf                 returns size of internal storage buffers used by FFT
*/

struct fft_plan_3d *fft_3d_create_plan_cfunc(
       MPI_Comm comm, int nfast, int nmid, int nslow,
       int in_ilo, int in_ihi, int in_jlo, int in_jhi,
       int in_klo, int in_khi,
       int out_ilo, int out_ihi, int out_jlo, int out_jhi,
       int out_klo, int out_khi,
       int howmany,
       int scaled, int permute, int *nbuf)

{
  struct fft_plan_3d *plan;
  int me,nprocs;
  int i,num,flag,remapflag,fftflag;
  int first_ilo,first_ihi,first_jlo,first_jhi,first_klo,first_khi;
  int second_ilo,second_ihi,second_jlo,second_jhi,second_klo,second_khi;
  int third_ilo,third_ihi,third_jlo,third_jhi,third_klo,third_khi;
  int out_size,first_size,second_size,third_size,copy_size,scratch_size;
  int np1,np2,ip1,ip2;
  int list[50];

/* system specific variables */

#ifdef FFT_INTEL
  FFT_DATA dummy;
#endif
#ifdef FFT_T3E
  FFT_DATA dummy[5];
  int isign,isys;
  double scalef;
#endif

/* query MPI info */

  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

/* compute division of procs in 2 dimensions not on-processor */

  bifactor(nprocs,&np1,&np2);
  ip1 = me % np1;
  ip2 = me/np1;

/* allocate memory for plan data struct */

  plan = (struct fft_plan_3d *) malloc(sizeof(struct fft_plan_3d));
  if (plan == NULL) return NULL;

/* remap from initial distribution to layout needed for 1st set of 1d FFTs
   not needed if all procs own entire fast axis initially
   first indices = distribution after 1st set of FFTs */

  if (in_ilo == 0 && in_ihi == nfast-1)
    flag = 0;
  else
    flag = 1;

  MPI_Allreduce(&flag,&remapflag,1,MPI_INT,MPI_MAX,comm);

  if (remapflag == 0) {
    first_ilo = in_ilo;
    first_ihi = in_ihi;
    first_jlo = in_jlo;
    first_jhi = in_jhi;
    first_klo = in_klo;
    first_khi = in_khi;
    plan->pre_plan = NULL;
  }
  else {
    first_ilo = 0;
    first_ihi = nfast - 1;
    first_jlo = ip1*nmid/np1;
    first_jhi = (ip1+1)*nmid/np1 - 1;
    first_klo = ip2*nslow/np2;
    first_khi = (ip2+1)*nslow/np2 - 1;
    plan->pre_plan =
      remap_3d_create_plan_cfunc(comm,in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
			   first_ilo,first_ihi,first_jlo,first_jhi,
			   first_klo,first_khi,howmany,
			   FFT_PRECISION,0,0,2);
    if (plan->pre_plan == NULL) return NULL;
  }

/* 1d FFTs along fast axis */

  plan->length1 = nfast;
  plan->total1 = nfast * (first_jhi-first_jlo+1) * (first_khi-first_klo+1);

/* remap from 1st to 2nd FFT
   choose which axis is split over np1 vs np2 to minimize communication
   second indices = distribution after 2nd set of FFTs */

  second_ilo = ip1*nfast/np1;
  second_ihi = (ip1+1)*nfast/np1 - 1;
  second_jlo = 0;
  second_jhi = nmid - 1;
  second_klo = ip2*nslow/np2;
  second_khi = (ip2+1)*nslow/np2 - 1;
  plan->mid1_plan =
      remap_3d_create_plan_cfunc(comm,
			   first_ilo,first_ihi,first_jlo,first_jhi,
			   first_klo,first_khi,
			   second_ilo,second_ihi,second_jlo,second_jhi,
			   second_klo,second_khi,howmany,
			   FFT_PRECISION,1,0,2);

  if (plan->mid1_plan == NULL) return NULL;

/* 1d FFTs along mid axis */

  plan->length2 = nmid;
  plan->total2 = (second_ihi-second_ilo+1) * nmid * (second_khi-second_klo+1);

/* remap from 2nd to 3rd FFT
   if final distribution is permute=2 with all procs owning entire slow axis
     then this remapping goes directly to final distribution
   third indices = distribution after 3rd set of FFTs */

  if (permute == 2 && out_klo == 0 && out_khi == nslow-1)
    flag = 0;
  else
    flag = 1;

  MPI_Allreduce(&flag,&remapflag,1,MPI_INT,MPI_MAX,comm);

  if (remapflag == 0) {
    third_ilo = out_ilo;
    third_ihi = out_ihi;
    third_jlo = out_jlo;
    third_jhi = out_jhi;
    third_klo = out_klo;
    third_khi = out_khi;
  }
  else {
    third_ilo = ip1*nfast/np1;
    third_ihi = (ip1+1)*nfast/np1 - 1;
    third_jlo = ip2*nmid/np2;
    third_jhi = (ip2+1)*nmid/np2 - 1;
    third_klo = 0;
    third_khi = nslow - 1;
  }
  
  plan->mid2_plan =
    remap_3d_create_plan_cfunc(comm,
			 second_jlo,second_jhi,second_klo,second_khi,
			 second_ilo,second_ihi,
			 third_jlo,third_jhi,third_klo,third_khi,
			 third_ilo,third_ihi,howmany,
			 FFT_PRECISION,1,0,2);
  if (plan->mid2_plan == NULL) return NULL;

/* 1d FFTs along slow axis */

  plan->length3 = nslow;
  plan->total3 = (third_ihi-third_ilo+1) * (third_jhi-third_jlo+1) * nslow;

/* remap from 3rd FFT to final distribution
   not needed if permute = 2 and third indices = out indices on all procs */

  if (permute == 2 &&
      out_ilo == third_ilo && out_ihi == third_ihi &&
      out_jlo == third_jlo && out_jhi == third_jhi &&
      out_klo == third_klo && out_khi == third_khi)
    flag = 0;
  else
    flag = 1;

  MPI_Allreduce(&flag,&remapflag,1,MPI_INT,MPI_MAX,comm);

  if (remapflag == 0)
    plan->post_plan = NULL;
  else {
    plan->post_plan =
      remap_3d_create_plan_cfunc(comm,
			   third_klo,third_khi,third_ilo,third_ihi,
			   third_jlo,third_jhi,
			   out_klo,out_khi,out_ilo,out_ihi,
			   out_jlo,out_jhi,howmany,
			   FFT_PRECISION,(permute+1)%3,0,2);
    if (plan->post_plan == NULL) return NULL;
  }

/* configure plan memory pointers and allocate work space
   out_size = amount of memory given to FFT by user
   first/second/third_size = amount of memory needed after pre,mid1,mid2 remaps
   copy_size = amount needed internally for extra copy of data
   scratch_size = amount needed internally for remap scratch space
   for each remap:
     use out space for result if big enough, else require copy buffer
     accumulate largest required remap scratch space */

  plan->fftsize = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
      (out_khi-out_klo+1);
  
  out_size = plan->fftsize;
  first_size = (first_ihi-first_ilo+1) * (first_jhi-first_jlo+1) * 
    (first_khi-first_klo+1);
  second_size = (second_ihi-second_ilo+1) * (second_jhi-second_jlo+1) * 
    (second_khi-second_klo+1);
  third_size = (third_ihi-third_ilo+1) * (third_jhi-third_jlo+1) * 
    (third_khi-third_klo+1);

  copy_size = 0;
  scratch_size = 0;

  if (plan->pre_plan) {
    if (first_size <= out_size)
      plan->pre_target = 0;
    else {
      plan->pre_target = 1;
      copy_size = MAX(copy_size,first_size);
    }
    scratch_size = MAX(scratch_size,first_size);
  }

  if (plan->mid1_plan) {
    if (second_size <= out_size)
      plan->mid1_target = 0;
    else {
      plan->mid1_target = 1;
      copy_size = MAX(copy_size,second_size);
    }
    scratch_size = MAX(scratch_size,second_size);
  }

  if (plan->mid2_plan) {
    if (third_size <= out_size)
      plan->mid2_target = 0;
    else {
      plan->mid2_target = 1;
      copy_size = MAX(copy_size,third_size);
    }
    scratch_size = MAX(scratch_size,third_size);
  }

  if (plan->post_plan)
    scratch_size = MAX(scratch_size,out_size);

  *nbuf = copy_size + scratch_size;

  if (copy_size) {
    plan->copy = (FFT_DATA *) malloc(howmany*copy_size*sizeof(FFT_DATA));
    if (plan->copy == NULL) return NULL;
  }
  else plan->copy = NULL;

  if (scratch_size) {
    plan->scratch = (FFT_DATA *) malloc(howmany*scratch_size*sizeof(FFT_DATA));
    if (plan->scratch == NULL) return NULL;
  }
  else plan->scratch = NULL;

  plan->howmany = howmany;

/* system specific pre-computation of 1d FFT coeffs 
   and scaling normalization */

#ifdef FFT_SGI

  plan->coeff1 = (FFT_DATA *) malloc((nfast+15)*sizeof(FFT_DATA));
  plan->coeff2 = (FFT_DATA *) malloc((nmid+15)*sizeof(FFT_DATA));
  plan->coeff3 = (FFT_DATA *) malloc((nslow+15)*sizeof(FFT_DATA));

  if (plan->coeff1 == NULL || plan->coeff2 == NULL ||
      plan->coeff3 == NULL) return NULL;

  FFT_1D_INIT(nfast,plan->coeff1);
  FFT_1D_INIT(nmid,plan->coeff2);
  FFT_1D_INIT(nslow,plan->coeff3);

  if (scaled == 0) 
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nmid*nslow);
  }

#endif

#ifdef FFT_INTEL

  flag = 0;

  num = 0;
  factor(nfast,&num,list);
  for (i = 0; i < num; i++)
    if (list[i] != 2 && list[i] != 3 && list[i] != 5) flag = 1;
  num = 0;
  factor(nmid,&num,list);
  for (i = 0; i < num; i++)
    if (list[i] != 2 && list[i] != 3 && list[i] != 5) flag = 1;
  num = 0;
  factor(nslow,&num,list);
  for (i = 0; i < num; i++)
    if (list[i] != 2 && list[i] != 3 && list[i] != 5) flag = 1;

  MPI_Allreduce(&flag,&fftflag,1,MPI_INT,MPI_MAX,comm);
  if (fftflag) {
    if (me == 0) printf("ERROR: FFTs are not power of 2,3,5\n");
    return NULL;
  }

  plan->coeff1 = (FFT_DATA *) malloc((3*nfast/2+1)*sizeof(FFT_DATA));
  plan->coeff2 = (FFT_DATA *) malloc((3*nmid/2+1)*sizeof(FFT_DATA));
  plan->coeff3 = (FFT_DATA *) malloc((3*nslow/2+1)*sizeof(FFT_DATA));

  if (plan->coeff1 == NULL || plan->coeff2 == NULL || 
      plan->coeff3 == NULL) return NULL;

  flag = 0;
  FFT_1D_INIT(&dummy,&nfast,&flag,plan->coeff1);
  FFT_1D_INIT(&dummy,&nmid,&flag,plan->coeff2);
  FFT_1D_INIT(&dummy,&nslow,&flag,plan->coeff3);

  if (scaled == 0) {
    plan->scaled = 1;
    plan->norm = nfast*nmid*nslow;
  }
  else
    plan->scaled = 0;

#endif

#ifdef FFT_DEC

  if (scaled == 0) {
    plan->scaled = 1;
    plan->norm = nfast*nmid*nslow;
  }
  else
    plan->scaled = 0;

#endif

#ifdef FFT_T3E

  plan->coeff1 = (double *) malloc((12*nfast)*sizeof(double));
  plan->coeff2 = (double *) malloc((12*nmid)*sizeof(double));
  plan->coeff3 = (double *) malloc((12*nslow)*sizeof(double));

  if (plan->coeff1 == NULL || plan->coeff2 == NULL || 
      plan->coeff3 == NULL) return NULL;

  plan->work1 = (double *) malloc((8*nfast)*sizeof(double));
  plan->work2 = (double *) malloc((8*nmid)*sizeof(double));
  plan->work3 = (double *) malloc((8*nslow)*sizeof(double));

  if (plan->work1 == NULL || plan->work2 == NULL || 
      plan->work3 == NULL) return NULL;

  isign = 0;
  scalef = 1.0;
  isys = 0;

  FFT_1D_INIT(&isign,&nfast,&scalef,dummy,dummy,plan->coeff1,dummy,&isys);
  FFT_1D_INIT(&isign,&nmid,&scalef,dummy,dummy,plan->coeff2,dummy,&isys);
  FFT_1D_INIT(&isign,&nslow,&scalef,dummy,dummy,plan->coeff3,dummy,&isys);

  if (scaled == 0) 
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nmid*nslow);
  }

#endif

#ifdef FFT_FFTW

  plan->plan_fast_forward = 
    fftw_create_plan(nfast,FFTW_FORWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);
  plan->plan_fast_backward = 
    fftw_create_plan(nfast,FFTW_BACKWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);

  if (nmid == nfast) {
    plan->plan_mid_forward = plan->plan_fast_forward;
    plan->plan_mid_backward = plan->plan_fast_backward;
  }
  else {
    plan->plan_mid_forward = 
      fftw_create_plan(nmid,FFTW_FORWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);
    plan->plan_mid_backward = 
      fftw_create_plan(nmid,FFTW_BACKWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);
  }

  if (nslow == nfast) {
    plan->plan_slow_forward = plan->plan_fast_forward;
    plan->plan_slow_backward = plan->plan_fast_backward;
  }
  else if (nslow == nmid) {
    plan->plan_slow_forward = plan->plan_mid_forward;
    plan->plan_slow_backward = plan->plan_mid_backward;
  }
  else {
    plan->plan_slow_forward = 
      fftw_create_plan(nslow,FFTW_FORWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);
    plan->plan_slow_backward = 
      fftw_create_plan(nslow,FFTW_BACKWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);
  }

  if (scaled == 0)
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nmid*nslow);
  }

#endif

  return plan;
}

/* ------------------------------------------------------------------- */
/* Destroy a 3d fft plan */

void fft_3d_destroy_plan_cfunc(struct fft_plan_3d *plan)

{
  if (plan->pre_plan) remap_3d_destroy_plan_cfunc(plan->pre_plan);
  if (plan->mid1_plan) remap_3d_destroy_plan_cfunc(plan->mid1_plan);
  if (plan->mid2_plan) remap_3d_destroy_plan_cfunc(plan->mid2_plan);
  if (plan->post_plan) remap_3d_destroy_plan_cfunc(plan->post_plan);

  if (plan->copy) free(plan->copy);
  if (plan->scratch) free(plan->scratch);

#ifdef FFT_SGI
  free(plan->coeff1);
  free(plan->coeff2);
  free(plan->coeff3);
#endif
#ifdef FFT_INTEL
  free(plan->coeff1);
  free(plan->coeff2);
  free(plan->coeff3);
#endif
#ifdef FFT_T3E
  free(plan->coeff1);
  free(plan->coeff2);
  free(plan->coeff3);
  free(plan->work1);
  free(plan->work2);
  free(plan->work3);
#endif
#ifdef FFT_FFTW
#define fftw_destroy_plan fftw_destroy_planK
  if (plan->plan_slow_forward != plan->plan_fast_forward) {
    fftw_destroy_plan(plan->plan_slow_forward);
    fftw_destroy_plan(plan->plan_slow_backward);
  }
  fftw_destroy_plan(plan->plan_fast_forward);
  fftw_destroy_plan(plan->plan_fast_backward);
#endif

  free(plan);
}

/* parallel remap functions - 1998, 1999

   Steve Plimpton, MS 1111, Dept 9221, Sandia National Labs
   (505) 845-7873
   sjplimp@cs.sandia.gov
*/

/* ------------------------------------------------------------------- */
/* F77 wrapper on remap */

void remap_3d_(double *in, double *out, double *buf,
	       struct remap_plan_3d **plan)

{
  remap_3d_cfunc(in,out,buf,*plan);
}

/* ------------------------------------------------------------------- */
/* F77 wrapper on remap_create_plan */

void remap_3d_create_plan_(
       MPI_Comm *comm,
       int *in_ilo, int *in_ihi, int *in_jlo, int *in_jhi,
       int *in_klo, int *in_khi,
       int *out_ilo, int *out_ihi, int *out_jlo, int *out_jhi,
       int *out_klo, int *out_khi,int *howmany,
       int *nqty, int *permute, int *memory, int *precision,
       struct remap_plan_3d **plan)

{
  int me;

/* convert F77 indices to C */

  *plan = remap_3d_create_plan_cfunc(*comm,
			       *in_ilo-1,*in_ihi-1,*in_jlo-1,*in_jhi-1,
			       *in_klo-1,*in_khi-1,
			       *out_ilo-1,*out_ihi-1,*out_jlo-1,*out_jhi-1,
			       *out_klo-1,*out_khi-1,*howmany,
			       *nqty, *permute, *memory, *precision);
  if (*plan == NULL) {
    MPI_Comm_rank(*comm,&me);
    printf("ERROR: REMAP 3d plan is NULL on proc %d\n",me);
  }
}

/* ------------------------------------------------------------------- */
/* F77 wrapper on remap_destroy_plan */

void remap_3d_destroy_plan_(struct remap_plan_3d **plan)

{
  remap_3d_destroy_plan_cfunc(*plan);
}

/* parallel remap functions - 1998, 1999

   Steve Plimpton, MS 1111, Dept 9221, Sandia National Labs
   (505) 845-7873
   sjplimp@cs.sandia.gov
*/

/* ------------------------------------------------------------------- */
/* Data layout for 3d remaps:

   data set of Nfast x Nmid x Nslow elements is owned by P procs
   each element = nqty contiguous datums
   on input, each proc owns a subsection of the elements
   on output, each proc will own a (presumably different) subsection
   my subsection must not overlap with any other proc's subsection,
     i.e. the union of all proc's input (or output) subsections must
     exactly tile the global Nfast x Nmid x Nslow data set
   when called from C, all subsection indices are 
     C-style from 0 to N-1 where N = Nfast or Nmid or Nslow
   when called from F77, all subsection indices are 
     F77-style from 1 to N where N = Nfast or Nmid or Nslow
   a proc can own 0 elements on input or output
     by specifying hi index < lo index
   on both input and output, data is stored contiguously on a processor
     with a fast-varying, mid-varying, and slow-varying index
*/
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
/* Perform 3d remap */

/* Arguments:

   in           starting address of input data on this proc
   out          starting address of where output data for this proc
                  will be placed (can be same as in)
   buf          extra memory required for remap
                if memory=0 was used in call to remap_3d_create_plan
		  then buf must be big enough to hold output result
		  i.e. nqty * (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) * 
		              (out_khi-out_klo+1)
		if memory=1 was used in call to remap_3d_create_plan
		  then buf is not used, can just be a dummy pointer
   plan         plan returned by previous call to remap_3d_create_plan
*/

void remap_3d_cfunc(double *in, double *out, double *buf,
	      struct remap_plan_3d *plan)

{
  MPI_Status status;
  int i,isend,irecv;
  double *scratch;
  int block;

  if (plan->memory == 0)
    scratch = buf;
  else
    scratch = plan->scratch;

/* post all recvs into scratch space */
  for (irecv = 0; irecv < plan->nrecv; irecv++)
    MPI_Irecv(scratch+plan->recv_bufloc[irecv],
	      plan->recv_size[irecv],
	      MPI_DOUBLE_PRECISION,
	      plan->recv_proc[irecv],0,
	      plan->comm,
	      &plan->request[irecv]);

/* send all messages to other procs */

  for (isend = 0; isend < plan->nsend; isend++) {
    for (block = 0; block < plan->howmany; block++) 
      plan->pack(in+plan->send_offset[isend]+block*plan->in_stride,
		 plan->sendbuf+block*plan->send_stride[isend],
		 &plan->packplan[isend]);
    MPI_Send(plan->sendbuf,plan->send_size[isend],MPI_DOUBLE_PRECISION,
	     plan->send_proc[isend],0,plan->comm);
  }       

/* copy in -> scratch -> out for self data */

  if (plan->self) {
    isend = plan->nsend;
    irecv = plan->nrecv;
    for (block = 0; block < plan->howmany; block++) 
      plan->pack(in+plan->send_offset[isend]+block*plan->in_stride,
		 scratch+plan->recv_bufloc[irecv]
		 +block*plan->send_stride[isend],
		 &plan->packplan[isend]);
    for (block = 0; block < plan->howmany; block++)
      plan->unpack(scratch+plan->recv_bufloc[irecv]
		   +block*plan->recv_stride[irecv],
		   out+plan->recv_offset[irecv]+block*plan->out_stride,
		   &plan->unpackplan[irecv]);
  }

/* unpack all messages from scratch -> out */

  for (i = 0; i < plan->nrecv; i++) {
    MPI_Waitany(plan->nrecv,plan->request,&irecv,&status);
    for (block = 0; block < plan->howmany; block++)
      plan->unpack(scratch+plan->recv_bufloc[irecv]
		   +block*plan->recv_stride[irecv],
		   out+plan->recv_offset[irecv]+block*plan->out_stride,
		   &plan->unpackplan[irecv]);
  }
}

/* ------------------------------------------------------------------- */
/* Create plan for performing a 3d remap */

/* Arguments:

   comm                 MPI communicator for the P procs which own the data
   in_ilo,in_ihi        input bounds of data I own in fast index
   in_jlo,in_jhi        input bounds of data I own in mid index
   in_klo,in_khi        input bounds of data I own in slow index
   out_ilo,out_ihi      output bounds of data I own in fast index
   out_jlo,out_jhi      output bounds of data I own in mid index
   out_klo,out_khi      output bounds of data I own in slow index
   nqty                 # of datums per element
   permute              permutation in storage order of indices on output
                          0 = no permutation
			  1 = permute once = mid->fast, slow->mid, fast->slow
			  2 = permute twice = slow->fast, fast->mid, mid->slow
   memory               user provides buffer memory for remap or system does
                          0 = user provides memory
			  1 = system provides memory
   precision            precision of data
                          1 = single precision (4 bytes per datum)
			  2 = double precision (8 bytes per datum)
*/

struct remap_plan_3d *remap_3d_create_plan_cfunc(
       MPI_Comm comm,
       int in_ilo, int in_ihi, int in_jlo, int in_jhi,
       int in_klo, int in_khi,
       int out_ilo, int out_ihi, int out_jlo, int out_jhi,
       int out_klo, int out_khi,int howmany,
       int nqty, int permute, int memory, int precision)

{
  struct remap_plan_3d *plan;
  MPI_Comm newcomm;
  struct extent_3d *array;
  struct extent_3d in,out,overlap;
  int i,iproc,nsend,nrecv,ibuf,size,me,nprocs;

/* query MPI info */

  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

/* single precision not yet supported */

  if (precision == 1) {
    if (me == 0) printf("Single precision not supported\n");
    return NULL;
  }

/* allocate memory for plan data struct */

  plan = (struct remap_plan_3d *) malloc(sizeof(struct remap_plan_3d));
  if (plan == NULL) return NULL;

/* store parameters in local data structs */

  in.ilo = in_ilo;
  in.ihi = in_ihi;
  in.isize = in.ihi - in.ilo + 1;

  in.jlo = in_jlo;
  in.jhi = in_jhi;
  in.jsize = in.jhi - in.jlo + 1;

  in.klo = in_klo;
  in.khi = in_khi;
  in.ksize = in.khi - in.klo + 1;

  out.ilo = out_ilo;
  out.ihi = out_ihi;
  out.isize = out.ihi - out.ilo + 1;

  out.jlo = out_jlo;
  out.jhi = out_jhi;
  out.jsize = out.jhi - out.jlo + 1;

  out.klo = out_klo;
  out.khi = out_khi;
  out.ksize = out.khi - out.klo + 1;

  plan->in_stride = in.ksize*in.jsize*in.isize*nqty;
  plan->out_stride = out.ksize*out.jsize*out.isize*nqty;
  plan->howmany = howmany;

/* combine output extents across all procs */

  array = (struct extent_3d *) malloc(nprocs*sizeof(struct extent_3d));
  if (array == NULL) return NULL;

  MPI_Allgather(&out,sizeof(struct extent_3d),MPI_BYTE,
		array,sizeof(struct extent_3d),MPI_BYTE,comm);

/* count send collides, including self */

  nsend = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    nsend += remap_3d_collide(&in,&array[iproc],&overlap);
  }

/* malloc space for send info */
plan->send_mallocd = 0; // Initialize as false

  if (nsend) {
    if (precision == 1)
      plan->pack = NULL;
    else
      plan->pack = pack_3d;

    plan->send_mallocd = 1;

    plan->send_offset = (int *) malloc(nsend*sizeof(int));
    plan->send_size = (int *) malloc(nsend*sizeof(int));
    plan->send_stride = (int *) malloc(nsend*sizeof(int));
    plan->send_proc = (int *) malloc(nsend*sizeof(int));
    plan->packplan = (struct pack_plan_3d *) 
      malloc(nsend*sizeof(struct pack_plan_3d));

    if (plan->send_offset == NULL || plan->send_size == NULL || 
	plan->send_stride == NULL || 
	plan->send_proc == NULL || plan->packplan == NULL) return NULL;
  }

/* store send info, with self as last entry */

  nsend = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (remap_3d_collide(&in,&array[iproc],&overlap)) {
      plan->send_proc[nsend] = iproc;
      plan->send_offset[nsend] = nqty * 
	((overlap.klo-in.klo)*in.jsize*in.isize + 
	((overlap.jlo-in.jlo)*in.isize + overlap.ilo-in.ilo));
      plan->packplan[nsend].nfast = nqty*overlap.isize;
      plan->packplan[nsend].nmid = overlap.jsize;
      plan->packplan[nsend].nslow = overlap.ksize;
      plan->packplan[nsend].nstride_line = nqty*in.isize;
      plan->packplan[nsend].nstride_plane = nqty*in.jsize*in.isize;
      plan->packplan[nsend].nqty = nqty;
      plan->send_stride[nsend]=nqty*overlap.isize*overlap.jsize*overlap.ksize;
      plan->send_size[nsend] = plan->send_stride[nsend]*plan->howmany;
      nsend++;
    }
  }

/* plan->nsend = # of sends not including self */

  if (nsend && plan->send_proc[nsend-1] == me)
    plan->nsend = nsend - 1;
  else
    plan->nsend = nsend;

/* combine input extents across all procs */

  MPI_Allgather(&in,sizeof(struct extent_3d),MPI_BYTE,
		array,sizeof(struct extent_3d),MPI_BYTE,comm);

/* count recv collides, including self */

  nrecv = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    nrecv += remap_3d_collide(&out,&array[iproc],&overlap);
  }
  
/* malloc space for recv info */
  plan->recv_mallocd = 0;
  if (nrecv) {
    if (precision == 1) {
      if (permute == 0)
	plan->unpack = NULL;
      else if (permute == 1) {
	if (nqty == 1)
	  plan->unpack = NULL;
	else if (nqty == 2)
	  plan->unpack = NULL;
	else
	  plan->unpack = NULL;
      }
      else if (permute == 2) {
	if (nqty == 1)
	  plan->unpack = NULL;
	else if (nqty == 2)
	  plan->unpack = NULL;
	else
	  plan->unpack = NULL;
      }
    }
    else if (precision == 2) {
      if (permute == 0)
	plan->unpack = unpack_3d;
      else if (permute == 1) {
	if (nqty == 1)
	  plan->unpack = unpack_3d_permute1_1;
	else if (nqty == 2)
	  plan->unpack = unpack_3d_permute1_2;
	else
	  plan->unpack = unpack_3d_permute1_n;
      }
      else if (permute == 2) {
	if (nqty == 1)
	  plan->unpack = unpack_3d_permute2_1;
	else if (nqty == 2)
	  plan->unpack = unpack_3d_permute2_2;
	else
	  plan->unpack = unpack_3d_permute2_n;
      }
    }

    plan->recv_mallocd = 1;

    plan->recv_offset = (int *) malloc(nrecv*sizeof(int));
    plan->recv_size = (int *) malloc(nrecv*sizeof(int));
    plan->recv_stride = (int *) malloc(nrecv*sizeof(int));
    plan->recv_proc = (int *) malloc(nrecv*sizeof(int));
    plan->recv_bufloc = (int *) malloc(nrecv*sizeof(int));
    plan->request = (MPI_Request *) malloc(nrecv*sizeof(MPI_Request));
    plan->unpackplan = (struct pack_plan_3d *) 
      malloc(nrecv*sizeof(struct pack_plan_3d));

    if (plan->recv_offset == NULL || plan->recv_size == NULL || 
	plan->recv_stride == NULL || 
	plan->recv_proc == NULL || plan->recv_bufloc == NULL ||
	plan->request == NULL || plan->unpackplan == NULL) return NULL;
  }

/* store recv info, with self as last entry */

  ibuf = 0;
  nrecv = 0;
  iproc = me;

  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (remap_3d_collide(&out,&array[iproc],&overlap)) {
      plan->recv_proc[nrecv] = iproc;
      plan->recv_bufloc[nrecv] = ibuf;

      if (permute == 0) {
	plan->recv_offset[nrecv] = nqty *
	  ((overlap.klo-out.klo)*out.jsize*out.isize +
	   (overlap.jlo-out.jlo)*out.isize + (overlap.ilo-out.ilo));
	plan->unpackplan[nrecv].nfast = nqty*overlap.isize;
	plan->unpackplan[nrecv].nmid = overlap.jsize;
	plan->unpackplan[nrecv].nslow = overlap.ksize;
	plan->unpackplan[nrecv].nstride_line = nqty*out.isize;
	plan->unpackplan[nrecv].nstride_plane = nqty*out.jsize*out.isize;
	plan->unpackplan[nrecv].nqty = nqty;
      }
      else if (permute == 1) {
	plan->recv_offset[nrecv] = nqty *
	  ((overlap.ilo-out.ilo)*out.ksize*out.jsize +
	   (overlap.klo-out.klo)*out.jsize + (overlap.jlo-out.jlo));
	plan->unpackplan[nrecv].nfast = overlap.isize;
	plan->unpackplan[nrecv].nmid = overlap.jsize;
	plan->unpackplan[nrecv].nslow = overlap.ksize;
	plan->unpackplan[nrecv].nstride_line = nqty*out.jsize;
	plan->unpackplan[nrecv].nstride_plane = nqty*out.ksize*out.jsize;
	plan->unpackplan[nrecv].nqty = nqty;
      }
      else {
	plan->recv_offset[nrecv] = nqty *
	  ((overlap.jlo-out.jlo)*out.isize*out.ksize +
	   (overlap.ilo-out.ilo)*out.ksize + (overlap.klo-out.klo));
	plan->unpackplan[nrecv].nfast = overlap.isize;
	plan->unpackplan[nrecv].nmid = overlap.jsize;
	plan->unpackplan[nrecv].nslow = overlap.ksize;
	plan->unpackplan[nrecv].nstride_line = nqty*out.ksize;
	plan->unpackplan[nrecv].nstride_plane = nqty*out.isize*out.ksize;
	plan->unpackplan[nrecv].nqty = nqty;
      }

      plan->recv_stride[nrecv]=nqty*overlap.isize*overlap.jsize*overlap.ksize;
      plan->recv_size[nrecv] = plan->recv_stride[nrecv]*plan->howmany;
      ibuf += plan->recv_size[nrecv];
      nrecv++;
    }
  }

/* plan->nrecv = # of recvs not including self */

  if (nrecv && plan->recv_proc[nrecv-1] == me)
    plan->nrecv = nrecv - 1;
  else
    plan->nrecv = nrecv;

/* init remaining fields in remap plan */

  plan->memory = memory;

  if (nrecv == plan->nrecv)
    plan->self = 0;
  else
    plan->self = 1;

/* free locally malloced space */

  free(array);

/* find biggest send message (not including self) and malloc space for it */

  size = 0;
  for (nsend = 0; nsend < plan->nsend; nsend++)
    size = MAX(size,plan->send_size[nsend]);

  if (size) {
    if (precision == 1)
      plan->sendbuf = NULL;
    else
      plan->sendbuf = (double *) malloc(size*sizeof(double));
    if (plan->sendbuf == NULL) return NULL;
  }
  else
    plan->sendbuf = NULL;

/* if requested, allocate internal scratch space for recvs,
   only need it if I will receive any data (including self) */

  if (memory == 1) {
    if (nrecv > 0) {
      if (precision == 1)
	plan->scratch = NULL;
      else
	plan->scratch =
	  (double *) malloc(plan->howmany*plan->out_stride*sizeof(double));
      if (plan->scratch == NULL) return NULL;
    }
  }

/* create new MPI communicator for remap */

  MPI_Comm_dup(comm,&plan->comm);

/* return pointer to plan */

  return plan;
}

/* ---------------------------------------------------------------------
   Destroy a 3d remap plan
--------------------------------------------------------------------- */

void remap_3d_destroy_plan_cfunc(struct remap_plan_3d *plan)
{
  /* free MPI communicator */

  MPI_Comm_free(&plan->comm);

  /* free internal arrays */

//  if (plan->nsend) {
  if(plan->send_mallocd) {
    free(plan->send_offset);
    free(plan->send_size);
    free(plan->send_stride);
    free(plan->send_proc);
    free(plan->packplan);

    if (plan->sendbuf != NULL)
      free(plan->sendbuf);

  }

//  if (plan->nrecv) {
  if(plan->recv_mallocd) {
    free(plan->recv_offset);
    free(plan->recv_size);
    free(plan->recv_stride);
    free(plan->recv_proc);
    free(plan->recv_bufloc);
    free(plan->request);
    free(plan->unpackplan);
    if (plan->memory) free(plan->scratch);
  }

  /* free plan itself */

  free(plan);
}

/* ----------------------------------------------------------------------
   collide 2 sets of indices to determine overlap
   compare bounds of block1 with block2 to see if they overlap
   return 1 if they do and put bounds of overlapping section in overlap
   return 0 if they do not overlap
------------------------------------------------------------------------- */

int remap_3d_collide(struct extent_3d *block1, struct extent_3d *block2, struct extent_3d *overlap)
{
  overlap->ilo = MAX(block1->ilo,block2->ilo);
  overlap->ihi = MIN(block1->ihi,block2->ihi);
  overlap->jlo = MAX(block1->jlo,block2->jlo);
  overlap->jhi = MIN(block1->jhi,block2->jhi);
  overlap->klo = MAX(block1->klo,block2->klo);
  overlap->khi = MIN(block1->khi,block2->khi);

  if (overlap->ilo > overlap->ihi ||
      overlap->jlo > overlap->jhi ||
      overlap->klo > overlap->khi) return 0;

  overlap->isize = overlap->ihi - overlap->ilo + 1;
  overlap->jsize = overlap->jhi - overlap->jlo + 1;
  overlap->ksize = overlap->khi - overlap->klo + 1;

  return 1;
}

/* ---------------------------------------------------------------------
   Pack and unpack functions:

   pack routines copy strided values from data into contiguous locs in buf
   unpack routines copy contiguous values from buf into strided locs in data
   different versions of unpack depending on permutation and # of values/element
      PACK_ARRAY    methods work via array indices (default)
      PACK_POINTER  methods work via pointers
      PACK_MEMCPY   methods work via pointers and memcpy function
                    no memcpy version of unpack_permute methods, just use POINTER versions
--------------------------------------------------------------------- */

#ifdef PACK_ARRAY

/* ---------------------------------------------------------------------
   pack from data -> buf
--------------------------------------------------------------------- */

void pack_3d(double *data, double *buf, struct pack_plan_3d *plan)
{
  register int in,out,fast,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  in = 0;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_plane;
    for (mid = 0; mid < nmid; mid++) {
      out = plane + mid*nstride_line;
      for (fast = 0; fast < nfast; fast++) {
        buf[in++] = data[out++];
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data
--------------------------------------------------------------------- */

void unpack_3d(double *buf, double *data, struct pack_plan_3d *plan)
{
  register int in,out,fast,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_plane;
    for (mid = 0; mid < nmid; mid++) {
      in = plane + mid*nstride_line;
      for (fast = 0; fast < nfast; fast++) {
        data[in++] = buf[out++];
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, one axis permutation, 1 value/element
--------------------------------------------------------------------- */

void unpack_3d_permute1_1(double *buf, double *data, struct pack_plan_3d *plan)
{
  register int in,out,fast,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      in = plane + mid;
      for (fast = 0; fast < nfast; fast++, in += nstride_plane) {
        data[in] = buf[out++];
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, one axis permutation, 2 values/element
--------------------------------------------------------------------- */

void unpack_3d_permute1_2(double *buf, double *data, struct pack_plan_3d *plan)
{
  register int in,out,fast,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      in = plane + 2*mid;
      for (fast = 0; fast < nfast; fast++, in += nstride_plane) {
        data[in] = buf[out++];
        data[in+1] = buf[out++];
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, one axis permutation, nqty values/element
--------------------------------------------------------------------- */

void unpack_3d_permute1_n(double *buf, double *data, struct pack_plan_3d *plan)
{
  register int in,out,iqty,instart,fast,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane,nqty;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      instart = plane + nqty*mid;
      for (fast = 0; fast < nfast; fast++, instart += nstride_plane) {
        in = instart;
        for (iqty = 0; iqty < nqty; iqty++) {
          data[in++] = buf[out++];
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, two axis permutation, 1 value/element
--------------------------------------------------------------------- */

void unpack_3d_permute2_1(double *buf, double *data, struct pack_plan_3d *plan)
{
  register int in,out,fast,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      in = slow + mid*nstride_plane;
      for (fast = 0; fast < nfast; fast++, in += nstride_line) {
        data[in] = buf[out++];
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, two axis permutation, 2 values/element
--------------------------------------------------------------------- */

void unpack_3d_permute2_2(double *buf, double *data, struct pack_plan_3d *plan)
{
  register int in,out,fast,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      in = 2*slow + mid*nstride_plane;
      for (fast = 0; fast < nfast; fast++, in += nstride_line) {
        data[in] = buf[out++];
        data[in+1] = buf[out++];
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, two axis permutation, nqty values/element
--------------------------------------------------------------------- */

void unpack_3d_permute2_n(double *buf, double *data, struct pack_plan_3d *plan)
{
  register int in,out,iqty,instart,fast,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,nqty;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      instart = nqty*slow + mid*nstride_plane;
      for (fast = 0; fast < nfast; fast++, instart += nstride_line) {
        in = instart;
        for (iqty = 0; iqty < nqty; iqty++) {
          data[in++] = buf[out++];
        }
      }
    }
  }
}

#endif
#ifdef PACK_POINTER

/* --------------------------------------------------------------------
   pack from data -> buf
--------------------------------------------------------------------- */

void pack_3d(double *data, double *buf, struct pack_plan_3d *plan)
{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  in = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_plane;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane+mid*nstride_line]);
      end = begin + nfast;
      for (out = begin; out < end; out++) {
        *(in++) = *out;
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data
--------------------------------------------------------------------- */

void unpack_3d(double *buf, double *data, struct pack_plan_3d *plan)
{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_plane;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane+mid*nstride_line]);
      end = begin + nfast;
      for (in = begin; in < end; in++) {
        *in = *(out++);
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, one axis permutation, 1 value/element
--------------------------------------------------------------------- */

void unpack_3d_permute1_1(double *buf, double *data, struct pack_plan_3d *plan)
{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane+mid]);
      end = begin + nfast*nstride_plane;
      for (in = begin; in < end; in += nstride_plane) {
        *in = *(out++);
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, one axis permutation, 2 values/element
--------------------------------------------------------------------- */

void unpack_3d_permute1_2(double *buf, double *data, struct pack_plan_3d *plan)
{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane+2*mid]);
      end = begin + nfast*nstride_plane;
      for (in = begin; in < end; in += nstride_plane) {
        *in = *(out++);
        *(in+1) = *(out++);
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, one axis permutation, nqty values/element
--------------------------------------------------------------------- */

void unpack_3d_permute1_n(double *buf, double *data, struct pack_plan_3d *plan)
{
  register double *in,*out,*instart,*begin,*end;
  register int iqty,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane,nqty;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane+nqty*mid]);
      end = begin + nfast*nstride_plane;
      for (instart = begin; instart < end; instart += nstride_plane) {
        in = instart;
        for (iqty = 0; iqty < nqty; iqty++) {
          *(in++) = *(out++);
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, two axis permutation, 1 value/element
--------------------------------------------------------------------- */

void unpack_3d_permute2_1(double *buf, double *data, struct pack_plan_3d *plan)
{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[slow+mid*nstride_plane]);
      end = begin + nfast*nstride_line;
      for (in = begin; in < end; in += nstride_line) {
        *in = *(out++);
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, two axis permutation, 2 values/element
--------------------------------------------------------------------- */

void unpack_3d_permute2_2(double *buf, double *data, struct pack_plan_3d *plan)
{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[2*slow+mid*nstride_plane]);
      end = begin + nfast*nstride_line;
      for (in = begin; in < end; in += nstride_line) {
        *in = *(out++);
        *(in+1) = *(out++);
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, two axis permutation, nqty values/element
--------------------------------------------------------------------- */

void unpack_3d_permute2_n(double *buf, double *data, struct pack_plan_3d *plan)
{
  register double *in,*out,*instart,*begin,*end;
  register int iqty,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,nqty;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[nqty*slow+mid*nstride_plane]);
      end = begin + nfast*nstride_line;
      for (instart = begin; instart < end; instart += nstride_line) {
        in = instart;
        for (iqty = 0; iqty < nqty; iqty++) {
          *(in++) = *(out++);
        }
      }
    }
  }
}

#endif
#ifdef PACK_MEMCPY

/* ---------------------------------------------------------------------
   pack from data -> buf
--------------------------------------------------------------------- */

void pack_3d(double *data, double *buf, struct pack_plan_3d *plan)
{
  register double *in,*out;
  register int mid,slow,size;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane,upto;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  size = nfast*sizeof(double);
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_plane;
    upto = slow*nmid*nfast;
    for (mid = 0; mid < nmid; mid++) {
      in = &(buf[upto+mid*nfast]);
      out = &(data[plane+mid*nstride_line]);
      memcpy(in,out,size);
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data
--------------------------------------------------------------------- */

void unpack_3d(double *buf, double *data, struct pack_plan_3d *plan)
{
  register double *in,*out;
  register int mid,slow,size;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane,upto;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  size = nfast*sizeof(double);
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_plane;
    upto = slow*nmid*nfast;
    for (mid = 0; mid < nmid; mid++) {
      in = &(data[plane+mid*nstride_line]);
      out = &(buf[upto+mid*nfast]);
      memcpy(in,out,size);
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, one axis permutation, 1 value/element
--------------------------------------------------------------------- */

void unpack_3d_permute1_1(double *buf, double *data, struct pack_plan_3d *plan)
{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane+mid]);
      end = begin + nfast*nstride_plane;
      for (in = begin; in < end; in += nstride_plane) {
        *in = *(out++);
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, one axis permutation, 2 values/element
--------------------------------------------------------------------- */

void unpack_3d_permute1_2(double *buf, double *data, struct pack_plan_3d *plan)
{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane+2*mid]);
      end = begin + nfast*nstride_plane;
      for (in = begin; in < end; in += nstride_plane) {
        *in = *(out++);
        *(in+1) = *(out++);
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, one axis permutation, nqty values/element
--------------------------------------------------------------------- */

void unpack_3d_permute1_n(double *buf, double *data, struct pack_plan_3d *plan)
{
  register double *in,*out,*instart,*begin,*end;
  register int iqty,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,plane,nqty;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    plane = slow*nstride_line;
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[plane+nqty*mid]);
      end = begin + nfast*nstride_plane;
      for (instart = begin; instart < end; instart += nstride_plane) {
        in = instart;
        for (iqty = 0; iqty < nqty; iqty++) {
          *(in++) = *(out++);
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, two axis permutation, 1 value/element
--------------------------------------------------------------------- */

void unpack_3d_permute2_1(double *buf, double *data, struct pack_plan_3d *plan)
{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[slow+mid*nstride_plane]);
      end = begin + nfast*nstride_line;
      for (in = begin; in < end; in += nstride_line) {
        *in = *(out++);
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, two axis permutation, 2 values/element
--------------------------------------------------------------------- */

void unpack_3d_permute2_2(double *buf, double *data, struct pack_plan_3d *plan)
{
  register double *in,*out,*begin,*end;
  register int mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[2*slow+mid*nstride_plane]);
      end = begin + nfast*nstride_line;
      for (in = begin; in < end; in += nstride_line) {
        *in = *(out++);
        *(in+1) = *(out++);
      }
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, two axis permutation, nqty values/element
--------------------------------------------------------------------- */

void unpack_3d_permute2_n(double *buf, double *data, struct pack_plan_3d *plan)
{
  register double *in,*out,*instart,*begin,*end;
  register int iqty,mid,slow;
  register int nfast,nmid,nslow,nstride_line,nstride_plane,nqty;

  nfast = plan->nfast;
  nmid = plan->nmid;
  nslow = plan->nslow;
  nstride_line = plan->nstride_line;
  nstride_plane = plan->nstride_plane;
  nqty = plan->nqty;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    for (mid = 0; mid < nmid; mid++) {
      begin = &(data[nqty*slow+mid*nstride_plane]);
      end = begin + nfast*nstride_line;
      for (instart = begin; instart < end; instart += nstride_line) {
        in = instart;
        for (iqty = 0; iqty < nqty; iqty++) {
          *(in++) = *(out++);
        }
      }
    }
  }
}

#endif

/* ---------------------------------------------------------------------
   recursively divide n into small factors, return them in list
--------------------------------------------------------------------- */

void factor(int n, int *num, int *list)
{
  if (n == 1) {
    return;
  }
  else if (n % 2 == 0) {
    *list = 2;
    (*num)++;
    factor(n/2,num,list+1);
  }
  else if (n % 3 == 0) {
    *list = 3;
    (*num)++;
    factor(n/3,num,list+1);
  }
  else if (n % 5 == 0) {
    *list = 5;
    (*num)++;
    factor(n/5,num,list+1);
  }
  else if (n % 7 == 0) {
    *list = 7;
    (*num)++;
    factor(n/7,num,list+1);
  }
  else if (n % 11 == 0) {
    *list = 11;
    (*num)++;
    factor(n/11,num,list+1);
  }
  else if (n % 13 == 0) {
    *list = 13;
    (*num)++;
    factor(n/13,num,list+1);
  }
  else {
    *list = n;
    (*num)++;
    return;
  }
}

/* ---------------------------------------------------------------------
   divide n into 2 factors of as equal size as possible
--------------------------------------------------------------------- */

void bifactor(int n, int *factor1, int *factor2)
{
  int n1,n2,facmax;

  facmax = sqrt((double) n);

  for (n1 = facmax; n1 > 0; n1--) {
    n2 = n/n1;
    if (n1*n2 == n) {
      *factor1 = n1;
      *factor2 = n2;
      return;
    }
  }
}

/* ------------------------------------------------------------------ */
