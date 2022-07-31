/* ---------------------------------------------------------------------------------------------------------------------------------
   Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.
   See the README file in the top-level directory.

   Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).
   This software is distributed uner the modified Berkeley Software Distribution (BSD) License.
   Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.
--------------------------------------------------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------
   2d data remapping functions

   Original author: Steve Plimpton, Sandia National Laboratories
                    Parallel FFT Package - 1998, 1999

   Data layout for 2d remaps:

   data set of Nfast x Nslow elements is owned by P procs
   each element = nqty contiguous datums
   on input, each proc owns a subsection of the elements
   on output, each proc will own a (presumably different) subsection
   my subsection must not overlap with any other proc's subsection,
      i.e. the union of all proc's input (or output) subsections must
      exactly tile the global Nfast x Nslow data set
   when called from C, all subsection indices are
      C-style from 0 to N-1 where N = Nfast or Nslow
   when called from F77, all subsection indices are
      F77-style from 1 to N where N = Nfast or Nslow
   a proc can own 0 elements on input or output
      by specifying hi index < lo index
   on both input and output, data is stored contiguously on a processor
      with a fast-varying and slow-varying index
--------------------------------------------------------------------- */

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

/* ------------------------------------------------------------------ */

// Timing info

double recv_time, send_time, waitany_time;

// Collision between 2 regions

struct extent_2d {
  int ilo,ihi,isize;
  int jlo,jhi,jsize;
};

// Loop counters for performing a pack/unpack

struct pack_plan_2d {
  int nfast;                         /* # of elements in fast index */
  int nslow;                         /* # of elements in slow index */
  int nstride;                       /* stride between succesive slow indices */
  int nqty;                          /* # of values/element */
};

// Function prototypes for the pack/unpack routines

void pack_2d(double *, double *, struct pack_plan_2d *);
void unpack_2d(double *, double *, struct pack_plan_2d *);
void unpack_2d_permute_1(double *, double *, struct pack_plan_2d *);
void unpack_2d_permute_2(double *, double *, struct pack_plan_2d *);
void unpack_2d_permute_n(double *, double *, struct pack_plan_2d *);

// Details of how to perform a 2d remap

struct remap_plan_2d {
  double *sendbuf;                   /* buffer for MPI sends */
  double *scratch;                   /* scratch buffer for MPI recvs */
  void (*pack)();                    /* which pack function to use */
  void (*unpack)();                  /* which unpack function to use */
  int *send_offset;                  /* extraction loc for each send */
  int *send_size;                    /* size of each send message */
  int *send_proc;                    /* proc to send each message to */
  struct pack_plan_2d *packplan;     /* pack plan for each send message */
  int *recv_offset;                  /* insertion loc for each recv */
  int *recv_size;                    /* size of each recv message */
  int *recv_proc;                    /* proc to recv each message from */
  int *recv_bufloc;                  /* offset in scratch buf for each recv */
  MPI_Request *request;              /* MPI request for each posted recv */
  struct pack_plan_2d *unpackplan;   /* unpack plan for each recv message */
  int nrecv;                         /* # of recvs from other procs */
  int nsend;                         /* # of sends to other procs */
  int self;                          /* whether I send/recv with myself */
  int memory;                        /* user provides scratch space or not */
  MPI_Comm comm;                     /* group of procs performing remap */
};

// Function prototypes for the remap routines

struct remap_plan_2d *remap_2d_create_plan_cfunc(MPI_Comm,
  int, int, int, int, int, int, int, int, int, int, int, int);
void remap_2d_cfunc(double *, double *, double *, struct remap_plan_2d *);
void remap_2d_destroy_plan_cfunc(struct remap_plan_2d *);
int remap_2d_collide(struct extent_2d *, struct extent_2d *, struct extent_2d *);

/* ---------------------------------------------------------------------
   Fortran wrapper on remap_2d_create_plan_cfunc
--------------------------------------------------------------------- */

void FTN(remap_2d_create_plan)(FOR_MPI_COMM *for_comm,
         int *in_ilo, int *in_ihi, int *in_jlo, int *in_jhi,
         int *out_ilo, int *out_ihi, int *out_jlo, int *out_jhi,
         int *nqty, int *permute, int *memory, int *precision,
         struct remap_plan_2d **plan)
{
  int me;
  MPI_Comm comm;

#if (MPI_VERSION == 1)
  comm = *for_comm;
#else
  comm = MPI_Comm_f2c(*for_comm);
#endif

  // Create plan for performing a 2d remap
  // Note: Convert Fortran indices to C

  *plan = remap_2d_create_plan_cfunc(comm,
          *in_ilo-1,*in_ihi-1,*in_jlo-1,*in_jhi-1,
          *out_ilo-1,*out_ihi-1,*out_jlo-1,*out_jhi-1,
          *nqty,*permute,*memory,*precision);

  if (*plan == NULL) {
    MPI_Comm_rank(comm,&me);
    printf("ERROR: REMAP 2d plan is NULL on proc %d\n",me);
  }
}

/* ---------------------------------------------------------------------
   Fortran wrapper on remap_2d_cfunc
--------------------------------------------------------------------- */

void FTN(remap_2d)(double *in, double *out, double *buf, struct remap_plan_2d **plan)
{
  remap_2d_cfunc(in,out,buf,*plan);
}

/* ---------------------------------------------------------------------
   Fortran wrapper on remap_2d_destroy_plan_cfunc
--------------------------------------------------------------------- */

void FTN(remap_2d_destroy_plan)(struct remap_plan_2d **plan)
{
  remap_2d_destroy_plan_cfunc(*plan);
}

/* ---------------------------------------------------------------------
   Perform a 2d remap

   arguments:
   in           starting address of input data on this proc
   out          starting address of where output data for this proc
                   will be placed (can be same as in)
   buf          extra memory required for remap
                if memory=0 was used in call to remap_2d_create_plan
                   then buf must be big enough to hold output result
                   i.e. nqty * (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1)
                if memory=1 was used in call to remap_2d_create_plan
                   then buf is not used, can just be a dummy pointer
   plan         plan returned by previous call to remap_2d_create_plan
--------------------------------------------------------------------- */

void remap_2d_cfunc(double *in, double *out, double *buf, struct remap_plan_2d *plan)
{
  MPI_Status status;
  int i,isend,irecv;
  double *scratch;
  double st, et;

  if (plan->memory == 0) {
    scratch = buf;
  }
  else {
    scratch = plan->scratch;
  }

  // post all recvs into scratch space

  st = MPI_Wtime();
  for (irecv = 0; irecv < plan->nrecv; irecv++) {
    MPI_Irecv(&scratch[plan->recv_bufloc[irecv]],plan->recv_size[irecv],
              MPI_DOUBLE_PRECISION,plan->recv_proc[irecv],0,plan->comm,&plan->request[irecv]);
  }
  recv_time = recv_time + MPI_Wtime() - st;

  // send all messages to other procs

  st = MPI_Wtime();
  for (isend = 0; isend < plan->nsend; isend++) {
    plan->pack(&in[plan->send_offset[isend]],plan->sendbuf,&plan->packplan[isend]);
    MPI_Send(plan->sendbuf,plan->send_size[isend],MPI_DOUBLE_PRECISION,plan->send_proc[isend],0,plan->comm);
  }
  send_time = send_time + MPI_Wtime() - st;

  // copy in -> scratch -> out for self data

  if (plan->self) {
    isend = plan->nsend;
    irecv = plan->nrecv;
    plan->pack(&in[plan->send_offset[isend]],&scratch[plan->recv_bufloc[irecv]],&plan->packplan[isend]);
    plan->unpack(&scratch[plan->recv_bufloc[irecv]],&out[plan->recv_offset[irecv]],&plan->unpackplan[irecv]);
  }

  // unpack all messages from scratch -> out

  st = MPI_Wtime();
  for (i = 0; i < plan->nrecv; i++) {
    MPI_Waitany(plan->nrecv,plan->request,&irecv,&status);
    plan->unpack(&scratch[plan->recv_bufloc[irecv]],&out[plan->recv_offset[irecv]],&plan->unpackplan[irecv]);
  }

  waitany_time = waitany_time + MPI_Wtime() - st;
}

/* ----------------------------------------------------------------------
   Create a plan for performing a 2d remap

   arguments:
   comm              MPI communicator for the P procs which own the data
   in_ilo,in_ihi     input bounds of data I own in fast index
   in_jlo,in_jhi     input bounds of data I own in slow index
   out_ilo,out_ihi   output bounds of data I own in fast index
   out_jlo,out_jhi   output bounds of data I own in slow index
   nqty              # of datums per element
   permute           permutation in storage order of indices on output
                        0 = no permutation
                        1 = permute = slow->fast, fast->slow
   memory            user provides buffer memory for remap or system does
                        0 = user provides memory
                        1 = system provides memory internally
   precision         precision of data
                        1 = single precision (4 bytes per datum)
                        2 = double precision (8 bytes per datum)
------------------------------------------------------------------------- */

struct remap_plan_2d *remap_2d_create_plan_cfunc(MPI_Comm comm,
       int in_ilo, int in_ihi, int in_jlo, int in_jhi,
       int out_ilo, int out_ihi, int out_jlo, int out_jhi,
       int nqty, int permute, int memory, int precision)
{
  MPI_Comm newcomm;
  struct remap_plan_2d *plan;
  struct extent_2d *array;
  struct extent_2d in,out,overlap;
  int i,iproc,nsend,nrecv,ibuf,size,me,nprocs;

  recv_time = 0;  send_time = 0; waitany_time = 0;

  // query MPI info

  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

  // single precision not yet supported

  if (precision == 1) {
    if (me == 0) printf("Single precision not supported\n");
    return NULL;
  }

  // allocate memory for plan data struct

  plan = (struct remap_plan_2d *) malloc(sizeof(struct remap_plan_2d));
  if (plan == NULL) return NULL;

  // store parameters in local data structs

  in.ilo = in_ilo;
  in.ihi = in_ihi;
  in.isize = in.ihi - in.ilo + 1;

  in.jlo = in_jlo;
  in.jhi = in_jhi;
  in.jsize = in.jhi - in.jlo + 1;

  out.ilo = out_ilo;
  out.ihi = out_ihi;
  out.isize = out.ihi - out.ilo + 1;

  out.jlo = out_jlo;
  out.jhi = out_jhi;
  out.jsize = out.jhi - out.jlo + 1;

  // combine output extents across all procs

  array = (struct extent_2d *) malloc(nprocs*sizeof(struct extent_2d));
  if (array == NULL) return NULL;

  MPI_Allgather(&out,sizeof(struct extent_2d),MPI_BYTE,array,sizeof(struct extent_2d),MPI_BYTE,comm);

  // count send collides, including self

  nsend = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    nsend += remap_2d_collide(&in,&array[iproc],&overlap);
  }

  // malloc space for send info

  if (nsend) {
    if (precision == 1) {
      plan->pack = NULL;
    }
    else {
      plan->pack = pack_2d;
    }

    plan->send_offset = (int *) malloc(nsend*sizeof(int));
    plan->send_size = (int *) malloc(nsend*sizeof(int));
    plan->send_proc = (int *) malloc(nsend*sizeof(int));
    plan->packplan = (struct pack_plan_2d *) malloc(nsend*sizeof(struct pack_plan_2d));

    if (plan->send_offset == NULL || plan->send_size == NULL ||
        plan->send_proc == NULL || plan->packplan == NULL) return NULL;
  }

  // store send info, with self as last entry

  nsend = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (remap_2d_collide(&in,&array[iproc],&overlap)) {
      plan->send_proc[nsend] = iproc;
      plan->send_offset[nsend] = nqty*((overlap.jlo-in.jlo)*in.isize + (overlap.ilo-in.ilo));
      plan->packplan[nsend].nfast = nqty*overlap.isize;
      plan->packplan[nsend].nslow = overlap.jsize;
      plan->packplan[nsend].nstride = nqty*in.isize;
      plan->packplan[nsend].nqty = nqty;
      plan->send_size[nsend] = nqty*overlap.isize*overlap.jsize;
      nsend++;
    }
  }

  // plan->nsend = # of sends not including self

  if (nsend && plan->send_proc[nsend-1] == me) {
    plan->nsend = nsend - 1;
  }
  else {
    plan->nsend = nsend;
  }

  // combine input extents across all procs

  MPI_Allgather(&in,sizeof(struct extent_2d),MPI_BYTE,array,sizeof(struct extent_2d),MPI_BYTE,comm);

  // count recv collides, including self

  nrecv = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    nrecv += remap_2d_collide(&out,&array[iproc],&overlap);
  }

  // malloc space for recv info

  if (nrecv) {
    if (precision == 1) {
      if (permute == 0)
	plan->unpack = NULL;
      else if (nqty == 1)
	plan->unpack = NULL;
      else if (nqty == 2)
	plan->unpack = NULL;
      else
	plan->unpack = NULL;
    }
    else if (precision == 2) {
      if (permute == 0)
	plan->unpack = unpack_2d;
      else if (nqty == 1)
	plan->unpack = unpack_2d_permute_1;
      else if (nqty == 2)
	plan->unpack = unpack_2d_permute_2;
      else
	plan->unpack = unpack_2d_permute_n;
    }

    plan->recv_offset = (int *) malloc(nrecv*sizeof(int));
    plan->recv_size = (int *) malloc(nrecv*sizeof(int));
    plan->recv_proc = (int *) malloc(nrecv*sizeof(int));
    plan->recv_bufloc = (int *) malloc(nrecv*sizeof(int));
    plan->request = (MPI_Request *) malloc(nrecv*sizeof(MPI_Request));
    plan->unpackplan = (struct pack_plan_2d *) malloc(nrecv*sizeof(struct pack_plan_2d));

    if (plan->recv_offset == NULL || plan->recv_size == NULL ||
	plan->recv_proc == NULL || plan->recv_bufloc == NULL ||
	plan->request == NULL || plan->unpackplan == NULL) return NULL;
  }

  // store recv info, with self as last entry

  ibuf = 0;
  nrecv = 0;
  iproc = me;

  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (remap_2d_collide(&out,&array[iproc],&overlap)) {
      plan->recv_proc[nrecv] = iproc;
      plan->recv_bufloc[nrecv] = ibuf;

      if (permute == 0) {
	plan->recv_offset[nrecv] = nqty*((overlap.jlo-out.jlo)*out.isize + (overlap.ilo-out.ilo));
	plan->unpackplan[nrecv].nfast = nqty*overlap.isize;
	plan->unpackplan[nrecv].nslow = overlap.jsize;
	plan->unpackplan[nrecv].nstride = nqty*out.isize;
	plan->unpackplan[nrecv].nqty = nqty;
      }
      else {
	plan->recv_offset[nrecv] = nqty * ((overlap.ilo-out.ilo)*out.jsize + (overlap.jlo-out.jlo));
	plan->unpackplan[nrecv].nfast = overlap.isize;
	plan->unpackplan[nrecv].nslow = overlap.jsize;
	plan->unpackplan[nrecv].nstride = nqty*out.jsize;
	plan->unpackplan[nrecv].nqty = nqty;
      }

      plan->recv_size[nrecv] = nqty*overlap.isize*overlap.jsize;
      ibuf += plan->recv_size[nrecv];
      nrecv++;
    }
  }

  // plan->nrecv = # of recvs not including self

  if (nrecv && plan->recv_proc[nrecv-1] == me) {
    plan->nrecv = nrecv - 1;
  }
  else {
    plan->nrecv = nrecv;
  }

  // init remaining fields in remap plan

  plan->memory = memory;

  if (nrecv == plan->nrecv) {
    plan->self = 0;
  }
  else {
    plan->self = 1;
  }

  // free locally malloced space

  free(array);

  // find biggest send message (not including self) and malloc space for it

  size = 0;
  for (nsend = 0; nsend < plan->nsend; nsend++) {
    size = MAX(size,plan->send_size[nsend]);
  }

  if (size) {
    if (precision == 1)
      plan->sendbuf = NULL;
    else
      plan->sendbuf = (double *) malloc(size*sizeof(double));
    if (plan->sendbuf == NULL) return NULL;
  }

  // if requested, allocate internal scratch space for recvs,
  // only need it if I will receive any data (including self)

  if (memory == 1) {
    if (nrecv > 0) {
      if (precision == 1)
        plan->scratch = NULL;
      else
        plan->scratch = (double *) malloc(nqty*out.isize*out.jsize*sizeof(double));
      if (plan->scratch == NULL) return NULL;
    }
  }

  // create new MPI communicator for remap

  MPI_Comm_dup(comm,&plan->comm);

  // return pointer to plan

  return plan;
}

/* ---------------------------------------------------------------------
   Destroy a 2d remap plan
--------------------------------------------------------------------- */

void remap_2d_destroy_plan_cfunc(struct remap_plan_2d *plan)
{
  MPI_Comm_free(&plan->comm);

  if (plan->nsend) {
    free(plan->send_offset);
    free(plan->send_size);
    free(plan->send_proc);
    free(plan->packplan);
    free(plan->sendbuf);
  }

  if (plan->nrecv) {
    free(plan->recv_offset);
    free(plan->recv_size);
    free(plan->recv_proc);
    free(plan->recv_bufloc);
    free(plan->request);
    free(plan->unpackplan);
    if (plan->memory) free(plan->scratch);
  }

  free(plan);
}

/* ---------------------------------------------------------------------
   collide 2 sets of indices to determine overlap
   compare bounds of block1 with block2 to see if they overlap
   return 1 if they do and put bounds of overlapping section in overlap
   return 0 if they do not overlap
--------------------------------------------------------------------- */

int remap_2d_collide(struct extent_2d *block1, struct extent_2d *block2, struct extent_2d *overlap)
{
  overlap->ilo = MAX(block1->ilo,block2->ilo);
  overlap->ihi = MIN(block1->ihi,block2->ihi);
  overlap->jlo = MAX(block1->jlo,block2->jlo);
  overlap->jhi = MIN(block1->jhi,block2->jhi);

  if (overlap->ilo > overlap->ihi || overlap->jlo > overlap->jhi) return 0;

  overlap->isize = overlap->ihi - overlap->ilo + 1;
  overlap->jsize = overlap->jhi - overlap->jlo + 1;

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

void pack_2d(double *data, double *buf, struct pack_plan_2d *plan)
{
  register int in,out,fast,slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  in = 0;
  for (slow = 0; slow < nslow; slow++) {
    out = slow*nstride;
    for (fast = 0; fast < nfast; fast++) {
      buf[in++] = data[out++];
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data
--------------------------------------------------------------------- */

void unpack_2d(double *buf, double *data, struct pack_plan_2d *plan)
{
  register int in,out,fast,slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    in = slow*nstride;
    for (fast = 0; fast < nfast; fast++) {
      data[in++] = buf[out++];
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, axis permutation, 1 value/element
--------------------------------------------------------------------- */

void unpack_2d_permute_1(double *buf, double *data, struct pack_plan_2d *plan)
{
  register int in,out,fast,slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    in = slow;
    for (fast = 0; fast < nfast; fast++, in += nstride) {
      data[in] = buf[out++];
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, axis permutation, 2 values/element
--------------------------------------------------------------------- */

void unpack_2d_permute_2(double *buf, double *data, struct pack_plan_2d *plan)
{
  register int in,out,fast,slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    in = 2*slow;
    for (fast = 0; fast < nfast; fast++, in += nstride) {
      data[in] = buf[out++];
      data[in+1] = buf[out++];
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, axis permutation, nqty values/element
--------------------------------------------------------------------- */

void unpack_2d_permute_n(double *buf, double *data, struct pack_plan_2d *plan)
{
  register int in,out,iqty,instart,fast,slow;
  register int nfast,nslow,nstride,nqty;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;
  nqty = plan->nqty;

  out = 0;
  for (slow = 0; slow < nslow; slow++) {
    instart = nqty*slow;
    for (fast = 0; fast < nfast; fast++, instart += nstride) {
      in = instart;
      for (iqty = 0; iqty < nqty; iqty++) {
        data[in++] = buf[out++];
      }
    }
  }
}

#endif
#ifdef PACK_POINTER

/* ---------------------------------------------------------------------
   pack from data -> buf
--------------------------------------------------------------------- */

void pack_2d(double *data, double *buf, struct pack_plan_2d *plan)
{
  register double *in,*out,*begin,*end;
  register int slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  in = buf;
  for (slow = 0; slow < nslow; slow++) {
    begin = &(data[slow*nstride]);
    end = begin + nfast;
    for (out = begin; out < end; out++) {
      *(in++) = *out;
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data
--------------------------------------------------------------------- */

void unpack_2d(double *buf, double *data, struct pack_plan_2d *plan)
{
  register double *in,*out,*begin,*end;
  register int slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    begin = &(data[slow*nstride]);
    end = begin + nfast;
    for (in = begin; in < end; in++) {
      *in = *(out++);
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, axis permutation, 1 value/element
--------------------------------------------------------------------- */

void unpack_2d_permute_1(double *buf, double *data, struct pack_plan_2d *plan)
{
  register double *in,*out,*begin,*end;
  register int slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    begin = &(data[slow]);
    end = begin + nfast*nstride;
    for (in = begin; in < end; in += nstride) {
      *in = *(out++);
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, axis permutation, 2 values/element
--------------------------------------------------------------------- */

void unpack_2d_permute_2(double *buf, double *data, struct pack_plan_2d *plan)
{
  register double *in,*out,*begin,*end;
  register int slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    begin = &(data[2*slow]);
    end = begin + nfast*nstride;
    for (in = begin; in < end; in += nstride) {
      *in = *(out++);
      *(in+1) = *(out++);
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, axis permutation, nqty values/element
--------------------------------------------------------------------- */

void unpack_2d_permute_n(double *buf, double *data, struct pack_plan_2d *plan)
{
  register double *in,*out,*instart,*begin,*end;
  register int iqty,slow;
  register int nfast,nslow,nstride,nqty;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;
  nqty = plan->nqty;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    begin = &(data[nqty*slow]);
    end = begin + nfast*nstride;
    for (instart = begin; instart < end; instart += nstride) {
      in = instart;
      for (iqty = 0; iqty < nqty; iqty++) {
        *(in++) = *(out++);
      }
    }
  }
}

#endif
#ifdef PACK_MEMCPY

/* ---------------------------------------------------------------------
   pack from data -> buf
--------------------------------------------------------------------- */

void pack_2d(double *data, double *buf, struct pack_plan_2d *plan)
{
  register double *in,*out;
  register int slow,size;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  size = nfast*sizeof(double);
  for (slow = 0; slow < nslow; slow++) {
    in = &(buf[slow*nfast]);
    out = &(data[slow*nstride]);
    memcpy(in,out,size);
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data
--------------------------------------------------------------------- */

void unpack_2d(double *buf, double *data, struct pack_plan_2d *plan)
{
  register double *in,*out;
  register int slow,size;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  size = nfast*sizeof(double);
  for (slow = 0; slow < nslow; slow++) {
    in = &(data[slow*nstride]);
    out = &(buf[slow*nfast]);
    memcpy(in,out,size);
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, one axis permutation, 1 value/element
--------------------------------------------------------------------- */

void unpack_2d_permute_1(double *buf, double *data, struct pack_plan_2d *plan)
{
  register double *in,*out,*begin,*end;
  register int slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    begin = &(data[slow]);
    end = begin + nfast*nstride;
    for (in = begin; in < end; in += nstride) {
      *in = *(out++);
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, axis permutation, 2 values/element
--------------------------------------------------------------------- */

void unpack_2d_permute_2(double *buf, double *data, struct pack_plan_2d *plan)
{
  register double *in,*out,*begin,*end;
  register int slow;
  register int nfast,nslow,nstride;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    begin = &(data[2*slow]);
    end = begin + nfast*nstride;
    for (in = begin; in < end; in += nstride) {
      *in = *(out++);
      *(in+1) = *(out++);
    }
  }
}

/* ---------------------------------------------------------------------
   unpack from buf -> data, axis permutation, nqty values/element
--------------------------------------------------------------------- */

void unpack_2d_permute_n(double *buf, double *data, struct pack_plan_2d *plan)
{
  register double *in,*out,*instart,*begin,*end;
  register int iqty,slow;
  register int nfast,nslow,nstride,nqty;

  nfast = plan->nfast;
  nslow = plan->nslow;
  nstride = plan->nstride;
  nqty = plan->nqty;

  out = buf;
  for (slow = 0; slow < nslow; slow++) {
    begin = &(data[nqty*slow]);
    end = begin + nfast*nstride;
    for (instart = begin; instart < end; instart += nstride) {
      in = instart;
      for (iqty = 0; iqty < nqty; iqty++) {
        *(in++) = *(out++);
      }
    }
  }
}

#endif

/* ------------------------------------------------------------------- */
