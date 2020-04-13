/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil -*- */
/*
 * Copyright (c) 2004-2005 The Trustees of Indiana University and Indiana
 *                         University Research and Technology
 *                         Corporation.  All rights reserved.
 * Copyright (c) 2004-2005 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2004-2005 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 * Copyright (c) 2004-2005 The Regents of the University of California.
 *                         All rights reserved.
 * Copyright (c) 2011-2012 Cisco Systems, Inc.  All rights reserved.
 * Copyright (c) 2014      Los Alamos National Security, LLC. All rights
 *                         reserved.
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

#include "ompi_config.h"

#include "ompi/mpi/fortran/mpif-h/bindings.h"

#if OPAL_HAVE_WEAK_SYMBOLS && OMPI_PROFILE_LAYER
#pragma weak PMPI_WIN_FLUSH_LOCAL = ompi_win_flush_local_f
#pragma weak pmpi_win_flush_local = ompi_win_flush_local_f
#pragma weak pmpi_win_flush_local_ = ompi_win_flush_local_f
#pragma weak pmpi_win_flush_local__ = ompi_win_flush_local_f

#pragma weak PMPI_Win_flush_local_f = ompi_win_flush_local_f
#pragma weak PMPI_Win_flush_local_f08 = ompi_win_flush_local_f
#elif OMPI_PROFILE_LAYER
OMPI_GENERATE_F77_BINDINGS (PMPI_WIN_FLUSH_LOCAL,
                           pmpi_win_flush_local,
                           pmpi_win_flush_local_,
                           pmpi_win_flush_local__,
                           pompi_win_flush_local_f,
                           (MPI_Fint *rank, MPI_Fint *win, MPI_Fint *ierr),
                           (rank, win, ierr) )
#endif

#if OPAL_HAVE_WEAK_SYMBOLS
#pragma weak MPI_WIN_FLUSH_LOCAL = ompi_win_flush_local_f
#pragma weak mpi_win_flush_local = ompi_win_flush_local_f
#pragma weak mpi_win_flush_local_ = ompi_win_flush_local_f
#pragma weak mpi_win_flush_local__ = ompi_win_flush_local_f

#pragma weak MPI_Win_flush_local_f = ompi_win_flush_local_f
#pragma weak MPI_Win_flush_local_f08 = ompi_win_flush_local_f
#endif

#if ! OPAL_HAVE_WEAK_SYMBOLS && ! OMPI_PROFILE_LAYER
OMPI_GENERATE_F77_BINDINGS (MPI_WIN_FLUSH_LOCAL,
                           mpi_win_flush_local,
                           mpi_win_flush_local_,
                           mpi_win_flush_local__,
                           ompi_win_flush_local_f,
                           (MPI_Fint *rank, MPI_Fint *win, MPI_Fint *ierr),
                           (rank, win, ierr) )
#endif


#if OMPI_PROFILE_LAYER && ! OPAL_HAVE_WEAK_SYMBOLS
#include "ompi/mpi/fortran/mpif-h/profile/defines.h"
#endif

void ompi_win_flush_local_f(MPI_Fint *rank, MPI_Fint *win, MPI_Fint *ierr)
{
    int c_ierr;
    MPI_Win c_win = MPI_Win_f2c(*win);

    c_ierr = MPI_Win_flush_local(OMPI_FINT_2_INT(*rank), c_win);
    if (NULL != ierr) *ierr = OMPI_INT_2_FINT(c_ierr);
}
