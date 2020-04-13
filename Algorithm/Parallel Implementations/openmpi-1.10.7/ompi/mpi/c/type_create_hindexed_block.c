/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil -*- */
/*
 * Copyright (c) 2012      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2013      Los Alamos National Security, LLC.  All rights
 *                         reserved.
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

#include "ompi_config.h"

#include "ompi/mpi/c/bindings.h"
#include "ompi/runtime/params.h"
#include "ompi/communicator/communicator.h"
#include "ompi/errhandler/errhandler.h"
#include "ompi/datatype/ompi_datatype.h"
#include "ompi/memchecker.h"

#if OPAL_HAVE_WEAK_SYMBOLS && OMPI_PROFILING_DEFINES
#pragma weak MPI_Type_create_hindexed_block = PMPI_Type_create_hindexed_block
#endif

#if OMPI_PROFILING_DEFINES
#include "ompi/mpi/c/profile/defines.h"
#endif

static const char FUNC_NAME[] = "MPI_Type_create_hindexed_block";


int MPI_Type_create_hindexed_block(int count,
                                   int blocklength,
                                   const MPI_Aint array_of_displacements[],
                                   MPI_Datatype oldtype,
                                   MPI_Datatype *newtype)
{
    int rc;

    MEMCHECKER(
        memchecker_datatype(oldtype);
        );

    if( MPI_PARAM_CHECK ) {
        OMPI_ERR_INIT_FINALIZE(FUNC_NAME);
        if( count < 0 ) {
            return OMPI_ERRHANDLER_INVOKE(MPI_COMM_WORLD, MPI_ERR_COUNT,
                                          FUNC_NAME);
        } else if( (count > 0) && (blocklength < 0 || NULL == array_of_displacements) ) {
            return OMPI_ERRHANDLER_INVOKE(MPI_COMM_WORLD, MPI_ERR_ARG,
                                          FUNC_NAME );
        } else if (NULL == oldtype || MPI_DATATYPE_NULL == oldtype ||
                   NULL == newtype) {
            return OMPI_ERRHANDLER_INVOKE(MPI_COMM_WORLD, MPI_ERR_TYPE,
                                          FUNC_NAME );
        }
    }

    OPAL_CR_ENTER_LIBRARY();

    rc = ompi_datatype_create_hindexed_block( count, blocklength, array_of_displacements,
                                              oldtype, newtype );
    if( rc != MPI_SUCCESS ) {
        ompi_datatype_destroy( newtype );
        OMPI_ERRHANDLER_RETURN( rc, MPI_COMM_WORLD, rc, FUNC_NAME );
    }
    {
        const int* a_i[2] = {&count, &blocklength};
        ompi_datatype_set_args( *newtype, 2, a_i, count, array_of_displacements, 1, &oldtype,
                                MPI_COMBINER_HINDEXED_BLOCK );
    }

    OPAL_CR_EXIT_LIBRARY();
    return MPI_SUCCESS;
}
