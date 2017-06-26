/*
 *	This file is part of qpPresolver.
 *
 *	qpPresolver -- An implementation of presolving (= preprocessing) techniques
 *  for Quadratic Programming.
 *	Copyright (C) 2017 by Dominik Cebulla et al. All rights reserved.
 *
 *	qpPresolver is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpPresolver is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpPresolver; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file include/qpPresolver/stack.h
 *	\author Dominik Cebulla
 *	\version 1.0 Beta
 *	\date 2017
 *
 *  Declaration of a stack (and corresponding data structures) which stores the
 *  order of the applied preprocessing techniques.
 */


#ifndef QPPRESOLVER_STACK_H
#define QPPRESOLVER_STACK_H

#include <string.h>
#include <stdlib.h>

#include <qpPresolver/return_values.h>
#include <qpPresolver/types.h>


/** \brief Definition of stack element types. */
typedef enum stackelement_type
{
    QPP_SE_UNKNOWN = 0,             /**< Unknown stack element. */
    QPP_SE_FIXED_VARIABLE = 10,     /**< Fixing a variable procedure. */
    QPP_SE_STRICT_Y,                /**< (Not used) Identification of strictly positive / negative
                                         multipliers corresponding to the linear constraints. */
    QPP_SE_EMPTY_ROW,               /**< Empty rows procedure. */
    QPP_SE_SINGLETON_ROW,           /**< Singleton rows procedure. */
    QPP_SE_DOUBLETON_ROW,           /**< Doubleton rows procedure. */
    QPP_SE_FORCING_PC,              /**< Forcing primal constraints procedure. */
    QPP_SE_REDUNDANT_PC,            /**< Redundant primal constraints procedure. */
    QPP_SE_NEWBOUNDS_PC,            /**< Tightening variables bounds based on primal
                                         constraints procedure. */
    QPP_SE_SINGLETON_COLUMN,        /**< Singleton columns procedure. */
    QPP_SE_DUPLICATE_COLUMN,        /**< Duplicate columns procedure. */
    QPP_SE_SPARSIFICATION,          /**< Sparsification procedure. */
    QPP_SE_NEWBOUNDS_DC,            /**< (Not used) Tightening primal variable bounds
                                         based on dual constraints procedure. */
    QPP_SE_SCALING                  /**< Scaling procedure. */
} qpp_stackelement_type_t;


/** \brief Stack element containing data of the scaling procedure. */
typedef struct se_scaling
{
    qpp_real_t *C;      /**< Column scaling factors. */
    qpp_real_t *R;      /**< Row scaling factors. */
} qpp_se_scaling_t;


/** \brief Stack element containing data for the "fixing a variable" procedure. */
typedef struct se_fixed_variable
{
    qpp_int_t j;                    /**< Index of the fixed variable. */
    qpp_real_t value;               /**< Value at which the variable is fixed. */
    qpp_active_type_t actv_type;    /**< Determines if the variable is active at a bound. */
} qpp_se_fixed_variable_t;


/** \brief Stack element containing data for the empty rows procedure. */
typedef struct se_empty_rows
{
    qpp_int_t i;            /**< Index of the empty row / constraint. */
} qpp_se_empty_rows_t;


/** \brief Stack element containing data for the singleton rows procedure. */
typedef struct se_singleton_row
{
    qpp_int_t i;                /**< Index of the singleton row. */
    qpp_int_t j;                /**< Index of the variable occurring in row \p i. */
    qpp_real_t a_ij;            /**< Value of the (only) nonzero element in row \p i. */
    /*qpp_real_t xl_tight_old;
    qpp_real_t xu_tight_old;
    qpp_real_t xl_medium_old;
    qpp_real_t xu_medium_old;*/
    qpp_tighter_bounds_type_t tbt_t; /**< Indicates which tightest bounds got tighter. */
    qpp_tighter_bounds_type_t tbt_m; /**< Indicates which medium bounds got tighter. */
} qpp_se_singleton_row_t;


/** \brief Stack element containing data for the forcing primal constraints procedure. */
typedef struct se_forcing_pc
{
    qpp_int_t i;            /**< Index of the forcing row. */
    qpp_int_t num_neg;      /**< Number of negative nonzero elements in row \p i. */
    qpp_int_t num_pos;      /**< Number of positive nonzero elements in row \p i. */
    qpp_int_t *indices;     /**< Subscripts of the nonzero elements of row \p i. */
    qpp_real_t *values;     /**< Values of the nonzero elements of row \p i. */
    qpp_int_t sign;         /**< Determines if constraint \p i is active at its lower
                                 (\p sign = 1) or at its upper (\p sign = -1) bound. */
} qpp_se_forcing_pc_t;


/** \brief Stack element containing data for the redundant primal constraints procedure. */
typedef struct se_redundant_pc
{
    qpp_int_t i;            /**< Index of the redundant primal constraint. */
} qpp_se_redundant_pc_t;


/** \brief Stack element containing data for the bound tightening procedure based
 *      on the primal constraints. */
typedef struct se_new_bounds_pc
{
    qpp_int_t i;                /**< Index of the constraint from which the new bounds are derived. */
    qpp_int_t *lb_indices;      /**< Subscripts of the variables with new tighter lower bounds. */
    qpp_int_t *ub_indices;      /**< Subscripts of the variables with new tighter upper bounds. */
    qpp_real_t *lb_values;      /**< Values of the respective old (!) lower bounds. */
    qpp_real_t *ub_values;      /**< Values of the respective old(!) upper bounds. */
    qpp_int_t lb_num;           /**< Length of arrays \p lb_indices and \p lb_values. */
    qpp_int_t ub_num;           /**< Length of arrays \p ub_indices and \p ub_values. */
} qpp_se_new_bounds_pc_t;


/** \brief Stack element containing data for the singleton columns procedure. */
typedef struct se_singleton_column
{
    qpp_int_t i;                /**< Index of the constraint in which the only nonzero
                                     element of the singleton column occurs. */
    qpp_int_t j;                /**< Index of the singleton column. */
    qpp_int_t split;            /**< Indicates if row \p i has been split in order to
                                     free variable \p j. */
    qpp_real_t al;              /**< Old value of the lower bound on constraint \p i. */
    qpp_real_t au;              /**< Old value of the upper bound on constraint \p i. */
    qpp_real_t a_ij;            /**< Value of the (only) nonzero element in column \p j. */
} qpp_se_singleton_column_t;


/** \brief Stack element containing data for the duplicate columns procedure. */
typedef struct se_duplicate_column
{
    qpp_int_t j;                /**< Index of column 1. */
    qpp_int_t k;                /**< Index of column 2. */
    qpp_real_t alpha;           /**< It holds: Column \p j = \p alpha * column \p k. */
    qpp_real_t xl_tight_old;    /**< Old value of the tightest lower bound on variable \p k. */
    qpp_real_t xu_tight_old;    /**< Old value of the tightest upper bound on variable \p k. */
    qpp_real_t xl_medium_old;   /**< Old value of the medium lower bound on variable \p k. */
    qpp_real_t xu_medium_old;   /**< Old value of the medium upper bound on variable \p k. */
} qpp_se_duplicate_column_t;


/** \brief Stack element containing data for the sparsification procedure. */
typedef struct se_sparsification
{
    qpp_int_t i;            /**< Index of the \b pivot constraint. */
    qpp_int_t k;            /**< Index of the altered constraint. */
    qpp_int_t j;            /**< Index of the pivot column. */
    qpp_real_t alpha;       /**< We compute: Row \p k := Row \p k - \p alpha * row \p i. */
    qpp_real_t al;          /**< Old lower bound on constraint \p k. */
    qpp_real_t au;          /**< Old upper bound on constraint \p k. */
} qpp_se_sparsification_t;


/** \brief Stack element of the presolve stack. */
typedef struct stackelement
{
    qpp_stackelement_type_t type;       /**< Indicates, to which preprocessing procedure
                                             the stack element corresponds. */
    union
    {
        qpp_se_fixed_variable_t *fv;
        qpp_se_empty_rows_t *er;
        qpp_se_singleton_row_t *sr;
        qpp_se_forcing_pc_t *fpc;
        qpp_se_redundant_pc_t *rpc;
        qpp_se_new_bounds_pc_t *nbpc;
        qpp_se_singleton_column_t *sc;
        qpp_se_duplicate_column_t *dc;
        qpp_se_sparsification_t *spf;
        qpp_se_scaling_t *scal;
    } elements;
} qpp_stackelement_t;


/** \brief Presolve stack that stores the order of all applied preprocessing techniques. */
typedef struct stack
{
    qpp_stackelement_t *elem;   /**< Stack elements stored in an array. */
    qpp_int_t size;             /**< Actual number of elements in the stack. */
    qpp_int_t max_size;         /**< Maximum size of the stack (can be reallocated!). */
    qpp_int_t init_size;        /**< Memory will be (re-)allocated in \p init_size steps. */
} qpp_stack_t;


/*  ====================================================================================
    Interface functions.
    ====================================================================================*/

/** \brief Allocates memory for a new and empty \p stack.
 *
 *  \param size Sets the initial \p max_size of the \p stack. Also \p init_size is set
 *      to this value. Note that the maximum allowed size of the \p stack can be
 *      changed via \p qppStackReallocIf().
 *
 *  \return Pointer to a newly created \p stack on success, NULL on failure.
 */
qpp_stack_t *qppStackAlloc(const qpp_int_t size);


/** \brief Reallocates memory if there is not enough space left in the given \p stack.
 *
 *  It is checked if the stack is able to hold \p n_additional_elements new elements.
 *  If yes, then the stack is not altered, otherwise we try to enlarge the capacity of
 *  the stack by reallocating memory.
 *
 *  \param stack Valid pointer to a \p stack.
 *  \param num_additional_elements Number of elements that shall be added to the \p stack.
 *
 *  \return QPP_OK \n QPP_OUT_OF_MEMORY
 */
qpp_return_value_t qppStackReallocIf(qpp_stack_t *stack,
                                     const qpp_int_t num_additional_elements);


/** \brief Deallocates memory allocated by the given \p stack.
 *
 *  Note that this function also deallocates the memory allocated by the stack elements,
 *  i.e. if a stack element itself contains pointers, then the corresponding memory is
 *  also deallocated.
 *
 *  \param stack Pointer to a \p stack whose memory will be deallocated.
 */
void qppStackFree(qpp_stack_t *stack);


/** \brief Pushes a new element to the given \p stack.
 *
 *  \param stack Valid pointer to a \p stack.
 *  \param element Valid pointer to an arbitrary stack element that shall be added to the
 *      \p stack. All such elements have the name pattern \p se_*.
 *  \param type Type of the new stack element, see \p stackelement_type. This is necessary
 *      in order to identify which information the new \p element holds.
 */
void qppStackPush(qpp_stack_t *stack,
                  void *element,
                  const qpp_stackelement_type_t type);


/** Pops an element from the given \p stack.
 *
 *  \param stack Valid pointer to a \p stack.
 *
 *  \return Pointer to the top-most stack element.
 */
qpp_stackelement_t *qppStackPop(qpp_stack_t *stack);


/** Checks if the given \p stack is empty.
 *
 *  \param Valid pointer to a \p stack.
 *
 *  \return QPP_BT_TRUE if the stack is empty, QPP_BT_FALSE otherwise.
 */
qpp_bool_type_t qppStackIsEmpty(const qpp_stack_t *stack);

#endif /* QPPRESOLVER_STACK_H */
