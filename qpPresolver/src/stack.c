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
 *	\file src/stack.c
 *	\author Dominik Cebulla
 *	\version 1.0 Beta
 *	\date 2017
 *
 *  Implementation of a stack (and corresponding data structures) which stores the
 *  order of the applied preprocessing techniques.
 */


#include <qpPresolver/stack.h>


/*  ======================================================================================
    Implementation of interface functions.
    ====================================================================================*/

qpp_stack_t *qppStackAlloc(const qpp_int_t size)
{
    qpp_stack_t *stack;

    if (size <= 0)
    {
        return NULL;
    }

    stack = (qpp_stack_t*) malloc(sizeof(qpp_stack_t));
    if (stack == NULL)
    {
        return NULL;
    }

    stack->size = 0;
    stack->elem = (qpp_stackelement_t*) calloc(size,  sizeof(qpp_stackelement_t));

    if (stack->elem == NULL)
    {
        free(stack);
        return NULL;
    }
    stack->max_size = size;
    stack->init_size = size;
    return stack;
}


void qppStackFree(qpp_stack_t *stack)
{
    qpp_stackelement_t *elem;

    if (stack != NULL)
    {
        while (!qppStackIsEmpty(stack))
        {
            elem = qppStackPop(stack);
            switch (elem->type)
            {
            case QPP_SE_FORCING_PC:
                free(elem->elements.fpc->indices);
                free(elem->elements.fpc->values);
                break;

            /*case QPP_SET_NEWBOUNDS_DC:
                free(elem->elements.nbdc->lb_indices);
                free(elem->elements.nbdc->ub_indices);
                free(elem->elements.nbdc->lb_values);
                free(elem->elements.nbdc->ub_values);
                break;*/

            case QPP_SE_NEWBOUNDS_PC:
                free(elem->elements.nbpc->lb_indices);
                free(elem->elements.nbpc->ub_indices);
                free(elem->elements.nbpc->lb_values);
                free(elem->elements.nbpc->ub_values);
                break;

            case QPP_SE_SCALING:
                free(elem->elements.scal->C);
                free(elem->elements.scal->R);

            default:
                break;
            }
            free((void*) elem->elements.fv);
        }
        free(stack->elem);
        free(stack);
    }
}


qpp_return_value_t qppStackReallocIf(qpp_stack_t *stack,
                                     const qpp_int_t num_additional_elements)
{
    qpp_int_t factor;
    qpp_stackelement_t *elem;

    if (stack->size + num_additional_elements < stack->max_size)
    {
        return QPP_OK;
    }

    factor = (qpp_int_t) ((stack->size + num_additional_elements - stack->max_size)
                          / stack->init_size) + 1;

    elem = (qpp_stackelement_t*) realloc(stack->elem, (sizeof(qpp_stackelement_t) *
                                         (stack->max_size + factor*stack->init_size)));

    if (elem == NULL)
    {
        return QPP_OUT_OF_MEMORY;
    }

    /*  stack->elem is deallocated by realloc */
    stack->elem = elem;
    stack->max_size += factor * stack->init_size;

    return QPP_OK;
}


inline void qppStackPush(qpp_stack_t *stack,
                         void *element,
                         const qpp_stackelement_type_t type)
{
    stack->elem[stack->size].elements.fv = element;
    stack->elem[stack->size].type = type;
    ++stack->size;
}


inline qpp_stackelement_t *qppStackPop(qpp_stack_t *stack)
{
    if (stack->size <= 0)
    {
        return NULL;
    }
    --stack->size;
    return &stack->elem[stack->size];
}


inline qpp_bool_type_t qppStackIsEmpty(const qpp_stack_t *stack)
{
    if (stack->size <= 0)
    {
        return QPP_BT_TRUE;
    }
    return QPP_BT_FALSE;
}
