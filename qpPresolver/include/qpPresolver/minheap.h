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
 *	\file include/qpPresolver/minheap.h
 *	\author Dominik Cebulla
 *	\version 1.0 Beta
 *	\date 2017
 *
 *  Declaration of a min-heap containing a key-value pair. Will be used as a priority
 *  queue in the presolver.
 */


#ifndef QPPRESOLVER_MINHEAP_H
#define QPPRESOLVER_MINHEAP_H

#include <stdlib.h>
#include <string.h>

#include <qpPresolver/constants.h>
#include <qpPresolver/types.h>



/** \brief Struct containing necessary data for the implementation of a min-heap.
 *
 *  The min-heap contains key-value pairs and is sorted with respect to its keys.
 */
typedef struct minheap
{
    qpp_int_t *key;          /**< The min-heap is sorted with respect to its keys. */
    qpp_int_t *value;        /**< Values corresponding to the respective keys. */
    qpp_int_t max_size;      /**< Length of the arrays \p key and \p value. Can be
                                  changed by \p mhRealloc(). */
    qpp_int_t size;          /**< Actual number of entries in the min-heap
                                  (i.e. \p size <= \p max_size). */
} qpp_minheap_t;



/*  ======================================================================================
    Interface functions.
    ====================================================================================*/

/** \brief Allocates memory for a new \p minheap.
 *
 *  \param size Maximum number of elements in the \p minheap, i.e. length of arrays
 *      \p key and \p value. Must be >= 0.
 *
 *  \return Valid pointer to a newly created \p minheap, or null pointer on failure.
 */
qpp_minheap_t *mhAlloc(const qpp_int_t size);


/** \brief Reallocates memory for the given heap.
 *
 *  Changes the size of the underlying arrays key and value. The new size can be
 *  either smaller (may lead to loss of data!) or larger than the original size.
 *
 *  \param heap Valid pointer to a \p minheap.
 *  \param new_size New size of the \p minheap.
 *
 *  \return QPP_OK \n QPP_FATAL_ERROR \n QPP_OUT_OF_MEMORY \n
 *      QPP_INVALID_ARGUMENT \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t mhRealloc(qpp_minheap_t *heap,
                             const qpp_int_t new_size);


/** \brief Deallocates memory allocated by the given heap.
 *
 *  \param heap Pointer to a \p minheap whose memory will be deallocated.
 */
void mhFree(qpp_minheap_t *heap);


/** \brief Places an existing element in the given heap at a position such that
 *      the (min-)heap property is satisfied.
 *
 *  \param heap Valid pointer to a \p minheap.
 *  \param i Array subscript of an element in the heap which will be placed in
 *      a correct position.
 *
 *  \return QPP_OK \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t mhHeapify(qpp_minheap_t *heap,
                             qpp_int_t i);


/** \brief Modifies a given heap such that the (min-)heap property is satisfied.
 *
 *  \param heap Valid pointer to a \p minheap.
 *
 *  \return QPP_OK \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t mhBuild(qpp_minheap_t *heap);


/** \brief Decreases the value of a key element and sorts it into the given heap
 *  such that the (min-)heap property is satisfied.
 *
 *  \param heap Valid pointer to a \p minheap.
 *  \param i Array subscript of the element whose \p key value shall be decreased.
 *  \param new_key New (smaller!) key value.
 *
 *  \return QPP_OK \n QPP_INVALID_ARGUMENT \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t mhDecreaseKey(qpp_minheap_t *heap,
                                 qpp_int_t i,
                                 const qpp_int_t new_key);


/** \brief Inserts a new element into the given heap (if possible). The min-heap
 *      property will be satisfied afterwards.
 *
 *  \param heap Valid pointer to a \p minheap.
 *  \param key Key of the new element.
 *  \param value Value of the new element.
 *
 *  \return QPP_OK \n QPP_OUT_OF_MEMORY \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t mhInsert(qpp_minheap_t *heap,
                            const qpp_int_t key,
                            const qpp_int_t value);


/** \brief Deletes an element from the given heap while maintaining the (min-)heap property.
 *
 *  \param heap Valid pointer to a \p minheap.
 *  \param i Array subscript of the element that will be removed from the heap.
 *
 *  \return QPP_OK \n QPP_INVALID_ARGUMENT \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t mhDelete(qpp_minheap_t *heap,
                            const qpp_int_t i);


/** \brief Extracts the value of the given heap with the smallest key.
 *
 *  Returns the value of the element with the smallest key, i.e. heap->value[0].
 *  Afterwards, the element is removed from the heap. The heap is modified such that the
 *  (min-)heap property is satisfied.
 *
 *  \param heap Valid pointer to a \p minheap.
 *  \param min If the heap is not empty, then \p *min stores the value (!) of the element
 *      with the smallest key.
 *
 *  \return QPP_OK \n QPP_INVALID_ARGUMENT \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t mhExtractMin(qpp_minheap_t *heap,
                                qpp_int_t *min);


/** \brief Searches in the given heap for an element with a specific value.
 *
 *  \param heap Valid pointer to a \p minheap.
 *  \param value We search for the first occurrence of a heap element with value \p value.
 *  \param index If an element with value \p value exists, then \p *index stores the array
 *      subscript of this element.
 *
 *  \return QPP_OK \n QPP_HEAP_ELEMENT_NOT_FOUND \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t mhSearchForValue(const qpp_minheap_t *heap,
                                    const qpp_int_t value,
                                    qpp_int_t *index);


/** \brief Checks if the given heap is empty.
 *
 *  \param heap Valid pointer to a \p minheap.
 *
 *  \return QPP_BT_TRUE if the heap does not contain any element, QPP_BT_FALSE otherwise.
 */
qpp_bool_type_t mhIsEmpty(const qpp_minheap_t *heap);


/** \brief Checks if the given heap satisfies the (min-)heap property.
 *
 *  \param heap Valid pointer to a \p minheap.
 *
 *  \return QPP_BT_TRUE if the heap satisfies the (min-)heap property,
 *      QPP_BT_FALSE otherwise.
 */
qpp_bool_type_t mhIsHeap(const qpp_minheap_t *heap);

#endif /* QPPRESOLVER_MINHEAP_H */
