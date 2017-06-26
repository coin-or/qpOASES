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
 *	\file src/minheap.c
 *	\author Dominik Cebulla
 *	\version 1.0 Beta
 *	\date 2017
 *
 *  Implementation of a min-heap containing a key-value pair. Will be used as a priority
 *  queue in the presolver.
 */


#include <qpPresolver/minheap.h>

/*  ======================================================================================
    Static (= private) functions.
    ====================================================================================*/

static inline qpp_int_t mhParent(const qpp_int_t i)
{
    return (i-1)/2;
}

static inline qpp_int_t mhLeft(const qpp_int_t i)
{
    return 2*i + 1;
}

static inline qpp_int_t mhRight(const qpp_int_t i)
{
    return 2 * (i+1);
}

static inline void swap(qpp_int_t *x,
                        qpp_int_t *y)
{
	qpp_int_t t;
	t = *x; *x = *y; *y = t;
}

/*  ======================================================================================
    Implementation of interface functions.
    ====================================================================================*/

qpp_minheap_t *mhAlloc(const qpp_int_t size)
{
    qpp_minheap_t *heap;

    if (size < 0)
    {
        return NULL;
    }

    heap = (qpp_minheap_t*) malloc(sizeof(qpp_minheap_t));

    if (heap == NULL)
    {
        return NULL;
    }

    if (size == 0)
    {
        heap->key = NULL;
        heap->value = NULL;
        heap->max_size = 0;
        heap->size = 0;
    }
    else
    {
        heap->key = (qpp_int_t*) malloc(2 * size * sizeof(qpp_int_t));
        if (heap->key == NULL)
        {
            free(heap);
            return NULL;
        }

        heap->value = &heap->key[size];
        heap->max_size = size;
        heap->size = 0;
    }

    return heap;
}


qpp_return_value_t mhRealloc(qpp_minheap_t *heap,
                             const qpp_int_t new_size)
{
    qpp_int_t *iptr;

    if (heap == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    if (new_size <= 0)
    {
        return QPP_INVALID_ARGUMENT;
    }

    if (heap->max_size == new_size)
    {
        return QPP_OK;
    }

    if (new_size < heap->max_size)
    {
        memmove(&heap->key[new_size], heap->value, new_size * sizeof(qpp_int_t));
        if (new_size < heap->size)
        {
            heap->size = new_size;
        }
    }

    iptr = (qpp_int_t*) realloc(heap->key, 2 * new_size * sizeof(qpp_int_t));
    if ( (iptr == NULL) && (new_size < heap->max_size) )
    {
        return QPP_FATAL_ERROR;
    }
    if ( (iptr == NULL) && (new_size > heap->max_size) )
    {
        return QPP_OUT_OF_MEMORY;
    }

    heap->key = iptr; iptr += new_size;
    heap->value = iptr;

    if (new_size > heap->max_size)
    {
        memmove(heap->value, &heap->key[heap->max_size], heap->size);
    }

    heap->max_size = new_size;

    return QPP_OK;
}


void mhFree(qpp_minheap_t *heap)
{
    if (heap != NULL)
    {
        free(heap->key);
        free(heap);
    }
}


qpp_return_value_t mhHeapify(qpp_minheap_t *heap,
                             qpp_int_t i)
{
    qpp_int_t l, r, m, hsize;
    qpp_int_t *hkey, *hval;

    if (heap == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    hsize = heap->size;
    hkey = heap->key;
    hval = heap->value;

    while ( 1 )
    {
        l = mhLeft(i);
        r = mhRight(i);
        m = i;

        if ( (l < hsize) && (hkey[l] < hkey[m]) )
        {
            m = l;
        }

        if ( (r < hsize) && (hkey[r] < hkey[m]) )
        {
            m = r;
        }

        if (m != i)
        {
            swap(&hkey[i], &hkey[m]);
            swap(&hval[i], &hval[m]);
            i = m;
        }
        else
        {
            break;
        }
    }
    return QPP_OK;
}


qpp_return_value_t mhBuild(qpp_minheap_t *heap)
{
    qpp_int_t i;

    if (heap == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    for (i = heap->size/2-1; i >= 0; --i)
    {
        mhHeapify(heap, i);
    }
    return QPP_OK;
}


qpp_return_value_t mhDecreaseKey(qpp_minheap_t *heap,
                                 qpp_int_t i,
                                 const qpp_int_t new_key)
{
    qpp_int_t p;
    qpp_int_t *hkey, *hval;

    if (heap == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

	hkey = heap->key;
	hval = heap->value;

    if (new_key >= hkey[i])
    {
        return QPP_INVALID_ARGUMENT;
    }

    hkey[i] = new_key;
    p = mhParent(i);

    while ( (i > 0) && (hkey[p] > hkey[i]) )
    {
        swap(&hkey[p], &hkey[i]);
        swap(&hval[p], &hval[i]);
        i = p;
        p = mhParent(i);
    }
    return QPP_OK;
}


qpp_return_value_t mhInsert(qpp_minheap_t *heap,
                            const qpp_int_t key,
                            const qpp_int_t value)
{
    if (heap == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    if (heap->size >= heap->max_size)
    {
        /* ToDo: Maybe another error code? Or reallocating? */
        return QPP_OUT_OF_MEMORY;
    }

    heap->key[heap->size] = QPP_MAX_INT;
    heap->value[heap->size] = value;
    mhDecreaseKey(heap, heap->size, key);
    ++heap->size;

    return QPP_OK;
}


qpp_return_value_t mhDelete(qpp_minheap_t *heap,
                            const qpp_int_t i)
{
    qpp_int_t *hkey, *hvalue;

    if (heap == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    if ( (i >= heap->size) || (i < 0) )
    {
        return QPP_INVALID_ARGUMENT;
    }

    hkey = heap->key;
    hvalue = heap->value;

    --heap->size;
    hkey[i] = hkey[heap->size];
    hvalue[i] = hvalue[heap->size];
    mhHeapify(heap, i);

    return QPP_OK;
}


qpp_return_value_t mhExtractMin(qpp_minheap_t *heap,
                                qpp_int_t *min)
{
    if (heap == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }
    if (heap->size <= 0)
    {
        return QPP_INVALID_ARGUMENT;
    }

    *min = heap->value[0];
    heap->size--;
    heap->key[0] = heap->key[heap->size];
    heap->value[0] = heap->value[heap->size];
    mhHeapify(heap, 0);

    return QPP_OK;
}


qpp_return_value_t mhSearchForValue(const qpp_minheap_t *heap,
                                    const qpp_int_t value,
                                    qpp_int_t *index)
{
    qpp_int_t i;
    qpp_int_t *hval;

    if ( (heap == NULL) || (index == NULL) )
    {
        return QPP_NULL_ARGUMENT;
    }

	hval = heap->value;

    for (i = 0; i < heap->size; ++i)
    {
        if (hval[i] == value)
        {
            *index = i;
            return QPP_OK;
        }
    }
    return QPP_HEAP_ELEMENT_NOT_FOUND;
}


qpp_bool_type_t mhIsEmpty(const qpp_minheap_t *heap)
{
    if (heap == NULL)
    {
        return QPP_BT_TRUE;
    }

    if (heap->size <= 0)
    {
        return QPP_BT_TRUE;
    }
    return QPP_BT_FALSE;
}


qpp_bool_type_t mhIsHeap(const qpp_minheap_t *heap)
{
    qpp_int_t i, l, r, n;
    qpp_int_t *hkey;

    if (heap == NULL)
    {
        return QPP_BT_FALSE;
    }

    hkey = heap->key;
    n = heap->size;

    for (i = 0; i < n/2; ++i)
    {
        l = mhLeft(i);
        r = mhRight(i);
        if ( (l < n) && (hkey[i] > hkey[l]) )
        {
            return QPP_BT_FALSE;
        }
        if ( (r < n) && (hkey[i] > hkey[r]) )
        {
            return QPP_BT_FALSE;
        }
    }
    return QPP_BT_TRUE;
}
