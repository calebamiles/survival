#include <Python.h>
#include <stdio.h>
#include <time.h>

#include "math.h"

#include "structmember.h"

#include <numpy/arrayobject.h>
#include <numpy/ufuncobject.h>


/* 
 * A random number generator implementation to replace
 * the random module
 */

/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/


/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */



/*
 *
 * Represents a single pyEvent that can happen once.  The typical example is death.
 * time - (float) The time to death or censorship
 * censored - (bool) True if the pyEvent was right censored; False otherwise.
 *
 */
typedef struct {
    PyObject_HEAD
    double event_time;
    int censored;
} pyEvent;

typedef struct {
    PyObject_HEAD
    pyEvent e;
} _event;

static NPY_INLINE int
pyevent_nonzero(pyEvent pe) {
    return pe.event_time != 0.0;
}

static PyTypeObject _event_Type;

static NPY_INLINE int
_event_Check(PyObject *object) {
    return PyObject_IsInstance(object,(PyObject*)&_event_Type);
}

static PyObject*
_event_FrompyEvent(pyEvent *pe) {
    _event *e = (_event *)_event_Type.tp_alloc(&_event_Type,0);
    if (e) {
        e->e = *pe;
    }

    return (PyObject *)e;
}

static void
event_dealloc(_event *self)
{
    Py_XDECREF(&self->e);
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
event_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    pyEvent pe;

    pe.event_time = 0.0;
    pe.censored = 0;

    return _event_FrompyEvent(&pe);
}

static PyObject*
_event_get_time(_event *self, void* closure) {
    return PyFloat_FromDouble(self->e.event_time);
}

static PyObject*
_event_get_censored(_event *self, void *closure) {
    return PyInt_FromLong(self->e.censored);
}

static int
_event_set_time(_event *self, PyObject *value, void *closure)
{
  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the time attribute");
    return -1;
  }

  if (!PyFloat_Check(value)) {
    PyErr_SetString(PyExc_TypeError, "Time attribute must be float");
    return -1;
  }

  self->e.event_time = PyFloat_AS_DOUBLE(value);

  return 0;
}

static int
_event_set_censored(_event *self, PyObject *value, void *closure)
{
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the last attribute");
        return -1;
    }

    if (!PyInt_Check(value)) {
        PyErr_SetString(PyExc_TypeError, "The last attribute value must be an integer (boolean)");
        return -1;
    }

    // should force C style boolean semantics here....
    self->e.censored = PyInt_AsLong(value);

    return 0;
}

static PyGetSetDef event_getset[] = {
    {"time",(getter)_event_get_time,(setter)_event_set_time,"time", NULL},
    {"censored",(getter)_event_get_censored,(setter)_event_set_censored,"censored", NULL},
    {NULL} /* sentinel */
};

static int
event_init(_event *self, PyObject *args, PyObject *kwds)
{
    pyEvent pe;
    int censored = 0;

    static char *kwlist[] = {"time", "censored", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "d|i", kwlist, &pe.event_time,
            &censored)) {

        PyErr_SetString(PyExc_RuntimeError, "event_init: unable to parse input arguments");
        return -1;
    }

    if (censored < 0) {
        PyErr_SetString(PyExc_ValueError, "event_init: invalid censored value");
        return -1;
    }

    if (censored == 0)
        pe.censored = 0;
    else if (censored > 0)
        pe.censored = 1;

    self->e = pe;

    return 0;
}

static PyMemberDef event_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef event_methods[] = {
    {NULL}  /* Sentinel */
};

/* _event presentation methods */
static PyObject *
event_repr(_event *self)
{
    char *val = PyOS_double_to_string(self->e.event_time, 'f', 3, 0, NULL);

    if (val != NULL)
        return PyString_FromFormat("_event at %p, time: %s", self, val);
    else
        return PyString_FromFormat("_event at %p", self);
}

static PyObject *
event_print(_event *self)
{
    char *val = PyOS_double_to_string(self->e.event_time, 'f', 3, 0, NULL);

    if (val != NULL)
        return PyString_FromFormat("%s", val);
    else
        return PyString_FromString("");
}

static PyTypeObject _event_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "_Survival._event",        /*tp_name*/
    sizeof(_event),            /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)event_dealloc, /*tp_dealloc*/
    (printfunc)event_print,    /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    (reprfunc)event_repr,      /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Non-repeatable single event",       /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    event_methods,                     /* tp_methods */
    event_members,             /* tp_members */
    event_getset,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)event_init,      /* tp_init */
    0,                         /* tp_alloc */
    event_new,                 /* tp_new */
};

/* Numpy support */

static PyObject*
npy_event_getitem(void* data, void* arr) {
    pyEvent pe;
    memcpy(&pe, data, sizeof(pyEvent));
    return _event_FrompyEvent(&pe);
}

static int
npy_event_setitem(PyObject* item, void* data, void* arr) {
    pyEvent pe;
    if (_event_Check(item)) {
        pe = ((_event *)item)->e;
    } else {
        return -1; // sorry no graceful error handling or casting...
    }
    memcpy(data, &pe, sizeof(pyEvent));
    return 0;
}

static NPY_INLINE void
byteswapint(int* x) {
    char* p = (char*)x;
    size_t i;
    for (i = 0; i < sizeof(*x)/2; i++) {
        int j = sizeof(*x)-1-i;
        char t = p[i];
        p[i] = p[j];
        p[j] = t;
    }
}

static NPY_INLINE void
byteswapdouble(double* x) {
    char* p = (char*)x;
    size_t i;
    for (i = 0; i < sizeof(*x)/2; i++) {
        int j = sizeof(*x)-1-i;
        char t = p[i];
        p[i] = p[j];
        p[j] = t;
    }
}


static void
npy_event_copyswapn(void* dst_, npy_intp dstride, void* src_, npy_intp sstride, npy_intp n, int swap, void* arr) {
    char *dst = (char*)dst_, *src = (char*)src_;
    if (!src) {
        return;
    }
    npy_intp i;
    if (swap) {
        for (i = 0; i < n; i++) {
            pyEvent* pe = (pyEvent*)(dst+dstride*i);
            memcpy(pe,src+sstride*i,sizeof(pyEvent));
            byteswapdouble(&pe->event_time);
            byteswapint(&pe->censored);
        }
    } else if (dstride==sizeof(pyEvent) && sstride==sizeof(pyEvent)) {
        memcpy(dst,src,n*sizeof(pyEvent));
    } else {
        for (i = 0; i < n; i++) {
            memcpy(dst+dstride*i,src+sstride*i,sizeof(pyEvent));
        }
    }
}

static void
npy_event_copyswap(void* dst, void* src, int swap, void* arr) {
    if (!src) {
        return;
    }
    pyEvent* pe = (pyEvent*)dst;
    memcpy(pe,src,sizeof(pyEvent));
    if (swap) {
        byteswapdouble(&pe->event_time);
        byteswapint(&pe->censored);
    }
}

static npy_bool
npy_event_nonzero(void* data, void* arr) {
    pyEvent pe;
    memcpy(&pe,data,sizeof(pe));
    return pyevent_nonzero(pe)?NPY_TRUE:NPY_FALSE;
}

static PyArray_ArrFuncs npy_event_arrfuncs;

typedef struct { char c; pyEvent pe; } e_align_test;

PyArray_Descr npy_event_descr = {
    PyObject_HEAD_INIT(0)
    &_event_Type,       /* typeobj */
    'O',                    /* kind */
    'e',                    /* type */
    '=',                    /* byteorder */
    NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM, /* hasobject */
    0,                      /* type_num */
    sizeof(pyEvent),       /* elsize */
    offsetof(e_align_test,pe), /* alignment */
    0,                      /* subarray */
    0,                      /* fields */
    0,                      /* names */
    &npy_event_arrfuncs,  /* f */
};

/* multiEvent type */

/* This is kind of a wart, ideally I believe the best thing to do would be 
 * to subclass the ndarray because that's really what this object is, but
 * alas hopefully this first pass will suffice
 */

/*
 * Represents a repeatable event, such as a heart attack.
 * uncensored_events - (array) An array of Event objects, all of which have censored=False.
 *
 */
typedef struct {
    PyObject_HEAD
    PyArrayObject *events; // events should contain array of _event objects
} pyMultiEvent;

typedef struct {
    PyObject_HEAD
    pyMultiEvent m;
} _multiEvent;

static NPY_INLINE int
pyMultievent_nonzero(pyMultiEvent *pm) {
    // get a pointer to the events array
    PyArrayObject *events = (PyArrayObject *)pm->events;

    // assert that the number of non-zero events is greater than zero
    return PyArray_CountNonzero(events) != 0;
}

static PyTypeObject _multiEvent_Type;

static NPY_INLINE int
_multiEvent_Check(PyObject *object) {
    return PyObject_IsInstance(object,(PyObject*)&_multiEvent_Type);
}

static PyObject*
_multiEvent_From_pyMultiEvent(pyMultiEvent pm) {
    _multiEvent *m = (_multiEvent*)_multiEvent_Type.tp_alloc(&_multiEvent_Type,0);

   if (m) {
        m->m = pm;
    }

    return (PyObject *)m;
}

static void
_multiEvent_dealloc(_multiEvent *self)
{
    //Py_XDECREF(self->m.events);
    PyArray_XDECREF(self->m.events);
    Py_XDECREF(&self->m);
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
_multiEvent_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    // create the python multiEvent Object
    pyMultiEvent pm;

    // initialize empty events array
    pm.events = NULL; //(PyArrayObject *)PyArray_Empty(ndim, dims, &npy_event_descr, 0);

    return _multiEvent_From_pyMultiEvent(pm);
}

static PyArrayObject*
_multiEvent_get_events(_multiEvent *self, void* closure) {
    PyArrayObject *events = self->m.events;

    return (events);
}

static int
_multiEvent_set_events(_multiEvent *self, PyObject *value, void *closure)
{
  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the events attribute");
    return -1;
  }

  // could probably check specific array type
  if (!PyArray_Check(value)) {
    PyErr_SetString(PyExc_TypeError, "The events attribute must be an array");
    return -1;
  }

  self->m.events = (PyArrayObject *)(value);

  return 0;
}

static PyGetSetDef _multiEvent_getset[] = {
    {"events",(getter)_multiEvent_get_events,
              (setter)_multiEvent_set_events,
              "events",
              NULL},
    {NULL} /* sentinel */
};



static int
_multiEvent_init(_multiEvent *self, PyObject *args, PyObject *kwds)
{
    PyArrayObject *events = NULL;

    NpyIter *iter;
    NpyIter_IterNextFunc *iter_next;
    char **dataptrarray;

    double last_time = 0.0;
    int decumulate_times = 0;

    static char *kwlist[] = {"events", "decumulate_times", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|i", kwlist,
            &PyArray_Type, &events, &decumulate_times)) {

        PyErr_SetString(PyExc_ValueError, "multiEvent_init: unable to parse input arguments");
        return -1;
    }

    if (events == NULL) {
        PyErr_SetString(PyExc_ValueError, "invalid array of events");
        return -1;
    }

    npy_uint32 flags = NPY_ITER_C_INDEX | NPY_ITER_REFS_OK | NPY_ITER_READONLY;
    iter = NpyIter_New(events, flags, NPY_CORDER, NPY_SAFE_CASTING, &npy_event_descr);

    if (iter == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "unable to iterate over events"); // fix error code
        return -1; 
    }

    iter_next = NpyIter_GetIterNext(iter, NULL);
    dataptrarray = NpyIter_GetDataPtrArray(iter);

    do  {
        pyEvent *event = (pyEvent *)dataptrarray[0]; 

        // fail if any of the events are censored
        if (event->censored == 1) {
            PyErr_SetString(PyExc_ValueError, "right censored event passed"); // fix error code
            return -1; 
        }

        if (decumulate_times)
            event->event_time -= last_time;

    } while (iter_next(iter));

    NpyIter_Deallocate(iter);

    pyMultiEvent pm;
    pm.events = events;

    self->m = pm;

    return 0;
}

static PyObject*
_multiEvent_Time(_multiEvent *self)
{
    double event_times = 0.0;

    NpyIter *iter;
    NpyIter_IterNextFunc *iter_next;
    char **dataptrarray;

    npy_uint32 flags = NPY_ITER_C_INDEX | NPY_ITER_REFS_OK | NPY_ITER_READONLY;
    iter = NpyIter_New(self->m.events, flags, NPY_CORDER, NPY_SAFE_CASTING, &npy_event_descr);
    if (iter == NULL) {
        PyErr_SetString(PyExc_RuntimeError,
                "_multiEvent_Time: unable to create events iterator");

        goto fail;
    }

    iter_next = NpyIter_GetIterNext(iter, NULL);
    if (iter_next == NULL) {
        PyErr_SetString(PyExc_RuntimeError,
                "_multiEvent_Time: unable to iterate over events");

        goto fail;
    }

    dataptrarray = NpyIter_GetDataPtrArray(iter);

    do {
        pyEvent *event = (pyEvent *)dataptrarray[0];
        event_times += event->event_time;
    } while(iter_next(iter));

    Py_INCREF(self->m.events); // will SEGFAULT otherwise :(
    NpyIter_Deallocate(iter);

    return PyFloat_FromDouble(event_times);

fail:
    NpyIter_Deallocate(iter);
    return Py_BuildValue("");
}


/* multiEvent presentation functions */
static PyObject *
multiEvent_repr(_multiEvent *self)
{
    double event_times = 0.0;
    int num_events = 0;
    PyArrayIterObject *iter;

    iter = (PyArrayIterObject *)PyArray_IterNew((PyObject *)self->m.events);

    while (iter->index < iter->size)  {
        pyEvent *event = (pyEvent *)iter->dataptr;

       event_times += event->event_time;
       num_events++;

       PyArray_ITER_NEXT(iter);
    }

    char *val = PyOS_double_to_string(event_times, 'f', 3, 0, NULL);
    if (val != NULL)
        return PyString_FromFormat("_multiEvent at %p, %d event(s), time: %s", self, num_events, val);
    else
        return PyString_FromFormat("_multiEvent at %p", self);
}

static PyObject *
multiEvent_print(_multiEvent *self)
{
    double event_times = 0.0;
    int num_events = 0;
    PyArrayIterObject *iter;

    iter = (PyArrayIterObject *)PyArray_IterNew((PyObject *)self->m.events);

    while (iter->index < iter->size)  {
        pyEvent *event = (pyEvent *)iter->dataptr;

       event_times += event->event_time;
       num_events++;

       PyArray_ITER_NEXT(iter);
    }

    char *val = PyOS_double_to_string(event_times, 'f', 3, 0, NULL);
    if (val != NULL)
        return PyString_FromFormat("%d event(s), time: %s", num_events, val);
    else
        return PyString_FromFormat("%d event(s)", num_events);
}

static PyMemberDef _multiEvent_members[] = {
    {NULL}  /* Sentinel */
};

// modify so that arguments update the underlying fields
static PyMethodDef _multiEvent_methods[] = {
    {"time", (PyCFunction)_multiEvent_Time, METH_NOARGS,
     "Return the sum of the uncensored events"
    },
    {NULL}  /* Sentinel */
};


static PyTypeObject _multiEvent_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                               /*ob_size*/
    "_Survival._multiEvent",         /*tp_name*/
    sizeof(_multiEvent),             /*tp_basicsize*/
    0,                               /*tp_itemsize*/
    (destructor)_multiEvent_dealloc, /*tp_dealloc*/
    (printfunc)multiEvent_print,     /*tp_print*/
    0,                               /*tp_getattr*/
    0,                               /*tp_setattr*/
    0,                               /*tp_compare*/
    (reprfunc)multiEvent_repr,       /*tp_repr*/
    0,                               /*tp_as_number*/
    0,                               /*tp_as_sequence*/
    0,                               /*tp_as_mapping*/
    0,                               /*tp_hash */
    0,                               /*tp_call*/
    0,                               /*tp_str*/
    0,                               /*tp_getattro*/
    0,                               /*tp_setattro*/
    0,                               /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Represents a repeatable event, such as a heart attack.",       /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    _multiEvent_methods,                     /* tp_methods */
    _multiEvent_members,             /* tp_members */
    _multiEvent_getset,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)_multiEvent_init,      /* tp_init */
    0,                         /* tp_alloc */
    _multiEvent_new,                 /* tp_new */
};

/* Numpy support */

static PyObject*
npy_multiEvent_getitem(void* data, void* arr) {
    pyMultiEvent pm;
    memcpy(&pm ,data,sizeof(pyMultiEvent));
    return _multiEvent_From_pyMultiEvent(pm);
}

static int
npy_multiEvent_setitem(PyObject* item, void* data, void* arr) {
    pyMultiEvent pm;
    if (_multiEvent_Check(item)) {
        pm = ((_multiEvent *)item)->m;
    }
    else {
        // sorry no graceful error handling or casting...
        return -1;
    }
    memcpy(data,&pm,sizeof(pyMultiEvent));
    return 0;
}


static void
npy_multiEvent_copyswapn(void* dst_, npy_intp dstride, void* src_, npy_intp sstride, npy_intp n, int swap, void* arr) {
    char *dst = (char*)dst_, *src = (char*)src_;
    if (!src) {
        return;
    }
    npy_intp i;
    if (swap) {
        for (i = 0; i < n; i++) {
            pyMultiEvent* pm = (pyMultiEvent*)(dst+dstride*i);
            memcpy(pm,src+sstride*i,sizeof(pyMultiEvent));
            //byteswapdouble(&pm->time);
            //byteswapint(&pm->censored);
        }
    }
    else if (dstride==sizeof(pyMultiEvent) && sstride==sizeof(pyMultiEvent)) {
        memcpy(dst,src,n*sizeof(pyMultiEvent));
    }
    else {
        for (i = 0; i < n; i++) {
            memcpy(dst+dstride*i,src+sstride*i,sizeof(pyMultiEvent));
        }
    }
}

static void
npy_multiEvent_copyswap(void* dst, void* src, int swap, void* arr) {
    if (!src) {
        return;
    }
    pyMultiEvent* pm = (pyMultiEvent*)dst;
    memcpy(pm,src,sizeof(pyMultiEvent));
    if (swap) {
        //byteswapdouble(&pm->time);
        //byteswapint(&pm->censored);
    }
}

static npy_bool
npy_multiEvent_nonzero(void* data, void* arr) {
    pyMultiEvent pm;
    memcpy(&pm,data,sizeof(pm));
    return pyMultievent_nonzero(&pm)?NPY_TRUE:NPY_FALSE;
}

static PyArray_ArrFuncs npy_multiEvent_arrfuncs;

typedef struct { char c; pyMultiEvent pm; } m_align_test;

PyArray_Descr npy_multiEvent_descr = {
  PyObject_HEAD_INIT(0)
  &_multiEvent_Type,       /* typeobj */
  'O',                    /* kind */
  'e',                    /* type */
  '=',                    /* byteorder */
  NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM, /* hasobject */
  0,                      /* type_num */
  sizeof(pyMultiEvent),       /* elsize */
  offsetof(m_align_test,pm), /* alignment */
  0,                      /* subarray */
  0,                      /* fields */
  0,                      /* names */
  &npy_multiEvent_arrfuncs,  /* f */
};

/* ufunc definitions */

/*
 * The calling signature for _logp, that is this function maps an
 * array (value), a scalar (hazard), and a scalar (period) to a 
 * scalar
 */
char *_logp_signature = "(i), (), () -> ()";


static void
_logp_loop(char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(func))
{
  npy_intp dN = *dimensions++;
  npy_intp N_;                
  npy_intp s0 = *steps++;
  npy_intp s1 = *steps++;
  npy_intp s2 = *steps++;
  npy_intp s3 = *steps++;
  npy_intp di = dimensions[0];
  npy_intp i;

  npy_intp event_times_stride=steps[0];

  for (N_ = 0; N_ < dN; N_++, args[0] += s0, args[1] += s1, args[2] += s2, args[3] += s3) {
   char *event_times=args[0], *hazard=args[1], *period=args[2], *result=args[3];
    double remaining_time = (*(double *)period);
    double tmp = 0;
    for (i = 0; i < di; i++) {
      tmp += log((*(double *)hazard)) - ((*(double *)hazard) * (*(double *)event_times));
      remaining_time -= (*(double *)event_times);      
      event_times += event_times_stride;
    }
    
    if (remaining_time < 0) {
      *(double *)result = 0.0;
    } else {
       tmp += -1.0 * (*(double *)hazard) * remaining_time;
      *(double *)result = tmp;
    }
  }
}

static PyUFuncGenericFunction _logp_functions[] = { _logp_loop };
static void * _logp_data[] = { (void *)NULL};
char _logp_signatures[] = { NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE };

static PyObject *
_event_logp(PyObject *self, PyObject *args, PyObject *keywds)
{
    PyObject *events_obj = NULL, *hazards_obj = NULL, *periods_obj = NULL;
    PyArrayObject *events = NULL, *hazards = NULL, *periods = NULL;
    PyArrayObject *op[3];

    NpyIter *iter;
    NpyIter_IterNextFunc *iter_next;
    char **dataptrarray;


    static char *kwlist[] = {"value", "hazard", "period", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, keywds,  "OOO|", kwlist,
            &events_obj, &hazards_obj, &periods_obj))
        goto fail;

    // we can't assume that we were passed arrays
    int reqmts = NPY_IN_ARRAY;
    events =  (PyArrayObject *)PyArray_FROM_OTF(events_obj, NPY_NOTYPE, reqmts);
    hazards = (PyArrayObject *)PyArray_FROM_OTF(hazards_obj, NPY_DOUBLE, reqmts);
    periods = (PyArrayObject *)PyArray_FROM_OTF(periods_obj, NPY_DOUBLE, reqmts);

    if (events == NULL || hazards == NULL || periods == NULL) {
        goto fail;
    }
 
    npy_uint32 flags = NPY_ITER_REFS_OK |  NPY_ITER_MULTI_INDEX;
    npy_uint32 op_flags[3];
    op_flags[0] = NPY_ITER_READONLY;
    op_flags[1] = NPY_ITER_READONLY;
    op_flags[2] = NPY_ITER_READONLY;

    PyArray_Descr *dtypes[3];
    dtypes[0] = &npy_event_descr;
    dtypes[1] = PyArray_DescrFromType(NPY_DOUBLE);
    dtypes[2] = PyArray_DescrFromType(NPY_DOUBLE);
   
    op[0] = events;
    op[1] = hazards;
    op[2] = periods;
 
    iter = NpyIter_MultiNew(3, op, flags, NPY_KEEPORDER, 
            NPY_NO_CASTING, op_flags, dtypes);
    if (iter == NULL)
        goto fail;

    dataptrarray = NpyIter_GetDataPtrArray(iter);
    iter_next = NpyIter_GetIterNext(iter, NULL);
    if (iter_next == NULL) {
        NpyIter_Deallocate(iter);
        PyErr_SetString(PyExc_RuntimeError, "event_logp: unable to iterate over input arguments");
        goto fail;
    }

    double result = 0.0;

    do {
        pyEvent *event = (pyEvent *)dataptrarray[0];
        double hazard = *(double *)dataptrarray[1];
        double period = *(double *)dataptrarray[2];

        if (event->event_time > period)
            continue;

        result += ((1 - event->censored) * (hazard) - ((event->event_time) * (hazard)));
    } while (iter_next(iter));

    Py_DECREF(hazards);
    Py_DECREF(periods);
    Py_DECREF(events);

    NpyIter_Deallocate(iter);

    return PyFloat_FromDouble(result);

fail:
    Py_XDECREF(hazards);
    Py_XDECREF(periods);
    Py_XDECREF(events);
    Py_XDECREF(hazards_obj);
    Py_XDECREF(periods_obj);
    Py_XDECREF(events_obj);

    return Py_BuildValue(""); // essentailly a null object
}

static PyObject *
_event_random(PyObject *self, PyObject *args, PyObject *keywds)
{
    PyArrayObject *op[3];
    PyObject *hazards_obj = NULL, *periods_obj = NULL;
    PyArrayObject  *hazards = NULL, *periods = NULL, *out = NULL;
    NpyIter *iter;
    NpyIter_IterNextFunc *iter_next;
    char **dataptrarray;
 
    // seed the random number generator
    unsigned long seed = time(NULL);
    init_genrand(seed);

    static char *kwlist[] = {"hazard", "period", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, keywds,  "OO|", kwlist,
                &hazards_obj, &periods_obj)) {
        PyErr_SetString(PyExc_ValueError,
                "event_random: unable to parse input arguments");

        goto fail;
    }

    // we can't assume that we were passed arrays
    int reqmts = NPY_IN_ARRAY;
    hazards = (PyArrayObject *)PyArray_FROM_OTF(hazards_obj, NPY_DOUBLE, reqmts);
    periods = (PyArrayObject *)PyArray_FROM_OTF(periods_obj, NPY_DOUBLE, reqmts);

    if (hazards == NULL || periods == NULL) {
        goto fail;
    }

    npy_uint32 flags = NPY_ITER_REFS_OK |  NPY_ITER_MULTI_INDEX;
    npy_uint32 op_flags[3];
    op_flags[0] = NPY_ITER_READONLY;
    op_flags[1] = NPY_ITER_READONLY; 
    op_flags[2] = NPY_ITER_WRITEONLY | NPY_ITER_ALLOCATE;

    PyArray_Descr *dtypes[3];
    dtypes[0] = PyArray_DescrFromType(NPY_DOUBLE);
    dtypes[1] = PyArray_DescrFromType(NPY_DOUBLE);
    dtypes[2] = &npy_event_descr;

    op[0] = hazards;
    op[1] = periods;
    op[2] = NULL;

    iter = NpyIter_MultiNew(3, op, flags, NPY_KEEPORDER, 
            NPY_NO_CASTING, op_flags, dtypes);
    if (iter == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "event_random: unable to create iterator for input arguments");
        goto fail;
    }
    
    dataptrarray = NpyIter_GetDataPtrArray(iter);
    iter_next = NpyIter_GetIterNext(iter, NULL);

    if (iter_next == NULL) {
        NpyIter_Deallocate(iter);
        PyErr_SetString(PyExc_RuntimeError, "event_random: unable to iterate over input arguments");
        goto fail;
    }

    do {
        double hazard = *(double *)dataptrarray[0];
        double period = *(double *)dataptrarray[1];

        double f = genrand_real1();
        double event_time = -1.0 * log(1.0 - f) / hazard;
        int censored = (event_time > period);
        
        pyEvent event;
        event.event_time = event_time;
        event.censored = censored;

        memcpy(dataptrarray[2], &event, sizeof(pyEvent));
    } while(iter_next(iter));

    out = NpyIter_GetOperandArray(iter)[2];
    Py_INCREF(out);

    NpyIter_Deallocate(iter);

    Py_DECREF(hazards);
    Py_DECREF(periods);

    if (PyList_Check(hazards_obj))
        Py_DECREF(hazards_obj);

    if (PyList_Check(periods_obj))
        Py_DECREF(periods_obj);


    return (PyObject *)out; 

fail:
    Py_XDECREF(hazards);
    Py_XDECREF(periods);
    Py_XDECREF(hazards_obj);
    Py_XDECREF(periods_obj);

    return Py_BuildValue(""); // essentailly a null object
}

/* logp for arrays of MultiEvent type objects */

/*
 * We use the array iterator API to implement logp
 * a function maps an array of arrays (value), 
 * an array (hazard), and an array (period) to a 
 * scalar
 */

static PyObject * 
_multiEvent_logp(PyObject *self, PyObject *args, PyObject *keywds)
{
    PyArrayObject *op[3];
    PyObject *multiEvents_obj = NULL, *hazards_obj = NULL, *periods_obj = NULL;
    PyArrayObject *multiEvents = NULL, *hazards = NULL, *periods = NULL;

    NpyIter *iter;
    NpyIter_IterNextFunc *iter_next;
    char **dataptrarray;

    static char *kwlist[] = {"value", "hazard", "period", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds,  "OOO|", kwlist,
            &multiEvents_obj, &hazards_obj, &periods_obj)) {
        PyErr_SetString(PyExc_ValueError,
                "multiEvent_logp: unable to parse input arguments");

        goto fail;
    }

    // we can't assume that we were passed arrays
    int reqmts = NPY_IN_ARRAY;
    multiEvents =  (PyArrayObject *)PyArray_FROM_OTF(multiEvents_obj, NPY_NOTYPE, reqmts);
    hazards = (PyArrayObject *)PyArray_FROM_OTF(hazards_obj, NPY_DOUBLE, reqmts);
    periods = (PyArrayObject *)PyArray_FROM_OTF(periods_obj, NPY_DOUBLE, reqmts);

    if (multiEvents == NULL || hazards == NULL || periods == NULL) {
        goto fail;
    }

    npy_uint32 flags = NPY_ITER_REFS_OK |  NPY_ITER_MULTI_INDEX;
    npy_uint32 op_flags[3];
    op_flags[0] = NPY_ITER_READONLY;
    op_flags[1] = NPY_ITER_READONLY; 
    op_flags[2] = NPY_ITER_READONLY;

    PyArray_Descr *dtypes[3];
    dtypes[0] = &npy_multiEvent_descr;
    dtypes[1] = PyArray_DescrFromType(NPY_DOUBLE);
    dtypes[2] = PyArray_DescrFromType(NPY_DOUBLE);

    op[0] = multiEvents; 
    op[1] = hazards;
    op[2] = periods;

    iter = NpyIter_MultiNew(3, op, flags, NPY_KEEPORDER, 
            NPY_NO_CASTING, op_flags, dtypes);

    if (iter == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "multiEvent_logp: unable to create iterator for input arguments");
        goto fail;
    }

    dataptrarray = NpyIter_GetDataPtrArray(iter);
    iter_next = NpyIter_GetIterNext(iter, NULL);

    if (iter_next == NULL) {
        NpyIter_Deallocate(iter);
        PyErr_SetString(PyExc_RuntimeError, "multiEvent_logp: unable to iterate over input arguments");
        goto fail;
    }

    double result = 0.0;

    do {
        pyMultiEvent *multiEvent = (pyMultiEvent *)dataptrarray[0];
        double hazard = *(double *)dataptrarray[1];
        double remaining_time = *(double *)dataptrarray[2];
        double tmp = 0.0;

        PyArrayObject *events = multiEvent->events;
        if (events == NULL) {
            PyErr_SetString(PyExc_RuntimeError, "multiEvent_logp: unable to get _events array");
            goto fail;
        }

        PyArrayIterObject *events_iter;
        events_iter = (PyArrayIterObject *)PyArray_IterNew((PyObject *)multiEvent->events);
 
        if (events_iter == NULL) {
            PyErr_SetString(PyExc_RuntimeError, "multiEvent_logp: unable to create iterator for events array");
            Py_XDECREF(events);
            goto fail;
        }

        while (events_iter->index < events_iter->size) {
            pyEvent *event = (pyEvent *)PyArray_ITER_DATA(events_iter);

            tmp += log(hazard) - hazard * event->event_time;
            remaining_time -= event->event_time;
            PyArray_ITER_NEXT(events_iter);
        }

        Py_DECREF(events_iter);

        if (remaining_time < 0)
            continue;

        result += tmp;
        result += -1.0 * hazard * remaining_time;

    } while (iter_next(iter));

    NpyIter_Deallocate(iter);

    if (PyList_Check(multiEvents_obj))
        Py_DECREF(multiEvents_obj);

    if (PyList_Check(hazards_obj))
        Py_DECREF(hazards_obj);

    if (PyList_Check(periods_obj))
        Py_DECREF(periods_obj);

    Py_DECREF(multiEvents);
    Py_DECREF(hazards);
    Py_DECREF(periods);

    return PyFloat_FromDouble(result);

fail:
    return Py_BuildValue(""); // essentailly a null object
}

static PyObject *
_multiEvent_random(PyObject *self, PyObject *args, PyObject *keywds)
{
    PyArrayObject *op[3];
    PyObject *hazards_obj = NULL, *periods_obj = NULL;
    PyArrayObject *hazards = NULL, *periods = NULL, *out = NULL;

    NpyIter *iter;
    NpyIter_IterNextFunc *iter_next;
    char **dataptrarray;
 
    static char *kwlist[] = {"hazard", "period", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds,  "OO|", kwlist,
                &hazards_obj, &periods_obj)) {
        PyErr_SetString(PyExc_ValueError,
                "multiEvent_random: unable to parse input arguments");

        goto fail;
    }

    // we can't assume that we were passed arrays
    int reqmts = NPY_IN_ARRAY;
    hazards = (PyArrayObject *)PyArray_FROM_OTF(hazards_obj, NPY_DOUBLE, reqmts);
    periods = (PyArrayObject *)PyArray_FROM_OTF(periods_obj, NPY_DOUBLE, reqmts);

    if (hazards == NULL || periods == NULL) {
        goto fail;
    }

    // seed the random number generator
    unsigned long seed = time(NULL);
    init_genrand(seed);

    npy_uint32 flags = NPY_ITER_REFS_OK |  NPY_ITER_MULTI_INDEX;
    npy_uint32 op_flags[3];
    op_flags[0] = NPY_ITER_READONLY;
    op_flags[1] = NPY_ITER_READONLY;
    op_flags[2] = NPY_ITER_WRITEONLY | NPY_ITER_ALLOCATE;

    PyArray_Descr *dtypes[3];
    dtypes[0] = PyArray_DescrFromType(NPY_DOUBLE);
    dtypes[1] = PyArray_DescrFromType(NPY_DOUBLE);
    dtypes[2] = &npy_multiEvent_descr;

    op[0] = hazards;
    op[1] = periods;
    op[2] = NULL; // this will be allocated automagically by the iterator :)

    iter = NpyIter_MultiNew(3, op, flags, NPY_KEEPORDER,
            NPY_SAFE_CASTING, op_flags, dtypes);
    if (iter == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "multiEvent_random: unable to create iterator for input arguments");
        goto fail;
    }

    dataptrarray = NpyIter_GetDataPtrArray(iter);
    iter_next = NpyIter_GetIterNext(iter, NULL);
    if (iter_next == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "multiEvent_random: unable to iterate over input arguments");
        goto fail;
    }

    do {
        int size = 0;
        double *head = NULL, *tail;

        double hazard = *(double *)dataptrarray[0];
        double remaining_time = *(double *)dataptrarray[1];
 
        while (remaining_time > 0) {
 
            double f = genrand_real1();
            double event_time = -1.0 * log(1.0 - f) / hazard;

            if (event_time > remaining_time)
                break;

            tail = (double *)realloc(head, (sizeof(double) * (size + 1)));

            if (tail == NULL) {
                PyErr_SetString(PyExc_MemoryError, "multiEvent_random: unable to allocate memory");
                free(head);
                goto fail;
            }

            head = tail;
            head[size] = event_time;

            remaining_time -= event_time;
            size++;
        }

        // create a new array of _event objects
        npy_intp dims[1] = {size};

        pyMultiEvent *multiEvent = (pyMultiEvent *)dataptrarray[2];
        multiEvent->events = (PyArrayObject *)PyArray_SimpleNewFromDescr(1, dims, &npy_event_descr);


        if (multiEvent->events == NULL) {
            PyErr_SetString(PyExc_RuntimeError, "multiEvent_random: unable to create events array");
            goto fail;
        }
 
        //Py_INCREF(multiEvent->events);

        PyArrayIterObject *events_iter;
        events_iter = (PyArrayIterObject *)PyArray_IterNew((PyObject *)multiEvent->events);

        if (events_iter == NULL) {
            PyErr_SetString(PyExc_RuntimeError, "multiEvent_random: unable to create iterator for events array");
            goto fail;
        }

        // now get the numbers out of our character array

        while(events_iter->index < events_iter->size) {
            pyEvent *event = (pyEvent *)PyArray_ITER_DATA(events_iter);
            event->event_time = head[events_iter->index];
            event->censored = 0;
            PyArray_ITER_NEXT(events_iter);
        }

        free(head);

        Py_DECREF(events_iter);
    } while (iter_next(iter));

    // return the array of _multiEvent objects
    out = NpyIter_GetOperandArray(iter)[2];
    Py_INCREF(out);

    //NpyIter_Deallocate(iter);
    Py_DECREF(iter);

    Py_DECREF(hazards);
    Py_DECREF(periods);

    return (PyObject *)out;

fail:
    Py_XDECREF(out);
    Py_XDECREF(hazards_obj);
    Py_XDECREF(periods_obj);
    Py_XDECREF(hazards);
    Py_XDECREF(periods);

    return Py_BuildValue(""); // essentailly a null object
}


/* Casting functions */

static NPY_INLINE double
pyevent_double(pyEvent pe)
{
    return (double)pe.event_time;	
}

static NPY_INLINE pyEvent
make_pyevent_double(double event_time)
{
	pyEvent pe = {event_time, 0};
	
	return pe;
}

static void
npycast_event_to_double(void *from_, void *to_, npy_intp n,
        void *fromarr, void *toarr)
{
    const pyEvent *from = (pyEvent *)from_;
    double *to = (double *)to_;
    npy_intp i;

    for (i = 0; i < n; i++) {
        pyEvent x = from[i];
        double y = pyevent_double(x);
        to[i] = y;
    }
}

static void
npycast_double_to_event(void *from_, void *to_, npy_intp n,
        void *fromarr, void *toarr)
{
    const double *from = (double *)from_;
    pyEvent *to = (pyEvent *)to_;
    npy_intp i;

    for (i = 0; i < n; i++) {
        double x = from[i];
        pyEvent y = make_pyevent_double(x);
        to[i] = y;
    }
}

/* module definitions */

static PyMethodDef module_methods[] = {
    {"_event_logp", (PyCFunction)_event_logp, METH_KEYWORDS, 
     "Implementation of logp using the array iterator API"},
    {"_event_random", (PyCFunction)_event_random, METH_KEYWORDS, 
     "Implementation of random using the array iterator API"},
    {"_multiEvent_logp", (PyCFunction)_multiEvent_logp, METH_KEYWORDS, 
     "Implementation of multiEvent_logp using the array iterator API"},
    {"_multiEvent_random", (PyCFunction)_multiEvent_random, METH_KEYWORDS, 
     "Implementation of multiEvent_random using the array iterator API"},
    {NULL}  /* Sentinel */
};

PyMODINIT_FUNC
init_Survival(void) {
    /* Initialize numpy */
    import_array();
    if (PyErr_Occurred()) {
        return;
    }
    import_umath();
    if (PyErr_Occurred()) {
        return;
    }
    
    PyObject* numpy_str = PyString_FromString("numpy");
    if (!numpy_str) {
        return;
    }
    
    PyObject* numpy = PyImport_Import(numpy_str);
    Py_DECREF(numpy_str);
    if (!numpy) {
        return;
    }

    /* Can't set this until we import numpy */
    _event_Type.tp_base = &PyGenericArrType_Type;

    /* Initialize _event type object */
    if (PyType_Ready(&_event_Type) < 0) {
        return;
    }

    /* Initialize _event descriptor */
    PyArray_InitArrFuncs(&npy_event_arrfuncs);
    npy_event_arrfuncs.getitem = npy_event_getitem;
    npy_event_arrfuncs.setitem = npy_event_setitem;
    npy_event_arrfuncs.copyswapn = npy_event_copyswapn;
    npy_event_arrfuncs.copyswap = npy_event_copyswap;
    npy_event_arrfuncs.nonzero = npy_event_nonzero;

    npy_event_descr.ob_type = &PyArrayDescr_Type;
    int npy_event = PyArray_RegisterDataType(&npy_event_descr);
    if (npy_event < 0) {
        return;
    }

    /* Support dtype(_event) syntax */
    if (PyDict_SetItemString(_event_Type.tp_dict,"dtype",(PyObject*)&npy_event_descr)<0) {
        return;
    }

    _multiEvent_Type.tp_base = &PyGenericArrType_Type;

    /* Initialize _multiEvent type object */
    if (PyType_Ready(&_multiEvent_Type) < 0) {
        return;
    }

    /* Initialize _multiEvent descriptor */
    PyArray_InitArrFuncs(&npy_multiEvent_arrfuncs);
    npy_multiEvent_arrfuncs.getitem = npy_multiEvent_getitem;
    npy_multiEvent_arrfuncs.setitem = npy_multiEvent_setitem;
    npy_multiEvent_arrfuncs.copyswapn = npy_multiEvent_copyswapn;
    npy_multiEvent_arrfuncs.copyswap = npy_multiEvent_copyswap;
    npy_multiEvent_arrfuncs.nonzero = npy_multiEvent_nonzero;


    npy_multiEvent_descr.ob_type = &PyArrayDescr_Type;
    int npy_multiEvent = PyArray_RegisterDataType(&npy_multiEvent_descr);
    if (npy_multiEvent < 0) {
        return;
    }

    /* Support dtype(_multiEvent) syntax */
    if (PyDict_SetItemString(_multiEvent_Type.tp_dict,"dtype",(PyObject*)&npy_multiEvent_descr)<0) {
        return;
    }

    /* Create module */
    PyObject* m = Py_InitModule3("_Survival", module_methods,
        "single and multiple event objects for survival analysis");
    if (!m) {
        return;
    }

    /* Add new types */
    Py_INCREF(&_event_Type);
    PyModule_AddObject(m,"_event",(PyObject*)&_event_Type);

    Py_INCREF(&_multiEvent_Type);
    PyModule_AddObject(m,"_multiEvent",(PyObject*)&_multiEvent_Type);

    /* Register casting functions */
    PyArray_Descr* pyevent_to_double_descr = &npy_event_descr;
    if (PyArray_RegisterCastFunc(pyevent_to_double_descr,(NPY_DOUBLE),npycast_event_to_double) < 0)
        return;

    PyArray_Descr* double_to_pyevent_descr = PyArray_DescrFromType(NPY_DOUBLE);
    if (PyArray_RegisterCastFunc(double_to_pyevent_descr,(npy_event),npycast_double_to_event) < 0)
        return;

    if (PyArray_RegisterCanCast(pyevent_to_double_descr,(NPY_DOUBLE),NPY_NOSCALAR) < 0)
        return;

    if (PyArray_RegisterCanCast(double_to_pyevent_descr,(npy_event),NPY_NOSCALAR) < 0)
        return;

    /* ufunc defninitions */
    PyObject *_logp;
    
    _logp = PyUFunc_FromFuncAndDataAndSignature(_logp_functions, _logp_data, _logp_signatures, 1,
                                    3, 1, PyUFunc_None, "_logp",
                                    "implentation of logp using the generalized ufunc API\n"\
                                    "     \"(i),(),()->()\" \n",
                                    0, _logp_signature);

    if (!_logp)
        return;

    // grab the dictionary of things defined in the module
    PyObject *d = PyModule_GetDict(m);

    // add _logp as a method of the module
    PyModule_AddObject(m,"_logp",(PyObject*)_logp);

    PyDict_SetItemString(d, "_logp", _logp);

}
