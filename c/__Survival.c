#include <Python.h>
#include "structmember.h"
#include <numpy/arrayobject.h>
#include <numpy/ufuncobject.h>

typedef long bool; // a true waste of space, but Python wants to do C long <==> Python int conversion 
enum { false, true };

/*
 *
 * Represents a single PyEvent that can happen once.  The typical example is death.
 *  time - (float) The time to death or censorship
 *  censored - (bool) True if the PyEvent was right censored; False otherwise.
 *
 */
typedef struct {
    PyObject_HEAD
    double time;
    bool censored; 
} PyEvent;

typedef struct {
    PyObject_HEAD
    PyEvent e;	
} __Event;

static NPY_INLINE int
PyEvent_nonzero(PyEvent pe) {
    return pe.time != 0.0;
}

static PyTypeObject __Event_Type; 

static NPY_INLINE int
__Event_Check(PyObject* object) {
    return PyObject_IsInstance(object,(PyObject*)&__Event_Type);
}

static PyObject*
__Event_FromPyEvent(PyEvent pe) {
    __Event* e = (__Event*)__Event_Type.tp_alloc(&__Event_Type,0);
    if (e) {
        e->e = pe;
    }
    return (PyObject*)e;
}

static void
__Event_dealloc(__Event* self)
{
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
__Event_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyEvent pe;
    
    pe.time = 0.0;
    pe.censored = false;
    
    return __Event_FromPyEvent(pe);
}



static int
__Event_nonzero(PyObject* self) {
    PyEvent pe = ((__Event*)self)->e;
    return PyEvent_nonzero(pe);
}


static PyObject*
__Event_get_time(__Event* self, void* closure) {
    return PyFloat_FromDouble(self->e.time);
}

static PyObject*
__Event_get_censored(__Event* self, void* closure) {
    return PyInt_FromLong(self->e.censored);
}

static int
__Event_set_time(__Event *self, PyObject *value, void *closure)
{
  if (value == NULL) {
    //PyErr_SetString(PyExc_TypeError, "Cannot delete the last attribute");
    return -1;
  }
  
  if (! PyFloat_Check(value)) {
    //PyErr_SetString(PyExc_TypeError, "The last attribute value must be a string");
    return -1;
  }
      
  self->e.time = PyFloat_AS_DOUBLE(value);    

  return 0;
}

static int
__Event_set_censored(__Event *self, PyObject *value, void *closure)
{
  if (value == NULL) {
    //PyErr_SetString(PyExc_TypeError, "Cannot delete the last attribute");
    return -1;
  }
  
  if (! PyInt_Check(value)) {
    //PyErr_SetString(PyExc_TypeError, "The last attribute value must be a string");
    return -1;
  }
      
  self->e.censored = PyInt_AsLong(value);    

  return 0;
}


static PyGetSetDef __Event_getset[] = {
    {"time",(getter)__Event_get_time,(setter)__Event_set_time,"time", NULL},
    {"censored",(getter)__Event_get_censored,(setter)__Event_set_censored,"censored", NULL},
    {NULL} /* sentinel */
};

static PyObject*
__Event_time(__Event* self, void* closure) {
    return PyFloat_FromDouble(self->e.time);
}

static PyObject*
__Event_is_censored(__Event* self, void* closure) {
    return PyInt_FromLong(self->e.censored);
}


static int
__Event_init(__Event *self, PyObject *args, PyObject *kwds)
{
    
    PyEvent pe;
    long censored;
    
    static char *kwlist[] = {"time", "censored", NULL};

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "|di", kwlist, 
                                      &pe.time, 
                                      &censored))
        return -1;         

    // enforce C++ style semantics for a boolean   
    if (censored < 0)
       pe.censored = 0;
    else if (censored > 1)
       pe.censored = 1;
    
    self->e = pe;


    return 0;
}

static PyMemberDef __Event_members[] = {
    {"event", T_OBJECT_EX, offsetof(__Event, e), 0,
     "event"},
    {NULL}  /* Sentinel */
};

// modify so that arguments update the underlying fields
static PyMethodDef __Event_methods[] = {
    {"time", (PyCFunction)__Event_time, METH_NOARGS,
     "Return the time"
    },
    {"is_censored", (PyCFunction)__Event_is_censored, METH_NOARGS,
     "Return the censored flag"
    },
    {NULL}  /* Sentinel */
};


static PyTypeObject __Event_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "__Survival.__Event",       /*tp_name*/
    sizeof(__Event),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)__Event_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
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
    __Event_methods,                     /* tp_methods */
    __Event_members,             /* tp_members */
    __Event_getset,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)__Event_init,      /* tp_init */
    0,                         /* tp_alloc */
    __Event_new,                 /* tp_new */
};

static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};


/* Numpy support */

static PyObject*
npy__Event_getitem(void* data, void* arr) {
    PyEvent pe;
    memcpy(&pe ,data,sizeof(PyEvent));
    return __Event_FromPyEvent(pe);
}

static int
npy__Event_setitem(PyObject* item, void* data, void* arr) {
    PyEvent pe;
    if (__Event_Check(item)) {
        pe = ((__Event*)item)->e;
    }
    else {
    	// sorry no graceful error handling or casting...
    	return -1;
    }
    memcpy(data,&pe,sizeof(PyEvent));
    return 0;
}

static NPY_INLINE void
byteswapint(long* x) {
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
npy__Event_copyswapn(void* dst_, npy_intp dstride, void* src_, npy_intp sstride, npy_intp n, int swap, void* arr) {
    char *dst = (char*)dst_, *src = (char*)src_;
    if (!src) {
        return;
    }
    npy_intp i;
    if (swap) {
        for (i = 0; i < n; i++) {
            PyEvent* pe = (PyEvent*)(dst+dstride*i);
            memcpy(pe,src+sstride*i,sizeof(PyEvent));
            byteswapdouble(&pe->time);
            byteswapint(&pe->censored);
        }
    }
    else if (dstride==sizeof(PyEvent) && sstride==sizeof(PyEvent)) {
        memcpy(dst,src,n*sizeof(PyEvent));
    }
    else {
        for (i = 0; i < n; i++) {
            memcpy(dst+dstride*i,src+sstride*i,sizeof(PyEvent));
        }
    }
}

static void
npy__Event_copyswap(void* dst, void* src, int swap, void* arr) {
    if (!src) {
        return;
    }
    PyEvent* pe = (PyEvent*)dst;
    memcpy(pe,src,sizeof(PyEvent));
    if (swap) {
        byteswapdouble(&pe->time);
        byteswapint(&pe->censored);
    }
}

static npy_bool
npy__Event_nonzero(void* data, void* arr) {
    PyEvent pe;
    memcpy(&pe,data,sizeof(pe));
    return PyEvent_nonzero(pe)?NPY_TRUE:NPY_FALSE;
}

static PyArray_ArrFuncs npyEvent_arrfuncs;

typedef struct { char c; PyEvent pe; } align_test;

PyArray_Descr npyEvent_descr = {
    PyObject_HEAD_INIT(0)
    &__Event_Type,       /* typeobj */
    'V',                    /* kind */
    'e',                    /* type */
    '=',                    /* byteorder */
    NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM, /* hasobject */
    0,                      /* type_num */
    sizeof(PyEvent),       /* elsize */
    offsetof(align_test,pe), /* alignment */
    0,                      /* subarray */
    0,                      /* fields */
    0,                      /* names */
    &npyEvent_arrfuncs,  /* f */
};


PyMODINIT_FUNC
init__Survival(void) {
    /* Initialize numpy */
    import_array();
    if (PyErr_Occurred()) {
        return;
    }
/*    import_umath();
    if (PyErr_Occurred()) {
        return;
    }*/
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
    __Event_Type.tp_base = &PyGenericArrType_Type;

    /* Initialize rational type object */
    if (PyType_Ready(&__Event_Type) < 0) {
        return;
    }

    /* Initialize rational descriptor */
    PyArray_InitArrFuncs(&npyEvent_arrfuncs);
    npyEvent_arrfuncs.getitem = npy__Event_getitem;
    npyEvent_arrfuncs.setitem = npy__Event_setitem;
    npyEvent_arrfuncs.copyswapn = npy__Event_copyswapn;
    npyEvent_arrfuncs.copyswap = npy__Event_copyswap;
    npyEvent_arrfuncs.nonzero = npy__Event_nonzero;


    npyEvent_descr.ob_type = &PyArrayDescr_Type;
    int npy_Event = PyArray_RegisterDataType(&npyEvent_descr);
    if (npy_Event<0) {
        return;
    }

    /* Support dtype(__Event) syntax */
    if (PyDict_SetItemString(__Event_Type.tp_dict,"dtype",(PyObject*)&npyEvent_descr)<0) {
        return;
    }


    /* Create module */
    PyObject* m = Py_InitModule3("__Survival", module_methods,
        "Fixed precision rational numbers, including numpy support");
    if (!m) {
        return;
    }

    /* Add rational type */
    Py_INCREF(&__Event_Type);
    PyModule_AddObject(m,"__Event",(PyObject*)&__Event_Type);

}
