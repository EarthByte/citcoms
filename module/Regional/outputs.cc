// -*- C++ -*-
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  <LicenseText>
// 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

#include <portinfo>
#include <Python.h>
#include <cstdio>

#include "output_bindings.h"

extern "C" {
#include "mpi.h"
#include "global_defs.h"
#include "citcom_init.h"
#include "output.h"

}

char pyRegional_output_init__doc__[] = "";
char pyRegional_output_init__name__[] = "output_init";

PyObject * pyRegional_output_init(PyObject *self, PyObject *args)
{
    char *filename;
    FILE *fp;
    PyObject *obj;

    if (!PyArg_ParseTuple(args, "s", &filename))
        return NULL;

    fp = output_init(filename);
    obj = PyCObject_FromVoidPtr((void *)fp, NULL);

    return Py_BuildValue("O", obj);
}


char pyRegional_output_close__doc__[] = "";
char pyRegional_output_close__name__[] = "output_close";

PyObject * pyRegional_output_close(PyObject *self, PyObject *args)
{
    PyObject *obj;
    FILE *fp;

    if (!PyArg_ParseTuple(args, "O", &obj))
        return NULL;

    fp = static_cast<FILE*> (PyCObject_AsVoidPtr(obj));
    output_close(fp);

    Py_INCREF(Py_None);
    return Py_None;
}


char pyRegional_output_coord_header__doc__[] = "";
char pyRegional_output_coord_header__name__[] = "output_coord_header";

PyObject * pyRegional_output_coord_header(PyObject *self, PyObject *args)
{
    PyObject *obj;
    FILE *fp;

    if (!PyArg_ParseTuple(args, "O", &obj))
        return NULL;

    fp = static_cast<FILE*> (PyCObject_AsVoidPtr(obj));
    output_coord_header(E, fp);

    Py_INCREF(Py_None);
    return Py_None;
}


char pyRegional_output_coord__doc__[] = "";
char pyRegional_output_coord__name__[] = "output_coord";

PyObject * pyRegional_output_coord(PyObject *self, PyObject *args)
{
    PyObject *obj;
    FILE *fp;

    if (!PyArg_ParseTuple(args, "O", &obj))
        return NULL;

    fp = static_cast<FILE*> (PyCObject_AsVoidPtr(obj));
    output_coord(E, fp);

    Py_INCREF(Py_None);
    return Py_None;
}


char pyRegional_output_velo_header__doc__[] = "";
char pyRegional_output_velo_header__name__[] = "output_velo_header";

PyObject * pyRegional_output_velo_header(PyObject *self, PyObject *args)
{
    PyObject *obj;
    FILE *fp;
    int step;

    if (!PyArg_ParseTuple(args, "Oi", &obj, &step))
        return NULL;

    fp = static_cast<FILE*> (PyCObject_AsVoidPtr(obj));
    output_velo_header(E, fp, step);

    Py_INCREF(Py_None);
    return Py_None;
}


char pyRegional_output_velo__doc__[] = "";
char pyRegional_output_velo__name__[] = "output_velo";

PyObject * pyRegional_output_velo(PyObject *self, PyObject *args)
{
    PyObject *obj;
    FILE *fp;

    if (!PyArg_ParseTuple(args, "O", &obj))
        return NULL;

    fp = static_cast<FILE*> (PyCObject_AsVoidPtr(obj));
    output_velo(E, fp);

    Py_INCREF(Py_None);
    return Py_None;
}



// version
// $Id: outputs.cc,v 1.1 2003/05/21 21:28:00 tan2 Exp $

// End of file
