#include <Python.h>
#include <iostream>
#include <fstream>
#include "protein.hh"
#include "folder.hh"
#include "decoy-contact-folder.hh"

#include "Python.h"

static PyObject *FolderErrorObject;

/* ----------------------------------------------------- */

static char decoyfolder_fold__doc__[] =
""
;

static Folder* folder = NULL;
static PyObject *
decoyfolder_fold(PyObject *self /* Not used */, PyObject *args)
{
    const char *protein_sequence;
	if (!folder) {
		PyErr_SetString(FolderErrorObject, "uninitialized decoyfolder: call 'decoyfolder.init()'");
		return NULL;
	}
    if (!PyArg_ParseTuple(args, "s", &protein_sequence)) {
        return NULL;
	}
	Protein p(protein_sequence);
	FoldInfo folding_data = folder->fold(p);
    return Py_BuildValue("if", folding_data.getStructure(), folding_data.getFreeEnergy());
}

static char decoyfolder_getEnergy__doc__[] =
""
;

static PyObject *
decoyfolder_getEnergy(PyObject *self /* Not used */, PyObject *args)
{
    const char *protein_sequence;
	int sid;
	if (!folder) {
		PyErr_SetString(FolderErrorObject, "uninitialized decoyfolder: call 'decoyfolder.init()'");
		return NULL;
	}
    if (!PyArg_ParseTuple(args, "si", &protein_sequence, &sid)) {
        return NULL;
	}
	Protein p(protein_sequence);
	double energy = folder->getEnergy(p, sid);
    return Py_BuildValue("f", energy);
}

static char decoyfolder_init__doc__[] =
""
;

static PyObject *
decoyfolder_init(PyObject *self /* Not used */, PyObject *args)
{
	const char *map_file;
	const char *map_dir;
	double nconf;
	int length;
    if (!PyArg_ParseTuple(args, "ifss", &length, &nconf, &map_file, &map_dir))
        return NULL;
	string path = map_file;
	ifstream fin(path.c_str());
	if (!fin.good()) {
		string err = string("couldn't initialize decoyfolder: bad file ") + path;
		PyErr_SetString(FolderErrorObject, err.c_str());
		return NULL;
	}		
	folder = new DecoyContactFolder(length, nconf, fin, string(map_dir));
	fin.close();
	if (!folder->good()) {
		string err = string("couldn't initialize decoyfolder: bad data from ") + path;
		PyErr_SetString(FolderErrorObject, err.c_str());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}

/* List of methods defined in the module */

static struct PyMethodDef decoyfolder_methods[] = {
	{"fold",	(PyCFunction)decoyfolder_fold,	METH_VARARGS,	decoyfolder_fold__doc__},
	{"init",	(PyCFunction)decoyfolder_init,	METH_VARARGS,	decoyfolder_init__doc__},
	{"getEnergy",	(PyCFunction)decoyfolder_getEnergy,	METH_VARARGS,	decoyfolder_getEnergy__doc__},
	{NULL,	 (PyCFunction)NULL, 0, NULL}		/* sentinel */
};


/* Initialization function for the module (*must* be called initdecoyfolder) */

static char decoyfolder_module_documentation[] = 
""
;

PyMODINIT_FUNC
initdecoyfolder()
{
	PyObject *m, *d;

	/* Create the module and add the functions */
	m = Py_InitModule4("decoyfolder", decoyfolder_methods,
		decoyfolder_module_documentation,
		(PyObject*)NULL,PYTHON_API_VERSION);

	/* Add some symbolic constants to the module */
	d = PyModule_GetDict(m);
	FolderErrorObject = PyString_FromString("decoyfolder.error");
	PyDict_SetItemString(d, "error", FolderErrorObject);

	/* XXXX Add constants here */
	
	/* Check for errors */
	if (PyErr_Occurred())
		Py_FatalError("can't initialize module decoyfolder");
}

