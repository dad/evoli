#include <Python.h>
#include <iostream>
#include "protein.hh"
#include "folder.hh"
#include "compact-lattice-folder.hh"

#include "Python.h"

static PyObject *FolderErrorObject;

/* ----------------------------------------------------- */

static char folder_foldProtein__doc__[] =
""
;

static CompactLatticeFolder* folder = NULL;
static PyObject *
folder_foldProtein(PyObject *self /* Not used */, PyObject *args)
{
    const char *protein_sequence;
	if (!folder) {
		PyErr_SetString(FolderErrorObject, "uninitialized folder: call 'folder.init(n)' with n = side-length of lattice protein");
		return NULL;
	}
    if (!PyArg_ParseTuple(args, "s", &protein_sequence)) {
        return NULL;
	}
	Protein p(protein_sequence);
	FoldInfo folding_data = folder->fold(p);
    return Py_BuildValue("if", folding_data.getStructure(), folding_data.getFreeEnergy());
}

static char folder_init__doc__[] =
""
;

static PyObject *
folder_init(PyObject *self /* Not used */, PyObject *args)
{
	int side_length;
    if (!PyArg_ParseTuple(args, "i", &side_length))
        return NULL;
	folder = new CompactLatticeFolder(side_length);
	folder->enumerateStructures();
	Py_INCREF(Py_None);
	return Py_None;
}

/* List of methods defined in the module */

static struct PyMethodDef folder_methods[] = {
	{"foldProtein",	(PyCFunction)folder_foldProtein,	METH_VARARGS,	folder_foldProtein__doc__},
 {"init",	(PyCFunction)folder_init,	METH_VARARGS,	folder_init__doc__},
 
	{NULL,	 (PyCFunction)NULL, 0, NULL}		/* sentinel */
};


/* Initialization function for the module (*must* be called initfolder) */

static char folder_module_documentation[] = 
""
;

PyMODINIT_FUNC
initfolder()
{
	PyObject *m, *d;

	/* Create the module and add the functions */
	m = Py_InitModule4("folder", folder_methods,
		folder_module_documentation,
		(PyObject*)NULL,PYTHON_API_VERSION);

	/* Add some symbolic constants to the module */
	d = PyModule_GetDict(m);
	FolderErrorObject = PyString_FromString("folder.error");
	PyDict_SetItemString(d, "error", FolderErrorObject);

	/* XXXX Add constants here */
	
	/* Check for errors */
	if (PyErr_Occurred())
		Py_FatalError("can't initialize module folder");
}

