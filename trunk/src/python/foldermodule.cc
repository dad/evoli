/*
This file is part of the evoli project.
Copyright (C) 2006 Allan Drummond <dadrummond@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1
*/


#include <Python.h> // needs to be first include
#include "protein.hh"
#include "folder.hh"
#include "compact-lattice-folder.hh"

static PyObject *FolderErrorObject;

/* ----------------------------------------------------- */

static char folder_fold__doc__[] =
""
;

static Folder* folder = NULL;
static PyObject *
folder_fold(PyObject *self /* Not used */, PyObject *args)
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
	auto_ptr<FoldInfo> folding_data( folder->fold(p) );
	return Py_BuildValue("if", folding_data->getStructure(), folding_data->getDeltaG());
}

static char folder_getEnergy__doc__[] =
""
;

static PyObject *
folder_getEnergy(PyObject *self /* Not used */, PyObject *args)
{
	const char *protein_sequence;
	int sid;
	if (!folder) {
		PyErr_SetString(FolderErrorObject, "uninitialized folder: call 'folder.init(...)'");
		return NULL;
	}
	if (!PyArg_ParseTuple(args, "si", &protein_sequence, &sid)) {
		return NULL;
	}
	Protein p(protein_sequence);
	double energy = folder->getEnergy(p, sid);
	return Py_BuildValue("f", energy);
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
	Py_INCREF(Py_None);
	return Py_None;
}

/* List of methods defined in the module */

static struct PyMethodDef folder_methods[] = {
	{"fold",	(PyCFunction)folder_fold,	METH_VARARGS, folder_fold__doc__},
	{"init",	(PyCFunction)folder_init,	METH_VARARGS, folder_init__doc__},
	{"getEnergy",	(PyCFunction)folder_getEnergy,	METH_VARARGS, folder_getEnergy__doc__}, 
	{NULL,	 (PyCFunction)NULL, 0, NULL}		/* sentinel */
};


/* Initialization function for the module (*must* be called initfolder) */

static char folder_module_documentation[] = 
""
;

PyMODINIT_FUNC
initfolder(void)
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

