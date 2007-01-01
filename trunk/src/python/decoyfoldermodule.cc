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
#include <iostream>
#include <fstream>
#include "random.hh"
#include "protein.hh"
#include "folder.hh"
#include "decoy-contact-folder.hh"
#include "gene-util.hh"


static PyObject *FolderErrorObject;

/* ----------------------------------------------------- */

static DecoyContactFolder* folder = NULL;

static char decoyfolder_fold__doc__[] =
""
;

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
	auto_ptr<FoldInfo> fd( folder->fold(p) );
	//cout << folding_data.getStructure() << " " << folding_data.getDeltaG() << endl << p << endl;
	return Py_BuildValue("if", fd->getStructure(), fd->getDeltaG());
}

static char decoyfolder_foldStats__doc__[] =
""
;

static PyObject *
decoyfolder_foldStats(PyObject *self /* Not used */, PyObject *args)
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
	auto_ptr<DecoyFoldInfo> folding_data( dynamic_cast<DecoyFoldInfo*>( folder->fold(p) ) );
	//cout << folding_data.getStructure() << " " << folding_data.getDeltaG() << endl << p << endl;
	return Py_BuildValue("iffff", folding_data->getStructure(), folding_data->getDeltaG(), folding_data->getUnfoldedFreeEnergyMean(), folding_data->getUnfoldedFreeEnergyVariance(), folding_data->getMinEnergy());
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

static char decoyfolder_getSequenceForStructure__doc__[] =
""
;

// finds a random sequence with folding energy smaller than cutoff and structure given by struct_id
static PyObject *
decoyfolder_getSequenceForStructure( PyObject *self, PyObject *args)
{
	int protein_length;
	double free_energy_cutoff;
	int struct_id;

	if (!folder) {
		PyErr_SetString(FolderErrorObject, "uninitialized decoyfolder: call 'decoyfolder.init()'");
		return NULL;
	}
	if (!PyArg_ParseTuple(args, "idi", &protein_length, &free_energy_cutoff, &struct_id )) {
		return NULL;
	}
	Gene g = GeneUtil::getSequenceForStructure(*folder, 3*protein_length, free_energy_cutoff, struct_id);
	//cout << 3*protein_length << " " << free_energy_cutoff << " " << struct_id << " " << g << endl;
	Protein p = g.translate();
	
	return Py_BuildValue("s", p.toString().c_str());
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
	length = Random::runif();
	if (!PyArg_ParseTuple(args, "idss", &length, &nconf, &map_file, &map_dir))
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
	{"foldStats",	(PyCFunction)decoyfolder_foldStats,	METH_VARARGS,	decoyfolder_foldStats__doc__},
	{"init",	(PyCFunction)decoyfolder_init,	METH_VARARGS,	decoyfolder_init__doc__},
	{"getEnergy",	(PyCFunction)decoyfolder_getEnergy,	METH_VARARGS,	decoyfolder_getEnergy__doc__},
	{"getSequenceForStructure",	(PyCFunction)decoyfolder_getSequenceForStructure,	METH_VARARGS,	decoyfolder_getSequenceForStructure__doc__},
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

