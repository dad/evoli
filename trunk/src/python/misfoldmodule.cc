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
#include "tools.hh"
#include "protein.hh"
#include "folder.hh"
#include "compact-lattice-folder.hh"
#include "fitness-evaluator.hh"

static PyObject *MisfoldErrorObject;
static ErrorproneTranslation* fitness_evaluator = NULL;
static Folder* folder = NULL;

#if 1

/* ----------------------------------------------------- */

static char misfold_calcOutcomes__doc__[] =
""
;

static PyObject *
misfold_calcOutcomes(PyObject *self /* Not used */, PyObject *args)
{
	const char *gene_sequence;
	if (!fitness_evaluator) {
		PyErr_SetString(MisfoldErrorObject, "uninitialized misfold object: call 'misfold.init(...)' first");
		return NULL;
	}
	if (!PyArg_ParseTuple(args, "s", &gene_sequence)) {
		return NULL;
	}
	CodingDNA g(gene_sequence);
	double cffold, cfacc, cfrob, cftrunc;
	fitness_evaluator->calcOutcomes(g, cfacc, cfrob, cftrunc, cffold);
	//cout << cfacc << "\t" << cfrob << "\t" << cftrunc << "\t" << cffold << endl;
	return Py_BuildValue("dddd", cfacc, cfrob, cftrunc, cffold);
}

static char misfold_countOutcomes__doc__[] =
""
;

static PyObject *
misfold_countOutcomes(PyObject *self /* Not used */, PyObject *args)
{
	const char *gene_sequence;
	int num_to_fold;
	if (!fitness_evaluator) {
		PyErr_SetString(MisfoldErrorObject, "uninitialized misfold object: call 'misfold.init(...)' first");
		return NULL;
	}
	if (!PyArg_ParseTuple(args, "si", &gene_sequence, &num_to_fold)) {
		return NULL;
	}
	CodingDNA g(gene_sequence);
	int cffold, cfacc, cfrob, cftrunc;
	fitness_evaluator->countOutcomes(g, num_to_fold, cfacc, cfrob, cftrunc, cffold);
	return Py_BuildValue("iiii", cfacc, cfrob, cftrunc, cffold);
}

static char misfold_countOutcomeFractions__doc__[] =
""
;

static PyObject *
misfold_countOutcomeFractions(PyObject *self /* Not used */, PyObject *args)
{
	const char *gene_sequence;
	int num_to_fold;
	if (!fitness_evaluator) {
		PyErr_SetString(MisfoldErrorObject, "uninitialized misfold object: call 'misfold.init(...)' first");
		return NULL;
	}
	if (!PyArg_ParseTuple(args, "si", &gene_sequence, &num_to_fold)) {
		return NULL;
	}
	CodingDNA g(gene_sequence);
	int nfold, nacc, nrob, ntrunc;
	fitness_evaluator->countOutcomes(g, num_to_fold, nacc, nrob, ntrunc, nfold);
	double ntf = static_cast<double>(num_to_fold);
	return Py_BuildValue("dddd", nacc/ntf, nrob/(ntf-nacc), ntrunc/ntf, nfold/ntf);
}

static char misfold_setErrorRate__doc__[] =
""
;

static PyObject *
misfold_setErrorRate(PyObject *self /* Not used */, PyObject *args)
{
	double error_rate;
	if (!fitness_evaluator) {
		PyErr_SetString(MisfoldErrorObject, "uninitialized misfold object: call 'misfold.init(...)' first");
		return NULL;
	}
	if (!PyArg_ParseTuple(args, "d", &error_rate)) {
		return NULL;
	}
	fitness_evaluator->setErrorRate(error_rate);
	return Py_None;
}

static char misfold_getErrorRate__doc__[] =
""
;

static PyObject *
misfold_getErrorRate(PyObject *self /* Not used */, PyObject *args)
{
	double error_rate;
	if (!fitness_evaluator) {
		PyErr_SetString(MisfoldErrorObject, "uninitialized misfold object: call 'misfold.init(...)' first");
		return NULL;
	}
	error_rate = fitness_evaluator->getErrorRate();
	return Py_BuildValue("d", error_rate);
}

#endif

static char misfold_init__doc__[] =
""
;

static PyObject *
misfold_init(PyObject *self /* Not used */, PyObject *args)
{
	PyObject *folder_module;
	int target_structure_id, random_seed, protein_length;
	double max_free_energy, ca_cost, target_fraction_accurate;

	if (!PyArg_ParseTuple(args, "Oiidddi", &folder_module, &protein_length, &target_structure_id, &max_free_energy, &ca_cost, &target_fraction_accurate, &random_seed)) {
		PyErr_SetString(MisfoldErrorObject, "bad arguments to decoyfolder.init");
		return NULL;
	}
	if (!PyModule_Check(folder_module)) {
		PyErr_SetString(MisfoldErrorObject, "first argument not a Python module");
		return NULL;
	}
	Random::seed(random_seed);
	//cout << side_length << " " << target_structure_id << " " << max_free_energy << " " << ca_cost << " " << target_fraction_accurate << endl;
	
	//folder = new CompactLatticeFolder(5);
	PyObject *c_folder = PyObject_GetAttrString(folder_module, "_C_FOLDER");
	//cout << "misfold c_folder = " << c_folder << endl;
	if (PyCObject_Check(c_folder)){
		folder = (Folder*)PyCObject_AsVoidPtr(c_folder);
		//cout << "misfold c_folder as void* = " << c_folder << endl;
	}
	
	//ErrorproneTranslation(Folder *protein_folder, const int protein_length, const StructureID protein_structure_ID, const double max_free_energy, const double tr_cost, const double ca_cost, const double target_fraction_accurate );

	if (fitness_evaluator) {
		//cout << "misfold: cleaning up old FE" << endl;
		delete fitness_evaluator;
	}
	fitness_evaluator = new ErrorproneTranslation(folder, protein_length, target_structure_id, max_free_energy, 1.0, ca_cost, target_fraction_accurate);
	//cout << "misfold: made new FE" << endl;
	Py_INCREF(Py_None);
	return Py_None;
}

/* List of methods defined in the module */

static struct PyMethodDef misfold_methods[] = {
	{"calcOutcomes",	(PyCFunction)misfold_calcOutcomes,	METH_VARARGS, misfold_calcOutcomes__doc__},
	{"countOutcomes",	(PyCFunction)misfold_countOutcomes,	METH_VARARGS, misfold_countOutcomes__doc__},
	{"countOutcomeFractions",	(PyCFunction)misfold_countOutcomeFractions,	METH_VARARGS, misfold_countOutcomeFractions__doc__},
	{"setErrorRate",	(PyCFunction)misfold_setErrorRate,	METH_VARARGS, misfold_setErrorRate__doc__},
	{"getErrorRate",	(PyCFunction)misfold_getErrorRate,	METH_VARARGS, misfold_getErrorRate__doc__},
	{"init",	(PyCFunction)misfold_init,	METH_VARARGS, misfold_init__doc__},
	{NULL,	 (PyCFunction)NULL, 0, NULL}		/* sentinel */
};


/* Initialization function for the module (*must* be called initmisfold) */

static char misfold_module_documentation[] =
""
;

PyMODINIT_FUNC
initmisfold(void)
{
	PyObject *m, *d;

	/* Create the module and add the functions */
	m = Py_InitModule4("misfold", misfold_methods,
		misfold_module_documentation,
		(PyObject*)NULL,PYTHON_API_VERSION);

	/* Add some symbolic constants to the module */
	d = PyModule_GetDict(m);
	MisfoldErrorObject = PyString_FromString("misfold.error");
	PyDict_SetItemString(d, "error", MisfoldErrorObject);

	/* XXXX Add constants here */

	/* Check for errors */
	if (PyErr_Occurred())
		Py_FatalError("can't initialize module misfold");
}
