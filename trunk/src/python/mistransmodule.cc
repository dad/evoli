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
#include "fitness-evaluator.hh"

static PyObject *MistransErrorObject;
static ErrorproneTranslation* fitness_evaluator = NULL;
static Folder* folder = NULL;

#if 1

/* ----------------------------------------------------- */

static char mistrans_getTranslatedProteins__doc__[] =
""
;

static PyObject *
mistrans_getTranslatedProteins(PyObject *self /* Not used */, PyObject *args)
{
	const char *gene_sequence;
	uint num_to_translate;
	bool only_mistranslated;

	if (!fitness_evaluator) {
		PyErr_SetString(MistransErrorObject, "uninitialized mistrans object: call 'mistrans.init(...)' first");
		return NULL;
	}
	if (!PyArg_ParseTuple(args, "sii", &gene_sequence, &num_to_translate, &only_mistranslated)) {
		return NULL;
	}
	vector<Protein> proteins;
	CodingDNA g(gene_sequence);
	fitness_evaluator->translationOutcomes(g, num_to_translate, only_mistranslated, proteins);
	//cout << cfacc << "\t" << cfrob << "\t" << cftrunc << "\t" << cffold << endl;
	PyObject* prot_list = PyList_New(0);
	for (vector<Protein>::iterator it=proteins.begin(); it!=proteins.end(); it++) {
		PyObject* str = PyString_FromString((*it).toString().c_str());
		PyList_Append(prot_list, str);
	}
	return Py_BuildValue("O", prot_list);
}


#endif

static char mistrans_init__doc__[] =
""
;

static PyObject *
mistrans_init(PyObject *self /* Not used */, PyObject *args)
{
	PyObject *folder_module;
	int target_structure_id, random_seed, protein_length;
	double max_free_energy, ca_cost, target_fraction_accurate;
	//mistrans.init(decoyfolder, target_sid, max_dg, ca_cost, target_frac_accurate, 111)
	if (!PyArg_ParseTuple(args, "Oiidddi", &folder_module, &protein_length, &target_structure_id, &max_free_energy, &ca_cost, &target_fraction_accurate, &random_seed)) {
		PyErr_SetString(MistransErrorObject, "bad arguments to mistrans.init");
		return NULL;
	}
	if (!PyModule_Check(folder_module)) {
		PyErr_SetString(MistransErrorObject, "first argument not a Python module");
		return NULL;
	}
	Random::seed(random_seed);
	// Get Folder from module that's been passed in.
	PyObject *c_folder = PyObject_GetAttrString(folder_module, "_C_FOLDER");
	if (PyCObject_Check(c_folder)){
		folder = (Folder*)PyCObject_AsVoidPtr(c_folder);
		//cout << "mistrans c_folder as void* = " << c_folder << endl;
	}
	else {
		PyErr_SetString(MistransErrorObject, "mistrans.init: module does not have a Folder object (_C_FOLDER)");
		return NULL;
	}

	//ErrorproneTranslation(Folder *protein_folder, const int protein_length, const StructureID protein_structure_ID, const double max_free_energy, const double tr_cost, const double ca_cost, const double target_fraction_accurate );

	if (fitness_evaluator) {
		//cout << "mistrans: cleaning up old FE" << endl;
		delete fitness_evaluator;
	}
	fitness_evaluator = new ErrorproneTranslation(folder, protein_length, target_structure_id, max_free_energy, 1.0, ca_cost, target_fraction_accurate);
	//cout << "mistrans: made new FE" << endl;
	Py_INCREF(Py_None);
	return Py_None;
}

/* List of methods defined in the module */

static struct PyMethodDef mistrans_methods[] = {
	{"getTranslatedProteins",	(PyCFunction)mistrans_getTranslatedProteins,	METH_VARARGS, mistrans_getTranslatedProteins__doc__},
	{"init",	(PyCFunction)mistrans_init,	METH_VARARGS, mistrans_init__doc__},
	{NULL,	 (PyCFunction)NULL, 0, NULL}		/* sentinel */
};


/* Initialization function for the module (*must* be called initmistrans) */

static char mistrans_module_documentation[] =
""
;

PyMODINIT_FUNC
initmistrans(void)
{
	PyObject *m, *d;

	/* Create the module and add the functions */
	m = Py_InitModule4("mistrans", mistrans_methods,
		mistrans_module_documentation,
		(PyObject*)NULL,PYTHON_API_VERSION);

	/* Add some symbolic constants to the module */
	d = PyModule_GetDict(m);
	MistransErrorObject = PyString_FromString("mistrans.error");
	PyDict_SetItemString(d, "error", MistransErrorObject);

	/* XXXX Add constants here */

	/* Check for errors */
	if (PyErr_Occurred())
		Py_FatalError("can't initialize module mistrans");
}
