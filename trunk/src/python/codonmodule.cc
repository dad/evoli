/*
This file is part of the evoli project.
Copyright (C) 2006 Allan Drummond <dadrummond@gmail.com>,
                               Claus O. Wilke <cwilke@mail.utexas.edu>

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

#include "genetic-code.hh"
#include "codon.hh"


static PyObject *CodonErrorObject;

/* ----------------------------------------------------- */


bool validCodon( const char* codon )
{
	bool valid = true;
	const char* accept_chars = "aAgGtTuUcC"; // we only allow these characters in codons
	for ( int i=0; i<3; i++ )
	{
		bool found = false;
		for ( int j=0; j<10; j++ )
		{
			if ( codon[i]==accept_chars[j] )
			{
				found = true;
				break;
			}
		}
		if ( !found )
		{
			valid = false;
			break;
		}
	}
	return valid;
}

static char codon_validCodon__doc__[] =
"Tests whether a codon is valid (i.e., contains exactly three letters, composed of the letters a, A, c, C, u, U, t, T, g, G) or not."
;

static PyObject *
codon_validCodon(PyObject *self /* Not used */, PyObject *args)
{
	const char *codon;
	int size;
	if ( !PyArg_ParseTuple( args, "s#", &codon, &size ) ) {
		return NULL;
	}

	int valid = 1;
	if ( size != 3 ) 
		valid = 0;
	else
		if ( !validCodon( codon ) )
			valid = 0;
	
	return Py_BuildValue( "i", valid );
}

static char codon_calcDnDs__doc__[] =
"Calculates the number of nonsynonymous and synonymous substitutions between two codons, and returns them in a list."
;

static PyObject *
codon_calcDnDs(PyObject *self /* Not used */, PyObject *args)
{
	const char *codon1;
	const char *codon2;
	int size1, size2;
	if ( !PyArg_ParseTuple( args, "s#s#", &codon1, &size1, &codon2, &size2 ) ) {
		return NULL;
	}

	if ( size1 != 3 || size2 != 3 ) 
	{
		PyErr_SetString(CodonErrorObject, "Each codon must be exactly three characters long.");
		return NULL;
	}

	if ( !validCodon( codon1 ) || !validCodon( codon2 ) )
	{
		PyErr_SetString(CodonErrorObject, "Invalid codon. Codons can contain only the letters a, A, u, U, t, T, g, G, c, C.");
		return NULL;
	}
	
	int c1 = CodonUtil::lettersToCodon( codon1[0], codon1[1], codon1[2] );
	int c2 = CodonUtil::lettersToCodon( codon2[0], codon2[1], codon2[2] );

	double dn, ds;
	// calculate dn and ds and store in variables
	GeneticCodeUtil::calcDnDs( dn, ds, c1, c2 );

	return Py_BuildValue( "ff", dn, ds );
}

static char codon_calcNS_doc__[] =
"Calculates the number of nonsynonymous and synonymous sites in a codon. Counts stop sites as nonsynonymous sites."
;

static PyObject *
codon_calcNS(PyObject *self /* Not used */, PyObject *args)
{
	const char *codon1;
	int size1;
	if ( !PyArg_ParseTuple( args, "s#", &codon1, &size1 ) ) {
		return NULL;
	}

	if ( size1 != 3 ) 
	{
		PyErr_SetString(CodonErrorObject, "Codon must be exactly three characters long.");
		return NULL;
	}

	if ( !validCodon( codon1 ) )
	{
		PyErr_SetString(CodonErrorObject, "Invalid codon. Codons can contain only the letters a, A, u, U, t, T, g, G, c, C.");
		return NULL;
	}
	
	int c1 = CodonUtil::lettersToCodon( codon1[0], codon1[1], codon1[2] );

	// calculate N and S
	double S = GeneticCodeUtil::calcSynonymousSites( c1 );
	double N = 3-S;

	return Py_BuildValue( "ff", N, S );
}

/* List of methods defined in the module */

static struct PyMethodDef codon_methods[] = {
	{"validCodon",	(PyCFunction)codon_validCodon,	METH_VARARGS, codon_validCodon__doc__},
	{"calcDnDs",	(PyCFunction)codon_calcDnDs,	METH_VARARGS, codon_calcDnDs__doc__},
	{"calcNS",	(PyCFunction)codon_calcNS,	METH_VARARGS, codon_calcNS_doc__},
	{NULL,	 (PyCFunction)NULL, 0, NULL}		/* sentinel */
};


/* Initialization function for the module (*must* be called initfolder) */

static char codon_module_documentation[] = 
""
;

PyMODINIT_FUNC
initcodon(void)
{
	PyObject *m, *d;

	/* Create the module and add the functions */
	m = Py_InitModule4("codon", codon_methods,
		codon_module_documentation,
		(PyObject*)NULL,PYTHON_API_VERSION);

	/* Add some symbolic constants to the module */
	d = PyModule_GetDict(m);
	CodonErrorObject = PyString_FromString("codon.error");
	PyDict_SetItemString(d, "error", CodonErrorObject);

	/* XXXX Add constants here */
	
	/* Check for errors */
	if (PyErr_Occurred())
		Py_FatalError("can't initialize module codon");
}

