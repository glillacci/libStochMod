/*
 *  iFF.c
 *  StochMod
 *
 *	Incoherent feed-forward loop
 *
 *  Created by Gabriele Lillacci in September 2012.
 *	Latest revision: September 2012.
 *
 *
 *	This free software is available under the Creative Commons Attribution Share Alike License.
 *	You are permitted to use, redistribute and adapt this software as long as appropriate credit
 *	is given to the original author, and all derivative works are distributed under the same
 *	license or a compatible one.
 *	For more information, visit http://creativecommons.org/licenses/by-sa/3.0/ or send a letter to
 *	Creative Commons, 171 2nd Street, Suite 300, San Francisco, California, 94105, USA.
 */

#include <stochmod.h>

// Number of species
#define N 4
// Number of reactions
#define R 9
// Number of parameters
#define L 6
// Number of inputs
#define Z 0
// Number of outputs
#define P 1


/**
 *
  === SPECIES ===
	X1	--->	A		Positive regulator
	X2	--->	B		Negative regulator
	X3	--->	M		Regulated species
	X4	--->	P		Reporter gene

 === REACTIONS ===
	reaction 1.		NULL --(k1)--> A		Constitutive production of positive regulator
	reaction 2.		A --(k2)--> NULL		Degradation of positive regulator
	reaction 3.		A --(k3)--> A + B		Production of negative regulator
	reaction 4.		B --(k4)--> NULL		Degradation of negative regulator
	reaction 5.		A --(k5)--> A + M		Activation of regulated species
	reaction 6.		B + M --(k6)--> B		Repression of regulated species
	reaction 7.		M --(k7)--> NULL		Degradation of regulated species
	reaction 8.		M --(k8)--> M + P		Production of reporter gene
	reaction 9.		P --(k9)--> NULL		Degradation of reporter gene
  */


/**
 Propensity evaluation function for iFF.
 */
int iff_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop)
{
	// Check sizes of vectors
	if ((X->size != N) || (params->size != L+Z) || (prop->size != R))
	{
		fprintf (stderr, "error in iff_propensity_eval: vector sizes are not correct\n");
		fprintf (stderr, "\tstate: %d - params: %d - propensities: %d\n", (int) X->size, (int) params->size, (int) prop->size);
		return GSL_EFAILED;
	}

	// Recover species from X vector
	double X1 = gsl_vector_get (X, 0);
	double X2 = gsl_vector_get (X, 1);
	double X3 = gsl_vector_get (X, 2);
	double X4 = gsl_vector_get (X, 3);

	// Parameter recovery statements
	double k1 = gsl_vector_get (params, 0);
	double k2 = gsl_vector_get (params, 1);
	double k3 = gsl_vector_get (params, 2);
	double k4 = gsl_vector_get (params, 3);
	double k5 = gsl_vector_get (params, 4);
	double k6 = gsl_vector_get (params, 5);
	double k7 = 1.0;
	double k8 = 1.0;
	double k9 = 1.0;

	// Propensity evaluation statements
	gsl_vector_set (prop, 0, (k1));
	gsl_vector_set (prop, 1, (k2)*X1);
	gsl_vector_set (prop, 2, (k3)*X1);
	gsl_vector_set (prop, 3, (k4)*X2);
	gsl_vector_set (prop, 4, (k5)*X1);
	gsl_vector_set (prop, 5, (k6)*X2*X3);
	gsl_vector_set (prop, 6, (k7)*X3);
	gsl_vector_set (prop, 7, (k8)*X3);
	gsl_vector_set (prop, 8, (k9)*X4);

	// Signal that computation was completed successfully
	return GSL_SUCCESS;
}


/**
 State update function for iFF.
 */
int iff_state_update (gsl_vector * X, size_t rxnid)
{
	// Check sizes of state vector
	if (X->size != N)
	{
		fprintf (stderr, "error in iff_state_update: state vector size is not correct\n");
		return GSL_EFAILED;
	}

	// Check that reaction id is correct
	if (rxnid >= R)
	{
		fprintf (stderr, "error in iff_state_update: reaction id is not correct\n");
		return GSL_EFAILED;
	}

	switch (rxnid)
	{
	case 0:
		gsl_vector_set (X, 0, gsl_vector_get (X, 0) + (1));
		break;

	case 1:
		gsl_vector_set (X, 0, gsl_vector_get (X, 0) + (-1));
		break;

	case 2:
		gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (1));
		break;

	case 3:
		gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (-1));
		break;

	case 4:
		gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (1));
		break;

	case 5:
		gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (-1));
		break;

	case 6:
		gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (-1));
		break;

	case 7:
		gsl_vector_set (X, 3, gsl_vector_get (X, 3) + (1));
		break;

	case 8:
		gsl_vector_set (X, 3, gsl_vector_get (X, 3) + (-1));
		break;
	}


	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Sample a new random initial state for iFF.
 */
int iff_initial_conditions (gsl_vector * X0, const gsl_rng * r)
{
	// Check sizes of state vector
	if (X0->size != N)
	{
		fprintf (stderr, "error in iff_initial_conditions: state vector size is not correct\n");
		return GSL_EFAILED;
	}

	// Sample a new initial state
	gsl_vector_set (X0, 0, 0);
	gsl_vector_set (X0, 1, 0);
	gsl_vector_set (X0, 2, 0);
	gsl_vector_set (X0, 3, 0);

	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Output function for iFF.
 */
int iff_output (gsl_matrix * out)
{
	if ((out->size1 != P) || (out->size2 != N))
	{
		fprintf (stderr, "error in iff_output: output matrix size is not correct\n");
		return GSL_EFAILED;
	}

	// Reset the output matrix
	gsl_matrix_set_all (out, 0.0);

	// Set the non-zero terms
	gsl_matrix_set (out, 0, 2, 1.0);

	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Model information function for iFF.
 */
void iff_mod_setup (stochmod * model)
{
	model->propensity = &iff_propensity_eval;
	model->update = &iff_state_update;
	model->initial = &iff_initial_conditions;
	model->output = &iff_output;
	model->nspecies = N;
	model->nrxns = R;
	model->nparams = L;
	model->nin = Z;
	model->nout = P;
	model->name = "Incoherent feed-forward loop (iFF)";
}

