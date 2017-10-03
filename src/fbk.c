/*
 *  fbk.c
 *  StochMod
 *
 *  This file is part of libStochMod.
 *  Copyright 2011-2017 Gabriele Lillacci.
 *
 *  libStochMod is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  libStochMod is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libStochMod.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../stochmod.h"

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
	reaction 6.		B + A --(k6)--> B		Repression of positive regulator
	reaction 7.		M --(k7)--> NULL		Degradation of regulated species
	reaction 8.		M --(k8)--> M + P		Production of reporter gene
	reaction 9.		P --(k9)--> NULL		Degradation of reporter gene

	*/


/**
 Propensity evaluation function for FBK.
 */
int fbk_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop)
{
	// Check sizes of vectors
	if ((X->size != N) || (params->size != L+Z) || (prop->size != R))
	{
		fprintf (stderr, "error in fbk_propensity_eval: vector sizes are not correct\n");
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
	gsl_vector_set (prop, 5, (k6)*X2*X1);
	gsl_vector_set (prop, 6, (k7)*X3);
	gsl_vector_set (prop, 7, (k8)*X3);
	gsl_vector_set (prop, 8, (k9)*X4);

	// Signal that computation was completed successfully
	return GSL_SUCCESS;
}


/**
 State update function for FBK.
 */
int fbk_state_update (gsl_vector * X, size_t rxnid)
{
	// Check sizes of state vector
	if (X->size != N)
	{
		fprintf (stderr, "error in fbk_state_update: state vector size is not correct\n");
		return GSL_EFAILED;
	}

	// Check that reaction id is correct
	if (rxnid >= R)
	{
		fprintf (stderr, "error in fbk_state_update: reaction id is not correct\n");
		return GSL_EFAILED;
	}

	// Update the state vector according to which reaction fired
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
		gsl_vector_set (X, 0, gsl_vector_get (X, 0) + (-1));
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
 Sample a new random initial state for FBK.
 */
int fbk_initial_conditions (gsl_vector * X0, const gsl_rng * r)
{
	// Check sizes of state vector
	if (X0->size != N)
	{
		fprintf (stderr, "error in fbk_initial_conditions: state vector size is not correct\n");
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
 Output function for FBK.
 */
int fbk_output (gsl_matrix * out)
{
	if ((out->size1 != P) || (out->size2 != N))
	{
		fprintf (stderr, "error in fbk_output: output matrix size is not correct\n");
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
 Model information function for FBK.
 */
void fbk_mod_setup (stochmod * model)
{
	model->propensity = &fbk_propensity_eval;
	model->update = &fbk_state_update;
	model->initial = &fbk_initial_conditions;
	model->output = &fbk_output;
	model->nspecies = N;
	model->nrxns = R;
	model->nparams = L;
	model->nin = Z;
	model->nout = P;
	model->name = "Feedback loop (FBK)";
}
