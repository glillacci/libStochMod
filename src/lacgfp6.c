/*
 *  lacgfp6.c
 *  StochMod
 *
 *	Lac-GFP construct model v6
 *	Model with GFP maturation
 *
 *  Created by Gabriele Lillacci in August 2012.
 *	Latest revision: August 2012.
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


/**
 *
  === SPECIES ===
	X(0)	--->	lacI		lacI mRNA
	X(1)	--->	LACI		LACI protein monomer
	X(2)	--->	LACI2		LACI dimer
	X(3)	--->	PLac		Unoccopied (active) Lac promoter
	X(4)	--->	O2Lac		Occupied promoter bound to LACI dimer
	X(5)	--->	O4Lac		Occupied promoter bound to LACI tetramer
	X(6)	--->	gfp			gfp mRNA
	X(7)	--->	GFP			GFP protein (dark)
	X(8)	--->	mGFP		GFP protein (mature)


 === REACTIONS ===
	reaction 1.		NULL --(k1)--> lacI					Transcription of lacI mRNA (constitutive)
	reaction 2.		lacI --(k2)--> NULL					Degradation of lacI mRNA (constitutive)
	reaction 3.		lacI --(k3)--> lacI + LACI			Translation of LACI protein
	reaction 4.		LACI --(k4+k5*u1)--> NULL			Degradation of LACI protein, increased by the input (IPTG)
	reaction 5.		LACI + LACI --(k6)--> LACI2			Dimerization of LACI protein
	reaction 6.		LACI2 --(k7)--> LACI + LACI			Dissociation of LACI dimer
	reaction 7.		LACI2 + PLac --(k8)--> O2Lac		Binding of LACI dimer to Lac operator sequence
	reaction 8.		O2Lac --(k9)--> LACI2 + PLac		Dissociation of LACI dimer from operator sequence
	reaction 9.		O2Lac + LACI2 --(k10)--> O4Lac		Binding of two LacI/operator complexes and tetramerization
	reaction 10.	O4Lac --(k11)--> LACI2 + O2Lac		Dissociation of tetramer structure
	reaction 11.	PLac --(k12)--> PLac + gfp			Transcription of gfp mRNA from active Lac promoter
	reaction 12.	O2Lac --(k13)--> O2Lac + gfp		Transcription of gfp mRNA from Lac promoter bound to LacI dimer
	reaction 13.	O4Lac --(k14)--> O4Lac + gfp		Transcription of gfp mRNA from Lac promoter bound to LacI tetramer
	reaction 14.	gfp --(k15)--> NULL					Degradation of gfp mRNA
	reaction 15.	gfp --(k16)--> gfp + GFP			Translation of dark GFP protein
	reaction 16.	GFP --(k17)--> NULL					Degradation of dark GFP protein
	reaction 17.	GFP --(k18)--> mGFP					Maturation of GFP
	reaction 18.	mGFP --(k17)--> NULL				Degradation of mature GFP protein

  */


/**
 Propensity evaluation function for Lacgfp6.
 */
int lacgfp6_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop)
{
	// Check sizes of vectors
	if ((X->size != 9) || (params->size != 19) || (prop->size != 18))
	{
		fprintf (stderr, "error in lacgfp6_propensity_eval: vector sizes are not correct\n");
		fprintf (stderr, "\tstate: %d - params: %d - propensities: %d\n", (int) X->size, (int) params->size, (int) prop->size);
		return GSL_EFAILED;
	}

	// Recover species from X vector
	double X1 = gsl_vector_get (X, 0);
	double X2 = gsl_vector_get (X, 1);
	double X3 = gsl_vector_get (X, 2);
	double X4 = gsl_vector_get (X, 3);
	double X5 = gsl_vector_get (X, 4);
	double X6 = gsl_vector_get (X, 5);
	double X7 = gsl_vector_get (X, 6);
	double X8 = gsl_vector_get (X, 7);
	double X9 = gsl_vector_get (X, 8);

	// Parameter recovery statements
	double k1 = gsl_vector_get (params, 0);
	double k2 = gsl_vector_get (params, 1);
	double k3 = gsl_vector_get (params, 2);
	double k4 = gsl_vector_get (params, 3);
	double k5 = gsl_vector_get (params, 4);
	double k6 = gsl_vector_get (params, 5);
	double k7 = gsl_vector_get (params, 6);
	double k8 = gsl_vector_get (params, 7);
	double k9 = gsl_vector_get (params, 8);
	double k10 = gsl_vector_get (params, 9);
	double k11 = gsl_vector_get (params, 10);
	double k12 = gsl_vector_get (params, 11);
	double k13 = gsl_vector_get (params, 12);
	double k14 = gsl_vector_get (params, 13);
	double k15 = gsl_vector_get (params, 14);
	double k16 = gsl_vector_get (params, 15);
	double k17 = gsl_vector_get (params, 16);
	double k18 = gsl_vector_get (params, 17);

	// Input recovery statements
	double u1 = gsl_vector_get (params, 18);

	// Propensity evaluation statements
	gsl_vector_set (prop, 0, (k1));
	gsl_vector_set (prop, 1, (k2)*X1);
	gsl_vector_set (prop, 2, (k3)*X1);
	gsl_vector_set (prop, 3, (k4+k5*u1)*X2);
	gsl_vector_set (prop, 4, (k6)*X2*(X2-1));
	gsl_vector_set (prop, 5, (k7)*X3);
	gsl_vector_set (prop, 6, (k8)*X3*X4);
	gsl_vector_set (prop, 7, (k9)*X5);
	gsl_vector_set (prop, 8, (k10)*X5*X3);
	gsl_vector_set (prop, 9, (k11)*X6);
	gsl_vector_set (prop, 10, (k12)*X4);
	gsl_vector_set (prop, 11, (k13)*X5);
	gsl_vector_set (prop, 12, (k14)*X6);
	gsl_vector_set (prop, 13, (k15)*X7);
	gsl_vector_set (prop, 14, (k16)*X7);
	gsl_vector_set (prop, 15, (k17)*X8);
	gsl_vector_set (prop, 16, (k18)*X8);
	gsl_vector_set (prop, 17, (k17)*X9);

	// Signal that computation was completed successfully
	return GSL_SUCCESS;
}


/**
 State update function for Lacgfp6.
 */
int lacgfp6_state_update (gsl_vector * X, size_t rxnid)
{
	// Check sizes of state vector
	if (X->size != 9)
	{
		fprintf (stderr, "error in lacgfp6_state_update: state vector size is not correct\n");
		return GSL_EFAILED;
	}

	// Check that reaction id is correct
	if (rxnid >= 18)
	{
		fprintf (stderr, "error in lacgfp6_state_update: reaction id is not correct\n");
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
		gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (-2));
		gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (1));
		break;

	case 5:
		gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (2));
		gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (-1));
		break;

	case 6:
		gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (-1));
		gsl_vector_set (X, 3, gsl_vector_get (X, 3) + (-1));
		gsl_vector_set (X, 4, gsl_vector_get (X, 4) + (1));
		break;

	case 7:
		gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (1));
		gsl_vector_set (X, 3, gsl_vector_get (X, 3) + (1));
		gsl_vector_set (X, 4, gsl_vector_get (X, 4) + (-1));
		break;

	case 8:
		gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (-1));
		gsl_vector_set (X, 4, gsl_vector_get (X, 4) + (-1));
		gsl_vector_set (X, 5, gsl_vector_get (X, 5) + (1));
		break;

	case 9:
		gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (1));
		gsl_vector_set (X, 4, gsl_vector_get (X, 4) + (1));
		gsl_vector_set (X, 5, gsl_vector_get (X, 5) + (-1));
		break;

	case 10:
		gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (1));
		break;

	case 11:
		gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (1));
		break;

	case 12:
		gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (1));
		break;

	case 13:
		gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (-1));
		break;

	case 14:
		gsl_vector_set (X, 7, gsl_vector_get (X, 7) + (1));
		break;

	case 15:
		gsl_vector_set (X, 7, gsl_vector_get (X, 7) + (-1));
		break;

	case 16:
		gsl_vector_set (X, 7, gsl_vector_get (X, 7) + (-1));
		gsl_vector_set (X, 8, gsl_vector_get (X, 8) + (1));
		break;

	case 17:
		gsl_vector_set (X, 8, gsl_vector_get (X, 8) + (-1));
		break;
	}

	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Sample a new random initial state for Lacgfp6.
 */
int lacgfp6_initial_conditions (gsl_vector * X0, const gsl_rng * r)
{
	// Check sizes of state vector
	if (X0->size != 9)
	{
		fprintf (stderr, "error in lacgfp6_initial_conditions: state vector size is not correct\n");
		return GSL_EFAILED;
	}

	// Sample a new initial state
	gsl_vector_set (X0, 0, gsl_rng_uniform_int (r, 6));
	gsl_vector_set (X0, 1, gsl_rng_uniform_int (r, 11));
	gsl_vector_set (X0, 2, 0);
	gsl_vector_set (X0, 3, 0);
	gsl_vector_set (X0, 4, 0);
	gsl_vector_set (X0, 5, 50 + gsl_rng_uniform_int (r, 21));
	gsl_vector_set (X0, 6, 0);
	gsl_vector_set (X0, 7, 0);
	gsl_vector_set (X0, 8, 0);

	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Output function for Lacgfp6.
 */
int lacgfp6_output (gsl_matrix * out)
{
	if ((out->size1 != 1) || (out->size2 != 9))
	{
		fprintf (stderr, "error in lacgfp6_output: output matrix size is not correct\n");
		return GSL_EFAILED;
	}

	// Reset the output matrix
	gsl_matrix_set_all (out, 0.0);

	// Set the non-zero terms
	gsl_matrix_set (out, 0, 8, 1.0);

	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Model information function for lacgfp6.
 */
void lacgfp6_mod_setup (stochmod * model)
{
	model->propensity = &lacgfp6_propensity_eval;
	model->update = &lacgfp6_state_update;
	model->initial = &lacgfp6_initial_conditions;
	model->output = &lacgfp6_output;
	model->nspecies = 9;
	model->nrxns = 18;
	model->nparams = 18;
	model->nin = 1;
	model->nout = 1;
	model->name = "Lac-GFP construct model v6 (LACGFP6)";
}

