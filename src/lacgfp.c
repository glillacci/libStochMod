/*
 *  lacgfp.c
 *  StochMod
 *
 *	Lac-GFP construct model
 *
 *  Created by Gabriele Lillacci in January 2012.
 *	Latest revision: May 2012.
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
 === SPECIES ===
 	 X(0)	->	lacI		lacI mRNA
 	 X(1)	->	LACI		LACI protein
 	 X(2)	->	PLac		Unoccupied (active) Lac promoter
 	 X(3)  	-> 	O1Lac		Occupied Lac promoter, 1 repressor molecule bound
 	 X(4)  	-> 	O2Lac		Occupied Lac promoter, 2 repressor molecules bound
 	 X(5)  	-> 	O3Lac		Occupied Lac promoter, 3 repressor molecules bound
 	 X(6)  	-> 	O4Lac		Occupied Lac promoter, 4 repressor molecules bound
 	 X(7)  	-> 	gfp			gfp mRNA
 	 X(8)  	-> 	GFP			GFP protein

 === REACTIONS ===
 	 REACTION 1		NULL -(k1)-> lacI					Constitutive transcription of lacI mRNA
 	 REACTION 2		lacI -(k2)-> NULL					Degradation of lacI mRNA
 	 REACTION 3		lacI -(k3)-> lacI + LACI			Translation of LACI protein
 	 REACTION 4		LACI -(k4+k21*u)-> NULL				Degradation of LACI
 	 REACTION 5		LACI + PLac -(k5)-> O1Lac			Repressor binding
 	 REACTION 6		LACI + O1Lac -(k6)-> O2Lac
 	 REACTION 7		LACI + O2Lac -(k7)-> O3Lac
 	 REACTION 8		LACI + O3Lac -(k8)-> O4Lac
 	 REACTION 9		O1Lac -(k9)-> LACI + PLac			Repressor dissociation
 	 REACTION 10	O2Lac -(k10)-> LACI + O1Lac
 	 REACTION 11	O3Lac -(k11)-> LACI + O2Lac
 	 REACTION 12	O4Lac -(k12)-> LACI + O3Lac
 	 REACTION 13	PLac -(k13)-> PLac + gfp			Transcription of gfp mRNA
 	 REACTION 14 	O1Lac -(k14)-> O1Lac + gfp
 	 REACTION 15	O2Lac -(k15)-> O2Lac + gfp
 	 REACTION 16	O3Lac -(k16)-> O3Lac + gfp
 	 REACTION 17	O4Lac -(k17)-> O4Lac + gfp
 	 REACTION 18	gfp -(k18)-> NULL					Degradation of gfp mRNA
 	 REACTION 19	gfp -(k19)-> gfp + GFP				Translation of GFP
 	 REACTION 20	GFP -(k20)-> NULL					GFP degradation
  */


/**
 Propensity evaluation function for Lacgfp.
 */
int lacgfp_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop)
{
	// Check sizes of vectors
	if ((X->size != 9) || (params->size != 22) || (prop->size != 20))
	{
		fprintf (stderr, "error in lacgfp_propensity_eval: vector sizes are not correct\n");
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

	// Recover parameters from params vector
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
	double k19 = gsl_vector_get (params, 18);
	double k20 = gsl_vector_get (params, 19);
	double k21 = gsl_vector_get (params, 19);

	// Recover the input
	double u = gsl_vector_get (params, 21);

	// Evaluate the propensities
	gsl_vector_set (prop, 0, k1);
	gsl_vector_set (prop, 1, k2*X1);
	gsl_vector_set (prop, 2, k3*X1);
	gsl_vector_set (prop, 3, (k4+k21*u)*X2);
	gsl_vector_set (prop, 4, k5*X2*X3);
	gsl_vector_set (prop, 5, k6*X2*X4);
	gsl_vector_set (prop, 6, k7*X2*X5);
	gsl_vector_set (prop, 7, k8*X2*X6);
	gsl_vector_set (prop, 8, k9*X4);
	gsl_vector_set (prop, 9, k10*X5);
	gsl_vector_set (prop, 10, k11*X6);
	gsl_vector_set (prop, 11, k12*X7);
	gsl_vector_set (prop, 12, k13*X3);
	gsl_vector_set (prop, 13, k14*X4);
	gsl_vector_set (prop, 14, k15*X5);
	gsl_vector_set (prop, 15, k16*X6);
	gsl_vector_set (prop, 16, k17*X7);
	gsl_vector_set (prop, 17, k18*X8);
	gsl_vector_set (prop, 18, k19*X8);
	gsl_vector_set (prop, 19, k20*X9);

	// Signal that computation was completed successfully
	return GSL_SUCCESS;
}


/**
 State update function for Lacgfp.
 */
int lacgfp_state_update (gsl_vector * X, size_t rxnid)
{
	// Check sizes of state vector
	if (X->size != 9)
	{
		fprintf (stderr, "error in lacgfp_state_update: state vector size is not correct\n");
		return GSL_EFAILED;
	}

	// Check that reaction id is correct
	if (rxnid >= 20)
	{
		fprintf (stderr, "error in lacgfp_state_update: reaction id is not correct\n");
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
		gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (-1));
		gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (-1));
		gsl_vector_set (X, 3, gsl_vector_get (X, 3) + (1));
		break;

	case 5:
		gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (-1));
		gsl_vector_set (X, 3, gsl_vector_get (X, 3) + (-1));
		gsl_vector_set (X, 4, gsl_vector_get (X, 4) + (1));
		break;

	case 6:
		gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (-1));
		gsl_vector_set (X, 4, gsl_vector_get (X, 4) + (-1));
		gsl_vector_set (X, 5, gsl_vector_get (X, 5) + (1));
		break;

	case 7:
		gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (-1));
		gsl_vector_set (X, 5, gsl_vector_get (X, 5) + (-1));
		gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (1));
		break;

	case 8:
		gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (1));
		gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (1));
		gsl_vector_set (X, 3, gsl_vector_get (X, 3) + (-1));
		break;

	case 9:
		gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (1));
		gsl_vector_set (X, 3, gsl_vector_get (X, 3) + (1));
		gsl_vector_set (X, 4, gsl_vector_get (X, 4) + (-1));
		break;

	case 10:
		gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (1));
		gsl_vector_set (X, 4, gsl_vector_get (X, 4) + (1));
		gsl_vector_set (X, 5, gsl_vector_get (X, 5) + (-1));
		break;

	case 11:
		gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (1));
		gsl_vector_set (X, 5, gsl_vector_get (X, 5) + (1));
		gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (-1));
		break;

	case 12:
		gsl_vector_set (X, 7, gsl_vector_get (X, 7) + (1));
		break;

	case 13:
		gsl_vector_set (X, 7, gsl_vector_get (X, 7) + (1));
		break;

	case 14:
		gsl_vector_set (X, 7, gsl_vector_get (X, 7) + (1));
		break;

	case 15:
		gsl_vector_set (X, 7, gsl_vector_get (X, 7) + (1));
		break;

	case 16:
		gsl_vector_set (X, 7, gsl_vector_get (X, 7) + (1));
		break;

	case 17:
		gsl_vector_set (X, 7, gsl_vector_get (X, 7) + (-1));
		break;

	case 18:
		gsl_vector_set (X, 8, gsl_vector_get (X, 8) + (1));
		break;

	case 19:
		gsl_vector_set (X, 8, gsl_vector_get (X, 8) + (-1));
		break;
	}

	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Sample a new random initial state for Lacgfp.
 */
int lacgfp_initial_conditions (gsl_vector * X0, const gsl_rng * r)
{
	// Check sizes of state vector
	if (X0->size != 9)
	{
		fprintf (stderr, "error in lacgfp_initial_conditions: state vector size is not correct\n");
		return GSL_EFAILED;
	}

	// Sample new initial state
	gsl_vector_set (X0, 0, gsl_rng_uniform_int (r, 6));
	gsl_vector_set (X0, 1, gsl_rng_uniform_int (r, 11));
	gsl_vector_set (X0, 2, 1);
	gsl_vector_set (X0, 3, 0);
	gsl_vector_set (X0, 4, 0);
	gsl_vector_set (X0, 5, 0);
	gsl_vector_set (X0, 6, 0);
	gsl_vector_set (X0, 7, 0);
	gsl_vector_set (X0, 8, 0);

	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Output function for Lacgfp.
 */
int lacgfp_output (gsl_matrix * out)
{
	if ((out->size1 != 1) || (out->size2 != 9))
	{
		fprintf (stderr, "error in lacgfp_output: output matrix size is not correct\n");
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
 Model information function for Lacgfp.
 */
void lacgfp_mod_setup (stochmod * model)
{
	model->propensity = &lacgfp_propensity_eval;
	model->update = &lacgfp_state_update;
	model->initial = &lacgfp_initial_conditions;
	model->output = &lacgfp_output;
	model->nspecies = 9;
	model->nrxns = 20;
	model->nparams = 21;
	model->nin = 1;
	model->nout = 1;
	model->name = "Lac-GFP construct model (LACGFP)";
}

