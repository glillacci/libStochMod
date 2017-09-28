/*
 *  birthdeath.c
 *  StochMod
 *
 *	Birth-death reaction of a single species
 *
 *  Created by Gabriele Lillacci in July 2012.
 *	Latest revision: July 2012.
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
#define N 1
// Number of reactions
#define R 2
// Number of parameters
#define L 2
// Number of inputs
#define Z 0
// Number of outputs
#define P 1


/**
 === SPECIES ===
 	 X(0)  -> A - The one and only species
 
 === REACTIONS ===
 	 NULL --(k1)--> A		The birth reaction
 	 A --(k2)--> NULL		The death reaction
  */



/**
 Propensity evaluation function for BirthDeath.
 */
int birthdeath_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop)
{
	// Check sizes of vectors
	if ((X->size != N) || (params->size != L+Z) || (prop->size != R))
	{
		printf("\n\n>> error in birthdeath_propensity_eval: vector sizes are not correct...\n");
		return GSL_EFAILED;
	}
	
	// Recover species from X vector
	double X1 = gsl_vector_get (X, 0);
	
	// Recover parameters from params vector
	double k1 = gsl_vector_get (params, 0);
	double k2 = gsl_vector_get (params, 1);
	
	// Evaluate the propensities
	gsl_vector_set (prop, 0, k1);
	gsl_vector_set (prop, 1, k2*X1);
	
	// Signal that computation was completed successfully
	return GSL_SUCCESS;
}


/**
 State update function for Syncirc.
 */
int birthdeath_state_update (gsl_vector * X, size_t rxnid)
{
	// Check sizes of state vector
	if (X->size != N)
	{
		printf("\n\n>> error in birthdeath_state_update: state vector size is not correct...\n");
		return GSL_EFAILED;
	}
	
	// Check that reaction id is correct
	if (rxnid >= R)
	{
		printf("\n\n>> error in birthdeath_state_update: reaction id is not correct...\n");
		return GSL_EFAILED;
	}
	
	// Update the state vector according to which reaction fired
	switch (rxnid) {
		case 0:
			gsl_vector_set (X, 0, gsl_vector_get (X, 0) + 1);
			break;
			
		case 1:
			gsl_vector_set (X, 0, gsl_vector_get (X, 0) - 1);
			break;
	}
	
	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Sample a new random initial state for BirthDeath.
 */
int birthdeath_initial_conditions (gsl_vector * X0, const gsl_rng * r)
{
	// Check sizes of state vector
	if (X0->size != N)
	{
		printf("\n\n>> error in birthdeath_state_update: initial state vector size is not correct...\n");
		return GSL_EFAILED;
	}
	
	// Sample new initial state
	gsl_vector_set (X0, 0, gsl_rng_uniform_int (r, 11));
	
	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Output function for BirthDeath.
 */
int birthdeath_output (gsl_matrix * out)
{
	if ((out->size1 != N) || (out->size2 != P))
	{
		fprintf (stderr, "error in birthdeath_output: output matrix size is not correct\n");
		return GSL_EFAILED;
	}

	// Reset the output matrix
	gsl_matrix_set_all (out, 0.0);

	// Set the non-zero terms
	gsl_matrix_set (out, 0, 0, 1.0);

	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Model information function for BirthDeath.
 */
void birthdeath_mod_setup (stochmod * model)
{
	model->propensity = &birthdeath_propensity_eval;
	model->update = &birthdeath_state_update;
	model->initial = &birthdeath_initial_conditions;
	model->output = &birthdeath_output;
	model->nspecies = N;
	model->nrxns = R;
	model->nparams = L;
	model->nin = Z;
	model->nout = P;
	model->name = "Birth-Death process of a single chemical species (BIRTHDEATH)";
}

