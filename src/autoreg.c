/*
 *  autoreg.c
 *  StochMod
 *
 *	Stochastic gene autoregulation model
 *
 *  Created by Gabriele Lillacci in June 2011.
 *	Latest revision: June 2011.
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
 % X(0)  -> A - Active (unoccupied) promoter
 % X(1)  -> O - Occupied promoter
 % X(2)  -> m - mRNA
 % X(3)  -> p - protein
 % X(4)  -> pp - phospho protein
 
 === REACTIONS ===
  */



/**
 Propensity evaluation function for Autoreg.
 */
int autoreg_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop)
{
	// Check sizes of vectors
	if ((X->size != 5) || (params->size != 9) || (prop->size != 9))
	{
		printf("\n\n>> error in autoreg_propensity_eval: vector sizes are not correct...\n");
		return GSL_EFAILED;
	}
	
	// Recover species from X vector
	double X1 = gsl_vector_get (X, 0);
	double X2 = gsl_vector_get (X, 1);
	double X3 = gsl_vector_get (X, 2);
	double X4 = gsl_vector_get (X, 3);
	double X5 = gsl_vector_get (X, 4);
	
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
	
	// Evaluate the propensities
	gsl_vector_set (prop, 0, k1*X1*X5);
	gsl_vector_set (prop, 1, k2*X2);
	gsl_vector_set (prop, 2, k3*X1);
	gsl_vector_set (prop, 3, k4*X2);
	gsl_vector_set (prop, 4, k5*X3);
	gsl_vector_set (prop, 5, k6*X3);
	gsl_vector_set (prop, 6, k7*X4);
	gsl_vector_set (prop, 7, k8*X4/(1+X4));
	gsl_vector_set (prop, 8, k9*X5);
	
	// Signal that computation was completed successfully
	return GSL_SUCCESS;
}


/**
 State update function for Syncirc.
 */
int autoreg_state_update (gsl_vector * X, size_t rxnid)
{
	// Check sizes of state vector
	if (X->size != 5)
	{
		printf("\n\n>> error in autoreg_state_update: state vector size is not correct...\n");
		return GSL_EFAILED;
	}
	
	// Check that reaction id is correct
	if (rxnid>=9)
	{
		printf("\n\n>> error in autoreg_state_update: reaction id is not correct...\n");
		return GSL_EFAILED;
	}
	
	// Update the state vector according to which reaction fired
	switch (rxnid) {
		case 0:
			gsl_vector_set (X, 0, gsl_vector_get (X, 0) + (-1));
			gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (1));
			gsl_vector_set (X, 4, gsl_vector_get (X, 4) + (-1));
			break;
			
		case 1:
			gsl_vector_set (X, 0, gsl_vector_get (X, 0) + (1));
			gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (-1));
			gsl_vector_set (X, 4, gsl_vector_get (X, 4) + (1));
			break;
			
		case 2:
			gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (1));
			break;
			
		case 3:
			gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (1));
			break;
			
		case 4:
			gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (-1));
			break;
			
		case 5:
			gsl_vector_set (X, 3, gsl_vector_get (X, 3) + (1));
			break;
			
		case 6:
			gsl_vector_set (X, 3, gsl_vector_get (X, 3) + (-1));
			break;
			
		case 7:
			gsl_vector_set (X, 3, gsl_vector_get (X, 3) + (-1));
			gsl_vector_set (X, 4, gsl_vector_get (X, 4) + (1));
			break;
			
		case 8:
			gsl_vector_set (X, 3, gsl_vector_get (X, 3) + (1));
			gsl_vector_set (X, 4, gsl_vector_get (X, 4) + (-1));
			break;
			
	}
	
	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Sample a new random initial state for Autoreg.
 */
int autoreg_initial_conditions (gsl_vector * X0, const gsl_rng * r)
{
	// Check sizes of state vector
	if (X0->size != 5)
	{
		printf("\n\n>> error in autoreg_state_update: initial state vector size is not correct...\n");
		return GSL_EFAILED;
	}
	
	// Sample new initial state
	gsl_vector_set (X0, 0, gsl_rng_uniform_int (r, 3));
	gsl_vector_set (X0, 1, 2 - gsl_vector_get (X0, 0));
	gsl_vector_set (X0, 2, gsl_rng_uniform_int (r, 21));
	gsl_vector_set (X0, 3, gsl_rng_uniform_int (r, 201));
	gsl_vector_set (X0, 4, gsl_rng_uniform_int (r, 21));
	
	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Model information function for Autoreg.
 */
void autoreg_mod_setup (stochmod * model)
{
	model->propensity = &autoreg_propensity_eval;
	model->update = &autoreg_state_update;
	model->initial = &autoreg_initial_conditions;
	model->nspecies = 5;
	model->nrxns = 9;
	model->nparams = 9;
	model->nin = 0;
	model->nout = 1;
	model->name = "Stochastic Gene Autoregulation Model (AUTOREG)";
}

