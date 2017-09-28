/*
 *  stochrep.c
 *  StochMod
 *
 *	Stochastic repressilator model
 *
 *  Created by Gabriele Lillacci in March 2011.
 *	Latest revision: March 2011.
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
 
 
 === REACTIONS ===
 
   */


/**
 Propensity evaluation function for Stochrep.
 */
int stochrep_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop)
{
	// Check sizes of vectors
	if ((X->size != 21) || (params->size != 48) || (prop->size != 48))
	{
		printf("\n\n>> error in stochrep_propensity_eval: vector sizes are not correct...\n");
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
	double X10 = gsl_vector_get (X, 9);
	double X11 = gsl_vector_get (X, 10);
	double X12 = gsl_vector_get (X, 11);
	double X13 = gsl_vector_get (X, 12);
	double X14 = gsl_vector_get (X, 13);
	double X15 = gsl_vector_get (X, 14);
	double X16 = gsl_vector_get (X, 15);
	double X17 = gsl_vector_get (X, 16);
	double X18 = gsl_vector_get (X, 17);
	double X19 = gsl_vector_get (X, 18);
	double X20 = gsl_vector_get (X, 19);
	double X21 = gsl_vector_get (X, 20);
	
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
	double k21 = gsl_vector_get (params, 20);
	double k22 = gsl_vector_get (params, 21);
	double k23 = gsl_vector_get (params, 22);
	double k24 = gsl_vector_get (params, 23);
	double k25 = gsl_vector_get (params, 24);
	double k26 = gsl_vector_get (params, 25);
	double k27 = gsl_vector_get (params, 26);
	double k28 = gsl_vector_get (params, 27);
	double k29 = gsl_vector_get (params, 28);
	double k30 = gsl_vector_get (params, 29);
	double k31 = gsl_vector_get (params, 30);
	double k32 = gsl_vector_get (params, 31);
	double k33 = gsl_vector_get (params, 32);
	double k34 = gsl_vector_get (params, 33);
	double k35 = gsl_vector_get (params, 34);
	double k36 = gsl_vector_get (params, 35);
	double k37 = gsl_vector_get (params, 36);
	double k38 = gsl_vector_get (params, 37);
	double k39 = gsl_vector_get (params, 38);
	double k40 = gsl_vector_get (params, 39);
	double k41 = gsl_vector_get (params, 40);
	double k42 = gsl_vector_get (params, 41);
	double k43 = gsl_vector_get (params, 42);
	double k44 = gsl_vector_get (params, 43);
	double k45 = gsl_vector_get (params, 44);
	double k46 = gsl_vector_get (params, 45);
	double k47 = gsl_vector_get (params, 46);
	double k48 = gsl_vector_get (params, 47);
	
	// Evaluate the propensities
	gsl_vector_set (prop, 0, k1*X1*X21);
	gsl_vector_set (prop, 1, k2*X2*X21);
	gsl_vector_set (prop, 2, k3*X3*X21);
	gsl_vector_set (prop, 3, k4*X4*X21);
	gsl_vector_set (prop, 4, k5*X5);
	gsl_vector_set (prop, 5, k6*X4);
	gsl_vector_set (prop, 6, k7*X3);
	gsl_vector_set (prop, 7, k8*X2);
	gsl_vector_set (prop, 8, k9*X1);
	gsl_vector_set (prop, 9, k10*X2);
	gsl_vector_set (prop, 10, k11*X3);
	gsl_vector_set (prop, 11, k12*X4);
	gsl_vector_set (prop, 12, k13*X5);
	gsl_vector_set (prop, 13, k14*X6);
	gsl_vector_set (prop, 14, k15*X6);
	gsl_vector_set (prop, 15, k16*X7);
	gsl_vector_set (prop, 16, k17*X8*X7);
	gsl_vector_set (prop, 17, k18*X9*X7);
	gsl_vector_set (prop, 18, k19*X10*X7);
	gsl_vector_set (prop, 19, k20*X11*X7);
	gsl_vector_set (prop, 20, k21*X12);
	gsl_vector_set (prop, 21, k22*X11);
	gsl_vector_set (prop, 22, k23*X10);
	gsl_vector_set (prop, 23, k24*X9);
	gsl_vector_set (prop, 24, k25*X8);
	gsl_vector_set (prop, 25, k26*X9);
	gsl_vector_set (prop, 26, k27*X10);
	gsl_vector_set (prop, 27, k28*X11);
	gsl_vector_set (prop, 28, k29*X12);
	gsl_vector_set (prop, 29, k30*X13);
	gsl_vector_set (prop, 30, k31*X13);
	gsl_vector_set (prop, 31, k32*X14);
	gsl_vector_set (prop, 32, k33*X15*X14);
	gsl_vector_set (prop, 33, k34*X16*X14);
	gsl_vector_set (prop, 34, k35*X17*X14);
	gsl_vector_set (prop, 35, k36*X18*X14);
	gsl_vector_set (prop, 36, k37*X19);
	gsl_vector_set (prop, 37, k38*X18);
	gsl_vector_set (prop, 38, k39*X17);
	gsl_vector_set (prop, 39, k40*X16);
	gsl_vector_set (prop, 40, k41*X15);
	gsl_vector_set (prop, 41, k42*X16);
	gsl_vector_set (prop, 42, k43*X17);
	gsl_vector_set (prop, 43, k44*X18);
	gsl_vector_set (prop, 44, k45*X19);
	gsl_vector_set (prop, 45, k46*X20);
	gsl_vector_set (prop, 46, k47*X20);
	gsl_vector_set (prop, 47, k48*X21);
	
	// Signal that computation was completed successfully
	return GSL_SUCCESS;
}


/**
 State update function for Stochrep.
 */
int stochrep_state_update (gsl_vector * X, size_t rxnid)
{
	// Check sizes of state vector
	if (X->size != 21)
	{
		printf("\n\n>> error in stochrep_state_update: state vector size is not correct...\n");
		return GSL_EFAILED;
	}
	
	// Check that reaction id is correct
	if (rxnid>=48)
	{
		printf("\n\n>> error in stochrep_state_update: reaction id is not correct...\n");
		return GSL_EFAILED;
	}
	
	// Update the state vector according to which reaction fired
	switch (rxnid) {
		case 0:
			gsl_vector_set (X, 0, gsl_vector_get (X, 0) + (-1));
			gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (1));
			gsl_vector_set (X, 20, gsl_vector_get (X, 20) + (-1));
			break;
			
		case 1:
			gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (-1));
			gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (1));
			gsl_vector_set (X, 20, gsl_vector_get (X, 20) + (-1));
			break;
			
		case 2:
			gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (-1));
			gsl_vector_set (X, 3, gsl_vector_get (X, 3) + (1));
			gsl_vector_set (X, 20, gsl_vector_get (X, 20) + (-1));
			break;
			
		case 3:
			gsl_vector_set (X, 3, gsl_vector_get (X, 3) + (-1));
			gsl_vector_set (X, 4, gsl_vector_get (X, 4) + (1));
			gsl_vector_set (X, 20, gsl_vector_get (X, 20) + (-1));
			break;
			
		case 4:
			gsl_vector_set (X, 3, gsl_vector_get (X, 3) + (1));
			gsl_vector_set (X, 4, gsl_vector_get (X, 4) + (-1));
			gsl_vector_set (X, 20, gsl_vector_get (X, 20) + (1));
			break;
			
		case 5:
			gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (1));
			gsl_vector_set (X, 3, gsl_vector_get (X, 3) + (-1));
			gsl_vector_set (X, 20, gsl_vector_get (X, 20) + (1));
			break;
			
		case 6:
			gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (1));
			gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (-1));
			gsl_vector_set (X, 20, gsl_vector_get (X, 20) + (1));
			break;
			
		case 7:
			gsl_vector_set (X, 0, gsl_vector_get (X, 0) + (1));
			gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (-1));
			gsl_vector_set (X, 20, gsl_vector_get (X, 20) + (1));
			break;
			
		case 8:
			gsl_vector_set (X, 5, gsl_vector_get (X, 5) + (1));
			break;
			
		case 9:
			gsl_vector_set (X, 5, gsl_vector_get (X, 5) + (1));
			break;
			
		case 10:
			gsl_vector_set (X, 5, gsl_vector_get (X, 5) + (1));
			break;
			
		case 11:
			gsl_vector_set (X, 5, gsl_vector_get (X, 5) + (1));
			break;
			
		case 12:
			gsl_vector_set (X, 5, gsl_vector_get (X, 5) + (1));
			break;
			
		case 13:
			gsl_vector_set (X, 5, gsl_vector_get (X, 5) + (-1));
			break;
			
		case 14:
			gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (1));
			break;
			
		case 15:
			gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (-1));
			break;
			
		case 16:
			gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (-1));
			gsl_vector_set (X, 7, gsl_vector_get (X, 7) + (-1));
			gsl_vector_set (X, 8, gsl_vector_get (X, 8) + (1));
			break;
			
		case 17:
			gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (-1));
			gsl_vector_set (X, 8, gsl_vector_get (X, 8) + (-1));
			gsl_vector_set (X, 9, gsl_vector_get (X, 9) + (1));
			break;
			
		case 18:
			gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (-1));
			gsl_vector_set (X, 9, gsl_vector_get (X, 9) + (-1));
			gsl_vector_set (X, 10, gsl_vector_get (X, 10) + (1));
			break;
			
		case 19:
			gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (-1));
			gsl_vector_set (X, 10, gsl_vector_get (X, 10) + (-1));
			gsl_vector_set (X, 11, gsl_vector_get (X, 11) + (1));
			break;
			
		case 20:
			gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (1));
			gsl_vector_set (X, 10, gsl_vector_get (X, 10) + (1));
			gsl_vector_set (X, 11, gsl_vector_get (X, 11) + (-1));
			break;
			
		case 21:
			gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (1));
			gsl_vector_set (X, 9, gsl_vector_get (X, 9) + (1));
			gsl_vector_set (X, 10, gsl_vector_get (X, 10) + (-1));
			break;
			
		case 22:
			gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (1));
			gsl_vector_set (X, 8, gsl_vector_get (X, 8) + (1));
			gsl_vector_set (X, 9, gsl_vector_get (X, 9) + (-1));
			break;
			
		case 23:
			gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (1));
			gsl_vector_set (X, 7, gsl_vector_get (X, 7) + (1));
			gsl_vector_set (X, 8, gsl_vector_get (X, 8) + (-1));
			break;
			
		case 24:
			gsl_vector_set (X, 12, gsl_vector_get (X, 12) + (1));
			break;
			
		case 25:
			gsl_vector_set (X, 12, gsl_vector_get (X, 12) + (1));
			break;
			
		case 26:
			gsl_vector_set (X, 12, gsl_vector_get (X, 12) + (1));
			break;
			
		case 27:
			gsl_vector_set (X, 12, gsl_vector_get (X, 12) + (1));
			break;
			
		case 28:
			gsl_vector_set (X, 12, gsl_vector_get (X, 12) + (1));
			break;
			
		case 29:
			gsl_vector_set (X, 12, gsl_vector_get (X, 12) + (-1));
			break;
			
		case 30:
			gsl_vector_set (X, 13, gsl_vector_get (X, 13) + (1));
			break;
			
		case 31:
			gsl_vector_set (X, 13, gsl_vector_get (X, 13) + (-1));
			break;
			
		case 32:
			gsl_vector_set (X, 13, gsl_vector_get (X, 13) + (-1));
			gsl_vector_set (X, 14, gsl_vector_get (X, 14) + (-1));
			gsl_vector_set (X, 15, gsl_vector_get (X, 15) + (1));
			break;
			
		case 33:
			gsl_vector_set (X, 13, gsl_vector_get (X, 13) + (-1));
			gsl_vector_set (X, 15, gsl_vector_get (X, 15) + (-1));
			gsl_vector_set (X, 16, gsl_vector_get (X, 16) + (1));
			break;
			
		case 34:
			gsl_vector_set (X, 13, gsl_vector_get (X, 13) + (-1));
			gsl_vector_set (X, 16, gsl_vector_get (X, 16) + (-1));
			gsl_vector_set (X, 17, gsl_vector_get (X, 17) + (1));
			break;
			
		case 35:
			gsl_vector_set (X, 13, gsl_vector_get (X, 13) + (-1));
			gsl_vector_set (X, 17, gsl_vector_get (X, 17) + (-1));
			gsl_vector_set (X, 18, gsl_vector_get (X, 18) + (1));
			break;
			
		case 36:
			gsl_vector_set (X, 13, gsl_vector_get (X, 13) + (1));
			gsl_vector_set (X, 17, gsl_vector_get (X, 17) + (1));
			gsl_vector_set (X, 18, gsl_vector_get (X, 18) + (-1));
			break;
			
		case 37:
			gsl_vector_set (X, 13, gsl_vector_get (X, 13) + (1));
			gsl_vector_set (X, 16, gsl_vector_get (X, 16) + (1));
			gsl_vector_set (X, 17, gsl_vector_get (X, 17) + (-1));
			break;
			
		case 38:
			gsl_vector_set (X, 13, gsl_vector_get (X, 13) + (1));
			gsl_vector_set (X, 15, gsl_vector_get (X, 15) + (1));
			gsl_vector_set (X, 16, gsl_vector_get (X, 16) + (-1));
			break;
			
		case 39:
			gsl_vector_set (X, 13, gsl_vector_get (X, 13) + (1));
			gsl_vector_set (X, 14, gsl_vector_get (X, 14) + (1));
			gsl_vector_set (X, 15, gsl_vector_get (X, 15) + (-1));
			break;
			
		case 40:
			gsl_vector_set (X, 19, gsl_vector_get (X, 19) + (1));
			break;
			
		case 41:
			gsl_vector_set (X, 19, gsl_vector_get (X, 19) + (1));
			break;
			
		case 42:
			gsl_vector_set (X, 19, gsl_vector_get (X, 19) + (1));
			break;
			
		case 43:
			gsl_vector_set (X, 19, gsl_vector_get (X, 19) + (1));
			break;
			
		case 44:
			gsl_vector_set (X, 19, gsl_vector_get (X, 19) + (1));
			break;
			
		case 45:
			gsl_vector_set (X, 19, gsl_vector_get (X, 19) + (-1));
			break;
			
		case 46:
			gsl_vector_set (X, 20, gsl_vector_get (X, 20) + (1));
			break;
			
		case 47:
			gsl_vector_set (X, 20, gsl_vector_get (X, 20) + (-1));
			break;
			
	}
	
	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Model information function for Stochrep.
 */
void stochrep_mod_setup (stochmod * model)
{
	model->propensity = &stochrep_propensity_eval;
	model->update = &stochrep_state_update;
	model->initial = NULL;
	model->nspecies = 21;
	model->nrxns = 48;
	model->nparams = 48;
	model->nin = 0;
	model->nout = 3;
	model->name = "Stochastic Repressilator (STOCHREP)";
}

