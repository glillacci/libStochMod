/*
 *  syncirc.c
 *  StochMod
 *
 *	Three-gene synthetic repression cascade model
 *
 *  Created by Gabriele Lillacci in July 2010.
 *	Based on software created by Patrick Sheppard in June 2010.
 *	Latest revision: July 2010.
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
	0		a		TetR-RFP mRNA
	1		b		LI-GFP mRNA
	2		c		CFP mRNA
	3		A		TetR-RFP protein
	4		B		LI-GFP protein
	5		C		CFP protein
	6		Pb		LI-GFP promoter
	7		Pc		CFP promoter
	8		PbA		LI-GFP promoter bound to TetR-RFP 
	9		PcB		CFP promoter bound to LI-GFP
 
 === REACTIONS ===
	0		NULL --(kappa_a)--> a		Transcription at TetR-RFP promoter
	1		a --(gamma_a)--> NULL		Degradation of TetR-RFP mRNA
	2		a --(alpha_A)--> a + A		Translation of TetR-RFP mRNA
	3		A --(mu_A)--> NULL			Degradation of TetR-RFP protein
	4		A + Pb --(kd_A)--> PbA		Binding of TetR-RFP to LI-GFP promoter
	5		PbA --(kr_A)--> A + Pb		Dissociation of TetR-RFP from LI-GFP promoter
	6		Pb --(kappa_b)--> Pb + b	Transcription at LI-GFP promoter
	7		b --(gamma_b)--> NULL		Degradation of LI-GFP mRNA
	8		b --(alpha_B)--> b + B		Translation of LI-GFP mRNA
	9		B --(mu_B) --> NULL			Degradation of LI-GFP protein
	10		B + Pc --(kd_B)--> PcB		Binding of LI-GFP to CFP promoter
	11		PcB --(kr_B)-->Pc + B		Dissociation of LI-GFP from CFP promoter
	12		Pc --(kappa_c)--> Pc + c	Transcription at CFP promoter
	13		c --(gamma_c)--> NULL		Degradation of CFP mRNA
	14		c --(alpha_C)--> c + C		Translation of CFP mRNA
	15		C --(mu_C)--> NULL			Degradation of CFP protein
 */



/**
 Propensity evaluation function for Syncirc.
 */
int syncirc_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop)
{
	// Check sizes of vectors
	if ((X->size != 10) || (params->size != 16) || (prop->size != 16))
	{
		printf("\n\n>> error in syncirc_propensity_eval: vector sizes are not correct...\n");
		return GSL_EFAILED;
	}
	
	// Recover species from X vector
	double a = gsl_vector_get (X, 0);
	double b = gsl_vector_get (X, 1);
	double c = gsl_vector_get (X, 2);
	double A = gsl_vector_get (X, 3);
	double B = gsl_vector_get (X, 4);
	double C = gsl_vector_get (X, 5);
	double Pb = gsl_vector_get (X, 6);
	double Pc = gsl_vector_get (X, 7);
	double PbA = gsl_vector_get (X, 8);
	double PcB = gsl_vector_get (X, 9);
	
	// Recover parameters from params vector
	double kappa_a = gsl_vector_get (params, 0);
	double gamma_a = gsl_vector_get (params, 1);
	double alpha_A = gsl_vector_get (params, 2);
	double mu_A = gsl_vector_get (params, 3);
	double kd_A = gsl_vector_get (params, 4);
	double kr_A = gsl_vector_get (params, 5);
	double kappa_b = gsl_vector_get (params, 6);
	double gamma_b = gsl_vector_get (params, 7);
	double alpha_B = gsl_vector_get (params, 8);
	double mu_B = gsl_vector_get (params, 9);
	double kd_B = gsl_vector_get (params, 10);
	double kr_B = gsl_vector_get (params, 11);
	double kappa_c = gsl_vector_get (params, 12);
	double gamma_c = gsl_vector_get (params, 13);
	double alpha_C = gsl_vector_get (params, 14);
	double mu_C = gsl_vector_get (params, 15);
	
	// Evaluate the propensities
	gsl_vector_set (prop, 0, kappa_a*2);
	gsl_vector_set (prop, 1, gamma_a*a);
	gsl_vector_set (prop, 2, alpha_A*a);
	gsl_vector_set (prop, 3, mu_A*A);
	gsl_vector_set (prop, 4, kd_A*A*Pb);
	gsl_vector_set (prop, 5, kr_A*PbA);
	gsl_vector_set (prop, 6, kappa_b*Pb);
	gsl_vector_set (prop, 7, gamma_b*b);
	gsl_vector_set (prop, 8, alpha_B*b);
	gsl_vector_set (prop, 9, mu_B*B);
	gsl_vector_set (prop, 10, kd_B*B*Pc);
	gsl_vector_set (prop, 11, kr_B*PcB);
	gsl_vector_set (prop, 12, kappa_c*Pc);
	gsl_vector_set (prop, 13, gamma_c*c);
	gsl_vector_set (prop, 14, alpha_C*c);
	gsl_vector_set (prop, 15, mu_C*C);
	
	// Signal that computation was completed successfully
	return GSL_SUCCESS;
}


/**
 State update function for Syncirc.
 */
int syncirc_state_update (gsl_vector * X, size_t rxnid)
{
	// Check sizes of state vector
	if (X->size != 10)
	{
		printf("\n\n>> error in syncirc_state_update: state vector size is not correct...\n");
		return GSL_EFAILED;
	}
	
	// Check that reaction id is correct
	if (rxnid>=16)
	{
		printf("\n\n>> error in syncirc_state_update: reaction id is not correct...\n");
		return GSL_EFAILED;
	}
	
	// Update the state vector according to which reaction fired
	switch (rxnid) {
		case 0:
			// NULL --(kappa_a)--> a
			// Reaction 0 creates 1 molecule of a (species 0)
			gsl_vector_set (X, 0, gsl_vector_get (X, 0) + 1);
			break;
		
		case 1:
			// a --(gamma_a)--> NULL
			// Reaction 1 destroys 1 molecule of a (species 0)
			gsl_vector_set (X, 0, gsl_vector_get (X, 0) - 1);
			break;
		
		case 2:
			// a --(alpha_A)--> a + A
			// Reaction 2 creates 1 molecule of A (species 3)
			gsl_vector_set (X, 3, gsl_vector_get (X, 3) + 1);
			break;
		
		case 3:
			// A --(mu_A)--> NULL
			// Reaction 3 destroys 1 molecule of A (species 3)
			gsl_vector_set (X, 3, gsl_vector_get (X, 3) - 1);
			break;
		
		case 4:
			// A + Pb --(kd_A)--> PbA
			// Reaction 4 creates 1 molecule of PbA (species 8)
			// consumes 1 molecule of A (species 3) and
			// consumes 1 free promoter Pb (species 6)
			gsl_vector_set (X, 3, gsl_vector_get (X, 3) - 1);
			gsl_vector_set (X, 6, gsl_vector_get (X, 6) - 1);
			gsl_vector_set (X, 8, gsl_vector_get (X, 8) + 1);
			break;
			
		case 5:
			// PbA --(kr_A)--> A + Pb
			// Reaction 5 destroys 1 molecule of PbA (species 8)
			// creates 1 molecule of A (species 3) and
			// creates 1 free promoter Pb (species 6)
			gsl_vector_set (X, 8, gsl_vector_get (X, 8) - 1);
			gsl_vector_set (X, 3, gsl_vector_get (X, 3) + 1);
			gsl_vector_set (X, 6, gsl_vector_get (X, 6) + 1);
			break;
			
		case 6:
			// Pb --(kappa_b)--> Pb + b
			// Reaction 6 creates 1 molecule of b (species 1)
			gsl_vector_set (X, 1, gsl_vector_get (X, 1) + 1);
			break;
			
		case 7:
			// b --(gamma_b)--> NULL
			// Reaction 7 destroys 1 molecule of b (species 1)
			gsl_vector_set (X, 1, gsl_vector_get (X, 1) - 1);
			break;

		case 8:
			// b --(alpha_B)--> b + B
			// Reaction 8 creates 1 molecule of B (species 4)
			gsl_vector_set (X, 4, gsl_vector_get (X, 4) + 1);
			break;	
		
		case 9:
			// B --(mu_B) --> NULL
			// Reaction 9 destroys 1 molecule of B (species 4)
			gsl_vector_set (X, 4, gsl_vector_get (X, 4) - 1);
			break;
			
		case 10:
			// B + Pc --(kd_B)--> PcB	
			// Reaction 10 creates 1 molecule of PcB (species 9)
			// consumes 1 molecule of B (species 4) and
			// consumes 1 free promoter Pc (species 7)
			gsl_vector_set (X, 4, gsl_vector_get (X, 4) - 1);
			gsl_vector_set (X, 7, gsl_vector_get (X, 7) - 1);
			gsl_vector_set (X, 9, gsl_vector_get (X, 9) + 1);
			break;
			
		case 11:
			// PcB --(kr_B)-->Pc + B
			// Reaction 11 destroys 1 molecule of PcB (species 9)
			// creates 1 molecule of B (species 4) and
			// creates 1 free promoter Pc (species 7)
			gsl_vector_set (X, 9, gsl_vector_get (X, 9) - 1);
			gsl_vector_set (X, 4, gsl_vector_get (X, 4) + 1);
			gsl_vector_set (X, 7, gsl_vector_get (X, 7) + 1);
			break;
			
		case 12:
			// Pc --(kappa_c)--> Pc + c	
			// Reaction 12 creates 1 molecule of c (species 2)
			gsl_vector_set (X, 2, gsl_vector_get (X, 2) + 1);
			break;
		
		case 13:
			// c --(gamma_c)--> NULL	
			// Reaction 13 destroys 1 molecule of c (species 2)
			gsl_vector_set (X, 2, gsl_vector_get (X, 2) - 1);
			break;
		
		case 14:
			// c --(alpha_C)--> c + C	
			// Reaction 14 creates 1 molecule of C (species 5)
			gsl_vector_set (X, 5, gsl_vector_get (X, 5) + 1);
			break;
		
		case 15:
			// C --(mu_C)--> NULL	
			// Reaction 15 destroys 1 molecule of C (species 5)
			gsl_vector_set (X, 5, gsl_vector_get (X, 5) - 1);
			break;
			
	}
	
	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Model information function for Syncirc.
 */
void syncirc_mod_setup (stochmod * model)
{
	model->propensity = &syncirc_propensity_eval;
	model->update = &syncirc_state_update;
	model->initial = NULL;
	model->nspecies = 10;
	model->nrxns = 16;
	model->nparams = 16;
	model->nin = 0;
	model->nout = 3;
	model->name = "Three-gene synthetic repression cascade (SYNCIRC)";
}

