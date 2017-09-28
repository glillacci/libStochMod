/*
 *  synpi1.c
 *  StochMod
 *
 *	Synthetic PI v1
 *	Proportional feedback only
 *
 *  Created by Gabriele Lillacci in November 2013.
 *	Latest revision: November 2013.
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


//Number of species
#define N 8
//Number of reactions
#define R 14
//Number of parameters
#define L 13
//Number of inputs
#define Z 1
//Number of outputs
#define P 1


/**
=== SPECIES ===
	X1	--->	P4						Free P4 promoter
	X2	--->	P4-A					P4 promoter bound to AraC::GFP
	X3	--->	AraC::GFP				AraC::GFP fusion
	X4	--->	Ptac					Unoccopied (active) tac promoter
	X5	--->	Ptac-O2					Occupied tac promoter bound to Cherry::LacI::LVA dimer
	X6	--->	Ptac-O4					Occupied tac promoter bound to Cherry::LacI::LVA tetramer
	X7	--->	Cherry::LacI::LVA		Cherry::LacI::LVA monomer
	X8	--->	Cherry::LacI::LVA-2		Cherry::LacI::LVA dimer


 === REACTIONS ===
	reaction 1.		P4 + AraC::GFP --(k1)--> P4-A											Binding of AraC::GFP to P4 promoter
	reaction 2.		P4-A --(k2)--> P4 + AraC::GFP											Dissociation of AraC::GFP from P4 promoter
	reaction 3.		P4  --(k3)--> P4  + Cherry::LacI::LVA									Production of Cherry::LacI fusion from uninduced P4 promoter
	reaction 4.		P4-A --(k4)--> P4-A + Cherry::LacI::LVA									Production of Cherry::LacI fusion from induced P4 promoter
	reaction 5.		Cherry::LacI::LVA --(k5+k6*u1)--> NULL									Degradation of Cherry::LacI fusion, increased by the input (IPTG)
	reaction 6.		Cherry::LacI::LVA + Cherry::LacI::LVA --(k7)--> Cherry::LacI::LVA-2		Dimerization of Cherry::LacI::LVA
	reaction 7.		Cherry::LacI::LVA-2 --(k8)--> Cherry::LacI::LVA + Cherry::LacI::LVA		Dissociation of Cherry::LacI::LVA dimer
	reaction 8.		Ptac + Cherry::LacI::LVA-2 --(k9)--> Ptac-O2							Binding of Cherry::LacI::LVA dimer to tac promoter
	reaction 9.		Ptac-O2 + Cherry::LacI::LVA-2 --(k9)--> Ptac-O4							Binding of second Cherry::LacI::LVA dimer to tac promoter and tetramerization
	reaction 10.	Ptac-O4 --(k10)--> Ptac-O2 + Cherry::LacI::LVA-2						Dissociation of tetramer structure
	reaction 11.	Ptac --(k11)--> Ptac + AraC::GFP										Production of AraC::GFP fusion from free tac free promoter
	reaction 12.	Ptac-O2 --(k12)--> Ptac-O2 + AraC::GFP									Production of AraC::GFP fusion from Ptac promoter bound to a Cherry::LacI::LVA dimer
	reaction 13.	Ptac-O4 --(k12)--> Ptac-O4 + AraC::GFP									Production of AraC::GFP fusion from Ptac promoter bound to a Cherry::LacI::LVA tetramer
	reaction 14.	AraC::GFP --(k13)--> NULL												Degradation of AraC::GFP
  */


/**
 Propensity evaluation function for SynPI1.
 */
int synpi1_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop)
{
	// Check sizes of vectors
	if ((X->size != N) || (params->size != L+Z) || (prop->size != R))
	{
		fprintf (stderr, "error in synpi1_propensity_eval: vector sizes are not correct\n");
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

	// Input recovery statements
	double u1 = gsl_vector_get (params, 13);

	// Propensity evaluation statements
	gsl_vector_set (prop, 0, (k1)*X1*X3);
	gsl_vector_set (prop, 1, (k2)*X2);
	gsl_vector_set (prop, 2, (k3)*X1);
	gsl_vector_set (prop, 3, (k4)*X2);
	gsl_vector_set (prop, 4, (k5+k6*u1)*X7);
	gsl_vector_set (prop, 5, (k7)*X7*(X7-1));
	gsl_vector_set (prop, 6, (k8)*X8);
	gsl_vector_set (prop, 7, (k9)*X4*X8);
	gsl_vector_set (prop, 8, (k9)*X5*X8);
	gsl_vector_set (prop, 9, (k10)*X6);
	gsl_vector_set (prop, 10, (k11)*X4);
	gsl_vector_set (prop, 11, (k12)*X5);
	gsl_vector_set (prop, 12, (k12)*X6);
	gsl_vector_set (prop, 13, (k13)*X3);

	// Signal that computation was completed successfully
	return GSL_SUCCESS;
}


/**
 State update function for SynPI1.
 */
int synpi1_state_update (gsl_vector * X, size_t rxnid)
{
	// Check sizes of state vector
	if (X->size != N)
	{
		fprintf (stderr, "error in synpi1_state_update: state vector size is not correct\n");
		return GSL_EFAILED;
	}

	// Check that reaction id is correct
	if (rxnid >= R)
	{
		fprintf (stderr, "error in synpi1_state_update: reaction id is not correct\n");
		return GSL_EFAILED;
	}

	// Update the state vector according to which reaction fired
	switch (rxnid)
	{
	case 0:
		gsl_vector_set (X, 0, gsl_vector_get (X, 0) + (-1));
		gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (1));
		gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (-1));
		break;

	case 1:
		gsl_vector_set (X, 0, gsl_vector_get (X, 0) + (1));
		gsl_vector_set (X, 1, gsl_vector_get (X, 1) + (-1));
		gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (1));
		break;

	case 2:
		gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (1));
		break;

	case 3:
		gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (1));
		break;

	case 4:
		gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (-1));
		break;

	case 5:
		gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (-2));
		gsl_vector_set (X, 7, gsl_vector_get (X, 7) + (1));
		break;

	case 6:
		gsl_vector_set (X, 6, gsl_vector_get (X, 6) + (2));
		gsl_vector_set (X, 7, gsl_vector_get (X, 7) + (-1));
		break;

	case 7:
		gsl_vector_set (X, 3, gsl_vector_get (X, 3) + (-1));
		gsl_vector_set (X, 4, gsl_vector_get (X, 4) + (1));
		gsl_vector_set (X, 7, gsl_vector_get (X, 7) + (-1));
		break;

	case 8:
		gsl_vector_set (X, 4, gsl_vector_get (X, 4) + (-1));
		gsl_vector_set (X, 5, gsl_vector_get (X, 5) + (1));
		gsl_vector_set (X, 7, gsl_vector_get (X, 7) + (-1));
		break;

	case 9:
		gsl_vector_set (X, 4, gsl_vector_get (X, 4) + (1));
		gsl_vector_set (X, 5, gsl_vector_get (X, 5) + (-1));
		gsl_vector_set (X, 7, gsl_vector_get (X, 7) + (1));
		break;

	case 10:
		gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (1));
		break;

	case 11:
		gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (1));
		break;

	case 12:
		gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (1));
		break;

	case 13:
		gsl_vector_set (X, 2, gsl_vector_get (X, 2) + (-1));
		break;
	}

	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Sample a new random initial state for SynPI1.
 */
int synpi1_initial_conditions (gsl_vector * X0, const gsl_rng * r)
{
	// Check sizes of state vector
	if (X0->size != N)
	{
		fprintf (stderr, "error in synpi1_initial_conditions: state vector size is not correct\n");
		return GSL_EFAILED;
	}

	// Sample a new initial state
	gsl_vector_set (X0, 0, 20 + gsl_rng_uniform_int (r, 6) + gsl_rng_uniform_int (r, 6));
	gsl_vector_set (X0, 1, 0);
	gsl_vector_set (X0, 2, 0);
	gsl_vector_set (X0, 3, 20 + gsl_rng_uniform_int (r, 6) + gsl_rng_uniform_int (r, 6));
	gsl_vector_set (X0, 4, 0);
	gsl_vector_set (X0, 5, 0);
	gsl_vector_set (X0, 6, 0);
	gsl_vector_set (X0, 7, 0);

	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Output function for SynPI1.
 */
int synpi1_output (gsl_matrix * out)
{
	if ((out->size1 != P) || (out->size2 != N))
	{
		fprintf (stderr, "error in synpi1_output: output matrix size is not correct\n");
		return GSL_EFAILED;
	}

	// Reset the output matrix
	gsl_matrix_set_all (out, 0.0);

	// Set the non-zero terms
	//gsl_matrix_set (out, 0, 4, 2.0);
	//gsl_matrix_set (out, 0, 5, 4.0);
	gsl_matrix_set (out, 0, 6, 1.0);
	gsl_matrix_set (out, 0, 7, 2.0);

	// Signal that computation was completed correctly
	return GSL_SUCCESS;
}


/**
 Model information function for SynPI1.
 */
void synpi1_mod_setup (stochmod * model)
{
	model->propensity = &synpi1_propensity_eval;
	model->update = &synpi1_state_update;
	model->initial = &synpi1_initial_conditions;
	model->output = &synpi1_output;
	model->nspecies = N;
	model->nrxns = R;
	model->nparams = L;
	model->nin = Z;
	model->nout = P;
	model->name = "Synthetic PI version 1 (SYNPI1)";
}

