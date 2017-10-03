/*
 *  stochmod.h
 *  StochMod
 *
 *  Library-wide Header File
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

#ifndef _STOCHMOD_H_
#define _STOCHMOD_H_


/*
 System libraries includes
 */

#include <stdio.h>
#include <math.h>


/*
 GSL library includes
 */

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>



/*
 New data types
 */

// Model struct
typedef struct {
	int (* propensity) (const gsl_vector *, const gsl_vector *, gsl_vector *);
	int (* update) (gsl_vector *, size_t);
	int (* initial) (gsl_vector *, const gsl_rng *);
	int (* output) (gsl_matrix *);
	size_t nspecies;
	size_t nrxns;
	size_t nparams;
	size_t nin;
	size_t nout;
	char * name;
} stochmod;

// Enumeration for the models contained in the library
typedef enum {
	MODEL_SYNCIRC = 0,
	MODEL_STOCHREP = 1,
	MODEL_AUTOREG = 2,
	MODEL_LACGFP = 3,
	MODEL_LACGFP2 = 4,
	MODEL_LACGFP3 = 5,
	MODEL_LACGFP4 = 6,
	MODEL_LACGFP5 = 7,
	MODEL_BIRTHDEATH = 8,
	MODEL_LACGFP6 = 9,
	MODEL_LACGFP7 = 10,
	MODEL_LACGFP8 = 11,
	MODEL_IFF = 12,
	MODEL_FBK = 13,
	MODEL_LACGFP9 = 14,
	MODEL_LACGFP10 = 15,
	MODEL_SYNPI1 = 16,
} STOCHASTIC_MODEL;


/*
 Exported functions prototype declarations == SYNCIRC.C
 */
int syncirc_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop);
int syncirc_state_update (gsl_vector * X, size_t rxnid);
void syncirc_mod_setup (stochmod * model);


/*
 Exported functions prototype declarations == STOCHREP.C
 */
int stochrep_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop);
int stochrep_state_update (gsl_vector * X, size_t rxnid);
void stochrep_mod_setup (stochmod * model);


/*
 Exported functions prototype declarations == AUTOREG.C
 */
int autoreg_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop);
int autoreg_state_update (gsl_vector * X, size_t rxnid);
int autoreg_initial_conditions (gsl_vector * X0, const gsl_rng * r);
void autoreg_mod_setup (stochmod * model);


/*
 Exported functions prototype declarations == LACGFP.C
 */
int lacgfp_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop);
int lacgfp_state_update (gsl_vector * X, size_t rxnid);
int lacgfp_initial_conditions (gsl_vector * X0, const gsl_rng * r);
int lacgfp_output (gsl_matrix * out);
void lacgfp_mod_setup (stochmod * model);


/*
 Exported functions prototype declarations == LACGFP2.C
 */
int lacgfp2_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop);
int lacgfp2_state_update (gsl_vector * X, size_t rxnid);
int lacgfp2_initial_conditions (gsl_vector * X0, const gsl_rng * r);
int lacgfp2_output (gsl_matrix * out);
void lacgfp2_mod_setup (stochmod * model);


/*
 Exported functions prototype declarations == LACGFP3.C
 */
int lacgfp3_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop);
int lacgfp3_state_update (gsl_vector * X, size_t rxnid);
int lacgfp3_initial_conditions (gsl_vector * X0, const gsl_rng * r);
int lacgfp3_output (gsl_matrix * out);
void lacgfp3_mod_setup (stochmod * model);


/*
 Exported functions prototype declarations == LACGFP4.C
 */
int lacgfp4_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop);
int lacgfp4_state_update (gsl_vector * X, size_t rxnid);
int lacgfp4_initial_conditions (gsl_vector * X0, const gsl_rng * r);
int lacgfp4_output (gsl_matrix * out);
void lacgfp4_mod_setup (stochmod * model);


/*
 Exported functions prototype declarations == LACGFP5.C
 */
int lacgfp5_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop);
int lacgfp5_state_update (gsl_vector * X, size_t rxnid);
int lacgfp5_initial_conditions (gsl_vector * X0, const gsl_rng * r);
int lacgfp5_output (gsl_matrix * out);
void lacgfp5_mod_setup (stochmod * model);


/*
 Exported functions prototype declarations == BIRTHDEATH.C
 */
int birthdeath_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop);
int birthdeath_state_update (gsl_vector * X, size_t rxnid);
int birthdeath_initial_conditions (gsl_vector * X0, const gsl_rng * r);
int birthdeath_output (gsl_matrix * out);
void birthdeath_mod_setup (stochmod * model);

\
/*
 Exported functions prototype declarations == LACGFP6.C
 */
int lacgfp6_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop);
int lacgfp6_state_update (gsl_vector * X, size_t rxnid);
int lacgfp6_initial_conditions (gsl_vector * X0, const gsl_rng * r);
int lacgfp6_output (gsl_matrix * out);
void lacgfp6_mod_setup (stochmod * model);


/*
 Exported functions prototype declarations == LACGFP7.C
 */
int lacgfp7_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop);
int lacgfp7_state_update (gsl_vector * X, size_t rxnid);
int lacgfp7_initial_conditions (gsl_vector * X0, const gsl_rng * r);
int lacgfp7_output (gsl_matrix * out);
void lacgfp7_mod_setup (stochmod * model);


/*
 Exported functions prototype declarations == LACGFP8.C
 */
int lacgfp8_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop);
int lacgfp8_state_update (gsl_vector * X, size_t rxnid);
int lacgfp8_initial_conditions (gsl_vector * X0, const gsl_rng * r);
int lacgfp8_output (gsl_matrix * out);
void lacgfp8_mod_setup (stochmod * model);


/*
 Exported functions prototype declarations == IFF.C
 */
int iff_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop);
int iff_state_update (gsl_vector * X, size_t rxnid);
int iff_initial_conditions (gsl_vector * X0, const gsl_rng * r);
int iff_output (gsl_matrix * out);
void iff_mod_setup (stochmod * model);


/*
 Exported functions prototype declarations == FBK.C
 */
int fbk_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop);
int fbk_state_update (gsl_vector * X, size_t rxnid);
int fbk_initial_conditions (gsl_vector * X0, const gsl_rng * r);
int fbk_output (gsl_matrix * out);
void fbk_mod_setup (stochmod * model);


/*
 Exported functions prototype declarations == LACGFP9.C
 */
int lacgfp9_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop);
int lacgfp9_state_update (gsl_vector * X, size_t rxnid);
int lacgfp9_initial_conditions (gsl_vector * X0, const gsl_rng * r);
int lacgfp9_output (gsl_matrix * out);
void lacgfp9_mod_setup (stochmod * model);


/*
 Exported functions prototype declarations == LACGFP10.C
 */
int lacgfp10_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop);
int lacgfp10_state_update (gsl_vector * X, size_t rxnid);
int lacgfp10_initial_conditions (gsl_vector * X0, const gsl_rng * r);
int lacgfp10_output (gsl_matrix * out);
void lacgfp10_mod_setup (stochmod * model);

/*
 Exported functions prototype declarations == SYNPI1.C
 */
int synpi1_propensity_eval (const gsl_vector * X, const gsl_vector * params, gsl_vector * prop);
int synpi1_state_update (gsl_vector * X, size_t rxnid);
int synpi1_initial_conditions (gsl_vector * X0, const gsl_rng * r);
int synpi1_output (gsl_matrix * out);
void synpi1_mod_setup (stochmod * model);


#endif
