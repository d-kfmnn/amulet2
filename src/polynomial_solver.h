/*------------------------------------------------------------------------*/
/*! \file polynomial_solver.h
    \brief contains the polynomial solving routine

  Part of AMulet2.1 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#ifndef AMULET2_SRC_POLYNOMIAL_SOLVER_H_
#define AMULET2_SRC_POLYNOMIAL_SOLVER_H_
/*------------------------------------------------------------------------*/
#include "elimination.h"
/*------------------------------------------------------------------------*/
// / If final remainder is not equal to zero a counter example is generated and
// / printed to file <input_name>.wit, default is true, can be turned of
//  using command line input
extern bool gen_witness;
/*------------------------------------------------------------------------*/

/**
    Initializes the internal gate structure, with necessary informations
    for verification.
    Uses the AIG module
*/
void init_gates_verify();

/**
    Calls the preprocessing, slicing, reduction routines
    If certify is true, files need to be provided to store the proof.

    @param inp_f name of input file
    @param out_f1 name of first output file
    @param out_f2 name of second output file
    @param out_f3 name of third output file
    @param certify true when modus -certify is used
*/
void verify(const char * inp_f = 0,
            const char * out_f1 = 0,
            const char * out_f2 = 0,
            const char * out_f3 = 0,
            bool certify = 0);


#endif  // AMULET2_SRC_POLYNOMIAL_SOLVER_H_
