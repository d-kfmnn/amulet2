/*------------------------------------------------------------------------*/
/*! \file substitution_engine.h
    \brief contains the substitution engine

  Part of AMulet2.0 : AIG Multiplier Verification Tool.
  Copyright (C) 2020 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#ifndef _substitution_engine_H
#define _substitution_engine_H
/*------------------------------------------------------------------------*/
#include "substitution.h"
/*------------------------------------------------------------------------*/
/**
    Initializes the internal gate structure, with necessary informations
    for substitution.
    Uses the AIG module
*/
void init_gate_substitution();

/**
    Calls the substitution routine.
    We first detect the GP adder, then generate the equivalent miter
    and print the CNF miter and the rewritten AIG to the output files.

    @param out_f1 name of first output file
    @param out_f2 name of second output file
*/
void substitution(const char * out_f1, const char * out_f2);


#endif
