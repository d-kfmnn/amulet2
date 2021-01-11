/*------------------------------------------------------------------------*/
/*! \file nss.h
    \brief contains functions necessary to generate Nullstellensatz proofs

  Part of AMulet2.0 : AIG Multiplier Verification Tool.
  Copyright (C) 2020 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#ifndef _nss_H
#define _nss_H
/*------------------------------------------------------------------------*/
#include "gate.h"
/*------------------------------------------------------------------------*/

/**
    Prints the circuit poly to the file (without indices)

    @param file output file
*/
void print_circuit_poly_nss (FILE * file);

/**
    Prints the cofactors of the circuit polynomials

    @param polysfile output file for gate constraints
    @param file output file for cofactors
*/
void print_cofactors_poly_nss (FILE* polysfile, FILE * file);


/**
    Adding an ancestor polynomial to the ancestors of n

    @param n      Gate* to which ancestor is added
    @param anc    Gate* to Ancestor gate
    @param fac    Polynomial* definng ancestor polynomial
    @param depth  true if function is called internally (default = 0)
*/
void add_ancestors (Gate * n, Gate * anc, const Polynomial * fac, bool internal = 0);

/**
    Updates the cofactor of the modulo polynomial

    @param fac Polynomial* which is added to the cofactor of the modpoly
*/
void add_fac_mod (const Polynomial * fac);

/**
    Updates the cofactor of the gate n

    @param n      Gate* for which the cofactor is updated
    @param fac    Polynomial* which is added to the cofactor of n
*/
void add_fac (Gate * n, const Polynomial * fac);



#endif
