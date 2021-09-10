/*------------------------------------------------------------------------*/
/*! \file nss.h
    \brief contains functions necessary to generate Nullstellensatz proofs

  Part of AMulet2.1 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#ifndef AMULET2_SRC_NSS_H_
#define AMULET2_SRC_NSS_H_
/*------------------------------------------------------------------------*/
#include "gate.h"
/*------------------------------------------------------------------------*/
/**
    Prints the specification polynomial to the file

    @param file output file
*/
void print_spec_poly(FILE * file);
/**
    Prints the circuit poly to the file(without indices)

    @param file output file
*/
void print_circuit_poly_nss(FILE * file);

/**
    Prints the cofactors of the circuit polynomials

    @param file output file for cofactors
*/
void print_cofactors_poly_nss(FILE * file);


/**
    Adding an ancestor polynomial to the ancestors of n

    @param n      Gate* to which ancestor is added
    @param anc    Gate* to Ancestor gate
    @param fac    Polynomial* definng ancestor polynomial
    @param depth  true if function is called internally(default = 0)
*/
void add_ancestors(
  Gate * n, Gate * anc, const Polynomial * fac, bool internal = 0);

/**
    Updates the cofactor of the modulo polynomial

    @param fac Polynomial* which is added to the cofactor of the modpoly
*/
void add_fac_mod(const Polynomial * fac);

/**
    Updates the cofactor of the gate n

    @param n      Gate* for which the cofactor is updated
    @param fac    Polynomial* which is added to the cofactor of n
*/
void add_fac(Gate * n, const Polynomial * fac);

#endif  // AMULET2_SRC_NSS_H_
