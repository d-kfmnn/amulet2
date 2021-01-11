/*------------------------------------------------------------------------*/
/*! \file elimination.h
    \brief contains functions used in the polynomial solver

  This file contains all functions used for preprocessing the Gröbner basis
  and for reducing the specification by the slices.

  Part of AMulet2.0 : AIG Multiplier Verification Tool.
  Copyright (C) 2020 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#ifndef _elimination_H
#define _elimination_H
/*------------------------------------------------------------------------*/
#include <string.h>

#include "nss.h"
#include "pac.h"
#include "slicing.h"
/*------------------------------------------------------------------------*/
/// True, when the PAC proof format is involved
extern bool pac;

/// True, when the Nullstellensatz proof format is involved
extern bool nss;
/*------------------------------------------------------------------------*/
// Functions used for PAC proofs

/**
    Adds up the computed factors of a slice to compute the sliced specification
    Prints PAC rules for the process. Used only when pac=1.

    @param file output file for PAC rules

    @return polynomial, which is the sum of all computed factors of a slice
*/
Polynomial * add_up_factors(FILE * file);


/**
    Merges computed factors of a slice of the same level.
    Prints PAC rules for the process. Used only when pac=1.

    @param file output file for PAC rules
    @param container vector containing the factors, ordered by level

    @return vector of polynomials such that each polynomial has a unique level
*/
std::vector<Polynomial*> merge_factors
    (FILE * file, std::vector<Polynomial*> container);


/**
    Adds up the computed sliced specifications.
    Prints PAC rules for the process. Used only when pac=1.

    @param file output file for PAC rules
*/
void add_up_spec_of_slice(FILE * file);
/*------------------------------------------------------------------------*/
// These functions are used to eliminate by one gate constraint

/**
    Reduces the gate constraint of n1 by the gate constraint of n2

    @param n1 Gate that is reduced by n2
    @param n2 Gate
    @param file output file for PAC rules
*/
void eliminate_by_one_gate(Gate * n1, Gate *n2, FILE * file);

/**
    Reduces p1 by the gate constraint of n

    @param p1 Polynomial to be reduced
    @param n  Gate
    @param file output file for PAC rules

    @return the computed polynomial remainder
*/
Polynomial * reduce_by_one_poly(const Polynomial * p1, Gate * n, FILE * file);
/*------------------------------------------------------------------------*/
// Variable elimination heuristics in XOR-based slicing

/**
    Remove all internal gates in identified XOR gates

    @param file output file for PAC rules
*/
void remove_internal_xor_gates(FILE * file);

/**
    Remove all single-parent gates (only one loop)

    @param file output file for PAC rules
*/
void remove_single_occs_gates(FILE * file);

/**
    Remove all gates that have not been assigned to a slice

    @param file output file for PAC rules
*/
void remove_slice_minus_one_gates(FILE * file);

/*------------------------------------------------------------------------*/
// Variable elimination heuristics taken from AMulet1.0

/**
    Remove repeatedly all single gates (until completion)

    @param file output file for PAC rules
*/
void decomposing(FILE * file);


/**
    Remove all gates that are identified as a Booth pattern

    @param file output file for PAC rules
*/
void eliminate_booth_pattern(FILE * file);

/*------------------------------------------------------------------------*/

/**
    Defines the slice specification for slice i: -s_i + PP_i

    @param i unsigned integer for the index of the slice

    @return Polynomial for slice specification
*/
Polynomial * inc_spec_poly(unsigned i);

/**
    Reduces p1 by the constant 2^NN

    @param p1 Polynomial to be reduced
    @param print_rule bool indicating whether we print the PAC rules
    @param file output file for PAC rules

    @return the computed polynomial remainder
*/
Polynomial * mod_poly(const Polynomial *p1, bool print_rule, FILE * file);

/**
    Reduces the computed unsigned specification by the modulo constant

    @param p Polynomial to be reduced
    @param file output file for PAC rules
*/
void correct_pp_unsigned(const Polynomial * p, FILE * file);

/**
    Reducing the computed signed specification by the modulo constant

    @param p Polynomial to be reduced
    @param file output file for PAC rules
*/
void correct_pp_signed(const Polynomial * p, FILE * file);

/**
    Reducing the computed specification by the modulo constant,
    by calling appropriate subfunctions.
    @see correct_pp_unsigned(const Polynomial * p, FILE * file)
    @see correct_pp_signed(const Polynomial * p, FILE * file)

    @param p Polynomial to be reduced
    @param file output file for PAC rules
*/
void correct_pp(const Polynomial * p, FILE * file);

/*------------------------------------------------------------------------*/
/**
    Incremental verification algorithm, where we reduce the slice-wise
    specifcations by the gate constraints of the slices

    @param file output file for PAC rules

    @return the computed polynomial remainder, 0 if the circuit is correct

*/
const Polynomial * reduce(FILE * file);
/*------------------------------------------------------------------------*/
// Functions used to find counter examples
/**
    Checks whether the polynomial contains only input variables

    @param p Polynomial* to be checked

    @return true if p contains only input variables
*/
bool check_inputs_only(const Polynomial *p);

/**
    Writes a vector to the file such that all input variables in the term t
    are set to 1, all other to 0.

    @param t Term
    @param file output file for the counter example
*/
void write_witness_vector(const Term * t, FILE * file);

/**
    Search for the smallest term and calls write_witness_vector on that

    @param p Polynomial* for which counter examples are generated
    @param file output file for the counter example
*/
void write_witnesses(const Polynomial * p, FILE * file);

/**
    Generates a witness for the remainder polynomial p, by identifying
    the smallest term in the polynomial and setting all its variables
    to 1.

    @param p Polynomial* for which counter examples are generated
    @param name prefix name for the output file, suffix is '.cex'
*/
void generate_witness(const Polynomial * p, const char * name);

#endif
