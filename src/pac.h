/*------------------------------------------------------------------------*/
/*! \file pac.h
    \brief contains functions necessary to generate PAC proofs

  Part of AMulet2.0 : AIG Multiplier Verification Tool.
  Copyright (C) 2020 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#ifndef _pac_H
#define _pac_H
/*------------------------------------------------------------------------*/

#include "gate.h"
/*------------------------------------------------------------------------*/
/**
    Prints the specification polynomial to the file

    @param file output file
*/
void print_spec_poly (FILE * file);

/**
    Prints all initial gate constraints to the file (with indices)

    @param file output file
*/
void print_circuit_poly (FILE * file);

/**
    Prints a deletion rule

    @param p1 Polynomial*
    @param file output file
*/
void print_pac_del_rule (FILE * file, const Polynomial *p1);

/**
    Prints the modulo rule

    @param file output file
    @param p1   Polynomial*, factor
    @param p    Polynomial*, conclusion
*/
void print_pac_mod_rule (FILE * file, const Polynomial *p1, Polynomial *p);

/**
    Prints an addition rule for pac

    @param file output file
    @param p1   Polynomial*, First summand
    @param p2   Polynomial*, Second summand
    @param p    Polynomial*, Conclusion
*/
void print_pac_add_rule (FILE * file, const Polynomial *p1, const Polynomial *p2, Polynomial *p);

/**
    Prints the multiplication rule of pac

    @param file output file
    @param p1   Polynomial*, First factor
    @param p2   Polynomial*, Second factor
    @param p    Polynomial*, Conclusion
*/
void print_pac_mul_rule (FILE * file, const Polynomial *p1, const Polynomial *p2, Polynomial *p);



#endif
