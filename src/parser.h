/*------------------------------------------------------------------------*/
/*! \file parser.h
    \brief contains functions necessary to parse the AIG

  Part of AMulet2.1 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#ifndef AMULET2_SRC_PARSER_H_
#define AMULET2_SRC_PARSER_H_
/*------------------------------------------------------------------------*/
#include "aig.h"
/*------------------------------------------------------------------------*/
/**
    Checks whether 'model' contains an aiger node with
    output lhs and inputs rhs0 and rhs1.

    @param lhs unsigned integer
    @param rhs0 unsigned integer
    @param rhs1 unsigned integer

    @return True if 'model' contains such an aiger node
*/
bool match_and(unsigned lhs, unsigned rhs0, unsigned rhs1);

/**
    Identifies whether the input vectors are separated or interleaved.
*/
void determine_input_order();

/**
    Checks whether the input AIG fullfills our requirements.
*/
void init_aiger_with_checks();

/**
    Reads the input aiger given in the file called input_name to the aiger 'model'
    using the parserer function of aiger.h

    @param input_name char * ame of input file
*/
void parse_aig(const char * input_name);


#endif  // AMULET2_SRC_PARSER_H_
