/*------------------------------------------------------------------------*/
/*! \file aig.h
    \brief contains functions to parse and manipulate the input AIG

  Part of AMulet2.1 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#ifndef AMULET2_SRC_AIG_H_
#define AMULET2_SRC_AIG_H_
/*------------------------------------------------------------------------*/
#include <assert.h>

#include "signal_statistics.h"

extern "C" {
  #include "../includes/aiger.h"
}
/*------------------------------------------------------------------------*/

extern unsigned M;   // /< stores the maximum variable num of the input AIG
extern unsigned NN;  // /< stores the number of inputs of the input AIG

extern unsigned a0;    // /< input value for the LSB of input vector A
extern unsigned al;    // /< input value for the MSB of input vector A
extern unsigned ainc;  // /< distance between aiger values of A

extern unsigned b0;    // /< input value for the LSB of input vector B
extern unsigned bl;    // /< input value for the MSB of input vector B
extern unsigned binc;  // /< distance between aiger values of B


// / aiger * storing the generated miter during '-substitution'
// / Will be transformed into CNF and printed to the provided output file.
extern aiger * miter;

// / aiger * storing the generated rewritten AIG during '-substitution'
// / Will be printed to the provided output file.
extern aiger * rewritten;

/*------------------------------------------------------------------------*/

/**
    Initializes the 'aiger* miter' and 'aiger* rewritten'
*/
void init_aig_substitution();

/**
    Initializes the 'aiger* model', which is local to aig.cpp
*/
void init_aig_parsing();

/*------------------------------------------------------------------------*/

/**
    Resets the 'aiger* miter' and 'aiger* rewritten'
*/
void reset_aig_substitution();

/**
    Resets the 'aiger* model', which is local to aig.cpp
*/
void reset_aig_parsing();
/*------------------------------------------------------------------------*/
// Functions that interfer with aiger* model, that is used to store the
// input AIG.

/**
    Opens the input file and reads the contents to 'aiger* model'

    @param input_name a const char* refering to the name of the input file

    @return const char*, defining a possible error message
            Equal to zero if everything went right.
   )
*/
const char * aiger_open_and_read_to_model(const char * input_name);

/**
    Checks whether the given value corrensponds to an input of 'aiger* model'

    @param val an unsigned integer

    @return true when val is an input in 'aiger* model'.
*/
bool is_model_input(unsigned val);

/**
    Searches for the AIG node with value 'val' in 'aiger* model'

    @param val an unsigned integer

    @return an aiger_and* object, with lefthand side val
*/
aiger_and * is_model_and(unsigned val);

/**
    Returns the number of latches in 'aiger* model'

    @return an unsigned integer
*/
unsigned get_model_num_latches();

/**
    Returns the number of inputs in 'aiger* model'

    @return an unsigned integer
*/
unsigned get_model_num_inputs();

/**
    Returns the number of AND nodes in 'aiger* model'

    @return an unsigned integer
*/
unsigned get_model_num_ands();

/**
    Returns the number of outputs in 'aiger* model'

    @return an unsigned integer
*/
unsigned get_model_num_outputs();

/**
    Returns the maximum variable number of 'aiger* model'

    @return an unsigned integer
*/
unsigned get_model_maxvar();

/**
    Returns the value of the i'th input in 'aiger* model'

    @param i an unsigned integer, has to be smaller than NN

    @return an unsigned integer corresponding to the value of the i'th input
*/
unsigned get_model_inputs_lit(unsigned i);

/**
    Returns the name of the i'th input in 'aiger* model'

    @param i an unsigned integer, has to be smaller than NN

    @return a const char * corresponding to the name of the i'th input
*/
const char* get_model_inputs_name(unsigned i);


/**
    Returns the aiger value of the i'th output in 'aiger* model'

    @param i an unsigned integer, has to be smaller than NN

    @return an unsigned integer corresponding to the value of the i'th output
*/
unsigned slit(unsigned i);


/**
    Writes the 'aiger* model' to the provided file.

    @param file output file

    @return an integer value, zero if everything went right.
*/
int write_model(FILE *file);

#endif  // AMULET2_SRC_AIG_H_
