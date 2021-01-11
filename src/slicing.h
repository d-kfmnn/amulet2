/*------------------------------------------------------------------------*/
/*! \file slicing.h
    \brief contains functions slice and order the gates

  Part of AMulet2.0 : AIG Multiplier Verification Tool.
  Copyright (C) 2020 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#ifndef _slicing_H
#define _slicing_H
/*------------------------------------------------------------------------*/
#include <vector>

#include "gate.h"
/*------------------------------------------------------------------------*/
/// vector-list Gate* matrix to store slices
extern std::vector<std::list<Gate*>> slices ;
/*------------------------------------------------------------------------*/

/**
    Allocates the memory for the slices
    and adds the output gate of each slice
*/
void init_slices();

/**
    Prints the polynomials in the slices to stdout
*/
void print_slices();


/**
    Returns the child that has the highest position in a slice

    @param n Gate* for parent

    @return Gate*
*/
const Gate * topological_largest_child(const Gate * n);

/**
    Move gate n to slice i and insert it before it's topological_largest_child.

    @param n Gate*
    @param i integer, for new slice index
*/
void fix_slice(Gate *n, int i);

/**
    Attach gates to slices by following xor-chains from the output gates.
*/
void slice_by_xor_chains();

/**
    Inserts the parents of n before pre in the slice of pre

    @param n Gate*
    @param pre Gate*
*/
void upwards_slicing(const Gate * n, const Gate * pre);

/**
    Slice gates that are between xor-chains, using upwards_slicing
*/
void slice_jut_gates();

/**
    Fix some xors, that are assigned to different slices
    (needed for 7,3 counter trees )
*/
int fix_xors();

/**
    If fix_xors() is applied also the jut gates need to be fixed.
*/
void fix_jut_gates();

/**
    Overall slicing routine based on xor gates, calls all other functions
    automatically
*/
void slicing_xor();

/*------------------------------------------------------------------------*/

/**
    Mark the children of gate n to slice num, if not assigned otherwise.

    @param n Gate*
    @param num integer for  slice index
*/
void input_cone(Gate * n, int num);

/**
    Identify those gates that are input in bigger slices
*/
void find_carries();

/**
    Identifies whether polynomials use a booth pattern (as in the aoki benchmarks)

    @return True if booth pattern are found
*/
bool search_for_booth_pattern();

/**
   We repeatedly move gates to smaller slices.
*/
void merge_all();

/**
    We repeatedly move gates to bigger slices.
*/
void promote_all();

/**
    Fills the slices by adding the gates that are assigned to slices
*/
void fill_slices();

/**
    Overall slicing routine, that does not depend on xor-chains
*/
void slicing_non_xor();

#endif
