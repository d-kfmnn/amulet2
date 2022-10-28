/*------------------------------------------------------------------------*/
/*! \file slicing.h
    \brief contains functions slice and order the gates

  Part of AMulet2 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#ifndef AMULET2_SRC_SLICING_H_
#define AMULET2_SRC_SLICING_H_
/*------------------------------------------------------------------------*/
#include <vector>
#include <list>

#include "gate.h"
/*------------------------------------------------------------------------*/
// / vector-list Gate* matrix to store slices
extern std::vector<std::list<Gate*>> slices;
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
    Overall slicing routine based on xor gates, calls all other functions
    automatically
*/
bool slicing_xor();

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
    Identifies whether polynomials use a booth pattern(as in the aoki benchmarks)

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
    Undoes slicing attempt
*/
void clean_slices();

/**
    Overall slicing routine, that does not depend on xor-chains
*/
void slicing_non_xor();

#endif  // AMULET2_SRC_SLICING_H_
