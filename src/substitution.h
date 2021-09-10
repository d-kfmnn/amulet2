/*------------------------------------------------------------------------*/
/* AMulet2 : AIG Multiplier Verification Tool.
 *
 * Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
 *
 * substitution.h: functions used for adder substitution
 */
/*------------------------------------------------------------------------*/
#ifndef AMULET2_SRC_SUBSTITUTION_H_
#define AMULET2_SRC_SUBSTITUTION_H_
/*------------------------------------------------------------------------*/
#include <vector>

#include "gate.h"
/*------------------------------------------------------------------------*/

/**
    Checks whether all AIG outputs are used only once

    @return True if all AIG outputs are single gates
*/
bool all_single_output();

/**
    Checks whether all AIG outputs are XOR gates

    @return True if all outputs in slices 1 to NN-2 are XOR gates
*/
bool all_outputs_are_xor();

/**
    If the output of slice 2 has more than 3 parents, than the carry-in
    of the final stage is the output of slice 0(in the aoki benchmarks.)

    @return False if output of slice 2 has more than 3 parents but output
    of slice 0 is single gate
*/
bool slice_two_needs_carry_in_slice_zero();

/**
    Checks whether the output of slice 0 is a possible carry-in for the
    final stage adder

    @return True if output of slice 0 has more than one parent.
*/
bool cin_in_slice_0();

/**
    Returns the aiger value of last element of the inputs vector

    @param flip indicates whether the input has to be negated

    @return unsigned
*/
unsigned get_input(bool flip = 0);

/**
    Adds the gate n to the inputs vector

    @param n Gate* to be pushed
*/
void push_to_inputs(Gate * n);

/**
    Adds the gate n to the outputs vector and sets it's slice to i

    @param n Gate* to be pushed
    @param i integer indicating the slice index
*/
void push_to_outputs(Gate * n, int i);

/**
    Adds the gate n to the cins vector and sets it's slice to i

    @param n Gate* to be pushed
    @param i integer indicating the slice index
*/
void push_to_cins(Gate * n, int i);

/**
    Sets the gate n to the carry_in of the final stage adder

    @param n Gate*
*/
void set_carry_in(Gate * n);

/**
    Identifies the carry out of the final stage adder
*/
void identify_carry_out();

/**
    Identifies the propagate and generate gates of the final stage adder

    @return True if all checks succeeded
*/
bool identify_propagate_and_generate_gates();

/**
    For some multipliers we have to adjust the inputs, because they are sometimes
    propagate or generate gate themselves.
*/
void fix_inputs();

/**
    Follow all paths from the gate n and check that we stop at
    the identified inputs and carry-in of the final stage adder.
    All gates that are seen on the way are marked to belong to the final stage
    adder.

    @param n  Gate* pointing to root gate, from which we start the checks
    @param init: bool indicating whether we started from the carry out

    @return True if all checks succeeded
*/
bool follow_path_and_mark_gates(Gate * n, bool init = 0);

/**
    Calls follow_path_and_mark_gates for all identified output gates of
    the final stage adder.

    @return True if all checks succeeded
*/
bool follow_all_output_paths_and_mark_gates();

/**
    For all input gates we count how often they are used as inputs in the
    final stage adder.
*/
void correctly_mark_inputs();

/**
    Routine for identifying the final stage adder

    @return True if all checks succeeded
*/
bool identify_final_stage_adder();

/*------------------------------------------------------------------------*/
/**
    Adds the identified final stage adder to the miter
*/
void add_original_adder();

/**
    A vector called original_outputs is filled by the outputs of the
    final stage adder. These gates are then later used, when the miter
    is generated.
*/
void fill_original_outputs();

/**
    Generates an AIG of a half-adder with inputs i1 and i2.
    These are added to the rewritten AIG as well as to the miter.

    @param i1 unsigned, input 1 of the half adder
    @param i2 unsigned, input 2 of the half adder
    @param carry: bool indicating whether we add the carry of the half-adder
                  to the rewritten_outputs

    @return the sum output if carry = 0, and the carry output if carry = 1
*/
unsigned btor_ha(unsigned i1, unsigned i2, bool carry = 1);

/**
    Generates an AIG of a half-adder with inputs i1, i2 and i3.
    These are added to the rewritten AIG as well as to the miter.

    @param i1 unsigned, input 1 of the full adder
    @param i2 unsigned, input 2 of the full adder
    @param i3 unsigned, input 3 of the full adder
    @param carry: bool indicating whether we add the carry of the full-adder
                  to the rewritten_outputs

    @return the sum output if carry = 0, and the carry output if carry = 1
*/
unsigned btor_fa(unsigned i1, unsigned i2, unsigned i3, bool carry = 1);

/**
    Generates a ripple-carry adder using the inputs, stored in the
    inputs vector. The outpus of this adder are added to the rewritten_outputs
*/
void add_btor_adder();

/**
    Flips the last bit of a to negate the unsigned integer as
    used in the aiger format

    @param a unsigned integer

    @return the negated aiger value of a
*/
unsigned not_(unsigned a);

/**
    Generates an AIG gate with inputs a and b

    @param a unsigned integer
    @param b unsigned integer

    @return the value of the output gate
*/
unsigned and_(unsigned a, unsigned b);

/**
    Generates AIG gates for a implies b, using and_ and not_

    @param a unsigned integer
    @param b unsigned integer

    @return the value of the output gate
*/
unsigned implies_(unsigned a, unsigned b);

/**
    Generates AIG gates for a xnor b, using implies_ and and_

    @param a unsigned integer
    @param b unsigned integer

    @return the value of the output gate
*/
unsigned xnor_(unsigned a, unsigned b);

/**
    Builds the miter, where the the original_outputs and the rewritten_outputs
    are combined elementwise using a sequence of xnor gates

    @return True if all checks succeeded
*/
bool build_miter();

/**
    Routine to generate the miter aig

    @return True if all checks succeeded
*/
bool build_adder_miter();

/**
    Translate the miter aig to CNF and prints it to the file

    @param file Output file

    @return True if all checks succeeded
*/
bool miter_to_file(FILE * file);

/**
    Prints the trivial CNF  "a and not a" to the file

    @param file Output file

    @return True if all checks succeeded
*/
bool trivial_miter_to_file(FILE * file);
/*------------------------------------------------------------------------*/


/**
    Generates the rewritten aig by adding all gates that are not marked as
    the final stage adder.
    The gates of the ripple carry adder have been already added in add_btor_adder.
*/
void generate_rewritten_aig();



#endif  // AMULET2_SRC_SUBSTITUTION_H_
