/*------------------------------------------------------------------------*/
/*! \file polynomial_solver.cpp
    \brief contains the polynomial solving routine

  Part of AMulet2 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include "polynomial_solver.h"
/*------------------------------------------------------------------------*/
// Global variable
bool gen_witness = 1;
/*------------------------------------------------------------------------*/
// ERROR CODES:
static int err_writing = 41; // cannot write to
static int err_rem_poly =42; // remainder poly no witness
/*------------------------------------------------------------------------*/
void init_gates_verify() {
  init_mpz(NN);
  allocate_gates();
  mark_aig_outputs();
  set_parents_and_children(1);
  set_xor();
}

/*------------------------------------------------------------------------*/

bool verify(const char * inp_f, const char * out_f1,
            const char * out_f2, const char * out_f3, bool certify) {
  assert(!certify || inp_f);
  assert(!certify || out_f1);
  assert(!certify || out_f2);
  assert(!certify || out_f3);

  FILE * f1 = 0, *f2 = 0, *f3 = 0;
  if (certify) {
    if (!(f1 = fopen(out_f1, "w")))
    die(err_writing, "can not write output to '%s'", out_f1);

    if (!(f2 = fopen(out_f2, "w")))
    die(err_writing, "can not write output to '%s'", out_f2);

    if (!(f3 = fopen(out_f3, "w")))
    die(err_writing, "can not write output to '%s'", out_f3);
  }

  if (certify) {
    print_circuit_poly(f1);
    print_spec_poly(f3);
  }

  init_slices();

  mark_xor_chain_in_last_slice();
  init_time = process_time();

  remove_internal_xor_gates(f2);
  bool non_xor_slice = upper_half_xor_output();
  bool xor_slicing_failure = 0;
  if(non_xor_slice){
    msg("slicing based on xor");
    remove_single_occs_gates(f2);

    xor_slicing_failure = slicing_xor();
    if(!xor_slicing_failure) remove_slice_minus_one_gates(f2);
    else clean_slices();
  }

  if (!non_xor_slice || xor_slicing_failure) {
    msg("slicing based on input cones");
    xor_chain = 1;
    slicing_non_xor();

    if (search_for_booth_pattern()) eliminate_booth_pattern(f2);
    decomposing(f2);

  }



  slicing_elim_time = process_time();

  const Polynomial * rem = reduce(f2);
  bool res;
  if (rem && !rem->is_constant_zero_poly())  {
    if (!check_inputs_only(rem)){
      msg("REMAINDER IS");
      fputs("[amulet2] ", stdout);
      rem->print(stdout);
      msg("");
      die(err_rem_poly, "slicing failure - remainder polynomial contains non-inputs");

    }

    res = 0;
    msg("INCORRECT MULTIPLIER");
    msg("");

    if (inp_f && gen_witness) {
      msg("REMAINDER IS");
      fputs("[amulet2] ", stdout);
      rem->print(stdout);
      msg("");
      generate_witness(rem, inp_f);
    }
  } else {
    res = 1;
    msg("");
    msg("CORRECT MULTIPLIER");
    if (certify) {
      msg("");
      msg("writing gate constraints to '%s' ", out_f1);
      msg("writing proof certificate to '%s'", out_f2);
      msg("writing specification to '%s'    ", out_f3);
    }
  }
  delete(rem);

  if (proof == 3) print_cofactors_poly_nss(f2);

  reduction_time = process_time();
  if (certify) {
    fclose(f1);
    fclose(f2);
    fclose(f3);
  }
  return res;
}
