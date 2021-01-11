/*------------------------------------------------------------------------*/
/*! \file polynomial_solver.cpp
    \brief contains the polynomial solving routine

  Part of AMulet2.0 : AIG Multiplier Verification Tool.
  Copyright (C) 2020 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include "polynomial_solver.h"
/*------------------------------------------------------------------------*/
// Global variable
bool gen_witness = 1;
/*------------------------------------------------------------------------*/

void init_gates_verify(){
  init_mpz(NN);
  allocate_gates();
  mark_aig_outputs();
  set_parents_and_children(1);
  init_gate_constraints();
  set_xor();
}

/*------------------------------------------------------------------------*/

void verify(const char * inp_f, const char * out_f1,
            const char * out_f2, const char * out_f3, bool certify){

  assert(!certify || inp_f);
  assert(!certify || out_f1);
  assert(!certify || out_f2);
  assert(!certify || out_f3);

  FILE * f1 = 0, *f2 = 0, *f3 = 0;
  if(certify){
    if (!(f1 = fopen (out_f1, "w")))
    die ("can not write output to '%s'", out_f1);

    if (!(f2 = fopen (out_f2, "w")))
    die ("can not write output to '%s'", out_f2);

    if (!(f3 = fopen (out_f3, "w")))
    die ("can not write output to '%s'", out_f3);
  }

  if(certify && pac){
    print_circuit_poly(f1);
    print_spec_poly(f3);
  }
  else if(certify && nss){
    print_circuit_poly_nss(f1);
    print_spec_poly(f3);
  }

  init_slices();
  mark_xor_chain_in_last_slice();
  init_time = process_time();

  remove_internal_xor_gates(f2);
  if(!upper_half_xor_output()){
    msg("slicing based on input cones");
    xor_chain = 1;
    slicing_non_xor();

    if(search_for_booth_pattern()) eliminate_booth_pattern(f2);
    decomposing(f2);

  } else {
    msg("slicing based on xor");
    remove_single_occs_gates(f2);
    slicing_xor();
    remove_slice_minus_one_gates(f2);
  }

  if(verbose >= 2) {
    msg("after preprocessing: ");
    print_slices();
  }
  slicing_elim_time = process_time();

  const Polynomial * rem = reduce(f2);
  if(rem && !rem->is_constant_zero_poly())  {
    msg ("INCORRECT MULTIPLIER");
    msg("");

    if(inp_f && gen_witness) {
      msg("REMAINDER IS");
      fputs("[amulet2] ", stdout);
      rem->print(stdout);
      msg("");
      generate_witness(rem, inp_f);
    }
  }
  else {
    msg("");
    msg("CORRECT MULTIPLIER");
    if(certify){
      msg("");
      msg("writing gate constraints to '%s' ", out_f1);
      msg("writing proof certificate to '%s'", out_f2);
      msg("writing specification to '%s'    ", out_f3);
    }

  }
  delete(rem);

  if(nss) print_cofactors_poly_nss (f1, f2);

  reduction_time = process_time();
  if(certify) {
    fclose(f1);
    fclose(f2);
    fclose(f3);
  }
}
