/*------------------------------------------------------------------------*/
/*! \file substitution_engine.cpp
    \brief contains the substitution engine

  Part of AMulet2 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include "substitution_engine.h"

/*------------------------------------------------------------------------*/
// ERROR CODES:
static int err_write_file  = 51; // error write to file
static int err_write_aig   = 52; // error write aig
/*------------------------------------------------------------------------*/
void init_gate_substitution() {
  allocate_gates();
  mark_aig_outputs();
  set_parents_and_children(0);
  set_xor();
}

/*------------------------------------------------------------------------*/

bool substitution(const char * out_f1, const char * out_f2) {
  assert(out_f1);
  assert(out_f2);

  FILE * f1;
  if (!(f1 = fopen(out_f1, "w")))
      die(err_write_file, "can not write output to '%s'", out_f1);

  FILE *f2;
  if (!(f2 = fopen(out_f2, "w")))
      die(err_write_aig, "can not write output to '%s'", out_f2);

  init_aig_substitution();
  init_time = process_time();

  bool res;

  if (identify_final_stage_adder() && build_adder_miter()) {
    res = 1;

    if (!miter_to_file(f1))
      die(err_write_file, "failed to write miter to '%s'", out_f1);
    else
      msg("writing miter to %s", out_f1);

    generate_rewritten_aig();

    if (!aiger_write_to_file(rewritten, aiger_binary_mode, f2))
      die(err_write_aig, "failed to write rewritten aig to '%s'", out_f2);
    else
      msg("writing rewritten aig to '%s'", out_f2);


  } else {
    res = 0;

    if (!trivial_miter_to_file(f1))
      die(err_write_file, "failed to write trivial miter to '%s'", out_f1);
    else
      msg("writing trivial miter to %s", out_f1);

    if (!write_model(f2))
      die(err_write_aig, "failed to write original aig to '%s'", out_f2);
    else
      msg("writing original aig to '%s'", out_f2);
  }

  substitution_time = process_time();

  fclose(f1);
  fclose(f2);
  reset_aig_substitution();

  return res;
}
