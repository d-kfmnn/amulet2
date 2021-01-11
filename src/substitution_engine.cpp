/*------------------------------------------------------------------------*/
/*! \file substitution_engine.cpp
    \brief contains the substitution engine

  Part of AMulet2.0 : AIG Multiplier Verification Tool.
  Copyright (C) 2020 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include "substitution_engine.h"
/*------------------------------------------------------------------------*/
void init_gate_substitution(){
  allocate_gates();
  mark_aig_outputs();
  set_parents_and_children(0);
  set_xor();
}

/*------------------------------------------------------------------------*/

void substitution(const char * out_f1, const char * out_f2){
  assert(inp_f);
  assert(out_f1);
  assert(out_f2);

  FILE * f1;
  if (!(f1 = fopen (out_f1, "w")))
      die ("can not write output to '%s'", out_f1);

  FILE *f2;
  if (!(f2 = fopen (out_f2, "w")))
      die ("can not write output to '%s'", out_f2);

  init_aig_substitution();
  init_time = process_time();

  if(identify_final_stage_adder() && build_adder_miter()){
    if(!miter_to_file(f1))
       die ("failed to write miter to '%s'", out_f1);
    else msg("writing miter to %s", out_f1);

    generate_rewritten_aig();

    if (!aiger_write_to_file (rewritten, aiger_ascii_mode, f2))
      die ("failed to write rewritten aig to '%s'", out_f2);
    else msg ("writing rewritten aig to '%s'", out_f2);


  } else {
    if(!trivial_miter_to_file(f1))
       die ("failed to write trivial miter to '%s'", out_f1);
    else msg("writing trivial miter to %s", out_f1);

    if (!write_model(f2))
      die ("failed to write original aig to '%s'", out_f2);
    else msg ("writing original aig to '%s'",out_f2);
  }

  substitution_time = process_time();

  fclose(f1);
  fclose(f2);
  reset_aig_substitution();
}
