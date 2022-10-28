/*------------------------------------------------------------------------*/
/*! \file substitution.cpp
    \brief contains function to apply adder substitution

  Part of AMulet2 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include <list>

#include "substitution.h"
/*------------------------------------------------------------------------*/
// ERROR CODES:
static int err_miter = 71; // error in miter
/*------------------------------------------------------------------------*/
// Local variables

unsigned aig_idx;

bool no_cin;
bool single_gen_gate;

static Gate * carry_out;
static Gate * carry_in;
static std::vector<Gate*> outputs;
static std::vector<unsigned> original_outputs;
static std::vector<unsigned> rewritten_outputs;
static std::vector<Gate*> inputs;
static std::list<unsigned>plain_inputs;
static std::vector<Gate*> c_ins;
/*------------------------------------------------------------------------*/

bool all_single_output() {
  for (unsigned i = 0; i < NN-1; i++) {
    Gate * n = gate(slit(i));
    if (n->parents_size() > 1) return 0;
  }
  return 1;
}

/*------------------------------------------------------------------------*/

bool all_outputs_are_xor() {
  unsigned lit = slit(NN-1);
  if(lit < 2) return 0;

  for (unsigned i = 1; i < NN-1; i++) {
    unsigned lit = slit(i);
    if(lit < 2) return 0;
    Gate * n = gate(slit(i));
    if (!n->get_xor_gate()) return 0;
  }
  return 1;
}

/*------------------------------------------------------------------------*/

bool slice_two_needs_carry_in_slice_zero() {
  unsigned lit2 = slit(2);
  if(lit2 < 2) return 0;
  Gate * out2 = gate(slit(2));
  unsigned lit = slit(0);
  if(lit < 2) return 0;
  Gate * out0 = gate(slit(0));
  if (out2->parents_size() > 3 && out0->parents_size() == 1) return 0;
  return 1;
}

/*------------------------------------------------------------------------*/

bool cin_in_slice_0() {
  unsigned lit = slit(0);
  if(lit < 2) return 0;
  Gate * n = gate(slit(0));
  if (n->parents_size() > 1) return 1;
  return 0;
}

/*------------------------------------------------------------------------*/

unsigned get_input(bool flip) {
  Gate * n = inputs.back();
  inputs.pop_back();
  if (!flip && n->get_neg()) return not_(n->get_var_num());
  else if (flip && !n->get_neg()) return not_(n->get_var_num());
  return n->get_var_num();
}

/*------------------------------------------------------------------------*/

void push_to_inputs(Gate * n) {
  inputs.push_back(n);
  n->inc_fsa_inp();
  n->mark_fsa();
}

/*------------------------------------------------------------------------*/

void push_to_outputs(Gate * n, int i) {
  outputs.push_back(n);
  n->set_slice(i);
}

/*------------------------------------------------------------------------*/

void push_to_cins(Gate * n, int i) {
  c_ins.push_back(n);
  carry_in = n;
  n->mark_fsa();
  if (verbose >= 2)  msg("found cin of slice %i %s", i, n->get_var_name());
}

/*------------------------------------------------------------------------*/

void set_carry_in(Gate * n) {
  carry_in = n;
  n->mark_fsa();
  if (verbose >= 3)  msg("identified carry in %s ", n->get_var_name());
}

/*------------------------------------------------------------------------*/

void identify_carry_out() {
  Gate * largest_aig_output = gate(slit(NN-1));

  if (largest_aig_output->get_xor_gate() != 1) {
    carry_out = largest_aig_output;
    push_to_outputs(carry_out, NN-1);
  } else {
    Gate * l = xor_left_child(largest_aig_output);
    Gate * r = xor_right_child(largest_aig_output);
    if (r->get_level() > l->get_level()) carry_out = r;
    else
      carry_out = l;
    push_to_outputs(carry_out, -1);
  }
  if (verbose >= 3)  msg("identified carry out %s", carry_out->get_var_name());
}

/*------------------------------------------------------------------------*/

bool identify_propagate_and_generate_gates() {
// slice NN-1 contains carry, slice 0 is no XOR
  for (int i = NN-2; i > 0; i--) {
    Gate * n = gate(slit(i));

    if (i == 2 && n->parents_size() > 3) {
      assert(gate(slit(0))->parents_size() > 1);

      push_to_outputs(n, 2);
      push_to_outputs(gate(slit(1)), 1);
      push_to_outputs(gate(slit(0)), 0);

      push_to_inputs(n);
      push_to_inputs(gate(slit(1)));
      set_carry_in(gate(slit(0)));
      return 1;
    }

    Gate * internal_xor, *l = 0, *r = 0;
    if (i == 1 && n->parents_size() > 1) { internal_xor = n;
    } else {
      l = xor_left_child(n);
      r = xor_right_child(n);
      internal_xor = l->get_xor_gate() ? l : r;
    }

    int cmp = NN-1;
    if (internal_xor->parents_size() < 3) break;

    if (internal_xor->parents_size() == 3 && i < 3*cmp/4
        && !cin_in_slice_0()) {
      if (all_single_output()) break;
      else if (!booth && !signed_mult) break;  // sp-ba-csv
    }


    internal_xor->mark_prop_gen_gate();
    if (verbose >= 2)
      msg("found propagate gate %s", internal_xor->get_var_name());
    Gate *g_0 = 0, *g_1 = 0;

    if (internal_xor->get_xor_gate() &&
        xor_left_child(internal_xor)->parents_size() != 2 &&
        xor_right_child(internal_xor)->parents_size() != 2 &&
        (i != 1 || !signed_mult || n->parents_size() == 1 || booth)) {
      Gate * internal_and = derive_ha_and_gate(internal_xor);
      internal_and->mark_prop_gen_gate();
      if (verbose >= 2)
        msg("found generate gate %s", internal_and->get_var_name());

      aiger_and * par = is_model_and(internal_and->get_var_num());
      g_0 = gate(par->rhs0);
      g_1 = gate(par->rhs1);
      g_0->set_neg(aiger_sign(par->rhs0));
      g_1->set_neg(aiger_sign(par->rhs1));
      push_to_inputs(g_0);
      push_to_inputs(g_1);
    } else if (booth) {
      push_to_inputs(internal_xor);

      if (verbose >= 3)  msg("pushed xor %s", internal_xor->get_var_name());
      single_gen_gate = 1;
    }

    push_to_outputs(n, i);
    if (i != 1 || n->parents_size() == 1) {
      if (l->get_xor_gate()) push_to_cins(r, i);
      else
        push_to_cins(l, i);
    } else  {
      Gate * c = gate(slit(0));
      if (c->parents_size() > 1) {
         push_to_cins(c, i);
         push_to_outputs(c, 0);
      } else if (booth &&(g_0->get_xor_gate() || g_1->get_xor_gate())) {
        // bp-ba-csv
        Gate * not_xor_cin = g_0->get_xor_gate() ? g_1 : g_0;
        push_to_cins(not_xor_cin, i);
        no_cin = 1;
      }
    }
  }
  return 1;
}

/*------------------------------------------------------------------------*/

void fix_inputs() {
  if (!cin_in_slice_0() && !signed_mult) return;
  if (!cin_in_slice_0() && all_single_output()) return;

  static std::vector<Gate*> inputs_cpy;
  for (std::vector<Gate*>::const_iterator it = inputs.begin();
      it != inputs.end(); ++it) {
    Gate * n = *it;
    if (!n->get_prop_gen_gate()) { inputs_cpy.push_back(n);
    } else {
      aiger_and * and1 = is_model_and(n->get_var_num());
      if (aiger_sign(and1->rhs0) != aiger_sign(and1->rhs1)) {
        if (aiger_sign(and1->rhs0)) inputs_cpy.push_back(gate(and1->rhs0));
        if (aiger_sign(and1->rhs1)) inputs_cpy.push_back(gate(and1->rhs1));
      } else if (signed_mult && !booth && aiger_sign(and1->rhs0) &&
          !n->get_aig_output()) {
        gate(and1->rhs0)->inc_fsa_inp();
        inputs_cpy.push_back(gate(and1->rhs0));
      } else if (signed_mult && !booth) {
        n->unmark_prop_gen_gate();
      }
    }
  }
  inputs = inputs_cpy;
}

/*------------------------------------------------------------------------*/

bool follow_path_and_mark_gates(Gate * n, bool init) {
  if (n->get_input() && !n->get_fsa_inp()) return 0;

  n->mark_fsa();
  if (n == carry_in) return 1;
  if (n->get_fsa_inp()) return 1;

  aiger_and * and1 = is_model_and(n->get_var_num());
  Gate * l = gate(and1->rhs0);
  Gate * r = gate(and1->rhs1);

  if (!r->get_prop_gen_gate() && carry_in == r && init  && !r->get_neg()) {
     r->set_neg(aiger_sign(and1->rhs1));
  }
  if (!follow_path_and_mark_gates(r, init)) return 0;

  if (!l->get_prop_gen_gate() && carry_in == l && init && !l->get_neg()) {
     l->set_neg(aiger_sign(and1->rhs0));
  }
  if (!follow_path_and_mark_gates(l, init)) return 0;


  return 1;
}

/*------------------------------------------------------------------------*/

bool follow_all_output_paths_and_mark_gates() {
  msg("checking last stage adder");
  for (std::vector<Gate*>::const_iterator it = outputs.begin();
      it != outputs.end(); ++it) {
    Gate * n = *it;
    if (verbose >= 3)  msg("follow path starting with %s", n->get_var_name());
    bool init =(it == outputs.begin()) ? 1 : 0;
    if (!follow_path_and_mark_gates(n, init)) return 0;
  }
  return 1;
}

/*------------------------------------------------------------------------*/

void correctly_mark_inputs() {
    for (unsigned i = 0; i< inputs.size(); i++) {
      if (inputs[i]->get_prop_gen_gate()) continue;
      if (!inputs[i]->get_aig_output()) inputs[i]->reset_fsa_inp();
    }

    for (unsigned i = M-1; i >= 1; i--) {
      Gate * n = gates[i];
      if (!n->get_prop_gen_gate()) continue;
      if (single_gen_gate && n->get_fsa_inp()) continue;

      n->reset_fsa_inp();

      aiger_and * and1 = is_model_and(n->get_var_num());
      unsigned l = aiger_strip(and1->rhs0), r = aiger_strip(and1->rhs1);

      if (!n->get_xor_gate()) {
        gate(l)->inc_fsa_inp();
        gate(r)->inc_fsa_inp();
      }
    }

    carry_in->inc_fsa_inp();

    if (single_gen_gate) {
      for (unsigned i = 0; i< inputs.size(); i++) {
        if (!inputs[i]->get_fsa_inp()) {
          inputs[i]->inc_fsa_inp();
        }
      }
    }
}

/*------------------------------------------------------------------------*/

bool identify_final_stage_adder() {
  if (!all_outputs_are_xor()) {
    msg("substitution not possible - not all outputs are XORs");
    return 0;
  }
  if (!slice_two_needs_carry_in_slice_zero()) {
    msg("substitution not possible - carry in slice 0 not found");
    return 0;
  }

  identify_carry_out();
  if (!identify_propagate_and_generate_gates()) {
    msg("substitution not possible - propagate and generate gates not found");
    return 0;
  }
  fix_inputs();

  if (!follow_all_output_paths_and_mark_gates()) {
    msg("substitution not possible - no clear boundaries");
    return 0;
  }

  correctly_mark_inputs();

  return 1;
}
/*----------------------------------------------------------------------------*/
void add_original_adder() {
  for (unsigned i = 0; i < M-1; i++) {
    Gate * n = gates[i];

    if (!n->get_fsa()) continue;

    if (n->get_fsa_inp() || n == carry_in) {
      aiger_add_input(miter, n->get_var_num(), n->get_var_name());
      if (verbose >= 3) msg("miter input %s", n->get_var_name());
    } else {
      aiger_and * and1 = is_model_and(n->get_var_num());
      aiger_add_and(miter, and1->lhs, and1->rhs0, and1->rhs1);
      if (verbose >= 4)
        msg("original adder and %i %i %i", and1->lhs, and1->rhs0, and1->rhs1);
    }
  }
}

/*----------------------------------------------------------------------------*/

void fill_original_outputs() {
  // carry out in pos 0 needs special treatment
  for (int i = outputs.size()-1; i >= 1; i--) {
    Gate * n = outputs[i];
    unsigned res = slit(n->get_slice());
    original_outputs.push_back(res);
    if (verbose >= 3) msg("%i is output ", res);
  }
  Gate * c_out = outputs[0];
  if (c_out->get_aig_output()) {
    unsigned res = slit(NN-1);
    original_outputs.push_back(res);
    if (verbose >= 3) msg("%i is output ", res);

  } else {
    if (aiger_sign(c_out->get_var_num()))
      original_outputs.push_back(c_out->get_var_num());
    else
      original_outputs.push_back(not_(c_out->get_var_num()));

    if (verbose >= 3) msg("%i is output ", original_outputs.back());
  }
}

/*----------------------------------------------------------------------------*/

unsigned btor_ha(unsigned i1, unsigned i2, bool carry) {
  unsigned one, two, three;
  one = aig_idx + 2;
  two = one + 2;
  three = two + 2;

  aiger_add_and(miter, one, not_(i1), not_(i2));
  aiger_add_and(miter, two, i1, i2);
  aiger_add_and(miter, three, not_(one), not_(two));

  aiger_add_and(rewritten, one, not_(i1), not_(i2));
  aiger_add_and(rewritten, two, i1, i2);
  aiger_add_and(rewritten, three, not_(one), not_(two));

  aig_idx = aig_idx + 6;
  if (carry) {
    rewritten_outputs.push_back(three);
    if (verbose >= 2)
      msg("ha with outputs %i, %i, inputs  %i, %i", two, three, i1, i2);
    return two;
  } else {
    if (verbose >= 2)
      msg("ha with sum output %i, inputs  %i, %i", three, i1, i2);
    return three;
  }
}

/*------------------------------------------------------------------------*/

unsigned btor_fa(unsigned i1, unsigned i2, unsigned i3, bool carry) {
  unsigned one, two, three, four, five, six;
  one = aig_idx + 2;
  two = one + 2;
  three = two + 2;
  four = three + 2;
  five = four + 2;
  six = five +2;

  aiger_add_and(miter, one, not_(i1), not_(i2));
  aiger_add_and(miter, two, i1, i2);
  aiger_add_and(miter, three, not_(one), not_(two));
  aiger_add_and(miter, four, not_(i3), not_(three));
  aiger_add_and(miter, five, i3, three);
  aiger_add_and(miter, six, not_(four), not_(five));


  aiger_add_and(rewritten, one, not_(i1), not_(i2));
  aiger_add_and(rewritten, two, i1, i2);
  aiger_add_and(rewritten, three, not_(one), not_(two));
  aiger_add_and(rewritten, four, not_(i3), not_(three));
  aiger_add_and(rewritten, five, i3, three);
  aiger_add_and(rewritten, six, not_(four), not_(five));

  aig_idx = aig_idx + 12;
  rewritten_outputs.push_back(six);

  if (carry) {
    unsigned seven = six + 2;
    aiger_add_and(miter, seven, not_(two), not_(five));
    aiger_add_and(rewritten, seven, not_(two), not_(five));
    aig_idx = aig_idx + 2;
    if (verbose >= 2)
      msg("fa with outputs %i, %i, inputs  %i, %i, %i", seven, six, i1, i2, i3);
    return seven;
  } else {
    if (verbose >= 2)
      msg("fa no carry with output %i, inputs %i, %i, %i", six, i1, i2, i3);
    return six;
  }
}

/*------------------------------------------------------------------------*/

void add_btor_adder() {
  aig_idx = 2*get_model_maxvar()+2;

  unsigned i2 = 0, i3 = 0;
  unsigned c = carry_in->get_var_num();
  if (carry_in->get_neg()) c = not_(c);

  if ((!signed_mult || booth) && cin_in_slice_0()) {
    if (gate(slit(2))->parents_size() > 1) {
      // sp-ba-bc
      inputs.pop_back();
      inputs.pop_back();

      rewritten_outputs.push_back(slit(0));
      rewritten_outputs.push_back(slit(1));
      rewritten_outputs.push_back(slit(2));
      if (verbose >= 2) msg("single output %i, inputs  %i", slit(0), slit(0));
      if (verbose >= 2) msg("single output %i, inputs  %i", slit(1), slit(1));
      if (verbose >= 2) msg("single output %i, inputs  %i", slit(2), slit(2));

      i2 = get_input();
      i3 = get_input();
      c = btor_fa(not_(slit(2)), i2, i3);
      c = not_(c);
    } else {
      rewritten_outputs.push_back(slit(0));
      if (verbose >= 2) msg("single output %i, inputs  %i", slit(0), slit(0));
      i2 = get_input();
      i3 = get_input();
      if (!booth) { c = btor_ha(i2, i3);
      } else {
        // bp-ba
        c = btor_ha(c, not_(c), 0);
        c = btor_fa(c, i2, i3);
        c = not_(c);
      }
    }
  } else if (signed_mult && !booth && gate(slit(1))->parents_size() > 1) {
    if (cin_in_slice_0()) {
      rewritten_outputs.push_back(slit(0));
      rewritten_outputs.push_back(slit(1));

      if (verbose >= 2) msg("single output %i, inputs  %i", slit(0), slit(0));
      if (verbose >= 2) msg("single output %i, inputs  %i", slit(1), slit(1));
    } else {
      rewritten_outputs.push_back(c);
      if (verbose >= 2) msg("single output %i, inputs  %i", c, c);
    }

    Gate * n = inputs.back();
    c = n->get_var_num();
    rewritten_outputs.push_back(c);
    if (verbose >= 2) msg("single output %i, inputs  %i", c, c);
    c = not_(c);
    inputs.pop_back();
    i2 = get_input();
    i3 = get_input();
    c = btor_fa(c, i2, i3);
    c = not_(c);

  } else if (single_gen_gate) {
    i2 = get_input(1);
    c = btor_ha(c, i2);
  }

  for (int i = inputs.size()-1; i >= 0; i--) {
    Gate * v = inputs[i--];
    Gate * w = inputs[i];
    i2 = v->get_neg() ? not_(v->get_var_num()) : v->get_var_num();
    i3 = w->get_neg() ? not_(w->get_var_num()) : w->get_var_num();

    if (!v->get_fsa_inp()) { c = btor_ha(c, i3);
    } else if (v == w) {
      c = btor_ha(c, v->get_var_num());
      i++;
    } else if (booth &&  v == carry_in && gate(slit(1))->parents_size() > 1) {
      // bp-ba-csv
      c = btor_ha(c, not_(c), 0);
      c = btor_fa(c, v->get_var_num(), i3);
      c = not_(c);
    } else if (signed_mult && i == 0 && v->get_fsa_inp() == 2) {
      if (w->get_neg()) i3 = not_(i3);
      c = not_(btor_fa(c, i2, i3));
      c = btor_fa(c, i2, i3, 0);
      return;
    } else if (signed_mult && !booth && w == inputs[i-1]) {
      if (w->get_neg()) i3 = not_(i3);
      c = not_(btor_fa(c, i2, i3));
    } else if (v->get_fsa_inp() > 1) { c = not_(btor_fa(c, not_(i2), i3));
    } else {
      c = not_(btor_fa(c, i2, i3));
    }
  }

  if (signed_mult && carry_out->get_aig_output()) c = not_(c);
  rewritten_outputs.push_back(c);
  if (verbose >= 2) msg("msb output %i", c);
}

/*----------------------------------------------------------------------------*/
unsigned not_(unsigned a) { return a^1; }
/*------------------------------------------------------------------------*/
unsigned and_(unsigned a, unsigned b) {
  unsigned res;
  if (!a || !b || a == not_(b)) return 0;
  if (a == 1 || a == b) return b;
  if (b == 1) return a;
  res = 2*(miter->maxvar + 1);
  assert(a < res), assert(b < res);
  aiger_add_and(miter, res, a, b);
  if (verbose >= 4) msg("miter and %i %i %i", res, a, b);
  return res;
}

/*------------------------------------------------------------------------*/

unsigned implies_(unsigned a, unsigned b) {
  return not_(and_(a, not_(b)));
}

/*------------------------------------------------------------------------*/

unsigned xnor_(unsigned a, unsigned b) {
  return and_(implies_(a, b), implies_(b, a));
}

/*------------------------------------------------------------------------*/

bool build_miter() {
  if (original_outputs.size() != rewritten_outputs.size()) {
    fprintf(stdout, "orig output contains: ");
    for (std::vector<unsigned>::const_iterator it = original_outputs.begin();
        it != original_outputs.end(); ++it) {
      fprintf(stdout, "%i ", *it);
    }

    fprintf(stdout, "\n");

    fprintf(stdout, "rewritten output contains: ");
    for (std::vector<unsigned>::const_iterator it = rewritten_outputs.begin();
        it != rewritten_outputs.end(); ++it) {
      fprintf(stdout, "%i ", *it);
    }
    fprintf(stdout, "\n");

    msg("missmatch in outputs -> abort rewriting");
    return 0;
  }

  unsigned out = 1;
  for (unsigned i = 0; i < outputs.size(); i++)
    out = and_(out, xnor_(original_outputs[i], rewritten_outputs[i]));

  aiger_add_output(miter, not_(out), "miter");
  return 1;
}

/*------------------------------------------------------------------------*/

bool build_adder_miter() {
  msg("build adder miter");

  add_original_adder();
  fill_original_outputs();
  add_btor_adder();
  if (!build_miter()) return 0;
  return 1;
}

/*----------------------------------------------------------------------------*/
bool miter_to_file(FILE * file) {
  if (!file) return 0;
  if (miter->num_outputs > 1) die(err_miter, "miter has more than one output");

  msg("transform aiger miter to cnf miter");
  int *map, m, n;
  unsigned *refs, i, lit;

  aiger_reencode(miter);

  refs = reinterpret_cast<unsigned*>(calloc(2*(miter->maxvar+1), sizeof *refs));
  assert(refs);
  lit = miter->outputs[0].lit;
  refs[lit]++;

  i = miter->num_ands;
  while (i--) {
    lit = miter->ands[i].lhs;
    if (refs[lit])  {
      refs[miter->ands[i].rhs0]++;
      refs[miter->ands[i].rhs1]++;
    }
    if (refs[aiger_not(lit)]) {
      refs[aiger_not(miter->ands[i].rhs0)]++;
      refs[aiger_not(miter->ands[i].rhs1)]++;
    }
  }

  map = reinterpret_cast<int*>(calloc(2*(miter->maxvar+1), sizeof *map));
  m = 0;
  n = 1;
  if (refs[0] || refs[1]) {
    map[0] = -1;
    map[1] = 1;
    m++;
    n++;
  }
  for (lit = 2; lit <= 2*miter->maxvar; lit += 2) {
    if (!refs[lit] && !refs[lit+1]) continue;
    map[lit] = ++m;
    map[lit+1] = -m;

    if (lit <= 2*miter->num_inputs+1) continue;
    if (refs[lit]) n += 2;
    if (refs[lit+1]) n += 1;
  }

  fprintf(file, "p cnf %i %i\n", m, n);
  if (verbose >= 2) msg("p cnf %i %i", m, n);

  if (refs[0] || refs[1]) fprintf(file, "%d 0\n", map[1]);

  for (i = 0; i < miter->num_ands; i++) {
    lit = miter->ands[i].lhs;
    if (refs[lit]) {
      fprintf(file, "%d %d 0\n",
      map[aiger_not(lit)], map[miter->ands[i].rhs1]);
      fprintf(file, "%d %d 0\n",
      map[aiger_not(lit)], map[miter->ands[i].rhs0]);
    }
    if (refs[lit+1]) {
      fprintf(file, "%d %d %d 0\n", map[lit],
      map[aiger_not(miter->ands[i].rhs1)],
      map[aiger_not(miter->ands[i].rhs0)]);
    }
  }

  fprintf(file, "%d 0\n", map[miter->outputs[0].lit]);
  free(refs);
  free(map);
  return 1;
}

/*------------------------------------------------------------------------*/

bool trivial_miter_to_file(FILE * file) {
  if (!file) return 0;
  fprintf(file, "p cnf 1 2 \n");
  fprintf(file, "1 0\n");
  fprintf(file, "-1 0\n");
  return 1;
}

/*------------------------------------------------------------------------*/

void generate_rewritten_aig() {
  msg("generate rewritten aig");
  // adding orig inputs
  for (unsigned i = 0; i < get_model_num_inputs(); i++) {
    aiger_add_input(
      rewritten, get_model_inputs_lit(i), get_model_inputs_name(i));
  }

  // internal gates, that have not been replaced
  unsigned out, btor_carry = rewritten_outputs.back();

  for (unsigned i = NN; i < M-1; i++) {
    Gate * n = gates[i];
    if (n->get_fsa() && !n->get_fsa_inp() && n != carry_in) continue;
    aiger_and * and1 = is_model_and(n->get_var_num());

    if (static_cast<int>(aiger_strip(and1->rhs0)) ==
      carry_out->get_var_num()) {
      if (!aiger_sign(and1->rhs0)) out = not_(btor_carry);
      else
        out = btor_carry;
      aiger_add_and(rewritten, and1->lhs, out, and1->rhs1);
      if (verbose >= 4)
        msg("rewritten and %i %i %i", and1->lhs, out, and1->rhs1);
    } else if (static_cast<int>(aiger_strip(and1->rhs1)) ==
        carry_out->get_var_num()) {
      if (!aiger_sign(and1->rhs1)) out = not_(btor_carry);
      else
        out = btor_carry;
      aiger_add_and(rewritten, and1->lhs, out, and1->rhs0);
      if (verbose >= 4)
        msg("rewritten and %i %i %i", and1->lhs, out, and1->rhs0);
    } else {
      aiger_add_and(rewritten, and1->lhs, and1->rhs0, and1->rhs1);
      if (verbose >= 4)
        msg("rewritten and %i %i %i",  and1->lhs, and1->rhs0, and1->rhs1);
    }
  }

  // setting outputs
  #define BUFSIZE 80
  char buf[BUFSIZE + 1];
  unsigned j = 0;

  // first the ones, which are not included in fsa
  for (unsigned i = 0; i < NN; i++) {
    unsigned res = slit(i);
    if (!gate(res)->get_fsa()) {
      snprintf(buf, BUFSIZE, "O%d", j++);
      aiger_add_output(rewritten, res, buf);
      if (verbose >= 4) msg("rewritten output %i %s ", res, buf);
    } else {
      break;
    }
  }

  // then those that are replaced
  for (unsigned i = 0; i < rewritten_outputs.size()-1; i++) {
    snprintf(buf, BUFSIZE, "O%d", j++);
    aiger_add_output(rewritten, rewritten_outputs[i], buf);
    if (verbose >= 4) msg("rewritten output %i %s ", rewritten_outputs[i], buf);
  }

  // need to check whether last output was replaced or not
  snprintf(buf, BUFSIZE, "O%u", j);
  if (carry_out->get_aig_output()) {
    aiger_add_output(rewritten, btor_carry, buf);
    if (verbose >= 4) msg("rewritten output %i %s ", btor_carry, buf);
  } else {
    unsigned res = slit(NN-1);
    aiger_add_output(rewritten, res, buf);
    if (verbose >= 4) msg("rewritten output %i %s ", res, buf);
  }
  aiger_reencode(rewritten);
}
