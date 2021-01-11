/*------------------------------------------------------------------------*/
/*! \file elimination.cpp
    \brief contains functions used in the polynomial solver

  This file contains all functions used for preprocessing the Gröbner basis
  and for reducing the specification by the slices.

  Part of AMulet2.0 : AIG Multiplier Verification Tool.
  Copyright (C) 2020 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include "elimination.h"
/*------------------------------------------------------------------------*/
// Global variables
/// defines whether pac printing is enabled
bool pac = 0;

/// defines whether printing nullstellensatz proofs is enabled
bool nss = 0;
/*------------------------------------------------------------------------*/
// Local variables

/// used to collect the factors of each slice for PAC proofs
std::vector<Polynomial*> factors_per_slice;

/// used to collect the specification of each slice for PAC proofs
std::vector<Polynomial*> spec_of_slice;

/*-------------------------------------------------------------------------*/

Polynomial * add_up_factors(FILE * file){
  while(factors_per_slice.size()>1){
    Polynomial * p = factors_per_slice.back();
    factors_per_slice.pop_back();
    Polynomial * q = factors_per_slice.back();
    factors_per_slice.pop_back();

    Polynomial * add = add_poly(p,q);
    print_pac_add_rule(file, p, q, add);
    if(xor_chain) add = mod_poly(add, 1, file);

    print_pac_del_rule(file, p);
    print_pac_del_rule(file, q);

    delete(p);
    delete(q);
    factors_per_slice.push_back(add);
  }
  Polynomial * res = factors_per_slice.back();
  factors_per_slice.pop_back();
  return res;
}
/*------------------------------------------------------------------------*/

std::vector<Polynomial*> merge_factors
  (FILE * file, std::vector<Polynomial*> container){

  unsigned i = container.size();
  if(i == 1) return container;
  Polynomial * p = container[i-1];
  Polynomial * q = container[i-2];

  if(p->get_level() == q->get_level()){
    Polynomial * add = add_poly(p,q);
    print_pac_add_rule(file, p, q, add);
    if(xor_chain) add = mod_poly(add, 1, file);

    print_pac_del_rule(file, p);
    print_pac_del_rule(file, q);
    add->set_level(p->get_level()+1);
    delete(p);
    delete(q);
    container.pop_back();
    container.pop_back();
    container.push_back(add);
    container = merge_factors(file, container);
  }
  return container;
}

/*------------------------------------------------------------------------*/

void add_up_spec_of_slice(FILE * file){

  while(spec_of_slice.size()>1){
    Polynomial * p = spec_of_slice.back();
    spec_of_slice.pop_back();
    Polynomial * q = spec_of_slice.back();
    spec_of_slice.pop_back();

    Polynomial * add = add_poly(p,q);
    print_pac_add_rule(file, p, q, add);
    if(xor_chain) add = mod_poly(add, 1, file);

    print_pac_del_rule(file, p);
    print_pac_del_rule(file, q);

    delete(p);
    delete(q);
    spec_of_slice.push_back(add);
  }
  Polynomial * res = spec_of_slice.back();
  spec_of_slice.pop_back();
  if (xor_chain) correct_pp(res, file);
  delete(res);
}

/*------------------------------------------------------------------------*/

void eliminate_by_one_gate(Gate * n1, Gate *n2, FILE * file){
  Polynomial * p1 = n1->get_gate_constraint();
  Polynomial * p2 = n2->get_gate_constraint();

  const Polynomial * negfactor = divide_by_term(p1, p2->get_lt());
  if(negfactor->is_constant_zero_poly()) return;

  Polynomial * mult   = multiply_poly(negfactor, p2);
  Polynomial * rem    = add_poly(p1, mult);

  if(pac){
    assert(file);
    if(!negfactor->is_constant_one_poly())
      print_pac_mul_rule(file, p2, negfactor, mult);
    else mult->set_idx(p2->get_idx());
    print_pac_add_rule(file, p1, mult, rem);
    print_pac_del_rule(file, p1);
    print_pac_del_rule(file, mult);

  } else if(nss) {
    if(rem) add_ancestors(n1, n2, negfactor);
  }
  delete(mult);
  delete(negfactor);
  delete(p1);
  n1->set_gate_constraint(rem);
}

/*------------------------------------------------------------------------*/

Polynomial * reduce_by_one_poly(const Polynomial * p1, Gate * n, FILE * file){
  Polynomial * p2 = n->get_gate_constraint();

  const Polynomial * negfactor = divide_by_term(p1, p2->get_lt());

  if(negfactor->is_constant_zero_poly()) return p1->copy();

  Polynomial * mult   = multiply_poly(negfactor, p2);
  Polynomial * rem    = add_poly(p1, mult);

  if(pac){
    assert(file);
    if(!negfactor->is_constant_one_poly())
      print_pac_mul_rule(file, p2, negfactor, mult);
    else mult->set_idx(p2->get_idx());

    if (mult) {
      factors_per_slice.push_back(mult);
      factors_per_slice = merge_factors(file, factors_per_slice);

      if(!negfactor->is_constant_one_poly()) {
        print_pac_del_rule(file, p2);
      }
    }
    else delete(mult);
  } else delete(mult);

  if(nss) add_fac(n, negfactor);

  delete(negfactor);
  return rem;
}

/*------------------------------------------------------------------------*/

void remove_internal_xor_gates(FILE * file){
  msg("remove internal xor gates");
  int counter=0;
  for(unsigned i = NN; i < M-1; i++){
    Gate * n = gates[i];
    if(n->get_xor_gate() != 1) continue;
    if(n->get_elim()) continue;
    assert(n->children_size() == 2);

    Gate * l_gate = n->children_front();
    Gate * r_gate = n->children_back();
    if(l_gate->get_xor_gate() != 2) continue;
    if(r_gate->get_xor_gate() != 2) continue;
    assert(l_gate->children_size() == 2);
    assert(r_gate->children_size() == 2);
    Gate * ll_gate = l_gate->children_front();
    Gate * lr_gate = l_gate->children_back();

    //set xor children to ll and lr by overwriting l and r
    n->set_children_front(ll_gate);
    n->set_children_back(lr_gate);

    lr_gate->parents_push_back(n);
    ll_gate->parents_push_back(n);

    eliminate_by_one_gate(n, l_gate, file);

    if(l_gate->parents_size() == 1) {
      l_gate->mark_elim();
      if(pac) {
        assert(file);
        print_pac_del_rule(file, l_gate->get_gate_constraint());
      }
      if(verbose >= 3) msg("removed %s", l_gate->get_var_name());
      counter++;
      ll_gate->parents_remove(l_gate);
      lr_gate->parents_remove(l_gate);
    }

    eliminate_by_one_gate(n, r_gate, file);

    if(r_gate->parents_size() == 1) {
      r_gate->mark_elim();
      if(pac) {
        assert(file);
        print_pac_del_rule(file, r_gate->get_gate_constraint());
      }
      if(verbose >= 3) msg("removed %s", r_gate->get_var_name());
      counter++;
      ll_gate->parents_remove(r_gate);
      lr_gate->parents_remove(r_gate);
    }
  }
  if (verbose >= 1) msg("removed %i internal xor gates", counter);
}

/*------------------------------------------------------------------------*/

void remove_single_occs_gates(FILE * file){
  msg("remove single occurence gates");
  int counter = 0;
  for(unsigned i = NN; i<M-1; i++){

    Gate * n = gates[i];

    if(n->get_elim()) continue;
    if(n->parents_size() > 1) continue;

    Gate * parent = n->parents_front();
    if (parent->get_output()) continue;
    if (!n->get_xor_gate() && parent->get_xor_gate() == 1) continue;
    if (n->get_xor_chain()) continue;

    parent->children_remove(n);

    eliminate_by_one_gate(parent, n, file);

    for (std::list<Gate*>::const_iterator it = n->children_begin();
      it != n->children_end(); ++it){

      Gate * n_child = *it;
      if(!parent->is_child(n_child)) parent->children_push_back(n_child);

      n_child->parents_remove(n);
      if(!n_child->is_in_parents(parent)) n_child->parents_push_back(parent);
    }

    n->mark_elim();

    if(pac) {
      assert(file);
      print_pac_del_rule(file, n->get_gate_constraint());
    }
    if(verbose >= 3) msg("removed %s", n->get_var_name());
    counter++;
  }
  if (verbose >= 1) msg("removed %i single occurence gates", counter);
}

/*------------------------------------------------------------------------*/

void remove_slice_minus_one_gates(FILE * file){
  msg("remove gates that are not assigned to slices");
  int counter = 0;
  for(unsigned i = NN; i < M-1; i++){
    Gate * n = gates[i];
    if(n->get_elim()) continue;
    if(n->get_slice() > -1) continue;
    assert (!n->get_input());


    for (std::list<Gate*>::const_iterator it_c = n->children_begin();
      it_c != n->children_end(); ++it_c){

      Gate * n_child = *it_c;
      n_child->parents_remove(n);
    }

    for (std::list<Gate*>::const_iterator it_p = n->parents_begin();
      it_p != n->parents_end(); ++it_p){

      Gate * n_parent = *it_p;
      n_parent->children_remove(n);
      for (std::list<Gate*>::const_iterator it_c = n->children_begin();
        it_c != n->children_end(); ++it_c){

        Gate * n_child = *it_c;

        if(!n_parent->is_child(n_child))
          n_parent->children_push_back(n_child);
        if(!n_child->is_in_parents(n_parent))
          n_child->parents_push_back(n_parent);
      }
      eliminate_by_one_gate(n_parent, n, file);
    }

    n->mark_elim();

    if(pac) {
      assert(file);
      print_pac_del_rule(file, n->get_gate_constraint());
    }
    counter++;
    if(verbose >= 3) {
      msg("removed %s",n->get_var_name());
    }
  }
  if (verbose >= 1) msg("removed %i gates that are not assigned to slices", counter);
}

/*------------------------------------------------------------------------*/

void decomposing(FILE * file){
  msg("eliminate single occs");
  int counter = 0;
  int begin = xor_chain ? NN-2 : NN-1;
  for (int i = begin; i >= 0; i--) {
    bool change = 1;
    while(change){
      change = 0;
      for (std::list<Gate*>::const_iterator it =slices[i].begin();
        it != slices[i].end(); ++it){

        Gate * n = *it;
        if(n->parents_size() == 1 && !n->get_carry_gate()){
          Gate * parent = n->parents_front();
          parent->children_remove(n);

          eliminate_by_one_gate(parent, n, file);

          for (std::list<Gate*>::const_iterator it=n->children_begin();
            it != n->children_end(); ++it){

            Gate * n_child = *it;
            if(!parent->is_child(n_child))
              parent->children_push_back(n_child);

            n_child->parents_remove(n);

            if(!n_child->is_in_parents(parent))
              n_child->parents_push_back(parent);
          }
          n->mark_elim();
          if(pac) {
            assert(file);
            print_pac_del_rule(file, n->get_gate_constraint());
          }
          counter++;

          ++it;
          slices[i].remove(n);
          change = 1;

          if(verbose>=3) {
            msg("decomposed %s",n->get_var_name());
          }
        }
      }
    }
  }
  msg("decomposed %i variables", counter);
}

/*------------------------------------------------------------------------*/

void eliminate_booth_pattern(FILE * file){
  msg("eliminate booth pattern");
  int counter = 0;
  for(unsigned i = NN; i < M-1; i++){
    Gate * n = gates[i];
    if(!n->get_bo()) continue;
    if(n->get_elim()) continue;

    for (std::list<Gate*>::const_iterator it_c=n->children_begin();
      it_c != n->children_end(); ++it_c){

      Gate * n_child = *it_c;
      n_child->parents_remove(n);
    }

    for (std::list<Gate*>::const_iterator it_p=n->parents_begin();
      it_p != n->parents_end(); ++it_p){

      Gate * n_parent = *it_p;
      n_parent->children_remove(n);
      for (std::list<Gate*>::const_iterator it_c=n->children_begin();
        it_c != n->children_end(); ++it_c){

        Gate * n_child = *it_c;
        n_parent->children_push_back(n_child);
        n_child->parents_push_back(n_parent);
      }
      eliminate_by_one_gate(n_parent, n, file);
    }
    n->mark_elim();
    if(pac) {
      assert(file);
      print_pac_del_rule(file, n->get_gate_constraint());
    }
    counter++;
    if(verbose>=3) {
      msg("eliminated %s",n->get_var_name());
    }
  }
  msg("eliminated %i variables from booth pattern", counter);
}

/*------------------------------------------------------------------------*/

Polynomial * inc_spec_poly(unsigned i){
  mpz_t coeff;
  mpz_init(coeff);

  mpz_pow_ui(coeff, base, i);
  if(i == NN-1 && signed_mult) mpz_neg(coeff, coeff);

  const Var * s = gates[i+M-1]->get_var();
  Term * t1 = new_term(s, 0);
  Monomial * m1 = new Monomial(coeff, t1);
  push_mstack_end(m1);

  mpz_neg(coeff,coeff);

  int min_start = std::min(NN/2-1, i);
  for(int j = min_start; j >= 0; j--){
    if(mpz_sgn(coeff) == 1) mpz_neg(coeff,coeff);
    const Var * b = gates[b0+j*binc]->get_var();
    unsigned k = i - j;
    if (k > NN/2-1) break;
    if(k == NN/2-1 && signed_mult) mpz_neg(coeff,coeff);
    if(j == (int)NN/2-1 && signed_mult) mpz_neg(coeff,coeff);

    const Var * a = gates[a0+k*ainc]->get_var();
    add_to_vstack(b);
    add_to_vstack(a);
    Term * t1 = build_term_from_stack();
    Monomial * m1 = new Monomial(coeff, t1);
    push_mstack_end(m1);

  }
  mpz_clear(coeff);

  Polynomial * res = build_poly();
  return res;
}

/*------------------------------------------------------------------------*/

Polynomial * mod_poly(const Polynomial *p1, bool print_rule, FILE * file){

  bool pac_print = print_rule & pac;
  int exp = NN;

  mpz_t coeff;
  mpz_init(coeff);

  for(std::deque<Monomial*>::const_iterator it=p1->mon_begin();
    it!=p1->mon_end(); ++it){

    Monomial * m = *it;
    mpz_tdiv_r_2exp(coeff, m->coeff, exp);
    if(mpz_sgn(coeff) != 0){
      Monomial * tmp;
      if(m->get_term()) tmp =  new Monomial(coeff, m->get_term_copy());
      else tmp =  new Monomial(coeff, 0);
      push_mstack_end(tmp);
    }
  }
  mpz_clear(coeff);
  Polynomial * out = build_poly();
  out->set_idx(p1->get_idx());

  if(pac_print || nss){
    mpz_t quot;
    mpz_init(quot);
    for(std::deque<Monomial*>::const_iterator it=p1->mon_begin();
      it!=p1->mon_end(); ++it){

      Monomial * m = *it;
      if(pac_print || nss) {
        mpz_tdiv_q_2exp(quot, m->coeff, exp);
        if(mpz_sgn(quot) != 0){
          mpz_neg(quot, quot);
          Monomial * tmp;
          if(m->get_term()) tmp =  new Monomial(quot, m->get_term_copy());
          else tmp =  new Monomial(quot, 0);
          push_mstack_end(tmp);
        }
      }
    }
    mpz_clear(quot);
  }

  if(pac_print && !mstack_is_empty()){
    assert(file);
    Polynomial * p = build_poly();
    Polynomial * mod = multiply_poly_with_constant(p, mod_coeff);
    print_pac_mod_rule(file,  p, mod);
    print_pac_add_rule(file,  p1, mod, out);
    delete(p);
    delete(mod);
  }
  else if(nss && !mstack_is_empty()){
    Polynomial * p = build_poly();
    add_fac_mod(p);
    delete(p);
  }

  delete(p1);
  return out;
}

/*------------------------------------------------------------------------*/

void correct_pp_unsigned(const Polynomial * p, FILE * file){

  for (std::deque<Monomial *>::const_iterator it=p->mon_begin();
    it != p->mon_end(); ++it){

    Monomial * m = *it;

    if(mpz_sgn(m->coeff) == -1){
      Term * t = m->get_term();
      if(!t) continue;
      if(t->get_var_num() <= 0) continue;
      if(!gate(t->get_var_num())->get_input()) continue;
      push_mstack_end(new Monomial(one, t->copy()));
    }
  }

  if(!mstack_is_empty()){

    Polynomial * factor = build_poly();
    Polynomial * mod = multiply_poly_with_constant(factor, mod_coeff);
    print_pac_mod_rule(file,  factor, mod);
    Polynomial * add = add_poly(mod, p);
    print_pac_add_rule(file,  mod, p, add);
    delete(factor);
    delete(mod);
    delete(add);
  }
}

/*------------------------------------------------------------------------*/

void correct_pp_signed(const Polynomial * p, FILE * file){
  mpz_t half_mod;
  mpz_init(half_mod);
  mpz_pow_ui(half_mod, base, NN-1);

  mpz_t half_mod_neg;
  mpz_init(half_mod_neg);
  mpz_neg(half_mod_neg, half_mod);

  for (std::deque<Monomial *>::const_iterator it=p->mon_begin();
    it != p->mon_end(); ++it){

    Monomial * m = *it;
    Term * t = m->get_term();
    if(!t) continue;
    if(t->get_var_num() <= 0) continue;
    if(!gate(t->get_var_num())->get_input()) continue;

    if(mpz_cmp(m->coeff, half_mod)>0){
      Monomial * tmp = new Monomial(minus_one, m->get_term_copy());
      push_mstack_end(tmp);
    }
    else if(mpz_cmp(m->coeff, half_mod_neg)<0){
      Monomial * tmp = new Monomial(one,  m->get_term_copy());
      push_mstack_end(tmp);
    }
  }

  mpz_clear(half_mod);
  mpz_clear(half_mod_neg);

  Polynomial * factor = build_poly();
  if(!factor) return;

  Polynomial * mod = multiply_poly_with_constant(factor, mod_coeff);
  print_pac_mod_rule(file,  factor, mod);

  Polynomial * add = add_poly(mod, p);
  print_pac_add_rule(file,  mod, p, add);

  delete(mod);
  delete(factor);
  delete(add);
}

/*------------------------------------------------------------------------*/

void correct_pp(const Polynomial * p, FILE * file){
  if(signed_mult) correct_pp_signed(p, file);
  else correct_pp_unsigned(p, file);
}

/*------------------------------------------------------------------------*/

const Polynomial * reduce(FILE * file){
  msg("");
  msg("");
  msg("started reducing");
  Polynomial * rem = 0, * tmp;
  for (int i=NN-1; i>= 0; i--){
    if(verbose >= 1) msg("reducing by slice %i", i);
    Polynomial * inc_spec = inc_spec_poly(i);


    if(rem) {
      tmp = add_poly(inc_spec, rem);
      delete(inc_spec);
      delete(rem);
      rem = tmp;
    } else rem = inc_spec;


    std::list<Gate*> sl = slices[i];
    for (std::list<Gate*>::const_iterator it=sl.begin(); it != sl.end(); ++it){
      Gate * n = *it;
      if(n->get_elim()) continue;
      if(verbose >= 4){
        fputs("[amulet2] reducing by ", stdout);
        n->print_gate_constraint(stdout);

      }
      tmp = reduce_by_one_poly(rem, n, file);


      if(xor_chain) tmp = mod_poly(tmp, 0, file);

      delete(rem);
      rem = tmp;

      if(verbose >= 3){
        fputs("[amulet2] remainder is ", stdout);
        rem->print(stdout);
        msg(" ");
      }

    }

    if(verbose >= 2) {
      msg("after reducing by slice %i", i);
      fprintf(stdout, "[amulet2] remainder is ");
      rem->print(stdout);
      msg("");
    }

    if(pac) {
      Polynomial * pac_poly = add_up_factors(file);
      pac_poly->set_level(1);
      spec_of_slice.push_back(pac_poly);
      spec_of_slice = merge_factors(file, spec_of_slice);
    }
  }
  if(pac)  add_up_spec_of_slice(file);

  return rem;
}

/*------------------------------------------------------------------------*/

bool check_inputs_only(const Polynomial *p){
  for(std::deque<Monomial*>::const_iterator it=p->mon_begin();
    it!=p->mon_end(); ++it){

    Monomial * m = *it;
    if(!m->get_term()) continue;
    else{
      const Var * v = m->get_term()->get_var();
      Gate * n = gate(v->get_num());
      if (!n->get_input()) return 0;
    }
  }
  return 1;
}

/*------------------------------------------------------------------------*/

void write_witness_vector(const Term * t, FILE * file){


  fputs_unlocked ("[amulet2]   ", stdout);



  for(unsigned i = 0; i <= NN/2-1; i++){
    const Var * v = gates[a0+i*ainc]->get_var();

    if(t->contains(v)) {
      fprintf(file, "1");
      fprintf(stdout, "%s = ",v->get_name());
    }
    else fprintf(file, "0");

    const Var * w = gates[b0+i*binc]->get_var();

    if(t->contains(w)) {
      fprintf(file, "1");
      fprintf(stdout, "%s = ",w->get_name());
    }
    else fprintf(file, "0");
  }
  fprintf(stdout, "1, all other inputs = 0;\n");
  fprintf(file, "\n");
}

/*------------------------------------------------------------------------*/

void write_witnesses(const Polynomial * p, FILE * file){
  assert(check_inputs_only(p));

  int len = p->min_term_size();
  if(len == 0){
    msg ("  all inputs = 0;\n");
    for(unsigned i = 0; i <= NN/2-1; i++){
      fprintf(file, "00");
    }
    fprintf(file, "\n");

  } else {
    for(std::deque<Monomial*>::const_iterator it = p->mon_begin();
      it != p->mon_end(); ++it){

      Monomial * m = *it;
      if(m->get_term()) {
        int tlen = m->get_term_size();
        if(tlen == len) write_witness_vector(m->get_term(), file);
      }
    }
  }
  fprintf(file, ".");
}

/*------------------------------------------------------------------------*/

void generate_witness(const Polynomial * p, const char * name){
  if(!check_inputs_only(p))
  die ("cannot generate witness, as remainder polynomial contains non-inputs");

  char witness_name[100];
  memset(witness_name, '\0', sizeof(witness_name));
  for (int i = 0; name[i] != '.'; i++) {
    witness_name[i] = name[i];
  }
  strcat(witness_name, ".cex");


  FILE * witness_file;
  if (!(witness_file = fopen (witness_name, "w")))
  die ("cannot write output to '%s'", witness_name);

  msg ("");
  msg ("COUNTER EXAMPLES ARE: ");

  write_witnesses(p, witness_file);

  msg ("");
  msg ("");
  msg ("Counter examples are written to %s", witness_name);
  msg ("Call ");
  msg ("       'aiger/aigsim %s %s' ", name, witness_name);
  msg ("to simulate the provided counter example(s).");
  msg ("");
  msg ("Note: aiger/aigsim produces output in the form:");
  fputs_unlocked ("[amulet2] ", stdout);
  if (NN == 2)      fprintf(stdout, "  a[0]b[0]  s[0]\n");
  else if (NN == 4) fprintf(stdout, "  a[0]b[0]a[1]b[1]  s[0]s[1]s[2]s[3]\n");
  else
  fprintf(stdout, "  a[0]b[0]a[1]b[1]...a[%u]b[%u]  s[0]s[1]s[2]...s[%u]\n",
    NN/2-1, NN/2-1, NN-1);

  fclose(witness_file);
}
