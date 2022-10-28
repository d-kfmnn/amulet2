/*------------------------------------------------------------------------*/
/*! \file elimination.cpp
    \brief contains functions used in the polynomial solver

  This file contains all functions used for preprocessing the Gröbner basis
  and for reducing the specification by the slices.

  Part of AMulet2 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include <algorithm>
#include <list>

#include "elimination.h"
/*------------------------------------------------------------------------*/
// Global variables
int proof = 0;

/*------------------------------------------------------------------------*/
// ERROR CODES:
static int err_writing = 61; // cannot write to
static int err_witness = 62; // cannot write to
/*------------------------------------------------------------------------*/
// Local variables

// / used to collect the factors of each slice for PAC proofs
static std::vector<Polynomial*> factors_per_slice;

// / used to collect the cofactors of each slice for PAC proofs
static std::vector<const Polynomial*> co_factors;

// / used to collect the cofactors of each slice for PAC proofs
static std::vector<int> factor_indices;

// / used to collect the cofactors of each slice spec for PAC proofs
static std::vector<int> spec_indices;

// / used to collect the specification of each slice for PAC proofs
static std::vector<Polynomial*> spec_of_slice;

/*-------------------------------------------------------------------------*/

Polynomial * add_up_factors(FILE * file, bool print) {
  Polynomial * p = factors_per_slice.back();
  factors_per_slice.pop_back();

  while (!factors_per_slice.empty()) {
    Polynomial * q = factors_per_slice.back();
    factors_per_slice.pop_back();

    Polynomial * add = add_poly(p, q);
    if (print) {
      print_pac_add_rule(file, p, q, add);
      if (xor_chain) add = mod_poly(add, 1, file);
      print_pac_del_rule(file, p);
      print_pac_del_rule(file, q);
    }
    delete(p);
    delete(q);
    p = add;
  }

  return p;
}
/*------------------------------------------------------------------------*/
/**
    Merges computed factors of a slice of the same level.
    Prints PAC rules for the process. Used only when proof == 1 or proof == 2.

    @param file output file for PAC rules
    @param container vector containing the factors, ordered by level

    @return vector of polynomials such that each polynomial has a unique level
*/
static std::vector<Polynomial*> add_and_merge_factors( FILE * file,
    Polynomial * add_p, std::vector<Polynomial*> container, bool print) {
  if (container.empty()) {
    container.push_back(add_p);
    return container;
  }

  Polynomial * p = add_p;
  Polynomial * q = container.back();

  if (p->get_level() != q->get_level()) {
    container.push_back(add_p);
  } else {
    Polynomial * add = add_poly(p, q);

    if (print)  {
      print_pac_add_rule(file, p, q, add);
      if (xor_chain) add = mod_poly(add, 1, file);

      print_pac_del_rule(file, p);
      print_pac_del_rule(file, q);
    }
    add->set_level(p->get_level()+1);
    delete(p);
    delete(q);
    container.pop_back();
    container = add_and_merge_factors(file, add, container, print);
  }
  return container;
}

/*------------------------------------------------------------------------*/

Polynomial * add_up_spec_of_slice(FILE * file, bool print) {
  Polynomial * p = spec_of_slice.back();
  spec_of_slice.pop_back();

  while (!spec_of_slice.empty()) {
    Polynomial * q = spec_of_slice.back();
    spec_of_slice.pop_back();

    Polynomial * add = add_poly(p, q);
    if (print) {
      print_pac_add_rule(file, p, q, add);
      if (xor_chain) add = mod_poly(add, 1, file);
      print_pac_del_rule(file, p);
      print_pac_del_rule(file, q);
    }

    delete(p);
    delete(q);
    p = add;
  }
  return p;
}

/*------------------------------------------------------------------------*/

void eliminate_by_one_gate(Gate * n1, Gate *n2, FILE * file) {
  Polynomial * p1 = n1->get_gate_constraint();
  Polynomial * p2 = n2->get_gate_constraint();

  const Polynomial * negfactor = divide_by_term(p1, p2->get_lt());
  if (negfactor->is_constant_zero_poly()) return;

  Polynomial * mult   = multiply_poly(negfactor, p2);
  Polynomial * rem    = add_poly(p1, mult);

  if (proof == 1) {
    assert(file);

    if (!negfactor->is_constant_one_poly())
      print_pac_mul_rule(file, p2, negfactor, mult);
    else
      mult->set_idx(p2->get_idx());

    print_pac_add_rule(file, p1, mult, rem);
    print_pac_del_rule(file, p1);
    if (!negfactor->is_constant_one_poly()) print_pac_del_rule(file, mult);
  } else if (proof == 2) {
    assert(file);
    if (!negfactor->is_constant_one_poly())
      print_pac_combi_rule(file, p2, negfactor, p1, 0, rem);
    else
      print_pac_add_rule(file, p1, p2, rem);

    print_pac_del_rule(file, p1);
  } else if (proof == 3) {
    if (rem) add_ancestors(n1, n2, negfactor);
  }

  delete(mult);
  delete(negfactor);
  delete(p1);
  n1->set_gate_constraint(rem);
}

/*------------------------------------------------------------------------*/
const Polynomial * cofac_tmp, * base_tmp, * mul_tmp;

Polynomial * reduce_by_one_poly(
    const Polynomial * p1, Gate * n, FILE * file) {
  Polynomial * p2 = n->get_gate_constraint();
  const Polynomial * negfactor = divide_by_term(p1, p2->get_lt());
  if (negfactor->is_constant_zero_poly()) return p1->copy();

  Polynomial * mult   = multiply_poly(negfactor, p2);
  Polynomial * rem    = add_poly(p1, mult);

  if (!proof) {
    delete(mult);
    delete(negfactor);
  } else if (proof == 1) {
    assert(file);
    if (!negfactor->is_constant_one_poly())
      print_pac_mul_rule(file, p2, negfactor, mult);
    else
      mult->set_idx(p2->get_idx());

    if (mult) {
      factors_per_slice = add_and_merge_factors(file, mult, factors_per_slice, 1);

      if (!negfactor->is_constant_one_poly()) print_pac_del_rule(file, p2);
    } else {
      delete(mult);
    }

    delete(negfactor);
  } else if (proof == 2) {
    assert(file);
    if (mult) {
      co_factors.push_back(negfactor);
      factor_indices.push_back(p2->get_idx());
      factors_per_slice = add_and_merge_factors(file, mult, factors_per_slice, 0);
    } else {
      delete(mult);
      delete(negfactor);
    }
  } else if (proof == 3) {
    add_fac(n, negfactor);
    delete(mult);
    delete(negfactor);
  }

  return rem;
}

/*------------------------------------------------------------------------*/

void remove_internal_xor_gates(FILE * file) {
  msg("remove internal xor gates");
  int counter = 0;
  for (unsigned i = NN; i < M-1; i++) {
    Gate * n = gates[i];
    if (n->get_xor_gate() != 1) continue;
    if (n->get_elim()) continue;
    assert(n->children_size() == 2);

    Gate * l_gate = n->children_front();
    Gate * r_gate = n->children_back();
    if (l_gate->get_xor_gate() != 2) continue;
    if (r_gate->get_xor_gate() != 2) continue;
    assert(l_gate->children_size() == 2);
    assert(r_gate->children_size() == 2);
    if(l_gate->parents_size() != 1 && r_gate->parents_size() != 1) continue;
    Gate * ll_gate = l_gate->children_front();
    Gate * lr_gate = l_gate->children_back();

    // set xor children to ll and lr by overwriting l and r

    n->set_children_front(ll_gate);
    n->set_children_back(lr_gate);

    lr_gate->parents_push_back(n);
    ll_gate->parents_push_back(n);
    eliminate_by_one_gate(n, l_gate, file);


    if (l_gate->parents_size() == 1) {
      if (proof == 1 || proof == 2) {
        assert(file);
        print_pac_del_rule(file, l_gate->get_gate_constraint());
      }
      l_gate->mark_elim();
      delete(l_gate->get_gate_constraint());
      l_gate->set_gate_constraint(0);

      if (verbose >= 3) msg("removed %s", l_gate->get_var_name());
      counter++;

      ll_gate->parents_remove(l_gate);
      lr_gate->parents_remove(l_gate);
    } else {
      l_gate->parents_remove(n);
    }

    eliminate_by_one_gate(n, r_gate, file);

    if (r_gate->parents_size() == 1) {
      if (proof == 1 || proof == 2) {
        assert(file);
        print_pac_del_rule(file, r_gate->get_gate_constraint());
      }

      r_gate->mark_elim();
      delete(r_gate->get_gate_constraint());
      r_gate->set_gate_constraint(0);

      if (verbose >= 3) msg("removed %s", r_gate->get_var_name());
      counter++;

      ll_gate->parents_remove(r_gate);
      lr_gate->parents_remove(r_gate);
    } else {
      r_gate->parents_remove(n);
    }
  }
  if (verbose >= 1) msg("removed %i internal xor gates", counter);
}

/*------------------------------------------------------------------------*/

void remove_single_occs_gates(FILE * file) {
  msg("remove single occurence gates");
  int counter = 0;
  for (unsigned i = NN; i < M-1; i++) {
    Gate * n = gates[i];
    if (n->get_elim()) continue;
    if (n->parents_size() > 1) continue;

    //consider unconnected nodes
    if (n->parents_size() == 0){
      n->mark_elim();
      for (std::list<Gate*>::const_iterator it = n->children_begin();
          it != n->children_end(); ++it) {
        Gate * n_child = *it;
        n_child->parents_remove(n);
      }
      continue;
    }

    Gate * parent = n->parents_front();
    if (parent->get_output()) continue;
    if (!n->get_xor_gate() && parent->get_xor_gate() == 1) continue;
    if (n->get_xor_chain()) continue;

    eliminate_by_one_gate(parent, n, file);
    parent->children_remove(n);

    for (std::list<Gate*>::const_iterator it = n->children_begin();
        it != n->children_end(); ++it) {
      Gate * n_child = *it;
      if (!parent->is_child(n_child)) parent->children_push_back(n_child);

      n_child->parents_remove(n);
      if (!n_child->is_in_parents(parent)) n_child->parents_push_back(parent);
    }

    if (proof == 1 || proof == 2) {
      assert(file);
      print_pac_del_rule(file, n->get_gate_constraint());
    }

    n->mark_elim();
    delete(n->get_gate_constraint());
    n->set_gate_constraint(0);

    if (verbose >= 3) msg("removed %s", n->get_var_name());
    counter++;
  }
  if (verbose >= 1) msg("removed %i single occurence gates", counter);
}

/*------------------------------------------------------------------------*/
int remove_not_assigned_gate(FILE * file, Gate * n, int count){
  if(n->get_input()) return count;
  if(n->get_elim()) return count;
//  if (n->get_slice() > -1) return count;



  for (std::list<Gate*>::const_iterator it_c = n->children_begin();
      it_c != n->children_end(); ++it_c) {
    Gate * n_child = *it_c;
    n_child->parents_remove(n);
  }

  for (std::list<Gate*>::const_iterator it_p = n->parents_begin();
      it_p != n->parents_end(); ++it_p) {
    Gate * n_parent = *it_p;
    if (n_parent->get_elim()) continue;
    eliminate_by_one_gate(n_parent, n, file);
    n_parent->children_remove(n);

    for (std::list<Gate*>::const_iterator it_c = n->children_begin();
        it_c != n->children_end(); ++it_c) {
      Gate * n_child = *it_c;

      if (!n_parent->is_child(n_child))
        n_parent->children_push_back(n_child);
      if (!n_child->is_in_parents(n_parent))
        n_child->parents_push_back(n_parent);
    }
  }

  if (proof == 1 || proof == 2) {
    assert(file);
    print_pac_del_rule(file, n->get_gate_constraint());
  }

  n->mark_elim();
  count++;
  delete(n->get_gate_constraint());
  n->set_gate_constraint(0);

  if (verbose >= 3) msg("removed %s", n->get_var_name());

  for (std::list<Gate*>::const_iterator it_c = n->children_begin();
      it_c != n->children_end(); ++it_c) {
    Gate * n_child = *it_c;
    count = remove_not_assigned_gate(file, n_child, count);

  }

  return count;

}


void remove_slice_minus_one_gates(FILE * file) {
  msg("remove gates that are not assigned to slices");
  int counter = 0;
  for (unsigned i = NN; i < M-1; i++) {
    Gate * n = gates[i];

    if (n->get_elim()) continue;
    if (n->get_slice() > -1) continue;
    assert(!n->get_input());
  
    counter = remove_not_assigned_gate(file, n, counter);
  }
  if (verbose >= 1)
    msg("removed %i gates that are not assigned to slices", counter);
}

/*------------------------------------------------------------------------*/

void decomposing(FILE * file) {
  msg("eliminate single occs");
  int counter = 0;
  int begin = xor_chain ? NN-2 : NN-1;
  for (int i = begin; i >= 0; i--) {
    bool change = 1;
    while (change) {
      change = 0;
      for (std::list<Gate*>::const_iterator it =slices[i].begin();
          it != slices[i].end(); ++it) {
        Gate * n = *it;
        if (n->get_elim()) continue;

        if (n->parents_size() == 1 && !n->get_carry_gate()) {
          Gate * parent = n->parents_front();

          eliminate_by_one_gate(parent, n, file);
          parent->children_remove(n);


          for (std::list<Gate*>::const_iterator it=n->children_begin();
              it != n->children_end(); ++it) {
            Gate * n_child = *it;
            if (!parent->is_child(n_child))
              parent->children_push_back(n_child);

            n_child->parents_remove(n);

            if (!n_child->is_in_parents(parent))
              n_child->parents_push_back(parent);
          }
          if (proof == 1 || proof == 2) {
            assert(file);
            print_pac_del_rule(file, n->get_gate_constraint());
          }

          n->mark_elim();
          delete(n->get_gate_constraint());
          n->set_gate_constraint(0);

          counter++;

          ++it;
          slices[i].remove(n);
          change = 1;

          if (verbose >= 3)
            msg("decomposed %s", n->get_var_name());
        }
      }
    }
  }
  msg("decomposed %i variables", counter);
}

/*------------------------------------------------------------------------*/

void eliminate_booth_pattern(FILE * file) {
  msg("eliminate booth pattern");
  int counter = 0;
  for (unsigned i = NN; i < M-1; i++) {
    Gate * n = gates[i];
    if (!n->get_bo()) continue;
    if (n->get_elim()) continue;

    for (std::list<Gate*>::const_iterator it_c=n->children_begin();
        it_c != n->children_end(); ++it_c) {
      Gate * n_child = *it_c;
      n_child->parents_remove(n);
    }

    for (std::list<Gate*>::const_iterator it_p=n->parents_begin();
        it_p != n->parents_end(); ++it_p) {
      Gate * n_parent = *it_p;
      eliminate_by_one_gate(n_parent, n, file);
      n_parent->children_remove(n);
      for (std::list<Gate*>::const_iterator it_c=n->children_begin();
          it_c != n->children_end(); ++it_c) {
        Gate * n_child = *it_c;
        n_parent->children_push_back(n_child);
        n_child->parents_push_back(n_parent);
      }
    }

    if (proof == 1 || proof == 2) {
      assert(file);
      print_pac_del_rule(file, n->get_gate_constraint());
    }

    n->mark_elim();
    delete(n->get_gate_constraint());
    n->set_gate_constraint(0);

    counter++;
    if (verbose >= 3) {
      msg("eliminated %s", n->get_var_name());
    }
  }
  msg("eliminated %i variables from booth pattern", counter);
}

/*------------------------------------------------------------------------*/

Polynomial * inc_spec_poly(unsigned i) {
  mpz_t coeff;
  mpz_init(coeff);

  mpz_pow_ui(coeff, base, i);
  if (i == NN-1 && signed_mult) mpz_neg(coeff, coeff);

  const Var * s = gates[i+M-1]->get_var();
  Term * t1 = new_term(s, 0);
  Monomial * m1 = new Monomial(coeff, t1);
  push_mstack_end(m1);

  mpz_neg(coeff, coeff);

  int min_start = std::min(NN/2-1, i);
  for (int j = min_start; j >= 0; j--) {
    if (mpz_sgn(coeff) == 1) mpz_neg(coeff, coeff);
    const Var * b = gates[b0+j*binc]->get_var();
    unsigned k = i - j;
    if (k > NN/2-1) break;
    if (k == NN/2-1 && signed_mult) mpz_neg(coeff, coeff);
    if (j == static_cast<int>(NN/2-1) && signed_mult) mpz_neg(coeff, coeff);

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

Polynomial * mod_poly(const Polynomial *p1, bool print_rule, FILE * file) {
  bool pac_print = print_rule &(proof == 1 || proof == 2);
  int exp = NN;

  mpz_t coeff;
  mpz_init(coeff);

  for (size_t i = 0 ; i < p1->size(); i++) {
    Monomial * m = p1->get_mon(i);
    mpz_tdiv_r_2exp(coeff, m->coeff, exp);
    if (mpz_sgn(coeff) != 0) {
      Monomial * tmp;
      if (m->get_term()) tmp =  new Monomial(coeff, m->get_term_copy());
      else
        tmp =  new Monomial(coeff, 0);
      push_mstack_end(tmp);
    }
  }
  mpz_clear(coeff);
  Polynomial * out = build_poly();
  out->set_idx(p1->get_idx());

  if (pac_print || proof == 3) {
    mpz_t quot;
    mpz_init(quot);
    for (size_t i = 0 ; i < p1->size(); i++) {
      Monomial * m = p1->get_mon(i);

      mpz_tdiv_q_2exp(quot, m->coeff, exp);
      if (mpz_sgn(quot) != 0) {
        mpz_neg(quot, quot);
        Monomial * tmp;
        if (m->get_term()) tmp =  new Monomial(quot, m->get_term_copy());
        else tmp =  new Monomial(quot, 0);
        push_mstack_end(tmp);
      }
    }
    mpz_clear(quot);
  }

  if (pac_print && !mstack_is_empty()) {
    assert(file);
    Polynomial * p = build_poly();

    Polynomial * mod = multiply_poly_with_constant(p, mod_coeff);
    print_pac_mod_rule(file,  p, mod);
    print_pac_add_rule(file,  p1, mod, out);

    delete(p);
    delete(mod);

  } else if (proof == 3 && !mstack_is_empty()) {
    Polynomial * p = build_poly();
    add_fac_mod(p);
    delete(p);
  }

  delete(p1);
  return out;
}

/*------------------------------------------------------------------------*/

void correct_pp_unsigned(const Polynomial * p, FILE * file) {
  for (size_t i = 0 ; i < p->size(); i++) {
    Monomial * m = p->get_mon(i);

    if (mpz_sgn(m->coeff) == -1) {
      Term * t = m->get_term();
      if (!t) continue;
      if (t->get_var_num() <= 0) continue;
      if (!gate(t->get_var_num())->get_input()) continue;
      push_mstack_end(new Monomial(one, t->copy()));
    }
  }

  if (!mstack_is_empty()) {
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

void correct_pp_signed(const Polynomial * p, FILE * file) {
  mpz_t half_mod;
  mpz_init(half_mod);
  mpz_pow_ui(half_mod, base, NN-1);

  mpz_t half_mod_neg;
  mpz_init(half_mod_neg);
  mpz_neg(half_mod_neg, half_mod);

  for (size_t i = 0 ; i < p->size(); i++) {
    Monomial * m = p->get_mon(i);
    Term * t = m->get_term();
    if (!t) continue;
    if (t->get_var_num() <= 0) continue;
    if (!gate(t->get_var_num())->get_input()) continue;

    if (mpz_cmp(m->coeff, half_mod) > 0) {
      Monomial * tmp = new Monomial(minus_one, m->get_term_copy());
      push_mstack_end(tmp);
    } else if (mpz_cmp(m->coeff, half_mod_neg) < 0) {
      Monomial * tmp = new Monomial(one,  m->get_term_copy());
      push_mstack_end(tmp);
    }
  }

  mpz_clear(half_mod);
  mpz_clear(half_mod_neg);

  if (!mstack_is_empty()) {
    Polynomial * factor = build_poly();
    if (!factor) return;

    Polynomial * mod = multiply_poly_with_constant(factor, mod_coeff);
    print_pac_mod_rule(file,  factor, mod);

    Polynomial * add = add_poly(mod, p);
    print_pac_add_rule(file,  mod, p, add);

    delete(mod);
    delete(factor);
    delete(add);
  }
}

/*------------------------------------------------------------------------*/

void correct_pp(const Polynomial * p, FILE * file) {
  if (signed_mult) correct_pp_signed(p, file);
  else
    correct_pp_unsigned(p, file);
}

/*------------------------------------------------------------------------*/

const Polynomial * reduce(FILE * file) {
  msg("");
  msg("");
  msg("started reducing");
  Polynomial * rem = 0, * tmp;
  for (int i=NN-1; i>= 0; i--) {
    if (verbose >= 1) msg("reducing by slice %i", i);
    Polynomial * inc_spec = inc_spec_poly(i);

    if (rem) {
      tmp = add_poly(inc_spec, rem);
      delete(inc_spec);
      delete(rem);
      rem = tmp;
    } else {
      rem = inc_spec;
    }

    std::list<Gate*> sl = slices[i];
    for (std::list<Gate*>::const_iterator it=sl.begin(); it != sl.end(); ++it) {
      Gate * n = *it;
      if (n->get_elim()) continue;

      if (verbose >= 4  && n->get_gate_constraint()) {
        fputs("[amulet2] reducing by ", stdout);
        n->print_gate_constraint(stdout);
      }
      tmp = reduce_by_one_poly(rem, n, file);
      delete(n->get_gate_constraint());
      n->set_gate_constraint(0);

      if (xor_chain) tmp = mod_poly(tmp, 0, file);

      delete(rem);
      rem = tmp;

      if (verbose >= 3) {
        fputs("[amulet2] remainder is ", stdout);
        rem->print(stdout);
        msg(" ");
      }
    }

    if (verbose >= 2) {
      msg("after reducing by slice %i", i);
      fprintf(stdout, "[amulet2] remainder is ");
      rem->print(stdout);
      msg("");
    }

    if (proof == 1 || proof == 2) {
      Polynomial * pac_poly = add_up_factors(file, proof == 1);
      if (proof == 2){

        print_pac_vector_combi_rule(file, factor_indices, co_factors, pac_poly);
        factor_indices.clear();
        co_factors.clear();

        if (xor_chain) pac_poly = mod_poly(pac_poly, 1, file);
        spec_indices.push_back(pac_poly->get_idx());
      }
      pac_poly->set_level(1);
      spec_of_slice = add_and_merge_factors(file, pac_poly, spec_of_slice, proof == 1);

    }
  }

  if (proof == 1)  {
    Polynomial * res = add_up_spec_of_slice(file, 1);
    if (xor_chain) correct_pp(res, file);
    delete(res);
  } else if (proof == 2) {
    Polynomial * res = add_up_spec_of_slice(file, 0);
    print_pac_vector_add_rule(file, spec_indices, res);
    spec_indices.clear();

    if (xor_chain) res = mod_poly(res, 1, file);
    if (xor_chain) correct_pp(res, file);
    delete(res);
  }

  return rem;
}

/*------------------------------------------------------------------------*/

bool check_inputs_only(const Polynomial *p) {
  for (size_t i = 0 ; i < p->size(); i++) {
    Monomial * m = p->get_mon(i);
    if (!m->get_term()) { continue;
    } else {
      const Var * v = m->get_term()->get_var();
      Gate * n = gate(v->get_num());
      if (!n->get_input()) return 0;
    }
  }
  return 1;
}

/*------------------------------------------------------------------------*/

void write_witness_vector(const Term * t, FILE * file) {
  fputs_unlocked("[amulet2]   ", stdout);

  if(ainc == 2){
    for (unsigned i = 0; i <= NN/2-1; i++) {
      const Var * v = gates[a0+i*ainc]->get_var();

      if (t->contains(v)) {
        fprintf(file, "1");
        fprintf(stdout, "%s = ", v->get_name());
      } else {
        fprintf(file, "0");
      }

      const Var * w = gates[b0+i*binc]->get_var();

      if (t->contains(w)) {
        fprintf(file, "1");
        fprintf(stdout, "%s = ", w->get_name());
      } else {
        fprintf(file, "0");
      }
    }

  } else if (ainc == 1){
    for (unsigned i = 0; i <= NN/2-1; i++) {
      const Var * v = gates[a0+i*ainc]->get_var();

      if (t->contains(v)) {
        fprintf(file, "1");
        fprintf(stdout, "%s = ", v->get_name());
      } else {
        fprintf(file, "0");
      }
    }

    for (unsigned i = 0; i <= NN/2-1; i++) {
      const Var * w = gates[b0+i*binc]->get_var();

      if (t->contains(w)) {
        fprintf(file, "1");
        fprintf(stdout, "%s = ", w->get_name());
      } else {
        fprintf(file, "0");
      }
    }
  }


  fprintf(stdout, "1, all other inputs = 0;\n");
  fprintf(file, "\n");
}

/*------------------------------------------------------------------------*/

void write_witnesses(const Polynomial * p, FILE * file) {
  assert(check_inputs_only(p));

  int len = p->min_term_size();
  if (len == 0) {
    msg("  all inputs = 0;\n");
    for (unsigned i = 0; i <= NN/2-1; i++)
      fprintf(file, "00");

    fprintf(file, "\n");
  } else {
    for (size_t i = 0 ; i < p->size(); i++) {
      Monomial * m = p->get_mon(i);
      if (m->get_term()) {
        int tlen = m->get_term_size();
        if (tlen == len) write_witness_vector(m->get_term(), file);
      }
    }
  }
  fprintf(file, ".");
}

/*------------------------------------------------------------------------*/

void generate_witness(const Polynomial * p, const char * name) {
  if (!check_inputs_only(p))
  die(err_witness, "cannot generate witness, as remainder polynomial contains non-inputs");

  #define WITNESSLEN 100
  char witness_name[WITNESSLEN];
  memset(witness_name, '\0', sizeof(witness_name));
  for (int i = 0; name[i] != '.'; i++) {
    witness_name[i] = name[i];
  }
  snprintf(witness_name + strlen(witness_name),
           WITNESSLEN - strlen(witness_name), "%s", ".cex");

  FILE * witness_file;
  if (!(witness_file = fopen(witness_name, "w")))
  die(err_writing, "cannot write output to '%s'", witness_name);

  msg("");
  msg("COUNTER EXAMPLES ARE: ");

  write_witnesses(p, witness_file);

  msg("");
  msg("");
  msg("Counter examples are written to %s", witness_name);
  msg("You can run 'aigsim' from the AIGER library (http://fmv.jku.at/aiger/)");
  msg("to simulate the provided counter example(s).");
  msg("");
  msg("Note: 'aiger/aigsim %s %s' produces output in the form:",
  name, witness_name);
  if(ainc == 2){
    fputs_unlocked("[amulet2] ", stdout);
    if (NN == 2)      fprintf(stdout, "  a[0]b[0]  s[0]\n");
    else if (NN == 4) fprintf(stdout, "  a[0]b[0]a[1]b[1]  s[0]s[1]s[2]s[3]\n");
    else
    fprintf(stdout, "  a[0]b[0]a[1]b[1]...a[%u]b[%u]  s[0]s[1]s[2]...s[%u]\n",
      NN/2-1, NN/2-1, NN-1);
  } else {
    fputs_unlocked("[amulet2] ", stdout);
    if (NN == 2)      fprintf(stdout, "  a[0]b[0]  s[0]\n");
    else if (NN == 4) fprintf(stdout, "  a[0]a[1]b[0]b[1]  s[0]s[1]s[2]s[3]\n");
    else
    fprintf(stdout, "  a[0]a[1]...a[%u]b[0]b[1]...b[%u]  s[0]s[1]s[2]...s[%u]\n",
      NN/2-1, NN/2-1, NN-1);
  }

  fclose(witness_file);
}
