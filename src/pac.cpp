/*------------------------------------------------------------------------*/
/*! \file pac.cpp
    \brief contains functions necessary to generate PAC proofs

  Part of AMulet2.0 : AIG Multiplier Verification Tool.
  Copyright (C) 2020 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include "pac.h"
/*------------------------------------------------------------------------*/
static int poly_idx;
/*------------------------------------------------------------------------*/
void print_spec_poly (FILE * file){
  mpz_t coeff;
  mpz_init(coeff);

  // outputs
  for (int i = NN-1; i >= 0; i--){
    const Var * v = gates[i+M-1]->get_var();

    mpz_pow_ui(coeff, base, i);
    mpz_neg(coeff, coeff);
    if(i == (int)NN-1 && signed_mult) mpz_neg(coeff, coeff);

    Term * t1 = new_term(v, 0);
    Monomial * m1 = new Monomial(coeff, t1);
    push_mstack_end(m1);
  }

  // partial products
  for (int i = NN/2-1; i >= 0; i--) {
    const Var * a = gates[a0+i*ainc]->get_var();

    for (int j = NN/2-1; j >= 0; j--) {
      const Var * b = gates[b0+j*binc]->get_var();
      mpz_pow_ui(coeff, base, i+j);
      if (i == (int)NN/2-1 && signed_mult) mpz_neg(coeff, coeff);
      if (j == (int)NN/2-1 && signed_mult) mpz_neg(coeff, coeff);
      add_to_vstack(a);
      add_to_vstack(b);
      Term * t1 = build_term_from_stack();
      Monomial * m1 = new Monomial(coeff, t1);
      push_mstack_end(m1);
    }
  }
  mpz_clear(coeff);

  Polynomial * p = build_poly();
  p->print(file);
  delete(p);
}

/*------------------------------------------------------------------------*/

void print_circuit_poly (FILE * file){
  fputs("1 ", file);
  mpz_out_str(file, 10, mod_coeff);
  fputs(";\n", file);

  for (unsigned i = NN; i < num_gates; i++) {
    Polynomial * p = gates[i]->get_gate_constraint();
    assert(p);

    fprintf(file, "%i ", p->get_idx());
    p->print(file);
  }
  poly_idx = gates[num_gates-1]->get_gate_constraint()->get_idx();
  ++poly_idx;
}

/*------------------------------------------------------------------------*/

void print_pac_del_rule (FILE * file, const Polynomial *p1){
  assert(p1);

  fprintf(file, "%i d;\n", p1->get_idx());
}

/*------------------------------------------------------------------------*/

void print_pac_mod_rule (FILE * file, const Polynomial *p1, Polynomial *p){
  assert(p1 && !p1->is_constant_zero_poly());
  assert(p  && !p->is_constant_zero_poly());

  fprintf(file, "%u ", poly_idx);
  fputs("* 1, ", file);
  p1->print(file, 0); fputs(", ", file);
  p->print(file);

  p->set_idx(poly_idx++);
}

/*------------------------------------------------------------------------*/

void print_pac_add_rule
  (FILE * file, const Polynomial *p1, const Polynomial *p2, Polynomial *p){

  assert(p1 && !p1->is_constant_zero_poly());
  assert(p2 && !p2->is_constant_zero_poly());
  assert(p  && !p->is_constant_zero_poly());

  fprintf(file, "%u + %i, %i, ", poly_idx, p1->get_idx(), p2->get_idx());
  p->print(file);

  p->set_idx(poly_idx++);
}

/*------------------------------------------------------------------------*/

void print_pac_mul_rule
  (FILE * file, const Polynomial *p1, const Polynomial *p2, Polynomial *p){

  assert(p1 && !p1->is_constant_zero_poly());
  assert(p2 && !p2->is_constant_zero_poly());
  assert(p  && !p->is_constant_zero_poly());

  fprintf(file, "%u * %i, ", poly_idx, p1->get_idx());
  p2->print(file, 0); fputs(", ", file);
  p->print(file);

  p->set_idx(poly_idx++);
}

/*------------------------------------------------------------------------*/
