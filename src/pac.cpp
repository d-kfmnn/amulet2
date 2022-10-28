/*------------------------------------------------------------------------*/
/*! \file pac.cpp
    \brief contains functions necessary to generate PAC proofs

  Part of AMulet2 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include "pac.h"
/*------------------------------------------------------------------------*/
static int poly_idx;
/*------------------------------------------------------------------------*/
// ERROR CODES:
static int err_rule       = 81; // error in proof rule
/*------------------------------------------------------------------------*/

void print_circuit_poly(FILE * file) {
  fputs("1 ", file);
  mpz_out_str(file, 10, mod_coeff);
  fputs(";\n", file);

  for (unsigned i = NN; i < num_gates; i++) {
    Polynomial * p = gen_gate_constraint(i);
    assert(p);

    fprintf(file, "%i ", p->get_idx());
    p->print(file);
    delete(p);
  }
  poly_idx = M+1;
}

/*------------------------------------------------------------------------*/

void print_pac_del_rule(FILE * file, const Polynomial *p1) {
  assert(p1);

  fprintf(file, "%i d;\n", p1->get_idx());
}

/*------------------------------------------------------------------------*/

void print_pac_mod_rule(FILE * file, const Polynomial *p1, Polynomial *p) {
  assert(p1 && !p1->is_constant_zero_poly());
  assert(p  && !p->is_constant_zero_poly());

  fprintf(file, "%u ", poly_idx);
  fputs("% 1 *(", file);
  p1->print(file, 0); fputs("), ", file);
  p->print(file);

  p->set_idx(poly_idx++);
}

/*------------------------------------------------------------------------*/



void print_pac_add_rule(
    FILE * file, const Polynomial *p1, const Polynomial *p2, Polynomial *p) {
  assert(p1 && !p1->is_constant_zero_poly());
  assert(p2 && !p2->is_constant_zero_poly());
  assert(p  && !p->is_constant_zero_poly());

  fprintf(file, "%u %% %i + %i, ", poly_idx, p1->get_idx(), p2->get_idx());
  p->print(file);

  p->set_idx(poly_idx++);
}

/*------------------------------------------------------------------------*/

void print_pac_vector_add_rule(
  FILE * file, std::vector<int> indices, Polynomial * p){

  fprintf(file, "%u %% ", poly_idx);

  int ind;
  while (!indices.empty()) {
    ind = indices.back();
    indices.pop_back();

    fprintf(file, "%i",ind);
    if (indices.size() > 0) fputs(" + ", file);
  }
  fprintf(file, ", ");
  p->print(file);
  p->set_idx(poly_idx++);
}
/*------------------------------------------------------------------------*/

void print_pac_combi_rule(
    FILE * file, const Polynomial *p1, const Polynomial *p2,
    const Polynomial *p3, const Polynomial *p4, Polynomial *p) {
  assert(p1 && !p1->is_constant_zero_poly());
  assert(p2 && !p2->is_constant_zero_poly());
  assert(p3 && !p3->is_constant_zero_poly());
  assert(p  && !p->is_constant_zero_poly());

  fprintf(file, "%u %% %i", poly_idx, p1->get_idx());
  if (p2) {
    fputs(" *(", file);
    p2->print(file, 0); fputs(") ", file);
  }

  fprintf(file, "+ %i", p3->get_idx());
  if (p4) {
    fputs(" *(", file);
    p4->print(file, 0); fputs(") ", file);
  }

  fprintf(file, ", ");
  p->print(file);
  p->set_idx(poly_idx++);
}

/*------------------------------------------------------------------------*/

void print_pac_vector_combi_rule(
  FILE * file, std::vector<int> indices,
  std::vector<const Polynomial*> co_factors, Polynomial * p){

  if (co_factors.size() != indices.size()) die(err_rule, "combination rule receives invalid arguments;");

  fprintf(file, "%u %% ", poly_idx);

  const Polynomial * tmp;
  int ind;
  while (!co_factors.empty()) {
    tmp = co_factors.back();
    co_factors.pop_back();
    ind = indices.back();
    indices.pop_back();

    fprintf(file, "%i",ind);
    if (!tmp->is_constant_one_poly()){
      fputs(" *(", file);
      tmp->print(file, 0); fputs(")", file);
    }

    if (co_factors.size() > 0) fputs(" + ", file);
    delete(tmp);
  }
  fprintf(file, ", ");
  p->print(file);
  p->set_idx(poly_idx++);
}
/*------------------------------------------------------------------------*/

void print_pac_mul_rule(
    FILE * file, const Polynomial *p1, const Polynomial *p2, Polynomial *p) {
  assert(p1 && !p1->is_constant_zero_poly());
  assert(p2 && !p2->is_constant_zero_poly());
  assert(p  && !p->is_constant_zero_poly());

  fprintf(file, "%u %% %i *(", poly_idx, p1->get_idx());
  p2->print(file, 0); fputs("), ", file);
  p->print(file);

  p->set_idx(poly_idx++);
}

/*------------------------------------------------------------------------*/
