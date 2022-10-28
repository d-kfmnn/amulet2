/*------------------------------------------------------------------------*/
/*! \file nss.cpp
    \brief contains functions necessary to generate Nullstellensatz proofs

  Part of AMulet2 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include <map>

#include "nss.h"
/*------------------------------------------------------------------------*/
static Polynomial * mod_factor;
/*------------------------------------------------------------------------*/
void print_spec_poly(FILE * file) {
  mpz_t coeff;
  mpz_init(coeff);

  // outputs
  for (int i = NN-1; i >= 0; i--) {
    const Var * v = gates[i+M-1]->get_var();

    mpz_pow_ui(coeff, base, i);
    mpz_neg(coeff, coeff);
    if (i == static_cast<int>(NN-1) && signed_mult) mpz_neg(coeff, coeff);

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
      if (i == static_cast<int>(NN/2-1) && signed_mult) mpz_neg(coeff, coeff);
      if (j == static_cast<int>(NN/2-1) && signed_mult) mpz_neg(coeff, coeff);
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

void print_cofactors_poly_nss(FILE * file) {
  bool first = 1;
  fprintf(file, "%i %% ", M+1);
  for (unsigned i = num_gates-1; i >= NN ; i--) {
    Polynomial * p = gates[i]->get_cofactor();

    if (p && !p->is_constant_zero_poly()) {
      if (!first) fputs(" + ", file);
      fprintf(file, "%i *(", 2+i-NN);
      p->print(file, 0);
      fputs(")\n", file);
      first = 0;
    }
  }

  if (mod_factor && !mod_factor->is_constant_zero_poly()) {
    fprintf(file, " + 1 *(");
    mod_factor->print(file, 0);
    fputs(")\n", file);
    delete(mod_factor);
  }

  fputs(" , ", file);
  print_spec_poly(file);
}
/*------------------------------------------------------------------------*/

void add_ancestors(
    Gate * n, Gate * anc, const Polynomial * fac, bool internal) {
  assert(n);
  assert(anc);

  if (!fac) return;

  /* if function is called from extern on a gate 'n' we need to add the tuple
   *(n,1) to the ancestors, to indicate that the original polynomial of n
   * is used in the ancestors.
   */
  if (!internal) {
    std::map<Gate*, Polynomial*>::const_iterator it = n->search_in_anc(n);
    if (it == n->anc_end()) {

      Monomial * m1 = new Monomial(one, 0);
      push_mstack_end(m1);
      Polynomial * p = build_poly();
      n->set_ancestor(n, p);
    }
  }

  /* if 'anc' is original or add_ancestors is called internally
   * we check whether 'anc' is already contained in the ancestors.
   * If so, we add 'fac' to the corresponding co-factor of 'anc' in the ancestors.
   * Otherwise, the tuple(anc,fac) is added to the ancestors.
   */
  if (anc->orig() || internal) {
    std::map<Gate*, Polynomial*>::iterator it = n->search_in_anc(anc);
    if (it != n->anc_end()) {
      Polynomial * old_fac = it->second;
      Polynomial * new_fac = add_poly(fac, old_fac);
      it->second = new_fac;
      delete(old_fac);
    } else {
      n->set_ancestor(anc, fac->copy());
    }
  } else {
    /* if 'anc' is not original we traverse through its ancestors anc call
     * add_ancestors from internal
     */
    for (std::map<Gate*, Polynomial*>::const_iterator it=anc->anc_begin();
        it != anc->anc_end(); ++it) {
      Polynomial * new_fac = multiply_poly(fac, it->second);
      add_ancestors(n, it->first, new_fac, 1);
      delete(new_fac);
    }
  }
}
/*------------------------------------------------------------------------*/

void add_fac_mod(const Polynomial * fac) {
  if (!fac) return;

  if (!mod_factor) { mod_factor = fac->copy();
  } else {
    Polynomial * add = add_poly(mod_factor, fac);
    delete(mod_factor);
    mod_factor = add;
  }
}
/*------------------------------------------------------------------------*/

void add_fac(Gate * n, const Polynomial * fac) {
  assert(n);
  if (!fac) return;

  if (n->orig()) {
    /* If 'n' is original we add 'fac' as co-factor. */
    if (!n->get_cofactor()) { n->set_cofactor(fac->copy());
    } else {
      Polynomial * add = add_poly(n->get_cofactor(), fac);
      delete(n->get_cofactor());
      n->set_cofactor(add);
    }
  } else {
    /* Otherwise, we traverse through the ancestors of 'n' and add the product
     * of 'fac' and ancestor cofactor to the cofactors of the ancestor gates.
     */
    for (std::map<Gate*, Polynomial*>::const_iterator it = n->anc_begin();
        it != n->anc_end(); ++it) {
      Gate       * n_anc = it->first;
      Polynomial * p_anc = it->second;

      Polynomial * mul = multiply_poly(p_anc, fac);

      if (!n_anc->get_cofactor()) { n_anc->set_cofactor(mul);
      } else {
        Polynomial * add = add_poly(n_anc->get_cofactor(), mul);
        delete(mul);
        delete(n_anc->get_cofactor());
        n_anc->set_cofactor(add);
      }
    }
  }
}
