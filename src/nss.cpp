/*------------------------------------------------------------------------*/
/*! \file nss.cpp
    \brief contains functions necessary to generate Nullstellensatz proofs

  Part of AMulet2.0 : AIG Multiplier Verification Tool.
  Copyright (C) 2020 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include "nss.h"
/*------------------------------------------------------------------------*/
static Polynomial * mod_factor;
/*------------------------------------------------------------------------*/

void print_circuit_poly_nss (FILE * file){
  for (unsigned i = num_gates-1; i >= NN ; i--) {
    Polynomial * p = gates[i]->get_gate_constraint();
    p->print(file);
  }
}
/*------------------------------------------------------------------------*/

void print_cofactors_poly_nss (FILE * polysfile, FILE * file){

  for (unsigned i = num_gates-1; i >= NN ; i--) {
    Polynomial * p = gates[i]->get_cofactor();
    if(p) p->print(file);
    else  fputs("0;\n", file);
  }

  if(mod_factor){
    mpz_out_str(polysfile, 10, mod_coeff);
    fprintf(polysfile, ";\n");
    mod_factor->print(file);
    delete(mod_factor);
  }
}
/*------------------------------------------------------------------------*/

void add_ancestors (Gate * n, Gate * anc, const Polynomial * fac, bool internal) {
  assert(n);
  assert(anc);

  if(!fac) return;

  /* if function is called from extern on a gate 'n' we need to add the tuple
   * (n,1) to the ancestors, to indicate that the original polynomial of n
   * is used in the ancestors.
   */
  if(!internal){
      std::map<Gate*, Polynomial*>::const_iterator it = n->search_in_anc(n);
      if (it == n->anc_end()) {
        Polynomial * p = new Polynomial();
        Monomial * m1 = new Monomial(one, 0);
        p->mon_push_back(m1);
        n->set_ancestor(n, p);
      }
  }

  /* if 'anc' is original or add_ancestors is called internally
   * we check whether 'anc' is already contained in the ancestors.
   * If so, we add 'fac' to the corresponding co-factor of 'anc' in the ancestors.
   * Otherwise, the tuple (anc,fac) is added to the ancestors.
   */
  if(anc->orig() || internal){
    std::map<Gate*, Polynomial*>::iterator it = n->search_in_anc(anc);
    if (it != n->anc_end()){
      Polynomial * old_fac = it->second;
      Polynomial * new_fac = add_poly(fac, old_fac);
      it->second = new_fac;
      delete(old_fac);
    }
    else n->set_ancestor(anc, fac->copy());
  }

  /* if 'anc' is not original we traverse through its ancestors anc call
   * add_ancestors from internal
   */
  else {
    for (std::map<Gate*, Polynomial*>::const_iterator it=anc->anc_begin();
         it != anc->anc_end(); ++it)
    {
      Polynomial * new_fac = multiply_poly(fac, it->second);
      add_ancestors(n, it->first, new_fac, 1);
      delete(new_fac);
    }
  }
}
/*------------------------------------------------------------------------*/

void add_fac_mod (const Polynomial * fac){
  if(!fac) return;
  if(!mod_factor) mod_factor = fac->copy();
  else {
    Polynomial * add = add_poly(mod_factor, fac);
    delete(mod_factor);
    mod_factor = add;
  }
}
/*------------------------------------------------------------------------*/

void add_fac (Gate * n, const Polynomial * fac){
  assert(n);
  if(!fac) return;

  /* If 'n' is original we add 'fac' as co-factor. */
  if (n->orig()) {
    if(!n->get_cofactor()) n->set_cofactor(fac->copy());
    else {
      Polynomial * add = add_poly(n->get_cofactor(), fac);
      delete(n->get_cofactor());
      n->set_cofactor(add);
    }
  }
  /* Otherwise, we traverse through the ancestors of 'n'  and add the product
   * of 'fac' and the ancestor cofactor to the cofactors of the ancestor gates.
   */
  else {
    for (std::map<Gate*, Polynomial*>::const_iterator it = n->anc_begin();
    it != n->anc_end(); ++it)
    {
      Gate       * n_anc = it->first;
      Polynomial * p_anc = it->second;

      Polynomial * mul = multiply_poly(p_anc, fac);

      if (!n_anc->get_cofactor()) n_anc->set_cofactor(mul);
      else {
        Polynomial * add = add_poly(n_anc->get_cofactor(), mul);
        delete(mul);
        delete(n_anc->get_cofactor());
        n_anc->set_cofactor(add);
      }
    }
  }
}
