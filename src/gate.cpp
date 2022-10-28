/*------------------------------------------------------------------------*/
/*! \file gate.cpp
    \brief contains the class Gate and further functions to
    organize the gate structure, such as initializing the gate constraints

  Part of AMulet2 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include <string>
#include <list>

#include "gate.h"
/*------------------------------------------------------------------------*/
// Global variables
bool xor_chain = 0;
bool booth = 0;
bool signed_mult = 0;
/*------------------------------------------------------------------------*/
// ERROR CODES:
static int err_allocate       = 91; // failed to allocate gates
/*------------------------------------------------------------------------*/

Gate::Gate(int n_, std::string name_, int level_, bool input_, bool output_):
  v(new Var(name_, level_, n_)), input(input_), output(output_)  {
}

/*------------------------------------------------------------------------*/
bool Gate::orig() const {
  return ancestors.empty();
}

/*------------------------------------------------------------------------*/

Gate::~Gate() {
  if(v) delete(v);
  delete(co_factor);
  if (gate_constraint) delete(gate_constraint);
  for (std::map<Gate*, Polynomial*>::const_iterator it=ancestors.begin();
      it != ancestors.end(); ++it) {
    delete(it->second);
  }
}

/*------------------------------------------------------------------------*/
bool Gate::all_parents_are_sliced() const {
  for (auto it=parents_begin(); it != parents_end(); ++it) {
    Gate * parents = *it;
    if (parents->get_slice() == -1) return 0;
  }
  return 1;
}

/*------------------------------------------------------------------------*/

bool Gate::is_in_parents(const Gate * n) const {
  for (auto it=parents_begin(); it != parents_end(); ++it) {
    Gate * parents = *it;
    if (parents == n) return 1;
  }
  return 0;
}

/*------------------------------------------------------------------------*/

bool Gate::is_child(const Gate * n) const {
  for (std::list<Gate*>::const_iterator it = children_begin();
      it != children_end(); ++it) {
    Gate * child = *it;
    if (child == n) return 1;
  }
  return 0;
}

/*------------------------------------------------------------------------*/
Polynomial * Gate::get_gate_constraint() const {
  if (!gate_constraint && !elim) {
    // output aig are 0, -1, ...-NN+2
    if (output) init_gate_constraint(-1*get_var_num()+M-1);
    // gates are numbered 2,4,6,8,..
    else
      init_gate_constraint(get_var_num()/2-1);
  }
  return gate_constraint;
}

Gate ** gates;
unsigned num_gates;

/*------------------------------------------------------------------------*/

Gate * gate(unsigned lit) {
  assert(lit < 2*M);
  if(lit < 2) return 0;
  return gates[lit/2-1];
}

/*------------------------------------------------------------------------*/

void allocate_gates(bool assert) {
  unsigned aiger;
  num_gates = M + NN - 1;

  msg("allocating %i gates", num_gates);
  gates = new Gate*[num_gates];

  if (!gates) die(err_allocate, "failed to allocate gates");
  int level = 0;

  // inputs a
  for (unsigned i = a0; i <= al; i+=ainc) {
    aiger = 2*(i+1);
    if (assert) assert(is_model_input(aiger));

    std::string name = "a" + std::to_string((i-a0)/ainc);
    gates[i] = new Gate(aiger, name, ++level, 1);
  }

  // inputs b
  for (unsigned i = b0; i <= bl; i+=binc) {
    aiger = 2*(i+1);
    if (assert) assert(is_model_input(aiger));

    std::string name = "b" + std::to_string((i-b0)/binc);
    gates[i] = new Gate(aiger, name, ++level, 1);
  }

  // internal gates
  for (unsigned i = NN; i < M-1; i++) {
    aiger = 2*(i+1);
    if (assert) assert(is_model_and(aiger));

    std::string name = "l" + std::to_string(aiger);
    gates[i] = new Gate(aiger, name, ++level);
  }

  // output s
  for (unsigned i = M-1; i < num_gates; i++) {
    aiger = i-M+1;

    std::string name = "s" + std::to_string(aiger);
    gates[i] = new Gate(M-i-1, name, ++level, 0, 1);
    gates[i]->set_slice(aiger);
  }
}

/*------------------------------------------------------------------------*/

void mark_aig_outputs() {
  for (unsigned i = 0; i < NN; i++) {
    unsigned lit = slit(i);
    if(lit < 2) continue;
    Gate * n = gate(lit);
    n->mark_aig_output();
  }
}

/*------------------------------------------------------------------------*/

Gate * derive_ha_and_gate(const Gate *n) {
  Gate * ll = xor_left_child(n);
  Gate * rr = xor_right_child(n);
  for (std::list<Gate*>::const_iterator it_l=ll->parents_begin();
       it_l != ll->parents_end(); ++it_l) {
    Gate * ll_parent = *it_l;
    if (ll_parent -> get_xor_gate()) continue;

    aiger_and * and1 = is_model_and(ll_parent->get_var_num());
    Gate * c1 = gate(and1->rhs0);
    Gate * c2 = gate(and1->rhs1);

    if (c1 == ll && c2 == rr) return ll_parent;
    if (c1 == rr && c2 == ll) return ll_parent;
  }
  return 0;
}

/*------------------------------------------------------------------------*/

void set_xor() {
  int found_xor = 0;
  for (unsigned i = 0; i < M; i++) {
    Gate * n = gates[i];
    if (n->get_input() > 0) continue;
    if (n->get_xor_gate() > 0) continue;

    aiger_and * and1 = is_model_and(n->get_var_num());
    if (!and1) continue;
    unsigned l = and1->rhs0, r = and1->rhs1;
    if (!aiger_sign(l)) continue;
    if (!aiger_sign(r)) continue;
    if (l == r || l == aiger_not(r)) continue;
    l = aiger_strip(l);
    r = aiger_strip(r);
    aiger_and * land = is_model_and(l);
    if (!land) continue;
    aiger_and * rand = is_model_and(r);
    if (!rand) continue;

    unsigned ll = land->rhs0, lr = land->rhs1;
    unsigned rl = rand->rhs0, rr = rand->rhs1;
    if ((ll == aiger_not(rl) && lr == aiger_not(rr)) ||
       (ll == aiger_not(rr) && lr == aiger_not(rl))) {
      gate(l)->set_xor_gate(2);
      gate(r)->set_xor_gate(2);
      n->set_xor_gate(1);
      found_xor++;
      if (verbose >= 4) msg("xor-gate %s", n->get_var_name());
    }
  }
  if (verbose >= 1) msg("found %i xor-gates", found_xor);
}

/*------------------------------------------------------------------------*/

bool upper_half_xor_output() {
  for (unsigned i = num_gates-2; i > M-1; i--) {
    Gate * n = gates[i];
    if(!n->children_size()) return 0;
    n = n->children_front();
    if (!n->get_xor_gate()) return 0;
  }

  if(!gates[M-1]->children_size()) return 0;

  return 1;
}

/*------------------------------------------------------------------------*/

Gate * xor_left_child(const Gate * n) {
  if (!n->get_xor_gate() ) return 0;

  aiger_and * and1 = is_model_and(n->get_var_num());
  if (!and1) return 0;
  unsigned l = and1->rhs0;
  if (!aiger_sign(l)) return 0;
  l = aiger_strip(l);

  aiger_and * land = is_model_and(l);
  if (!land) return 0;

  unsigned ll = land->rhs0;
  return gate(ll);
}

/*------------------------------------------------------------------------*/

Gate * xor_right_child(const Gate * n) {
  if (!n->get_xor_gate() ) return 0;

  aiger_and * and1 = is_model_and(n->get_var_num());
  if (!and1) return 0;
  unsigned l = and1->rhs0;
  if (!aiger_sign(l)) return 0;
  l = aiger_strip(l);

  aiger_and * land = is_model_and(l);
  if (!land) return 0;

  unsigned lr = land->rhs1;
  return gate(lr);
}

/*------------------------------------------------------------------------*/

void mark_xor_chain_in_last_slice() {
  msg("marking xor chain gates");
  int counter = 0;

  Gate * out = gates[num_gates-1];
  if(out->children_size() == 0) return;
  assert(out->children_size() == 1);
  Gate * child = out->children_front();
  if(child->get_input()) return;

  std::queue<Gate*> downwards_queue;
  if (child->get_xor_gate() == 1) downwards_queue.push(child);

  while (!downwards_queue.empty()) {
    Gate * n = downwards_queue.front();
    downwards_queue.pop();

    aiger_and * and1 = is_model_and(n->get_var_num());
    unsigned l = and1->rhs0;
    aiger_and * land = is_model_and(l);
    unsigned ll = land->rhs0, lr = land->rhs1;
    Gate * ll_gate = gate(ll);
    Gate * lr_gate = gate(lr);

    if (ll_gate->get_xor_gate()) downwards_queue.push(ll_gate);
    if (lr_gate->get_xor_gate()) downwards_queue.push(lr_gate);

    n->mark_xor_chain();
    if (verbose >= 4) msg("xor-chain %s", n->get_var_name());
    counter++;
  }
  if (verbose >= 1) msg("marked %i xor gates in last slice", counter);
  if (counter) xor_chain = 1;
}

/*------------------------------------------------------------------------*/

void set_parents_and_children(bool set_children) {
  unsigned pp = 0;

  for (unsigned i = NN; i < M; i++) {
    Gate * n = gates[i];
    assert(!n->get_input());

    aiger_and * and1 = is_model_and(n->get_var_num());

    if (!and1) continue;
    unsigned l = and1->rhs0, r = and1->rhs1;
    Gate * l_gate = gate(l), *r_gate = gate(r);

    if (set_children) {
      n->children_push_back(l_gate);
      n->children_push_back(r_gate);
    }

    if(l_gate && r_gate) {
      n->set_level(l_gate->get_level() > r_gate->get_level() ?
                 l_gate->get_level()+1 : r_gate->get_level()+1);

      if (l_gate->get_input() && r_gate->get_input()
            && !aiger_sign(l) && !aiger_sign(r)) {
        n->mark_pp();
        pp++;
        if (verbose >= 4) msg("partial product %s", n->get_var_name());
      }

    } else if (l_gate) {
      n->set_level(l_gate->get_level()+1);
    } else if (r_gate){
      n->set_level(r_gate->get_level()+1);
    }

    if(l_gate) l_gate->parents_push_back(n);
    if(r_gate) r_gate->parents_push_back(n);


  }

  // set children for extra outputs
  for (unsigned i = 0; i < NN; i++) {
    Gate * n = gates[i+M-1];
    assert(n->get_output());
    unsigned lit = slit(i);
    if(lit < 2) continue;
    Gate * model_output_gate = gate(lit);
    if (set_children) n->children_push_back(model_output_gate);
    model_output_gate->parents_push_back(n);
  }

  if (verbose >= 1) msg("found %i partial products", pp);
  if (pp == NN*NN/4) {
    msg("assuming simple pp generator");
  } else {
    booth = 1;
    msg("assuming booth recoding");
  }
}

/*------------------------------------------------------------------------*/

bool parents_are_in_equal_or_larger_slice(const Gate * n, int i) {
  for (std::list<Gate*>::const_iterator it_n = n->parents_begin();
       it_n != n->parents_end(); ++it_n) {
    Gate * n_parent = *it_n;
    if (n_parent->get_slice() != -1 && n_parent->get_slice() < i) return 0;
  }
  return 1;
}

/*------------------------------------------------------------------------*/

Polynomial * negative_poly(const Var * v) {

  Term * t1 = new_term(v, 0);
  Monomial * m1 = new Monomial(minus_one, t1);
  push_mstack_end(m1);

  Monomial * m2 = new Monomial(one, 0);
  push_mstack_end(m2);

  return build_poly();
}

/*------------------------------------------------------------------------*/

Polynomial * positive_poly(const Var * v) {
  Term * t = new_term(v, 0);
  Monomial * m = new Monomial(one, t);
  push_mstack_end(m);

  return build_poly();
}

/*------------------------------------------------------------------------*/
static Polynomial * get_node_constraint(Gate * g, unsigned sign){
  if(g) {
    const Var * v1 = g->get_var();
    if (sign) return negative_poly(v1);
    else return positive_poly(v1);
  } else {
    if(sign) {
      push_mstack_end(new Monomial(one, 0));
      return build_poly();
    } else return 0;
  }
}

/*------------------------------------------------------------------------*/

Polynomial * gen_gate_constraint(unsigned i) {
  assert(i >= NN && i < M + NN - 1);
  Polynomial * p = 0;
  // gate constraint
  if (i < M-1) {
    Gate * n = gates[i];
    assert(!n->get_input());

    aiger_and * and1 = is_model_and(n->get_var_num());
    assert(and1);

    unsigned l = and1->rhs0, r = and1->rhs1;

    Gate * l_gate = gate(l), *r_gate = gate(r);
    const Var * v = n->get_var();
    Term * t1 = new_term(v, 0);
    Monomial * m1 = new Monomial(minus_one, t1);
    push_mstack_end(m1);
    Polynomial * p_h = build_poly();

    Polynomial * p1 = get_node_constraint(l_gate, aiger_sign(and1->rhs0));
    Polynomial * p2 = get_node_constraint(r_gate, aiger_sign(and1->rhs1));
    Polynomial * p_tl = multiply_poly(p1, p2);

    p = add_poly(p_h, p_tl);
    delete(p_h);
    delete(p_tl);

    delete(p1);
    delete(p2);
  } else {  // output
    Gate * n = gates[i];
    assert(n->get_output());

    Gate * model_output_gate = n->children_front();
    if (!model_output_gate) {
      unsigned lit = slit(i-M+1);
      const Var * v = n->get_var();
      Term * t1 = new_term(v, 0);
      Monomial * m1 = new Monomial(minus_one, t1);
      push_mstack_end(m1);
      if(lit == 1){
        push_mstack_end(new Monomial(one, 0));
      }
      p = build_poly();

    } else {
      const Var * v = n->get_var();
      Term * t1 = new_term(v, 0);
      Monomial * m1 = new Monomial(minus_one, t1);

      Polynomial * p_tl;
      if (aiger_sign(slit(i-M+1)))
        p_tl = negative_poly(model_output_gate->get_var());
      else
        p_tl = positive_poly(model_output_gate->get_var());

      push_mstack_end(m1);
      Polynomial * p_h = build_poly();

      p = add_poly(p_h, p_tl);
      delete(p_h);
      delete(p_tl);
    }

  }
  p->set_idx(2+i-NN);
  return p;
}


void init_gate_constraint(unsigned i) {
  assert(i >= NN && i < M + NN - 1);
  Polynomial * p = gen_gate_constraint(i);
  Gate * n = gates[i];
  n->set_gate_constraint(p);
}


void init_gate_constraints() {
  for (unsigned i = NN; i < M-1; i++) {
    init_gate_constraint(i);
  }

  for (unsigned i = 0; i < NN; i++) {
    init_gate_constraint(i+M-1);
  }
}

/*------------------------------------------------------------------------*/

void delete_gates() {
  for (unsigned i = 0; i < num_gates; i++) {
    delete(gates[i]);
  }
  delete[] gates;
}
