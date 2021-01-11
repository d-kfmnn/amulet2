/*------------------------------------------------------------------------*/
/*! \file slicing.cpp
    \brief contains functions slice and order the gates

  Part of AMulet2.0 : AIG Multiplier Verification Tool.
  Copyright (C) 2020 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include "slicing.h"
/*------------------------------------------------------------------------*/
// Global var
std::vector<std::list<Gate*>> slices ;
/*------------------------------------------------------------------------*/

void init_slices(){
 for (unsigned i=0; i<NN; i++){
   std::list<Gate*> l;
   Gate * n = gates[i+M-1];
   assert(n->get_output());
   l.push_back(n);
   slices.push_back(l);
 }
}

/*------------------------------------------------------------------------*/

 void print_slices(){
  for (int i=NN-1; i>=0; i--){
    msg("");
    msg("Slice %i", i);
    msg("");
    std::list<Gate*> l = slices[i];
    for (std::list<Gate*>::const_iterator it=l.begin(); it != l.end(); ++it){
      Gate * n = *it;
      fprintf(stdout, "[amulet2] ");
      n->print_gate_constraint(stdout);

    }
  }
}

/*------------------------------------------------------------------------*/

const Gate * topological_largest_child(const Gate * n){
  int slice = n->children_front()->get_slice();
  for (std::list<Gate*>::const_iterator it_n = n->children_begin();
       it_n != n->children_end(); ++it_n){

    Gate * n_child = *it_n;
    if(n_child->get_slice() != slice)
      die("error in topological_largest_child with %s", n->get_var_name());
  }


  Gate * n1 = n->children_front();

  for (std::list<Gate*>::const_iterator it_s = slices[n1->get_slice()].begin();
      it_s != slices[n1->get_slice()].end(); ++it_s){

    Gate * n_slice = *it_s;
    for (std::list<Gate*>::const_iterator it_n = n->children_begin();
         it_n != n->children_end(); ++it_n){

      Gate * n_child = *it_n;
      if(n_slice == n_child) return n_child;
    }
  }
  return 0;
}

/*------------------------------------------------------------------------*/

void fix_slice(Gate * n, int i){
  const Gate * after = topological_largest_child(n);

 if(n->get_slice() != -1) slices[n->get_slice()].remove(n);
 if(n->get_elim()) return;
 n->set_slice(i);
 for (std::list<Gate*>::const_iterator it = slices[i].begin();
      it != slices[i].end(); ++it){

   Gate * n_it = *it;
   if(n_it == after) {
     slices[i].insert(it, n);
     return;
   }
 }
 slices[i].push_back(n);
}

/*------------------------------------------------------------------------*/

void slice_by_xor_chains(){
 for (unsigned i = 0; i < NN; i++){
   Gate * output = gates[i+M-1];
   output->set_slice(i);

   assert(output->children_size() == 1);
   Gate * child = output->children_front();

   std::queue<Gate*> downwards_queue;

   if(i != NN-1) {
     downwards_queue.push(child);
     child->set_slice(i);
     slices[i].push_back(child);
   } else if (child->get_xor_gate() == 1 ) {
      downwards_queue.push(child);
      child->set_slice(i);
      slices[i].push_back(child);
   }

   while(!downwards_queue.empty()){
     Gate * n = downwards_queue.front();
     downwards_queue.pop();
      for (std::list<Gate*>::const_iterator it = n->children_begin();
           it != n->children_end(); ++it){
        Gate * n_child = *it;
        if (n_child->get_slice() == -1 && !n_child->get_aig_output() &&
           (n_child->get_xor_gate() == 1 || n_child->get_pp())){

           if(!n_child->get_pp()) downwards_queue.push(n_child);
           n_child->set_slice(i);
           slices[i].push_back(n_child);

        }
        else if (!n_child->get_input() && !n_child->get_carry_gate()){
           n_child->set_carry_gate(i);
        }
      }
    }
  }
}

/*------------------------------------------------------------------------*/

void upwards_slicing(const Gate * n, const Gate * pre){
  for (std::list<Gate*>::const_iterator it = n->parents_begin();
       it != n->parents_end(); ++it){
    Gate * n_parent = *it;
    if(n_parent->get_slice() != -1) continue;
    if(n_parent->get_xor_gate() == 1) continue;
    if(n_parent->get_output()) continue;
    n_parent->set_slice(pre->get_slice());
    if(!pre->is_child(n_parent)) {
      for (std::list<Gate*>::const_iterator it = slices[pre->get_slice()].begin();
          it != slices[pre->get_slice()].end(); ++it){

        Gate * cmp = *it;
        if(cmp == pre){
          slices[pre->get_slice()].insert(it, n_parent);
          break;
        }
      }
    }
    else {
      for (std::list<Gate*>::const_iterator it = slices[pre->get_slice()].begin();
           it != slices[pre->get_slice()].end(); ++it){

        Gate * cmp = *it;
        if(cmp == pre){
          slices[pre->get_slice()].insert(++it, n_parent);
          break;
        }
      }
    }
    if(!n_parent->get_carry_gate()) upwards_slicing(n_parent, n_parent);
  }
}

/*------------------------------------------------------------------------*/

void slice_jut_gates(){
 for (int i = NN-1; i >= 0; i--){
   Gate * output = gates[i+M-1];
   output->set_slice(i);

   assert(output->children_size() == 1);
   Gate * child = output->children_front();

   std::queue<Gate*> downwards_queue;
   if(child->get_xor_gate() == 1 || i != (int) NN-1) downwards_queue.push(child);

   while(!downwards_queue.empty()){
     Gate * n = downwards_queue.front();
     downwards_queue.pop();


     for (std::list<Gate*>::const_iterator it = n->children_begin();
          it != n->children_end(); ++it){

         Gate * n_child = *it;
         if(n_child->get_slice() == (int) i) downwards_queue.push(n_child);

         if ((n_child->get_xor_gate() == 1 || n_child->get_pp())){
             if(n_child->get_slice() == (int) i)
               upwards_slicing(n_child,n_child);
         }
         else if (n_child->get_carry_gate() == (int)i && !n_child->get_input()){
           upwards_slicing(n_child,n);
         }
       }
    }
  }
}

/*------------------------------------------------------------------------*/

int fix_xors(){
 int counter = 0;
 for(unsigned i = NN; i < M-1; i++){
   Gate * n = gates[i];
   if(n->get_elim()) continue;
   if(n->get_xor_gate() != 1) continue;
   if(n->get_aig_output()) continue;

   aiger_and * and1 = is_model_and(n->get_var_num());
   unsigned l = and1->rhs0, r = and1->rhs1;
   Gate *l_gate = gate(l), *r_gate = gate(r);

   aiger_and * land = is_model_and(l);
   unsigned ll = land->rhs0, lr = land->rhs1;
   Gate *ll_gate = gate(ll), *lr_gate = gate(lr);
   if(ll_gate->get_slice() < lr_gate->get_slice()) std::swap(ll_gate,lr_gate);

   if(ll_gate->parents_size() + lr_gate->parents_size() <= 3) {
       if(ll_gate->get_slice() == n->get_slice()-1 &&
          lr_gate->get_slice() == ll_gate->get_slice()) {

       if(!r_gate->get_elim()) fix_slice(r_gate, ll_gate->get_slice());
       if(!l_gate->get_elim()) fix_slice(l_gate, ll_gate->get_slice());

       fix_slice(n, n->get_slice()-1);
       counter++;
       n->mark_moved();
       if (verbose >= 3)
         msg("moved gate %s to slice %i", n->get_var_name(), n->get_slice());

     }
   }

   else if(ll_gate->parents_size() == 2 && lr_gate->parents_size() == 2) {
       if(ll_gate->get_slice() == n->get_slice()-1 &&
          lr_gate->get_slice() == ll_gate->get_slice() &&
         (ll_gate->get_moved() || lr_gate->get_moved())) {

       if(!r_gate->get_elim()) fix_slice(r_gate, ll_gate->get_slice());
       if(!l_gate->get_elim()) fix_slice(l_gate, ll_gate->get_slice());

       fix_slice(n, n->get_slice()-1);
       ++counter;
       if (verbose >= 3)
         msg("moved gate %s to slice %i", n->get_var_name(), n->get_slice());

     }
   }
 }
 if (verbose >= 1) msg("moved %i gates to smaller slices", counter);
 return counter;
}

/*------------------------------------------------------------------------*/

void fix_jut_gates(){
 int counter = 0;
 for(unsigned i = NN; i < M-1; i++){
   Gate * n = gates[i];

   if(n->get_xor_gate()) continue;
   if(n->get_elim()) continue;
   if(n->get_pp()) continue;
   if(n->children_size() < 4) continue; //CHANGED THIS FROM 3 to 4

   int slice = n->get_slice() - 1;
   int flag = 0;

   for (std::list<Gate*>::const_iterator it=n->children_begin();
        it != n->children_end(); ++it){

     Gate * n_child = *it;
     if(n_child->get_input()) continue;
     if(n_child->get_slice() != slice) {
       flag=1;
       break;
     }

   }
   if(flag) continue;

   fix_slice(n, n->get_slice()-1);
   ++counter;
   if (verbose >= 3)
     msg("moved gate %s to slice %i", n->get_var_name(), n->get_slice());
 }
 if (verbose >= 1) msg("moved %i adjacent gates to smaller slices", counter);

}

/*------------------------------------------------------------------------*/

void slicing_xor(){
   slice_by_xor_chains();
   slice_jut_gates();
   if(fix_xors()) fix_jut_gates();
}

/*------------------------------------------------------------------------*/

void input_cone(Gate * n, int num){
    assert(num >= 0);
    if (n->get_input()) return;
    if (n->get_slice() >= 0) return;

    assert (n->get_slice() < 0);
    assert (is_model_and(n->get_var_num()));
    n->set_slice(num);
    for (std::list<Gate*>::const_iterator it = n->children_begin();
         it != n->children_end(); ++it){

      Gate * n_child = *it;
      input_cone(n_child, num);
    }
}

/*------------------------------------------------------------------------*/

void find_carries(){
  for(unsigned j = M-1; j > NN; j--){
    Gate * n = gates[j];
    if(n->get_elim()) continue;
    n->set_carry_gate(0);
    for (std::list<Gate*>::const_iterator it=n->parents_begin();
         it != n->parents_end(); ++it){

      Gate * n_parent = *it;
      if(n_parent->get_slice() > n->get_slice()) {
        n->inc_carry_gate();
      }
    }
  }
}

/*------------------------------------------------------------------------*/

bool search_for_booth_pattern(){
  bool found_booth = 0;
  for(unsigned i = NN; i < M-1; i++){
    Gate * n = gates[i];
    if(n->get_elim()) continue;

    //Special treatment for slice 1
    if(n->get_slice() == 1){
      if(!n->get_pp()) continue;

      aiger_and * and1 = is_model_and(n->get_var_num());
      unsigned l = aiger_strip(and1->rhs0), r = aiger_strip(and1->rhs1);
      if (!gate(l)->get_input() || !gate(r)->get_input()) continue;

      if(l-r != 2) continue;

      if (verbose >=4) msg("found booth pattern %s", n->get_var_name());
      n->mark_bo();
      found_booth = 1;

    //all other slices
    } else {

      if(n->get_pp()) continue;
      if(!n->get_xor_gate()) continue;
      if(n->parents_size() != 1) continue;

      aiger_and * and_init = is_model_and(n->get_var_num());

      unsigned l = and_init->rhs0;
      if (!aiger_sign (l)) continue;

      aiger_and * land = is_model_and(l);
      unsigned ll = aiger_strip(land->rhs0), lr = aiger_strip(land->rhs1);
      if (!gate(ll)->get_input() || !gate(lr)->get_input()) continue;


      Gate * xor1 = n;
      aiger_and * and1 = is_model_and(xor1->get_var_num());

      Gate * vp = xor1->parents_front();
      unsigned p = vp->get_var_num();
      aiger_and * parent = is_model_and(p);
      unsigned pl = parent->rhs0, pr = parent->rhs1;
      Gate * xor2 = ((int)pl == xor1->get_var_num()) ? gate(pl) : gate(pr);

      if(xor2->parents_size() < NN/2+1) continue;

      if(xor2->get_slice() >= xor1->get_slice()) continue;

      aiger_and * and2 = is_model_and(xor2->get_var_num());
      if(!xor2->get_xor_gate()) continue;

      unsigned l2 = and2->rhs0;
      if (!aiger_sign (l2)) continue;
      aiger_and * land2 = is_model_and(l2);
      unsigned ll2 = aiger_strip(land2->rhs0), lr2 = aiger_strip(land2->rhs1);
      if (!gate(ll2)->get_input() || !gate(lr2)->get_input()) continue;
      if(ll != ll2 && ll != lr2 && lr != ll2 && lr != lr2) continue;

      found_booth = 1;
      xor1->mark_bo();
      gate(and1->rhs0)->mark_bo();
      gate(and1->rhs1)->mark_bo();
      xor2->mark_bo();
      gate(and2->rhs0)->mark_bo();
      gate(and2->rhs1)->mark_bo();
      vp->mark_bo();

      if (verbose >=4)
        msg("found booth pattern %s, %s, %s", xor1->get_var_name(),
                                              xor2->get_var_name(),
                                              vp->get_var_name());
    }
  }
  return found_booth;
}

/*------------------------------------------------------------------------*/

void merge_all(){
  int total_merge = 0;
  int merged = 1;
  while(merged){
    merged = 0;
    for(unsigned i = M-2; i > NN; i--){
      Gate * n = gates[i];

      if(n->get_slice() < 0) continue; // elim
      if(n->get_elim()) continue;
      if(is_model_input(n->get_var_num())) continue;

      if(n->get_xor_gate() == 2) {
        aiger_and * and1 = is_model_and(n->get_var_num());
        unsigned rhs0 = and1->rhs0;
        unsigned rhs1 = and1->rhs1;
        Gate * v0 = gate(rhs0);
        Gate * v1 = gate(rhs1);
        if(!(v0->get_slice() == v1->get_slice() &&
             v1->get_slice() < n->get_slice() &&
            !v0->get_pp() && !v1->get_pp()))
          continue;
      }

      if(n->get_xor_gate() == 1 && !n->get_aig_output()) continue;
      if(n->get_xor_gate() == 1 && n->parents_size() > 1 ) continue;

      bool flag = 0;
      for (std::list<Gate*>::const_iterator it=n->children_begin();
           it != n->children_end(); ++it){

        Gate * n_child = *it;
        if(n_child->get_input()) flag = 1;
        else if(n_child->get_slice() == n->get_slice()) flag = 1;
        else if(n_child->get_bo()) flag = 1;
        if (flag) break;
      }
      if (flag) continue;

      n->dec_slice();
      n->set_carry_gate(0);
      for (std::list<Gate*>::const_iterator it = n->parents_begin();
           it != n->parents_end(); ++it){

        Gate * n_parent = *it;
        if(n_parent->get_slice() > n->get_slice()) {
          n->inc_carry_gate();
        }
      }
      for (std::list<Gate*>::const_iterator it = n->children_begin();
           it != n->children_end(); ++it){

        Gate * n_child = *it;
        if(n->get_slice() == n_child->get_slice()){
          n_child->dec_carry_gate();
        }
      }

      merged = 1;
      if(verbose >= 3)
        msg("merged gate %s to slice %i", n->get_var_name(), n->get_slice());
      total_merge = total_merge +1;
    }
  }
  msg("totally merged %i variable(s)", total_merge);
}

/*------------------------------------------------------------------------*/

void promote_all(){
  int total_promote = 0;
  int promoted = 1;
  while(promoted){
    promoted = 0;
    for(unsigned i = NN; i < M; i++){
      Gate * n = gates[i];
      if(!n->get_carry_gate()) continue;
      if(n->get_pp()) continue;
      if((int) n->parents_size() != n->get_carry_gate()) continue;


      aiger_and * and1 = is_model_and(n->get_var_num());
      unsigned rhs0 = and1->rhs0;
      unsigned rhs1 = and1->rhs1;
      Gate * v0 = gate (rhs0);
      Gate * v1 = gate (rhs1);

      if(n->get_xor_gate()!=2 &&
         (!v0->get_carry_gate() || !v1->get_carry_gate()) &&
         (!v0->get_carry_gate() || !v1->get_input()) &&
         (!v1->get_carry_gate() || !v0->get_input()))  continue;

      bool flag = 0;
      for (std::list<Gate*>::const_iterator it = n->parents_begin();
        it != n->parents_end(); ++it){

        Gate * n_parent = *it;
        if(n_parent->get_slice() == n->get_slice()) flag = 1;
      }
      if (flag) continue;


      n->inc_slice();
      v0->inc_carry_gate();
      v1->inc_carry_gate();
      n->set_carry_gate(0);

      for (std::list<Gate*>::const_iterator it = n->parents_begin();
        it != n->parents_end(); ++it){

        Gate * n_parent = *it;
        if(n_parent->get_slice() > n->get_slice()) n->inc_carry_gate();
      }

      promoted = 1;
      total_promote++;

      if(verbose >= 3)
        msg("promoted gate %s to slice %i", n->get_var_name(), n->get_slice());
    }
  }
  msg("totally promoted %i variable(s)", total_promote);
}

/*------------------------------------------------------------------------*/

void fill_slices(){
  for (unsigned i = 0; i <= NN-1; i++) {
    for(unsigned j = M-2; j >= NN; j--){
      Gate * n = gates[j];
      if(n->get_slice() == (int) i) {
        slices[i].push_back(n);
      }
    }
  }
  msg("filled %i slices", NN);
}

/*------------------------------------------------------------------------*/

void slicing_non_xor(){
  for (unsigned i = 0; i < NN; i++)
    input_cone(gate(slit(i)), i);

  find_carries();
  merge_all();
  promote_all();
  fill_slices();
}
