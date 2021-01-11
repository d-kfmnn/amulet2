/*------------------------------------------------------------------------*/
/*! \file gate.h
    \brief contains the class Gate and further functions to
    organize the gate structure, such as initializing the gate constraints

  Part of AMulet2.0 : AIG Multiplier Verification Tool.
  Copyright (C) 2020 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#ifndef _gate_H
#define _gate_H
/*------------------------------------------------------------------------*/
#include <map>
#include <queue>

#include "aig.h"
#include "polynomial.h"
/*------------------------------------------------------------------------*/
/// set to true when the last slice contains an xor_chain
extern bool xor_chain;

/// set to true when the multiplier contains a booth pattern
extern bool booth;

/// set to true when a signed or unsigned multiplier is verified
extern bool signed_mult;
/*------------------------------------------------------------------------*/

/** \class Gate
  Internal structure to represent the AIG graph.
*/
class Gate {

  /// Variable of the gate, as used in the polynomials
  const Var * v;

  /// True if gate is an input
  bool input;

  /// True if the gate is an output s_i
  bool output;

  /// True if the gate is an output in the aig
  bool aig_output = 0;

  /// True if gate is identified as a partial product
  bool partial_product = 0;

  /// Distance to inputs
  int level = 0;

  /// True if gate belongs to xor_chain in last slice
  bool xor_chain = 0;

  /// is set to 1 for root node, 2 for internal nodes of XORs
  int  xor_gate = 0;

  /// Counts how often gate is used in bigger slice
  int carry_gate = 0;

  /// slice a gate is attached to
  int slice = -1;

  /// True if circuit is a prop_gen_gate (-substitute)
  bool prop_gen_gate = 0;

  /// True if gate is identified to belong to complex fsa (-substitute)
  bool fsa = 0;

  /// True if gate is input of complex fsa (-substitute)
  int fsa_inp = 0;

  /// True if gate occurs negative (-substitute)
  bool neg = 0;

  /// True if gate has been moved during fix_xors
  bool moved = 0;

  /// True if gate is identified to belong to booth pattern
  bool bo = 0;

  /// True if gate is eliminated during preprocessing
  bool elim = 0;

  /// Polynomial implied by the aig gate
  Polynomial * gate_constraint = 0;

  /// Polynomial generated as co-factor for nss proofs (-certify)
  Polynomial * co_factor = 0;

  /// Used to store all dependencies in nss proofs  (-certify)
  std::map<Gate*, Polynomial*> ancestors;

  /// list of gates that are parents
  std::list<Gate*> parents;

  /// list of gates that are children
  std::list<Gate*> children;

public:

  /**
      Constructor
      Calls constructor of Var

      @param n_ value, corresponding to aiger value
      @param name_ string name of the variable
      @param level position in order of the variable
  */
  Gate (int n_, std::string name_, int level_, bool input_ = 0, bool output_ = 0);

  /**
      Getter for v

      @return member v
  */
  const Var * get_var() const {return v;};

  /**
      Getter for number of v

      @return integer
  */
  int get_var_num() const {return v->get_num();};

  /**
      Getter for name of v

      @return char*
  */
  const char * get_var_name()const {return v->get_name();};

  /**
      Getter for input

      @return member input
  */
  bool get_input()  const {return input;};

  /**
      Getter for output

      @return member output
  */
  bool get_output() const {return output;};

  /**
      Getter for partial_product

      @return member partial_product
  */
  bool get_pp() const {return partial_product;};

  /**
      Sets partial_product to true
  */
  void mark_pp() {partial_product = 1;};

  /**
      Getter for  aig_output

      @return member  aig_output
  */
  bool get_aig_output() const {return aig_output;};

  /**
      Sets aig_output to true
  */
  void mark_aig_output() {aig_output = 1;};

  /**
      Getter for level

      @return member level
  */
  int get_level() const {return level;};

  /**
      Setter for level

      @param l integer
  */
  void set_level(int l) {level = l;};

  /**
      Getter for xor_chain

      @return member xor_chain
  */
  bool get_xor_chain() const {return xor_chain;};

  /**
      Sets xor_chain to true
  */
  void mark_xor_chain() {xor_chain = 1;};

  /**
      Getter for xor_gate

      @return member xor_gate
  */
  int get_xor_gate() const {return xor_gate;};

  /**
      Setter for xor_gate

      @param val integer
  */
  void set_xor_gate(int val) {xor_gate = val;};

  /**
      Getter for carry_gate

      @return member carry_gate
  */
  int get_carry_gate() const {return carry_gate;};

  /**
      Setter for carry_gate

      @param integer
  */
  void set_carry_gate(int val) {carry_gate = val;};

  /**
      Increases carry_gate
  */
  void inc_carry_gate() {carry_gate++;};

  /**
      Decreases carry_gate
  */
  void dec_carry_gate() {carry_gate--;};

  /**
      Getter for slice

      @return member slice
  */
  int get_slice() const {return slice;};

  /**
      Setter for slice

      @param val integer
  */
  void set_slice(int val) {slice = val;};

  /**
      Increases slice
  */
  void inc_slice() {slice++;};

  /**
      Decreases slice
  */
  void dec_slice() {slice--;};

  /**
      Getter for prop_gen_gate

      @return member prop_gen_gate
  */
  bool get_prop_gen_gate() const {return prop_gen_gate;};

  /**
      Sets prop_gen_gate to true
  */
  void mark_prop_gen_gate() {prop_gen_gate = 1;};

  /**
      Sets prop_gen_gate to false
  */
  void unmark_prop_gen_gate() {prop_gen_gate = 0;};

  /**
      Getter for fsa

      @return member fsa
  */
  bool get_fsa() const {return fsa;};

  /**
      Sets fsa to true
  */
  void mark_fsa() {fsa = 1;};

  /**
      Getter for fsa_inp

      @return member fsa_inp
  */
  int get_fsa_inp() const {return fsa_inp;};

  /**
      Increases fsa_inp
  */
  void inc_fsa_inp() {fsa_inp++;};

  /**
      Sets fsa_inp to 0
  */
  void reset_fsa_inp() {fsa_inp = 0;};

  /**
      Getter for neg

      @return member neg
  */
  bool get_neg() const {return neg;};

  /**
      Setter for neg

      @param val Boolean
  */
  void set_neg(bool val) {neg = val;};

  /**
      Getter for moved

      @return member moved
  */
  bool get_moved() const {return moved;};

  /**
      Sets moved to true
  */
  void mark_moved() {moved = 1;};

  /**
      Getter for bo

      @return member bo
  */
  bool get_bo() const {return bo;};

  /**
      Sets bo to true
  */
  void mark_bo() {bo = 1;};

  /**
      Getter for elim

      @return member elim
  */
  bool get_elim() const {return elim;};

  /**
      Sets elim to true
  */
  void mark_elim(){elim = 1;};

  /**
      Getter for gate_constraint

      @return member gate_constraint
  */
  Polynomial * get_gate_constraint() const {return gate_constraint;};

  /**
      Setter for gate_constraint

      @param p Polynomial *
  */
  void set_gate_constraint(Polynomial * p) {gate_constraint = p;};

  /**
      Prints the gate constraint

      @param file output file
  */
  void print_gate_constraint(FILE * file) const { gate_constraint->print(file);};

  /**
      Getter for co_factor

      @return member co_factor
  */
  Polynomial * get_cofactor() const {return co_factor;};

  /**
      Setter for co_factor

      @param p Polynomial*
  */
  void set_cofactor(Polynomial * p) {co_factor = p;};

  /**
      Getter for begin of ancestors

      @return std::map<Gate*, Polynomial*>::const_iterator
  */
  std::map<Gate*, Polynomial*>::const_iterator anc_begin() const {
    return ancestors.begin();
  };

  /**
      Getter for end of ancestors

      @return std::map<Gate*, Polynomial*>::const_iterator
  */
  std::map<Gate*, Polynomial*>::const_iterator anc_end()   const {
    return ancestors.end();
  };

  /**
      Searches for gate n in ancestors

      @param n Gate*
      @return std::map<Gate*, Polynomial*>::iterator
  */
  std::map<Gate*, Polynomial*>::iterator search_in_anc(Gate *n) {
    return ancestors.find(n);
  };

  /**
      Adds the tuple (n,p) to the ancestors

      @param n Gate*
      @param p Polynomial*
  */
  void set_ancestor (Gate *n, Polynomial * p) {ancestors[n] = p;};

  /**
      Getter for begin of parents

      @return std::list<Gate*>::const_iterator
  */
  std::list<Gate*>::const_iterator parents_begin() const {
    return parents.begin();
  };

  /**
      Getter for end of parents

      @return std::list<Gate*>::const_iterator
  */
  std::list<Gate*>::const_iterator parents_end()   const {
    return parents.end();
  };

  /**
      Getter for size of parents

      @return size_t
  */
  size_t parents_size() const {return parents.size();};

  /**
      Getter for front of parents

      @return Gate*
  */
  Gate * parents_front() const {return parents.front();};

  /**
      Appends gate n to the parents

      @param n Gate*
  */
  void parents_push_back(Gate * n) {parents.push_back(n);};

  /**
      Removes gate n from the parents

      @param n Gate*
  */
  void parents_remove(Gate * n) {parents.remove(n);};

  /**
      Getter for begin of children

      @return std::list<Gate*>::const_iterator
  */
  std::list<Gate*>::const_iterator children_begin() const {
    return children.begin();
  };

  /**
      Getter for end of children

      @return std::list<Gate*>::const_iterator
  */
  std::list<Gate*>::const_iterator children_end()  const {
    return children.end();
  };

  /**
      Getter for size of children

      @return size_t
  */
  size_t children_size() const {return children.size();};

  /**
      Getter for front of children

      @return Gate*
  */
  Gate * children_front() const {return children.front();};

  /**
      Getter for back of children

      @return Gate*
  */
  Gate * children_back() const {return children.back();};

  /**
      Setter for front of children

      @param n Gate*
  */
  void set_children_front(Gate * n) {children.front() = n;};

  /**
      Setter for back of children

      @param n Gate*
  */
  void set_children_back (Gate * n) {children.back() = n;};

  /**
      Appends gate n to the children

      @param n Gate*
  */
  void children_push_back(Gate * n) {children.push_back(n);};

  /**
      Removes gate n from the children

      @param n Gate*
  */
  void children_remove (Gate * n) {children.remove(n);};


  /**
      Determines whether gate constraint is original

      @return True if ancestors are empty, i.e., polynomial has not
      been reduced
  */
  bool orig() const;


  /**
      Destructor
  */
  ~Gate();

  /**
      Determines whether n is contained in the parents of this gate

      @param n Gate*

      @return True whether n is parent of this gate
  */
  bool is_in_parents(const Gate *n) const;

  /**
      Determines whether n is contained in the children of this gate

      @param n Gate*

      @return True whether n is child of this gate
  */
  bool is_child(const Gate *n) const;

};

/*------------------------------------------------------------------------*/
/// Gate ** where all gates are stored
extern Gate ** gates;

/// Counts the number of gates
extern unsigned num_gates;
/*------------------------------------------------------------------------*/

/**
    Returns the gate with aiger value 'lit'

    @param lit unsigned integer

    @returns Gate*
*/
Gate * gate (unsigned lit);

/**
    Allocate the Gate** gates and filling it
*/
void allocate_gates(bool assert = 1);

/**
    Mark those gates that are outputs in the aig
*/
void mark_aig_outputs();

/**
    For the xor output n, derive the corresponding carry output of an
    half-adder

    @param n Gate* that is an xor output

    @return Gate* that is the AND output of an half-adder
*/
Gate * derive_ha_and_gate(const Gate * n);

/**
    Identifies XOR gates in an AIG
*/
void set_xor ();

/**
    Identifies whether the output gates of slice N until NN-2 are
    XOR gates

    @return True when output gates of bigger slices are all xor gates
*/
bool upper_half_xor_output();

/**
    Returns the 'left' child of an XOR

    @param n Gate that is root of XOR

    @return Gate*
*/
Gate * xor_left_child (const Gate * n);

/**
    Returns the 'right' child of an XOR

    @param n Gate that is root of XOR

    @return Gate*
*/
Gate * xor_right_child (const Gate * n);

/**
    Mark all gates in the xor_chain in the last slice
*/
void mark_xor_chain_in_last_slice();

/**
    For all gates set the parents and children

    @param set_children indicate whether children will be set (0 in -substitute)
*/
void set_parents_and_children (bool set_children);


/**
    Defines the polynomial '1-v'

    @param v Var*

    @returns: Polynomial*
*/
Polynomial * negative_poly (const Var * v);

/**
    Defines the polynomial 'v'

    @param v Var*

    @return Polynomial*
*/
Polynomial * positive_poly (const Var * v);

/**
    Initializes all gate constraints in the gates[]
*/
void init_gate_constraints ();

/**
    Deletes Gate** gates by calling destructur of Gate
*/
void delete_gates();

#endif
