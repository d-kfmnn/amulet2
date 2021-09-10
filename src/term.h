/*------------------------------------------------------------------------*/
/*! \file term.h
    \brief contains the class Term and further functions to
    manipulate terms

  Part of AMulet2.1 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#ifndef AMULET2_SRC_TERM_H_
#define AMULET2_SRC_TERM_H_
/*------------------------------------------------------------------------*/
#include <stack>

#include "variable.h"
/*------------------------------------------------------------------------*/

/** \class Term
    This class is used to represent terms in a polynomial.
    Terms are represented as ordered linked lists of variables.
*/

class Term {
  // / head variable
  const Var * variable;

  // / tail in linked list
  Term * rest;

  // / reference counter
  uint64_t ref;

  // / hash value
  const uint64_t hash;

  // / hash collision chain link
  Term * next;

 public:
  /** Constructor

      @param _v Var*
      @param _r Term*
      @param _hash uint64_t
      @param _n Term*
  */
  Term(const Var * _v, Term * _r, uint64_t _hash, Term * _n);

  /** Getter for member variable

      @return Var*
  */
  const Var * get_var() const {return variable;}

  /** Getter for level of variable

      @return integer
  */
  int get_var_level() const {return variable->get_level();}

  /** Getter for num of variable

      @return integer
  */
  int get_var_num() const {return variable->get_num();}

  /** Getter for name of variable

      @return char*
  */
  const char * get_var_name() const {return variable->get_name();}

  /** Getter for member rest

      @return Term*
  */
  Term * get_rest() const {return rest;}

  /** Getter for member hash

      @return uint64_t
  */
  uint64_t get_hash() const {return hash;}

  /** Getter for member next

      @return Term*
  */
  Term * get_next() const {return next;}

  /** Setter for member term

      @param t Term*
  */
  void set_next(Term * t) {next = t;}

  /** Getter for member ref

      @return uint64_t
  */
  uint64_t get_ref() const {return ref;}

  /** Increases ref

      @return uint64_t
  */
  uint64_t inc_ref() {return ++ref;}

  /** Decreases ref

      @return uint64_t
  */
  uint64_t dec_ref() {return --ref;}

  /**
      Copy routine

      @return A copy of the current term
  */
  Term * copy();

  /**
      Printing routine

      @param file Output file
  */
  void print(FILE * file) const;

  /**
      Returns the number of variables in a term

      @return unsigned integer
  */
  unsigned size() const;

  /**
      Compares this term to term t using the levels of the variables

      @param t Term*

      @return +1 if this > t
              -1 if this < t
              0  if this = t
  */
  int cmp(const Term *t) const;

  /**
      Checks whether v is contained in Term

      @param v Var*

      @return true if v is contained in term
  */
  bool contains(const Var * v) const;
};

/*------------------------------------------------------------------------*/
// We organize terms in a hash table that is dynamically enlarged.
// Every time a new term is defined, we compute a hash value and insert
// the term. Terms are counted using a reference counter, which is incremented
// and decremented depending how often the term occurs in polynomials.

/**
    Compute hash_values

    @param variable Variable*
    @param rest Term*

    @return computed hash value for the term(variable, rest)
*/
uint64_t compute_hash_term(const Var * variable, const Term * rest);

/**
    Enlarges the hash table
*/
void enlarge_terms();

/**
    Builds a term, where variable is added at the front of rest

    @param variable Variable*
    @param rest     Term*

    @return Term*
*/
Term * new_term(const Var * variable, Term * rest);



/**
    Decrements the reference count of a term, and actually deletes a
    term if its reference count goes to zero.  In this case it also
    removes it from the hash table and applies the same procedure to the
    suffix 'rest'.

    @param t Term*
*/
void deallocate_term(Term * t);


/**
    Deallocates the hash table "term_table"
*/
void deallocate_terms();

/*------------------------------------------------------------------------*/
// Terms are generated using a stack "vstack"

/**
    Pushes a variable to the variable stack, thus variables have to be pushed
    in order

    @param v Var*
*/
void add_to_vstack(const Var* v);

/**
    Generates a term from the variable stack

    @return Term* generated from the variable stack
*/
Term * build_term_from_stack();

/**
    Multiplies two terms

    @param t1 Term*
    @param t2 Term*

    @return Term* which is product of t1*t2
*/
Term * multiply_term(Term * t1, Term * t2);

/**
    Remainder of t divided by v

    @param t Term*
    @param v Var*

    @return Term* where v is removed from t
*/
Term * remainder(const Term * t, const Var * v);

#endif  // AMULET2_SRC_TERM_H_
