/*------------------------------------------------------------------------*/
/*! \file monomial.h
    \brief contains the class Monomial and further functions to
    manipulate monomials

  Part of AMulet2.1 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#ifndef AMULET2_SRC_MONOMIAL_H_
#define AMULET2_SRC_MONOMIAL_H_
/*------------------------------------------------------------------------*/
#include <gmp.h>

#include "term.h"
/*------------------------------------------------------------------------*/

/** \class Monomial
    This class is used to represent monomials in a polynomial.
    A monomial consist of a coefficient and a term.
*/

class Monomial {
  // / Term*
  Term * term;

  // / reference counter
  unsigned ref;

 public:
  // / Coefficient
  mpz_t coeff;

  /** Constructor

      @param c mpz_t coefficient
      @param t Term*
  */
  Monomial(mpz_t _c, Term * _t);

  /** Getter for member term

      @return Term*
  */
  Term * get_term()       const {return term;}

  /** Getter for member term, calls copy routine of Term

      @return a copy of Term* term
  */
  Term * get_term_copy()  const {return term->copy();}

  /** Returns the size fo the term

      @return unsigned, the size of the term
  */
  unsigned get_term_size() const {return term->size();}

  /** Decreases the reference counter

      @return the decreases reference counter
  */
  unsigned dec_ref()             {return --ref;}

  /** Getter for the reference counter

      @return member ref
  */
  unsigned get_ref()       const {return ref;}

  /**
  Copy routine

  @return Same monomial with increased reference counter
  */
  Monomial * copy();

  /**
  Printing routine

  @param file Output file
  @param lm bool indicating whether we print a leading monomial
  */
  void print(FILE * file, bool lm = 0) const;

  /** Destructor */
  ~Monomial();
};
/*------------------------------------------------------------------------*/

/**
  Multiplies two monomials

  @param m1 Monomial*
  @param m2 Monomial*

  @return Product of m1*m2
*/
Monomial * multiply_monomial(const Monomial * m1, const Monomial *m2);

/**
  Wrapper for deconstructor, reduces the references of m until 0.
  Then deconstructor is called.

  @param m Monomial*
*/
void deallocate_monomial(Monomial * m);

#endif  // AMULET2_SRC_MONOMIAL_H_
