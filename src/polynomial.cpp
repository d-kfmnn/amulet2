/*------------------------------------------------------------------------*/
/*! \file polynomial.cpp
    \brief contains arithmetic operations for polynomials

  Part of AMulet2.0 : AIG Multiplier Verification Tool.
  Copyright (C) 2020 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include "polynomial.h"
/*------------------------------------------------------------------------*/

Polynomial::Polynomial() {}
Polynomial::Polynomial(int _idx):  idx(_idx) {}

/*------------------------------------------------------------------------*/

Polynomial * Polynomial::copy() const{
  Polynomial * out = new Polynomial(idx);
  for (std::deque<Monomial *>::const_iterator it = mon_begin();
       it != mon_end(); ++it){

    Monomial * m = *it;
    out->mon.push_back(m->copy());
  }
  return out;
}

/*------------------------------------------------------------------------*/

void Polynomial::print(FILE * file, bool end) const{
  if(mon.empty()) fputs_unlocked("0",file);

  else {
    for (std::deque<Monomial *>::const_iterator it = mon_begin();
         it != mon_end(); ++it){

      Monomial * m = *it;
      if(it == mon_begin()) m->print(file, 1);
      else m->print(file, 0);
    }
  }
  if(end) fputs(";\n",file);
}

/*------------------------------------------------------------------------*/

Polynomial::~Polynomial(){
  for (std::deque<Monomial *>::const_iterator it = mon_begin();
      it != mon_end(); ++it)

    deallocate_monomial(*it);
}

/*------------------------------------------------------------------------*/

bool Polynomial::is_constant_zero_poly() const{
  return mon.empty();
}

/*------------------------------------------------------------------------*/

bool Polynomial::is_constant_one_poly() const{
  if (mon.size() != 1) return 0;

  Monomial * m = mon.front();
  if(m->get_term()) return 0;
  if(mpz_cmp_si(m->coeff, 1) != 0) return 0;

  return 1;
}

/*------------------------------------------------------------------------*/

int Polynomial::min_term_size() const{
    int len = INT_MAX;

    for (std::deque<Monomial *>::const_iterator it = mon_begin();
      it != mon_end(); ++it){

      Monomial * m = *it;

      int tlen = 0;
      if(m->get_term()) tlen=m->get_term_size();
      if(tlen < len) len = tlen;
    }
    return len;
}

/*------------------------------------------------------------------------*/
// Local variables
static size_t size_mstack; ///< size of mstack
static size_t num_mstack = 0; ///< number of elements in mstack
static Monomial ** mstack; ///< Monomial** used for building poly
/*------------------------------------------------------------------------*/

void enlarge_mstack () {
  size_t new_size_mstack = size_mstack ? 2*size_mstack : 1;

  Monomial** newArr = new Monomial*[new_size_mstack];
  memcpy( newArr, mstack, size_mstack * sizeof(Monomial*) );
  delete [] mstack;
  mstack = newArr;
  size_mstack = new_size_mstack;
}

/*------------------------------------------------------------------------*/

void clear_mstack () { num_mstack = 0; }

/*------------------------------------------------------------------------*/

void deallocate_mstack () { delete[](mstack); }

/*------------------------------------------------------------------------*/

bool mstack_is_empty () { return num_mstack == 0;}
/*------------------------------------------------------------------------*/
void push_mstack_end(Monomial *m){

  if(size_mstack == num_mstack) enlarge_mstack();

  assert(m);
  if (mpz_sgn(m->coeff) == 0) {
    deallocate_monomial(m);
    return;
  }

  mstack[num_mstack++] = m;
}

/*------------------------------------------------------------------------*/

void push_mstack(Monomial *m){

  assert(m);
  if (mpz_sgn(m->coeff) == 0) {
    deallocate_monomial(m);
    return;
  }

  if (size_mstack == num_mstack) enlarge_mstack ();
  if(num_mstack == 0) {
    mstack[num_mstack++] = m;
    return;
  }

  if (!m->get_term()){
    Monomial * tmp = mstack[num_mstack-1];


    if(tmp->get_term()) mstack[num_mstack++] = m;
    else {
      mpz_t coeff;
      mpz_init(coeff);
      mpz_add(coeff, tmp->coeff, m->coeff);
      deallocate_monomial(m);
      deallocate_monomial(tmp);

      if(mpz_sgn(coeff) != 0){
        mstack[num_mstack-1] = new Monomial(coeff, 0);
      } else --num_mstack;
      mpz_clear(coeff);
    }
  }
  else {
    {
      assert(num_mstack>0);
      int i = num_mstack-1;
      int cmp = -1;
      Monomial * tmp = 0;

      while (i >= 0) {
        tmp = mstack[i];
        cmp = tmp->get_term()->cmp(m->get_term());

        if (cmp >= 0) break;
        i--;
      }

      if(cmp == 0){
        mpz_t coeff;
        mpz_init(coeff);
        mpz_add(coeff, tmp->coeff, m->coeff);


        if(mpz_sgn(coeff) == 0){
          for(unsigned j = i; j < num_mstack-1; j++)
          mstack[j] = mstack[j+1];

          num_mstack--;
        } else {
          mstack[i] = new Monomial(coeff, m->get_term_copy());
        }
        deallocate_monomial(m);
        deallocate_monomial(tmp);
        mpz_clear(coeff);
      }
      else{
        for(int j = num_mstack; j>i+1; j--)
        mstack[j] = mstack[j-1];


        mstack[i+1] = m;
        num_mstack++;
      }
    }
  }
}

/*------------------------------------------------------------------------*/

Polynomial * build_poly(){
   Polynomial * res = new Polynomial ();
   for(size_t i = 0; i< num_mstack; i++){
     res->mon_push_back(mstack[i]);
   }

   clear_mstack();
   return res;
}

/*------------------------------------------------------------------------*/

Polynomial * add_poly(const Polynomial * p1, const Polynomial * p2){
  assert(p1);
  assert(p2);

  std::deque<Monomial*>::const_iterator it1=p1->mon_begin()++;
  std::deque<Monomial*>::const_iterator it2=p2->mon_begin()++;
  Monomial * m1 = *it1;
  Monomial * m2 = *it2;
  mpz_t coeff;
  mpz_init(coeff);
  while(it1 != p1->mon_end() && it2 != p2->mon_end()){

    if(!m1->get_term() || !m2->get_term()){
      if(!m1->get_term() && !m2->get_term()) {

        mpz_add(coeff, m1->coeff, m2->coeff);
        if(mpz_sgn(coeff) != 0){
          Monomial * m = new Monomial(coeff, 0);
          push_mstack_end (m);
        }
        ++it1;
        ++it2;
        m1 = *it1;
        m2 = *it2;
      }
      else if (!m1->get_term()) {
        push_mstack_end(m2->copy());
        ++it2;
        m2 = *it2;
      }
      else {
        push_mstack_end(m1->copy());
        ++it1;
        m1 = *it1;
      }
    } else {
      int cmp = m1->get_term()->cmp(m2->get_term());
      if(cmp == 1){
        push_mstack_end(m1->copy());
        ++it1;
        m1 = *it1;
      }
      else if(cmp == -1){
        push_mstack_end(m2->copy());
        ++it2;
        m2 = *it2;
      }
      else{
        mpz_add(coeff, m1->coeff, m2->coeff);
        if(mpz_sgn(coeff) != 0){
          Monomial * m = new Monomial(coeff, m1->get_term_copy());
          push_mstack_end(m);
        }
        ++it1;
        ++it2;
        m1 = *it1;
        m2 = *it2;
      }
    }
  }
  mpz_clear(coeff);

  while(it1 != p1->mon_end()){
    push_mstack_end(m1->copy());
    ++it1;
    m1 = *it1;
  }
  while(it2 != p2->mon_end()){
    push_mstack_end(m2->copy());
    ++it2;
    m2 = *it2;
  }

  Polynomial * p = build_poly();
  return p;
}


/*------------------------------------------------------------------------*/

Polynomial * multiply_poly(const Polynomial * p1, const Polynomial * p2){
  mpz_t coeff;
  mpz_init(coeff);
  Term * t;
  for(std::deque<Monomial*>::const_iterator it1 = p1->mon_begin();
      it1 != p1->mon_end(); ++it1){

      Monomial * m1 = *it1;
    for(std::deque<Monomial*>::const_iterator it2 = p2->mon_begin();
        it2 != p2->mon_end(); ++it2){

      Monomial * m2 = *it2;
      if(mpz_cmp_si (m1->coeff, -1) == 0 ) mpz_neg(coeff, m2->coeff);
      else if (mpz_cmp_si (m2->coeff, -1) == 0 ) mpz_neg(coeff, m1->coeff);
      else mpz_mul(coeff, m1->coeff, m2->coeff);

      if(m1->get_term() && m2->get_term())
        t = multiply_term(m1->get_term(), m2->get_term());
      else if (m2->get_term())
        t = m2->get_term_copy();
      else if (m1->get_term())
        t = m1->get_term_copy();
      else t = 0;
      push_mstack(new Monomial(coeff, t));
    }
  }
  Polynomial * p = build_poly();
  mpz_clear(coeff);
  return p;
}

/*------------------------------------------------------------------------*/

Polynomial * multiply_poly_with_constant(const Polynomial *p1, mpz_t c){
  if (mpz_sgn(c) == 0) return 0;
  mpz_t coeff;
  mpz_init(coeff);

  for(std::deque<Monomial*>::const_iterator it = p1->mon_begin();
      it != p1->mon_end(); ++it){

    Monomial * m = *it;
    mpz_mul(coeff, m->coeff, c);
    if(m->get_term())
      push_mstack_end(new Monomial(coeff, m->get_term_copy()));
    else push_mstack_end(new Monomial(coeff, 0));

  }
  Polynomial * tmp = build_poly();
  mpz_clear(coeff);
  return tmp;
}

/*------------------------------------------------------------------------*/

Polynomial * divide_by_term(const Polynomial * p1, const Term * t){

  assert(t->size() == 1);
  const Var * v = t->get_var();


  for(std::deque<Monomial*>::const_iterator it = p1->mon_begin();
      it != p1->mon_end(); ++it){

    Monomial * lm_tmp = *it;

    if(!lm_tmp->get_term()) continue;
    if(lm_tmp->get_term()->cmp(t) == -1) break;

    if(lm_tmp->get_term()->contains(v)){
      Term * t_rem = remainder(lm_tmp->get_term(), v);
      if(t_rem) push_mstack_end(new Monomial(lm_tmp->coeff, t_rem->copy()));
      else{
        push_mstack_end(new Monomial(lm_tmp->coeff, 0));
        break;
      }
    }

  }
  Polynomial * p = build_poly();
  return p;
}

/*------------------------------------------------------------------------*/

mpz_t one;
mpz_t minus_one;
mpz_t base;
mpz_t mod_coeff;

/*------------------------------------------------------------------------*/

void init_mpz(unsigned exp){
  mpz_init_set_ui(one,1);
  mpz_init_set_si(minus_one, -1);
  mpz_init_set_ui(base,2);

  mpz_init(mod_coeff);
  mpz_pow_ui(mod_coeff, base, exp);
}

/*------------------------------------------------------------------------*/

void clear_mpz(){
  mpz_clear(one);
  mpz_clear(minus_one);
  mpz_clear(base);
  mpz_clear(mod_coeff);
}

/*------------------------------------------------------------------------*/
