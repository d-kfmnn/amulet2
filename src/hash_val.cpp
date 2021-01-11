/*------------------------------------------------------------------------*/
/*! \file hash_val.cpp
    \brief contains functions used to compute hash values for variables

  Part of AMulet2 : AIG Multiplier Verification Tool.
  Copyright (C) 2020 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/

#include "hash_val.h"
/*------------------------------------------------------------------------*/
// Local variables

/// array used to hold 64-bit random numbers
uint64_t nonces[32];

/// number of random numbers
const size_t num_nonces = sizeof nonces / sizeof nonces[0];
/*------------------------------------------------------------------------*/

uint16_t rand16 () {
  int tmp = rand ();
  assert (tmp >= 0);
  uint16_t res = tmp ^ (tmp >> 16);
  return res;
}

/*------------------------------------------------------------------------*/

uint64_t rand64 () {
  uint64_t res = 0;
  for (unsigned i = 0; i < 64; i += 16)
    res |= rand16 () << i;
  return res;
}

/*------------------------------------------------------------------------*/

void init_nonces () {
  srand (42);
  for (size_t i = 0; i < num_nonces; i++)
    nonces[i] = rand64 () | 1;
}

/*------------------------------------------------------------------------*/

uint64_t get_nonces_entry (size_t index){
  assert (index < 32);
  return nonces[index];
}

/*------------------------------------------------------------------------*/

uint64_t hash_string (const std::string& str) {
  uint64_t res = 0;
  size_t i = 0;
  for(std::string::const_iterator it = str.cbegin(); it != str.cend(); ++it){
    res += *it;
    res *= nonces[i++];
    if (i == num_nonces) i = 0;
  }
  return res;
}
