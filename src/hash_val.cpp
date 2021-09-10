/*------------------------------------------------------------------------*/
/*! \file hash_val.cpp
    \brief contains functions used to compute hash values for variables

  Part of AMulet2.1 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include <random>

#include "hash_val.h"
/*------------------------------------------------------------------------*/
// Local variables

// / array used to hold 64-bit random numbers
uint64_t nonces[32];

// / number of random numbers
const size_t num_nonces = sizeof nonces / sizeof nonces[0];

std::random_device rd;
std::mt19937_64 gen(42);
std::uniform_int_distribution<unsigned long long> dis;
/*------------------------------------------------------------------------*/


void init_nonces() {
  for (size_t i = 0; i < num_nonces; i++)
    nonces[i] = dis(gen) | 1;
}

/*------------------------------------------------------------------------*/

uint64_t get_nonces_entry(size_t index) {
  assert(index < 32);
  return nonces[index];
}

/*------------------------------------------------------------------------*/

uint64_t hash_string(const std::string& str) {
  uint64_t res = 0;
  size_t i = 0;
  for (std::string::const_iterator it = str.cbegin(); it != str.cend(); ++it) {
    res += *it;
    res *= nonces[i++];
    if (i == num_nonces) i = 0;
  }
  return res;
}
