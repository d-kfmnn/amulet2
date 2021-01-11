/*------------------------------------------------------------------------*/
/*! \file hash_val.h
    \brief contains functions used to compute hash values for variables

  Part of AMulet2.0 : AIG Multiplier Verification Tool.
  Copyright (C) 2020 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#ifndef _hash_val_H
#define _hash_val_H
/*------------------------------------------------------------------------*/
#include <assert.h>
#include <string>
/*------------------------------------------------------------------------*/
/**
    Generates a 16 bit random number

    @return a 16-bit random number
*/
uint16_t rand16 ();

/**
    Generates a 64 bit random number

    @return a 64-bit random number, based on rand16()
*/
uint64_t rand64 ();


/**
    Fills a 32-bit array with 64-bit random numbers
*/
void init_nonces ();

/**
    Returns the 64-bit random number in our array of random numbers

    @param index size_t

    @return returns the random number stored in nonces[index]
*/
uint64_t get_nonces_entry (size_t index);

/**
    Computes the hash value for the given string

    @param str a const std::string

    @return a uint64_t computed hash value for the input string
*/
uint64_t hash_string (const std::string& str);

#endif
