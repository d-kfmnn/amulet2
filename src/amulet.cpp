/*------------------------------------------------------------------------*/
/*! \file amulet.cpp
    \brief main file of our tool AMulet2

  Part of AMulet2.0 : AIG Multiplier Verification Tool.
  Copyright (C) 2020 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#define VERSION "2.0"
/*------------------------------------------------------------------------*/
/// Manual of AMulet2, will be printed with command line '-h'
static const char * USAGE =
"\n"
"### USAGE ###\n"
"usage : amulet <mode> <input.aig> <output files> [<option> ...] \n"
"\n"
"Depending on the <mode> the <output files> and <options> have to be set accordingly:"
"\n"
"\n"
"<mode> = -substitute:\n"
"    <output files> =  2 output files need to be passed in the following order \n"
"      <out.cnf>:        miter for checking the equivalence of the substituted adder \n"
"      <out.aig>:        rewritten aiger is stored in this file \n"
"\n"
"    <option> = the following options are available \n"
"      -h | --help       print this command line summary \n"
"      -v<1,2,3,4>       different levels of verbosity (default -v1) \n"
"      -signed           option for signed integer multipliers \n"
"\n"
"\n"
"<mode> = -verify:\n"
"    <output files> =  no output files are required \n"
"     \n "
"    <option> = the following options are available \n"
"       -h | --help           print this command line summary \n"
"       -v<1,2,3,4>           different levels of verbosity (default -v1) \n"
"       -signed               option for signed integer multipliers \n"
"       -no-counter-examples  do not generate and write counter examples\n"
"     \n"
"     \n"
"<mode> = -certify:\n"
"    <output files> =  3 output files need to be passed in the following order\n"
"      <out.polys>:      initial polynomial set \n"
"      <out.proof>:      proof rules (depending whether PAC proof or NS proof is generated) \n"
"      <out.spec> :      spec which should be checked \n"
"     \n "
"    <option> = the following options are available \n"
"       -h | --help           print this command line summary \n"
"       -v<1,2,3,4>           different levels of verbosity (default -v1) \n"
"       -signed               option for signed integer multipliers \n"
"       -no-counter-examples  do not generate and write counter examples\n"
"\n"
"       -pac                  produces proofs in condensed PAC format \n"
"       -nss                  produces proofs in Nullstellensatz format (default)\n"



;
/*------------------------------------------------------------------------*/
#include "parser.h"
#include "substitution_engine.h"
#include "polynomial_solver.h"
/*------------------------------------------------------------------------*/
/// Name of the input file
static const char * input_name = 0;

/// \brief
/// Name of first output file, which stores the CNF miter in '-substitute', and
/// the gate constraints in '-certify'.
static const char * output_name1 = 0;

/// \brief
/// Name of second output file, which stores the rewritten AIG in '-substitute',
/// and the core proof in '-certify'.
static const char * output_name2 = 0;

/// Name of third output file. Stores the specification in '-certify'.
static const char * output_name3 = 0;

/// Selected mode, '-substitute' = 1, '-verify' = 2, '-certify' = 3
static int mode;
/*------------------------------------------------------------------------*/
/**
    Calls the deallocaters of the involved data types
    @see reset_all_signal_handlers()
    @see delete_gates()
    @see deallocate_terms()
    @see deallocate_mstack()
    @see clear_mpz()
*/
static void reset_all(){
  reset_all_signal_handlers();
  delete_gates();
  deallocate_terms();
  deallocate_mstack ();
  clear_mpz();
}
/*------------------------------------------------------------------------*/
/**
    Main Function of AMulet2.
    Reads the given AIG and depending on the selected mode, either
    calls the substution engine or the polynomial solver.

    Prints statistics to stdout after finishing.
*/
int main (int argc, char ** argv) {
  msg ("AMulet Version " VERSION);
  msg ("Aiger multiplier examination tool");
  msg ("Copyright (C) 2020, Daniela Kaufmann, Johannes Kepler University Linz");

  for (int i = 1; i < argc; i++) {
    if (!strcmp (argv[i], "-h") ||
    !strcmp (argv[i], "--help")) {
      fputs (USAGE, stdout);
      fflush (stdout);
      exit (0);
    }
    else if (!strcmp (argv[i], "-v0")) verbose = 0;
    else if (!strcmp (argv[i], "-v1")) verbose = 1;
    else if (!strcmp (argv[i], "-v2")) verbose = 2;
    else if (!strcmp (argv[i], "-v3")) verbose = 3;
    else if (!strcmp (argv[i], "-v4")) verbose = 4;
    else if (!strcmp (argv[i], "-substitute")) {
      if(!mode) {
        msg("selected mode: adder substitution");
        mode = 1;
      } else die("mode has alreday been selected (try '-h')");
    } else if (!strcmp (argv[i], "-verify")) {
      if(!mode) {
        msg("selected mode: verification");
        mode = 2;
      } else die("mode has alreday been selected (try '-h')");
    } else if (!strcmp (argv[i], "-certify")) {
      if(!mode) {
        msg("selected mode: verification + certificates");
        mode = 3;
      } else die("mode has alreday been selected (try '-h')");
    } else if (!strcmp (argv[i], "-pac")) {
      pac = 1;
      if(nss) die("too many proof formats selected (try '-h')");
    } else if (!strcmp (argv[i], "-nss")) {
      nss = 1;
      if(pac) die("too many proof formats selected (try '-h')");
    } else if (!strcmp (argv[i], "-signed")) {
      signed_mult = 1;
    } else if (!strcmp (argv[i], "-no-counter-examples")) {
      gen_witness = 0;
    } else if (output_name3)
    die ("too many arguments '%s', '%s', '%s', '%s' and '%s' (try '-h')",
    input_name, output_name1, output_name2, output_name3, argv[i]);
    else if (output_name2) output_name3 = argv[i];
    else if (output_name1) output_name2 = argv[i];
    else if (input_name) output_name1 = argv[i];
    else input_name = argv[i];
  }

  if (!mode) die("select mode (try -h for more information)");
  if (!input_name)  die("no input file given (try '-h')");
  if (mode == 1){
    if (output_name3) die ("too many arguments '%s', '%s', '%s' and '%s' (try '-h')",
    input_name, output_name1, output_name2, output_name3);
    if (!output_name2) die ("too few arguments (try '-h')");
    if (pac) msg("option -pac is only possible in -certify (will be ignored)");
    if (nss) msg("option -nss is only possible in -certify (will be ignored)");
    pac = 0; nss = 0;
  } else if (mode == 2) {
    if (output_name1) die ("too many arguments (try '-h')");
    if (pac) msg("option -pac is only possible in -certify (will be ignored)");
    if (nss) msg("option -nss is only possible in -certify (will be ignored)");
    pac = 0; nss = 0;
  } else if (mode == 3){
    if (!output_name3) die ("too few arguments (try '-h')");
    if(!pac) nss = 1;
    if(nss) msg("proof format: Nullstellensatz");
    else msg("proof format: PAC");
  }


  init_all_signal_handers ();
  init_nonces ();

  parse_aig(input_name);

  if (mode == 1) {
    init_gate_substitution();
    substitution(output_name1, output_name2);
  }
  else {
    init_gates_verify();
    verify(input_name, output_name1, output_name2, output_name3, mode == 3);
  }
  reset_aig_parsing();
  reset_all();

  reset_time = process_time();
  print_statistics(mode);

  return 0;
}
