/*------------------------------------------------------------------------*/
/*! \file aig.cpp
    \brief contains functions to parse and manipulate the input AIG

  Part of AMulet2.1 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include "aig.h"
/*------------------------------------------------------------------------*/
// Global Variables
unsigned M, NN;
unsigned a0, al, ainc;
unsigned b0, bl, binc;

aiger * rewritten;
aiger * miter;
/*------------------------------------------------------------------------*/
// Local Variables

aiger * model;   // /< aiger* object, used for storing the given AIG graph
/*------------------------------------------------------------------------*/
void init_aig_substitution() {
  rewritten = aiger_init();
  miter = aiger_init();
}
/*------------------------------------------------------------------------*/
void reset_aig_substitution() {
  aiger_reset(rewritten);
  aiger_reset(miter);
}
/*------------------------------------------------------------------------*/
void init_aig_parsing() {
  model = aiger_init();
}
/*------------------------------------------------------------------------*/
void reset_aig_parsing() {
  aiger_reset(model);
}
/*------------------------------------------------------------------------*/
const char * aiger_open_and_read_to_model(const char * input_name) {
  return aiger_open_and_read_from_file(model, input_name);
}
/*------------------------------------------------------------------------*/

bool is_model_input(unsigned val) {
  return aiger_is_input(model, val);
}
/*------------------------------------------------------------------------*/
aiger_and * is_model_and(unsigned val) {
  return aiger_is_and(model, val);
}
/*------------------------------------------------------------------------*/
unsigned get_model_num_latches() {return model->num_latches;}
unsigned get_model_num_inputs() {return model->num_inputs;}
unsigned get_model_num_ands()   {return model->num_ands;}
unsigned get_model_num_outputs() {return model->num_outputs;}
unsigned get_model_maxvar()     {return model->maxvar;}
/*------------------------------------------------------------------------*/
unsigned get_model_inputs_lit(unsigned i) {
  assert(i < NN);
  return model->inputs[i].lit;
}
/*------------------------------------------------------------------------*/
const char* get_model_inputs_name(unsigned i) {
  assert(i < NN);
  return model->inputs[i].name;
}
/*------------------------------------------------------------------------*/
unsigned slit(unsigned i) {
  if (!model) return -1;
  assert(i < NN);
  return model->outputs[i].lit;
}
/*------------------------------------------------------------------------*/
int write_model(FILE *file) {
  return aiger_write_to_file(model, aiger_binary_mode, file);
}
/*------------------------------------------------------------------------*/
