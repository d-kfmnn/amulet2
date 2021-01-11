[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

AMulet 2.0 - AIG Multiplier Examination Tool
================================================================================

Our tool AMulet 2.0 is able to verify and certify unsigned and signed 
integer multipliers given as AIGs.

For further information we refer to the paper:

Daniela Kaufmann, Armin Biere. 
 [`AMulet 2.0 for Verifying Multiplier Circuits.`](https://danielakaufmann.at/kaufmannbiere-tacas21/)
To appear in Proc. 14th Intl. Conference on Tools and Algorithms for the Construction and Analysis of Systems (TACAS), 8 pages, 2021.



Dependencies: `libgmp` (https://gmplib.org/)

A doxygen documentation is provided in the folder `doc`, which can be seen
by opening `doc/index.html'

Use `./configure.sh && make` to configure and build `AMulet 2.0`.


usage : `amulet <mode> <input.aig> <output files> [<option> ...]`

Depending on the `<mode>` the `<output files>` and `<options>` have to be set accordingly:


    <mode> = -substitute:
      <output files> =  2 output files need to be passed in the following order 
        <out.cnf>:        miter for checking the equivalence of the substituted adder 
        <out.aig>:        rewritten aiger is stored in this file` 

      <option> = the following options are available 
        -h | --help       print this command line summary 
        -v<1,2,3,4>       different levels of verbosity  (default -v1)
        -signed           option for signed integer multipliers 


    <mode> = -verify:
      <output files> =  no output files are required 
      
      <option> = the following options are available 
         -h | --help           print this command line summary 
         -v<1,2,3,4>           different levels of verbosity (default -v1) 
         -signed               option for signed integer multipliers 
         -no-counter-examples  do not generate and write counter examples
     
     
    <mode> = -certify:
      <output files> =  3 output files need to be passed in the following order
        <out.polys>:      initial polynomial set 
        <out.proof>:      proof rules (depending whether PAC proof or NS proof is generated) 
        <out.spec> :      spec which should be checked 
      
      <option> = the following options are available 
         -h | --help           print this command line summary 
         -v<1,2,3,4>           different levels of verbosity (default -v1) 
         -signed               option for signed integer multipliers 
         -no-counter-examples  do not generate and write counter examples

         -pac                  produces proofs in condensed PAC format 
         -nss                  produces proofs in Nullstellensatz format (default)
