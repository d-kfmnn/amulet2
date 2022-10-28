[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

AMulet 2.2 - AIG Multiplier Examination Tool
================================================================================

Our tool AMulet 2.2 is able to verify and certify unsigned and signed 
integer multipliers given as AIGs.

For further information we refer to the paper

Daniela Kaufmann, Armin Biere. 
 [`AMulet 2.0 for Verifying Multiplier Circuits.`](https://danielakaufmann.at/kaufmannbiere-tacas21/)
In Proc. 14th Intl. Conference on Tools and Algorithms for the Construction and Analysis of Systems (TACAS), 8 pages, 2021.

and the corresponding website http://fmv.jku.at/amulet2  
  
----------------------------------------------------------------  
  
Dependencies: `libgmp` (https://gmplib.org/)

Use `./configure.sh && make` to configure and build `AMulet 2.2`.


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

         -p1          expanded proof (no linear combinations, only multiplication and addition)
         -p2          middle condensed proof(sequence of linear combinations, default)
         -p3          condensed proof(one single linear combination)

--------------------------------------------------
28.10.2022 AMulet 2.2:
  - Several bugfixes in the slicing routine, described in our TAP'22 paper [`Fuzzing and Delta Debugging And-Inverter Graph Verification Tools.`](https://danielakaufmann.at/wp-content/uploads/2022/07/TAP_Kaufmann.pdf)

10.09.2021 AMulet 2.1:
  - Reducing memory usage by changing data structure of polynomials

17.02.2021 AMulet 2.1:  
  - Instead of PAC and NSS we now support LPAC on different abstraction levels
    See https://github.com/d-kfmnn/pacheck2 for a corresponding proof checker.
           
  - Optimized polynomial generation 

