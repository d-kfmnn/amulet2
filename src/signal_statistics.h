/*------------------------------------------------------------------------*/
/*! \file signal_statistics.h
    \brief used to handle signals, messages and statistics

  Part of AMulet2.0 : AIG Multiplier Verification Tool.
  Copyright (C) 2020 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#ifndef _signal_statistics_H
#define _signal_statistics_H
/*------------------------------------------------------------------------*/
#include <stdarg.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <iostream>
/*------------------------------------------------------------------------*/
extern void (*original_SIGINT_handler)(int);
extern void (*original_SIGSEGV_handler)(int);
extern void (*original_SIGABRT_handler)(int);
extern void (*original_SIGTERM_handler)(int);

/**
    Returns name of signal

    @param sig integer

    @return char *
*/
const char * signal_name (int sig);

/**
    Initialize all signals
*/
void init_all_signal_handers ();


/**
    Catch signal and prints corresponding message

    @param sig integer
*/
void catch_signal (int sig);

/**
    Resets all signal handlers
*/
void reset_all_signal_handlers ();

/*------------------------------------------------------------------------*/
/// Level of output verbosity, ranges from 0 to 4
extern int verbose;

/**
    Prints an error message to stderr and exits the program

    @param char* fmt message
*/
void die (const char *fmt, ...);

/**
    Prints a message to stdout

    @param char* fmt message
*/
void msg (const char *fmt, ...);
/*------------------------------------------------------------------------*/

/// Time measures used for verify/certify modus
extern double init_time;         ///< measure for initializing time
extern double slicing_elim_time; ///< measure for time used to eliminate & slice
extern double reduction_time;    ///< measure for time used to reduce
extern double reset_time;        ///< measure for resetting time
extern double substitution_time; ///< measure for time used in substitution

/**
    Determines max used memory
*/
size_t maximum_resident_set_size ();

/**
    Determines the used process time
*/
double process_time ();


/**
    Print statistics of maximum memory and used process time depending on
    selected modus

    @param modus integer
*/
void print_statistics (int modus);




#endif
