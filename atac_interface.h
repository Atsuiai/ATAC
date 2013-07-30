#ifndef ATAC_INTERFACE_H
#define ATAC_INTERFACE_H

#include <armadillo>

using namespace std;

///////////////
// minimization input structure
///////////////

enum min_method_t
{
    NAMD,
    CT_OPT,
    CT_OPT_CUDA, // not implemented quite yet
    LAMMPS       // may or may not implement this
};

struct min_input_t
{
    min_method_t method;        // minimization method

    arma::mat positions;        // input positions
    arma::imat bondList;        // input bond list

    int interation_limit;       // hard iteration limit, 0 for no limit
    double max_force_cutoff;    // maximum force exit condition, 0 to disable

    // periodicity, 0 for non-periodic
    double x_period,
           y_period,
           z_period;
};

///////////////
// minimization output structure
///////////////

struct min_output_t
{
    // success boolean
    //   false if some atoms are not within cutoff/bond distance
    //   or if the linesearch fails
    bool success;

    arma::mat positions;        // the final position information either after convergence or before failure

    double total_potential,     // total system potential of all atoms
           average_potential,   // average potential per atom
           maximum_force,       // maximum force experienced by a single atom
           last_delta_e;        // the change in the potential between the last two iterations

    int iterations,             // iterations taken to converge
        function_evalutations;  // number of function evaluations
    double time_taken;          // seconds
};

///////////////
// minimization function
///////////////

min_output_t minimize(min_input_t input);


#endif // ATAC_INTERFACE_H
