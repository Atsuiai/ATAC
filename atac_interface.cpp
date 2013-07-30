#include "atac_interface.h"
//#include "../system_control_t.h"

min_output_t minimize(min_input_t input)
{
    min_output_t output;

    output.success = true;
    output.positions = input.positions;

    output.total_potential = 0;
    output.average_potential = 0;
    output.maximum_force = 0;
    output.last_delta_e = 0;

    output.iterations = 0;
    output.function_evalutations = 0;
    output.time_taken = 0;

    return output;
}
