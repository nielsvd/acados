/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "acados/ocp_nlp/ocp_nlp_sm_dopus.h"

#include <assert.h>
#include <stdlib.h>

#include "acados/sim/sim_common.h"
#include "acados/sim/sim_rk_common.h"
#include "acados/utils/math.h"
#include "acados/utils/types.h"

ocp_nlp_sm_dopus_args *ocp_nlp_sm_dopus_create_arguments() {

}

int_t ocp_nlp_sm_dopus_calculate_memory_size(const ocp_nlp_sm_in *sm_in,
                                             void *args_) {

}

char *ocp_nlp_sm_dopus_assign_memory(const ocp_nlp_sm_in *sm_in, void *args_,
                                     void **mem_, void *raw_memory) {

}

ocp_nlp_sm_dopus_memory *ocp_nlp_sm_dopus_create_memory(
    const ocp_nlp_sm_in *sm_in, void *args_) {

}

int_t ocp_nlp_sm_dopus_calculate_workspace_size(const ocp_nlp_sm_in *sm_in,
                                                void *args_) {

}

char *ocp_nlp_sm_dopus_assign_workspace(const ocp_nlp_sm_in *sm_in, void *args_,
                                        void **work_, void *raw_memory) {

}

ocp_nlp_sm_dopus_workspace *ocp_nlp_sm_dopus_create_workspace(
    const ocp_nlp_sm_in *sm_in, void *args_) {

}

int_t ocp_nlp_sm_dopus(const ocp_nlp_sm_in *sm_in, ocp_nlp_sm_out *sm_out,
                       void *args_, void *memory_, void *workspace_) {
    // If not initialized, compute sensitivities of all stages

    // Shift all sensitivities backwards, shift here or expect user to do that?

    // Shift all rimal and dual variables backwards, shift here or expect user to do that?

    // Gather indices of stages to recompute sensitivities "accurately"
    //  - Add last stage (always recompued)
    //  - Add M stages with largest value for norm-condition
    //    o Compute norm condition based on historical adjoints and current sensitivities

    // For all indices to recompute, call "accurate" SM method
    
    // For all remaining indicces, call "cheap" SM method

    // Compute adjoints (reuse workspace variables)

    // Compute adjoint gradient correction
}

void ocp_nlp_sm_dopus_initialize(const ocp_nlp_sm_in *sm_in, void *args_,
                                 void **mem, void **work) {

}

void ocp_nlp_sm_dopus_destroy(void *mem_, void *work_) {

}