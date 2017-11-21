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

#include "acados/ocp_nlp/ocp_nlp_sm_eh.h"

#include <assert.h>
#include <stdlib.h>

#include "acados/sim/sim_common.h"
#include "acados/sim/sim_rk_common.h"
#include "acados/utils/math.h"
#include "acados/utils/types.h"

ocp_nlp_sm_eh_args *ocp_nlp_sm_eh_create_arguments() {
    return 0;
}

int_t ocp_nlp_sm_eh_calculate_memory_size(const ocp_nlp_sm_in *sm_in,
                                          void *args_) {
    return 0;
}

char *ocp_nlp_sm_eh_assign_memory(const ocp_nlp_sm_in *sm_in, void *args_,
                                  void **mem_, void *raw_memory) {
    return 0;
}

ocp_nlp_sm_eh_memory *ocp_nlp_sm_eh_create_memory(const ocp_nlp_sm_in *sm_in,
                                                  void *args_) {
    return 0;
}

int_t ocp_nlp_sm_eh_calculate_workspace_size(const ocp_nlp_sm_in *sm_in,
                                             void *args_) {
    return 0;
}

char *ocp_nlp_sm_eh_assign_workspace(const ocp_nlp_sm_in *sm_in, void *args_,
                                     void **work_, void *raw_memory) {
    return 0;
}

ocp_nlp_sm_eh_workspace *ocp_nlp_sm_eh_create_workspace(
    const ocp_nlp_sm_in *sm_in, void *args_) {
    return 0;
}

int_t ocp_nlp_sm_eh(const ocp_nlp_sm_in *sm_in, ocp_nlp_sm_out *sm_out,
                    void *args_, void *memory_, void *workspace_) {
    // ocp_nlp_sm_eh_workspace *work = (ocp_nlp_sm_eh_workspace *) workspace_;
    // ocp_nlp_sm_eh_memory *mem = (ocp_nlp_sm_eh_memory *) memory_;

    // const int_t N = sm_in->N;
    // const int_t *nx = sm_in->nx;
    // const int_t *nu = sm_in->nu;
    // const int_t *ng = sm_in->ng;

    // real_t **hess_l = (real_t **)sm_out->hess_l;
    // real_t **grad_f = (real_t **)sm_out->grad_f;
    // real_t **jac_h = (real_t **)sm_out->jac_h;
    // real_t **jac_g = (real_t **)sm_out->jac_g;
    // real_t **h = (real_t **)sm_out->h;
    // real_t **g = (real_t **)sm_out->g;

    // sim_solver **sim = sm_in->sim;
    // ocp_nlp_function **cost = (ocp_nlp_function **) sm_in->cost;
    // ocp_nlp_function **path_constraints = sm_in->path_constraints;



    return 0;
}

void ocp_nlp_sm_eh_initialize(const ocp_nlp_sm_in *sm_in, void *args_,
                              void **mem, void **work) {
    
}

void ocp_nlp_sm_eh_destroy(void *mem_, void *work_) {

}