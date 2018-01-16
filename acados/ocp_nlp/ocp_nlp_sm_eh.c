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
    ocp_nlp_sm_eh_workspace *work = (ocp_nlp_sm_eh_workspace *) workspace_;
    // ocp_nlp_sm_eh_memory *mem = (ocp_nlp_sm_eh_memory *) memory_;

    const int_t N = sm_in->N;
    const int_t *nx = sm_in->nx;
    const int_t *nu = sm_in->nu;
    const int_t *ng = sm_in->ng;

    real_t **hess_l = sm_out->hess_l;
    real_t **grad_f = sm_out->grad_f;
    real_t **jac_h = sm_out->jac_h;
    real_t **jac_g = sm_out->jac_g;
    real_t **h = sm_out->h;
    real_t **g = sm_out->g;

    sim_solver **sim = sm_in->sim;
    ocp_nlp_function **cost = (ocp_nlp_function **) sm_in->cost;
    ocp_nlp_function **path_constraints = sm_in->path_constraints;

    for (int_t i = 0; i < N; i++) {
        // Pass state and control to integrator
        for (int_t j = 0; j < nx[i]; j++) sim[i]->in->x[j] = sm_in->x[i][j];
        for (int_t j = 0; j < nu[i]; j++) sim[i]->in->u[j] = sm_in->u[i][j];
        sim[i]->fun(sim[i]->in, sim[i]->out, sim[i]->args, sim[i]->mem,
                    sim[i]->work);

        // Sensitivities for the linearization of the system dynamics
        // TODO(rien): transition functions for changing dimensions not yet
        // implemented!
        for (int_t j = 0; j < nx[i]; j++) {
            h[i][j] = sim[i]->out->xn[j];
            for (int_t k = 0; k < nx[i] + nu[i]; k++)
                jac_h[i][k * nx[i] + j] = sim[i]->out->S_forw[k * nx[i] + j];
        }
    }

    for (int_t i = 0; i <= N; i++) {        
        // Gradient and Hessian of objective
        casadi_wrapper_in *f_in = cost[i]->in;
        casadi_wrapper_out *f_out = cost[i]->out;
        casadi_wrapper_args *f_args = cost[i]->args;
        casadi_wrapper_workspace *f_work = cost[i]->work;
        // Set input and output arguments
        f_in->x = sm_in->x[i];
        f_in->p = sm_in->u[i];
        f_in->x = NULL;
        f_in->lag = NULL;
        f_out->jac_y = sm_out->grad_f;
        f_out->hess_y = sm_out->hess_l;
        casadi_wrapper(f_in, f_out, f_args, f_work);

        // Add second order information of dynamics
        for (int_t j = 0; j < (nx[i] + nu[i]) * (nx[i] + nu[i]); j++) {
            hess_l[i][j] += work->hess_pi_f[i][j];
        }

        // Jacobian of inequality constraints
        if (sm_in->ng[i] > 0) {
            // Path constraints for shooting node i
            casadi_wrapper_in *pc_in = path_constraints[i]->in;
            casadi_wrapper_out *pc_out = path_constraints[i]->out;
            casadi_wrapper_args *pc_args = path_constraints[i]->args;
            casadi_wrapper_workspace *pc_work = path_constraints[i]->work;
            // Set input and output arguments
            f_in->x = sm_in->x[i];
            f_in->u = sm_in->u[i];
            f_in->p = NULL;
            f_in->lag = sm_in->lam[i];
            f_out->y = sm_out->g[i];
            f_out->jac_y = sm_out->jac_g[i];
            f_out->grad_adj = work->hess_lam_g[i];
            casadi_wrapper(pc_in, pc_out, pc_args, pc_work);

            // Add second order information of path constraints
            for (int_t j = 0; j < (nx[i] + nu[i]) * (nx[i] + nu[i]); j++) {
                hess_l[i][j] += work->hess_lam_g[i][j];
            }
        }

        // Regularize Hessian
    }

    return 0;
}

void ocp_nlp_sm_eh_initialize(const ocp_nlp_sm_in *sm_in, void *args_,
                              void **mem, void **work) {
    
}

void ocp_nlp_sm_eh_destroy(void *mem_, void *work_) {

}