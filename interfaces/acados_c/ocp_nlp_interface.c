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

#include "acados_c/ocp_nlp_interface.h"

// external
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "acados/ocp_nlp/ocp_nlp_cost_external.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"
#include "acados/ocp_nlp/ocp_nlp_cost_nls.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_disc.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/utils/mem.h"

int ocp_nlp_plan_calculate_size(int N)
{
    int bytes = sizeof(ocp_nlp_solver_plan);
    bytes += N * sizeof(sim_solver_plan);
    bytes += (N + 1) * sizeof(ocp_nlp_cost_t);
    bytes += N * sizeof(ocp_nlp_dynamics_t);
    return bytes;
}

ocp_nlp_solver_plan *ocp_nlp_plan_assign(int N, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_solver_plan *plan = (ocp_nlp_solver_plan *) c_ptr;
    c_ptr += sizeof(ocp_nlp_solver_plan);

    plan->sim_solver_plan = (sim_solver_plan *) c_ptr;
    c_ptr += N * sizeof(sim_solver_plan);

    plan->nlp_cost = (ocp_nlp_cost_t *) c_ptr;
    c_ptr += (N + 1) * sizeof(ocp_nlp_cost_t);

    plan->nlp_dynamics = (ocp_nlp_dynamics_t *) c_ptr;
    c_ptr += (N + 1) * sizeof(ocp_nlp_dynamics_t);

    // TODO(all): fix assert
    // assert( 0 == 0);

    return plan;
}

void ocp_nlp_plan_initialize_default(int N, ocp_nlp_solver_plan *plan)
{
    plan->nlp_solver = SQP_GN;
    for (int ii = 0; ii <= N; ii++)
    {
        plan->nlp_cost[ii] = NONLINEAR_LS;
        if (ii < N)
        {
            plan->sim_solver_plan[ii].sim_solver = ERK;
        }
    }
}

ocp_nlp_solver_plan *ocp_nlp_plan_create(int N)
{
    int bytes = ocp_nlp_plan_calculate_size(N);
    void *ptr = acados_malloc(bytes, 1);

    ocp_nlp_solver_plan *plan = ocp_nlp_plan_assign(N, ptr);

    ocp_nlp_plan_initialize_default(N, plan);

    return plan;
}

// TODO(dimitris): this leaks memory! Either provide free config or calculate size should be nested
ocp_nlp_solver_config *ocp_nlp_config_create(ocp_nlp_solver_plan plan, int N)
{
    int bytes = ocp_nlp_solver_config_calculate_size(N);
    void *config_mem = calloc(1, bytes);
    ocp_nlp_solver_config *config = ocp_nlp_solver_config_assign(N, config_mem);

    if (plan.nlp_solver == SQP_GN)
    {
        ocp_nlp_sqp_config_initialize_default(config);

        // QP solver
        config->qp_solver = ocp_qp_config_create(plan.ocp_qp_solver_plan);

        // cost
        for (int i = 0; i <= N; ++i)
        {
            switch (plan.nlp_cost[i])
            {
                case LINEAR_LS:
                    ocp_nlp_cost_ls_config_initialize_default(config->cost[i]);
                    break;
                case NONLINEAR_LS:
                    ocp_nlp_cost_nls_config_initialize_default(config->cost[i]);
                    break;
                case EXTERNALLY_PROVIDED:
                    ocp_nlp_cost_external_config_initialize_default(config->cost[i]);
                    break;
                default:
                    printf("Cost not available!\n");
                    exit(1);
            }
        }

        // Dynamics
        for (int i = 0; i < N; ++i)
        {
            switch (plan.nlp_dynamics[i])
            {
                case CONTINUOUS_MODEL:
                    ocp_nlp_dynamics_cont_config_initialize_default(config->dynamics[i]);
                    config->dynamics[i]->sim_solver = sim_config_create(plan.sim_solver_plan[i]);
                    break;
                case DISCRETE_MODEL:
                    ocp_nlp_dynamics_disc_config_initialize_default(config->dynamics[i]);
                    break;
                default:
                    printf("Dynamics not available!\n");
                    exit(1);
            }
        }

        // Constraints
        for (int i = 0; i <= N; ++i)
            ocp_nlp_constraints_config_initialize_default(config->constraints[i]);
    }
    else
    {
        printf("Solver not available!\n");
        exit(1);
    }
    return config;
}

ocp_nlp_dims *ocp_nlp_dims_create(void *config_)
{
    ocp_nlp_solver_config *config = config_;

    // int N = config->N;

    int bytes = ocp_nlp_dims_calculate_size(config);

    void *ptr = calloc(1, bytes);

    ocp_nlp_dims *dims = ocp_nlp_dims_assign(config, ptr);

    return dims;
}

ocp_nlp_in *ocp_nlp_in_create(ocp_nlp_solver_config *config, ocp_nlp_dims *dims)
{
    int bytes = ocp_nlp_in_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    ocp_nlp_in *nlp_in = ocp_nlp_in_assign(config, dims, ptr);

    return nlp_in;
}

int nlp_set_model_in_stage(ocp_nlp_solver_config *config, ocp_nlp_in *in, int stage,
                           const char *fun_type, void *fun_ptr)
{
    // NOTE(giaf) @dimitris, how do we do it with discrete model dynamics ?
    sim_solver_config *sim_config = config->dynamics[stage]->sim_solver;
    ocp_nlp_dynamics_cont_model *dynamics = in->dynamics[stage];

    int status = sim_set_model_internal(sim_config, dynamics->sim_model, fun_type, fun_ptr);

    return status;
}

ocp_nlp_out *ocp_nlp_out_create(ocp_nlp_solver_config *config, ocp_nlp_dims *dims)
{
    int bytes = ocp_nlp_out_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    ocp_nlp_out *nlp_out = ocp_nlp_out_assign(config, dims, ptr);

    // initialize to zeros
    for (int ii = 0; ii <= dims->N; ++ii)
        blasfeo_dvecse(dims->qp_solver->nu[ii] + dims->qp_solver->nx[ii], 0.0, nlp_out->ux + ii, 0);

    return nlp_out;
}

void *ocp_nlp_opts_create(ocp_nlp_solver_config *config, ocp_nlp_dims *dims)
{
    int bytes = config->opts_calculate_size(config, dims);

    void *ptr = calloc(1, bytes);

    void *opts = config->opts_assign(config, dims, ptr);

    config->opts_initialize_default(config, dims, opts);

    return opts;
}

int ocp_nlp_calculate_size(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, void *opts_)
{
    int bytes = sizeof(ocp_nlp_solver);

    bytes += config->memory_calculate_size(config, dims, opts_);
    bytes += config->workspace_calculate_size(config, dims, opts_);

    return bytes;
}

ocp_nlp_solver *ocp_nlp_assign(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, void *opts_,
                               void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_solver *solver = (ocp_nlp_solver *) c_ptr;
    c_ptr += sizeof(ocp_nlp_solver);

    solver->config = config;
    solver->dims = dims;
    solver->opts = opts_;

    solver->mem = config->memory_assign(config, dims, opts_, c_ptr);
    c_ptr += config->memory_calculate_size(config, dims, opts_);

    solver->work = (void *) c_ptr;
    c_ptr += config->workspace_calculate_size(config, dims, opts_);

    assert((char *) raw_memory + ocp_nlp_calculate_size(config, dims, opts_) == c_ptr);

    return solver;
}

ocp_nlp_solver *ocp_nlp_create(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, void *opts_)
{
    config->opts_update(config, dims, opts_);

    int bytes = ocp_nlp_calculate_size(config, dims, opts_);

    void *ptr = calloc(1, bytes);

    ocp_nlp_solver *solver = ocp_nlp_assign(config, dims, opts_, ptr);

    return solver;
}

int ocp_nlp_solve(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)
{
    return solver->config->evaluate(solver->config, solver->dims, nlp_in, nlp_out, solver->opts,
                                    solver->mem, solver->work);
}
