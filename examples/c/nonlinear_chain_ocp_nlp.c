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

// #include <iostream>
// #include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

// #include "catch/include/catch.hpp"
#include "acados/ocp_nlp/allocate_ocp_nlp.h"
#include "acados/ocp_nlp/ocp_nlp_sm_gn.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/sim/sim_casadi_wrapper.h"
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"
#include "acados/utils/casadi_wrapper.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#include "examples/c/chain_model/chain_model.h"
// #include "test/test_utils/eigen.h"
// #include "test/test_utils/read_matrix  .h"

real_t COMPARISON_TOLERANCE_IPOPT = 1e-6;

#define NN 15
#define TT 3.0
#define Ns 2

// using Eigen::MatrixXd;
// using Eigen::VectorXd;

int main() {
    // TODO(dimitris): fix for NMF > 1
    const int INEXACT = 0;
    const int d = 2;
    const int NMF = 1;
    if (INEXACT == 0) {
        printf(
            "\n----- NUMBER OF FREE MASSES = %d, d = %d (Exact Newton) -----\n",
            NMF, d);
    } else if (INEXACT == 1) {
        printf("\n----- NUMBER OF FREE MASSES = %d, d = %d (IN Scheme) -----\n",
               NMF, d);
    } else if (INEXACT == 2) {
        printf(
            "\n----- NUMBER OF FREE MASSES = %d, d = %d (INIS Scheme) -----\n",
            NMF, d);
    } else if (INEXACT == 3) {
        printf(
            "\n----- NUMBER OF FREE MASSES = %d, d = %d (FROZEN IN Scheme) "
            "-----\n",
            NMF, d);
    } else if (INEXACT == 4) {
        printf(
            "\n----- NUMBER OF FREE MASSES = %d, d = %d (FROZEN INIS Scheme) "
            "-----\n",
            NMF, d);
    }
    int_t NX = 6 * NMF;
    int_t NU = 3;
    int_t jj;

    real_t wall_pos = -0.01;
    int_t UMAX = 10;

    // Problem data
    int_t nx[NN + 1] = {0};
    int_t nu[NN + 1] = {0};
    int_t ny[NN + 1] = {0};
    int_t nb[NN + 1] = {0};
    int_t ng[NN + 1] = {0};
    for (int_t i = 0; i < NN; i++) {
        nx[i] = NX;
        nu[i] = NU;
        ny[i] = NX + NU;
    }
    nx[NN] = NX;
    ny[NN] = NX;

    /************************************************
     * cost
     ************************************************/
    real_t *W, *WN;
    real_t *uref;
    int_t max_sqp_iters = 20;
    real_t *x_end;
    real_t *u_end;

    d_zeros(&W, NX + NU, NX + NU);
    d_zeros(&WN, NX, NX);
    d_zeros(&uref, NU, 1);
    d_zeros(&x_end, NX, 1);
    d_zeros(&u_end, NU, 1);

    real_t x0[] = {0.0000000000000000e+00, 1.5000000000000000e+00,
                   5.0000000000000000e-01, 0.0000000000000000e+00,
                   0.0000000000000000e+00, 0.0000000000000000e+00};

    real_t ref[] = {
        1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00};

    for (int_t i = 0; i < NX; i++) W[i * (NX + NU + 1)] = 1e-2;
    for (int_t i = 0; i < NU; i++) W[(NX + i) * (NX + NU + 1)] = 1.0;
    for (int_t i = 0; i < NX; i++) WN[i * (NX + 1)] = 1e-2;

    ocp_nlp_ls_cost *ls_cost [NN+1];
    for (int_t k = 0; k < NN; k++) {
        allocate_ls_cost(nx[k], nu[k], nx[k]+nu[k], &ls_cost[k]);
        memcpy(ls_cost[k]->W, W, (ny[k] * ny[k]) * sizeof(real_t));
        memcpy(ls_cost[k]->y_ref, ref, ny[k] * sizeof(real_t));
        ls_cost[k]->fun->args->fun = &ls_cost_nm2;
        ls_cost[k]->fun->args->dims = &ls_cost_nm2_work;
        ls_cost[k]->fun->args->sparsity = &ls_cost_nm2_sparsity_out;
        casadi_wrapper_initialize(ls_cost[k]->fun->in, ls_cost[k]->fun->args,
                                  &ls_cost[k]->fun->work);
    }
    allocate_ls_cost(nx[NN], nu[NN], nx[NN] + nu[NN], &ls_cost[NN]);
    memcpy(ls_cost[NN]->W, WN, (ny[NN] * ny[NN]) * sizeof(real_t));
    memcpy(ls_cost[NN]->y_ref, ref, ny[NN] * sizeof(real_t));
    ls_cost[NN]->fun->args->fun = &ls_costN_nm2;
    ls_cost[NN]->fun->args->dims = &ls_costN_nm2_work;
    ls_cost[NN]->fun->args->sparsity = &ls_costN_nm2_sparsity_out;
    casadi_wrapper_initialize(ls_cost[NN]->fun->in, ls_cost[NN]->fun->args,
                              &ls_cost[NN]->fun->work);

    /************************************************
     * simulators
     ************************************************/
    real_t Ts = TT / NN;
    sim_in sim_in[NN];
    sim_out sim_out[NN];
    sim_info info[NN];
    sim_solver *integrators[NN];

    sim_RK_opts rk_opts[NN];
    sim_lifted_irk_memory irk_mem[NN];

    // TODO(rien): can I move this somewhere inside the integrator?
    // struct d_strmat str_mat[NN];
    // struct d_strmat str_sol[NN];

    for (jj = 0; jj < NN; jj++) {
        integrators[jj] = (sim_solver *)malloc(sizeof(sim_solver));
        integrators[jj]->in = &sim_in[jj];
        integrators[jj]->out = &sim_out[jj];
        integrators[jj]->args = &rk_opts[jj];
        if (d > 0) {
            integrators[jj]->fun = &sim_lifted_irk;
            integrators[jj]->mem = &irk_mem[jj];
        } else {
            integrators[jj]->fun = &sim_erk;
            integrators[jj]->mem = 0;
        }

        sim_in[jj].num_steps = Ns;
        sim_in[jj].step = Ts / sim_in[jj].num_steps;
        sim_in[jj].nx = NX;
        sim_in[jj].nu = NU;

        sim_in[jj].sens_forw = true;
        sim_in[jj].sens_adj = false;
        sim_in[jj].sens_hess = false;
        sim_in[jj].num_forw_sens = NX + NU;

        switch (NMF) {
            case 1:
                sim_in[jj].vde = &vde_chain_nm2;
                sim_in[jj].forward_vde_wrapper = &vde_fun;
                sim_in[jj].jac = &jac_chain_nm2;
                sim_in[jj].jacobian_wrapper = &jac_fun;
                break;
            case 2:
                sim_in[jj].vde = &vde_chain_nm3;
                sim_in[jj].forward_vde_wrapper = &vde_fun;
                sim_in[jj].jac = &jac_chain_nm3;
                sim_in[jj].jacobian_wrapper = &jac_fun;
                break;
            case 3:
                sim_in[jj].vde = &vde_chain_nm4;
                sim_in[jj].forward_vde_wrapper = &vde_fun;
                sim_in[jj].jac = &jac_chain_nm4;
                sim_in[jj].jacobian_wrapper = &jac_fun;
                break;
            default:
                break;
        }

        sim_in[jj].x = (real_t *)malloc(sizeof(*sim_in[jj].x) * (NX));
        sim_in[jj].u = (real_t *)malloc(sizeof(*sim_in[jj].u) * (NU));
        sim_in[jj].S_forw =
            (real_t *)malloc(sizeof(*sim_in[jj].S_forw) * (NX * (NX + NU)));
        for (int_t i = 0; i < NX * (NX + NU); i++) sim_in[jj].S_forw[i] = 0.0;
        for (int_t i = 0; i < NX; i++) sim_in[jj].S_forw[i * (NX + 1)] = 1.0;

        sim_in[jj].S_adj =
            (real_t *)malloc(sizeof(*sim_in[jj].S_adj) * (NX + NU));
        for (int_t i = 0; i < NX + NU; i++) sim_in[jj].S_adj[i] = 0.0;

        sim_in[jj].grad_K =
            (real_t *)malloc(sizeof(*sim_in[jj].grad_K) * (d * NX));
        for (int_t i = 0; i < d * NX; i++) sim_in[jj].grad_K[i] = 0.0;

        sim_out[jj].xn = (real_t *)malloc(sizeof(*sim_out[jj].xn) * (NX));
        sim_out[jj].S_forw =
            (real_t *)malloc(sizeof(*sim_out[jj].S_forw) * (NX * (NX + NU)));
        sim_out[jj].info = &info[jj];
        sim_out[jj].grad =
            (real_t *)malloc(sizeof(*sim_out[jj].grad) * (NX + NU));

        int_t workspace_size;
        if (d > 0) {
            sim_irk_create_arguments(&rk_opts[jj], d, "Gauss");
            if (INEXACT == 0) {
                sim_irk_create_Newton_scheme(&rk_opts[jj], d, "Gauss", exact);
            } else if (INEXACT == 1 || INEXACT == 3) {
                sim_irk_create_Newton_scheme(&rk_opts[jj], d, "Gauss",
                                             simplified_in);
            } else if (INEXACT == 2 || INEXACT == 4) {
                sim_irk_create_Newton_scheme(&rk_opts[jj], d, "Gauss",
                                             simplified_inis);
            }

            workspace_size = sim_lifted_irk_calculate_workspace_size(
                &sim_in[jj], &rk_opts[jj]);
            sim_lifted_irk_create_memory(&sim_in[jj], &rk_opts[jj],
                                         &irk_mem[jj]);
        } else {
            sim_erk_create_arguments(&rk_opts[jj], 4);
            workspace_size =
                sim_erk_calculate_workspace_size(&sim_in[jj], &rk_opts[jj]);
        }
        integrators[jj]->work = (void *)malloc(workspace_size);
    }

    /************************************************
     * box constraints
     ************************************************/

    int *idxb0;
    int_zeros(&idxb0, NX + NU, 1);
    real_t *lb0;
    d_zeros(&lb0, NX + NU, 1);
    real_t *ub0;
    d_zeros(&ub0, NX + NU, 1);
#ifdef FLIP_BOUNDS
    for (jj = 0; jj < NU; jj++) {
        lb0[jj] = -UMAX;  // umin
        ub0[jj] = UMAX;   // umax
        idxb0[jj] = jj;
    }
    for (jj = 0; jj < NX; jj++) {
        lb0[NU + jj] = x0[jj];  // xmin
        ub0[NU + jj] = x0[jj];  // xmax
        idxb0[NU + jj] = NU + jj;
    }
#else
    for (jj = 0; jj < NX; jj++) {
        lb0[jj] = x0[jj];  // xmin
        ub0[jj] = x0[jj];  // xmax
        idxb0[jj] = jj;
    }
    for (jj = 0; jj < NU; jj++) {
        lb0[NX + jj] = -UMAX;  // umin
        ub0[NX + jj] = UMAX;   // umax
        idxb0[NX + jj] = NX + jj;
    }
#endif

    nb[0] = NX + NU;

    int *idxb1;
    int_zeros(&idxb1, NMF + NU, 1);
    double *lb1[NN - 1];
    double *ub1[NN - 1];
    for (int_t i = 0; i < NN - 1; i++) {
        d_zeros(&lb1[i], NMF + NU, 1);
        d_zeros(&ub1[i], NMF + NU, 1);
#ifdef FLIP_BOUNDS
        for (jj = 0; jj < NU; jj++) {
            lb1[i][jj] = -UMAX;  // umin
            ub1[i][jj] = UMAX;   // umax
            idxb1[jj] = jj;
        }
        for (jj = 0; jj < NMF; jj++) {
            lb1[i][NU + jj] = wall_pos;  // wall position
            ub1[i][NU + jj] = 1e12;
            idxb1[NU + jj] = NU + 6 * jj + 1;
        }
#else
        for (jj = 0; jj < NMF; jj++) {
            lb1[i][jj] = wall_pos;  // wall position
            ub1[i][jj] = 1e12;
            idxb1[jj] = 6 * jj + 1;
        }
        for (jj = 0; jj < NU; jj++) {
            lb1[i][NMF + jj] = -UMAX;  // umin
            ub1[i][NMF + jj] = UMAX;   // umax
            idxb1[NMF + jj] = NX + jj;
        }
#endif
        nb[i + 1] = NMF + NU;
    }

    int *idxbN;
    int_zeros(&idxbN, NX, 1);
    real_t *lbN;
    d_zeros(&lbN, NX, 1);
    real_t *ubN;
    d_zeros(&ubN, NX, 1);
    for (jj = 0; jj < NX; jj++) {
        lbN[jj] = ref[jj];  // xmin
        ubN[jj] = ref[jj];  // xmax
        idxbN[jj] = jj;
    }
    nb[NN] = NX;

    real_t *hlb[NN + 1];
    real_t *hub[NN + 1];
    int *hidxb[NN + 1];

    hlb[0] = lb0;
    hub[0] = ub0;
    hidxb[0] = idxb0;
    for (int_t i = 1; i < NN; i++) {
        hlb[i] = lb1[i - 1];
        hub[i] = ub1[i - 1];
        hidxb[i] = idxb1;
    }
    hlb[NN] = lbN;
    hub[NN] = ubN;
    hidxb[NN] = idxbN;

    /************************************************
     * nonlinear path constraints
     ************************************************/
    ocp_nlp_function *path_constraints[NN + 1];
    for (int_t k = 0; k < NN; k++) {
        allocate_nlp_function(nx[k], nu[k], nx[k] + nu[k], &path_constraints[k]);
        path_constraints[k]->args->fun = &pathcon_nm2;
        path_constraints[k]->args->dims = &pathcon_nm2_work;
        path_constraints[k]->args->sparsity = &pathcon_nm2_sparsity_out;
        casadi_wrapper_initialize(path_constraints[k]->in,
                                  path_constraints[k]->args,
                                  &path_constraints[k]->work);
    }
    allocate_nlp_function(nx[NN], nu[NN], nx[NN] + nu[NN], &path_constraints[NN]);
    path_constraints[NN]->args->fun = &pathconN_nm2;
    path_constraints[NN]->args->dims = &pathconN_nm2_work;
    path_constraints[NN]->args->sparsity = &pathconN_nm2_sparsity_out;
    casadi_wrapper_initialize(path_constraints[NN]->in,
                              path_constraints[NN]->args,
                              &path_constraints[NN]->work);

    /************************************************
     * sensitivity method
     ************************************************/
    ocp_nlp_sm sensitivity_method;
    sensitivity_method.fun = &ocp_nlp_sm_gn;
    sensitivity_method.initialize = &ocp_nlp_sm_gn_initialize;
    sensitivity_method.destroy = &ocp_nlp_sm_gn_destroy;

    /************************************************
     * QP solver
     ************************************************/
    ocp_qp_solver qp_solver;
    qp_solver.fun = &ocp_qp_condensing_qpoases;
    qp_solver.initialize = &ocp_qp_condensing_qpoases_initialize;
    qp_solver.destroy = &ocp_qp_condensing_qpoases_destroy;
    qp_solver.qp_in = create_ocp_qp_in(NN, nx, nu, nb, ng);
    qp_solver.qp_out = create_ocp_qp_out(NN, nx, nu, nb, ng);
    // TODO(nielsvd): lines below should go
    int_t **idxb = (int_t **)qp_solver.qp_in->idxb;
    for (int_t i = 0; i <= NN; i++)
        for (int_t j = 0; j < nb[i]; j++) idxb[i][j] = hidxb[i][j];
    qp_solver.args =
        (void *)ocp_qp_condensing_qpoases_create_arguments(qp_solver.qp_in);

    /************************************************
     * SQP method
     ************************************************/
    ocp_nlp_in nlp_in;
    nlp_in.N = NN;
    nlp_in.nx = nx;
    nlp_in.nu = nu;
    nlp_in.nb = nb;
    nlp_in.ng = ng;
    nlp_in.idxb = (const int_t **) hidxb;
    nlp_in.lb = (const real_t **) hlb;
    nlp_in.ub = (const real_t **) hub;
    nlp_in.lg = NULL;
    nlp_in.ug = NULL;
    nlp_in.sim = (void **) integrators;
    nlp_in.cost = (void **) ls_cost;
    nlp_in.path_constraints = (void **) path_constraints;

    ocp_nlp_out nlp_out;
    nlp_out.x = (real_t **)malloc(sizeof(*nlp_out.x) * (NN + 1));
    nlp_out.u = (real_t **)malloc(sizeof(*nlp_out.u) * (NN + 1));
    nlp_out.pi = (real_t **)malloc(sizeof(*nlp_out.pi) * (NN + 1));
    nlp_out.lam = (real_t **)malloc(sizeof(*nlp_out.lam) * (NN + 1));
    // Allocate output variables
    for (int_t i = 0; i < NN; i++) {
        nlp_out.x[i] = (real_t *)malloc(sizeof(*nlp_out.x[i]) * (NX));
        nlp_out.u[i] = (real_t *)malloc(sizeof(*nlp_out.u[i]) * (NU));
        nlp_out.pi[i] = (real_t *)malloc(sizeof(*nlp_out.pi[i]) * (NX));
        nlp_out.lam[i] =
            (real_t *)malloc(sizeof(*nlp_out.lam[i]) * 2 * nb[i] + 2 * ng[i]);
    }
    nlp_out.x[NN] = (real_t *)malloc(sizeof(*nlp_out.x[NN]) * (NX));
    nlp_out.u[NN] = (real_t *)malloc(sizeof(*nlp_out.u[NN]) * 0);
    nlp_out.pi[NN] = (real_t *)malloc(sizeof(*nlp_out.pi[NN]) * (0));
    nlp_out.lam[NN] =
        (real_t *)malloc(sizeof(*nlp_out.lam[NN]) * 2 * nb[NN] + 2 * ng[NN]);

    ocp_nlp_sqp_args *nlp_args = ocp_nlp_sqp_create_arguments();
    nlp_args->maxIter = max_sqp_iters;
    nlp_args->sensitivity_method = &sensitivity_method;
    nlp_args->qp_solver = &qp_solver;

    // snprintf(nlp_args.qp_solver_name, sizeof(nlp_args.qp_solver_name), "%s",
    // "condensing_qpoases");

    ocp_nlp_sqp_memory *nlp_mem;
    ocp_nlp_sqp_workspace *nlp_work;
    ocp_nlp_sqp_initialize(&nlp_in, nlp_args, (void **)&nlp_mem,
                           (void **)&nlp_work);

    real_t **nlp_x_mem = (real_t **)nlp_mem->common->x;
    real_t **nlp_u_mem = (real_t **)nlp_mem->common->u;
    for (int_t i = 0; i < NN; i++) {
        for (int_t j = 0; j < NX; j++) nlp_x_mem[i][j] = ref[j];  // resX(j,i)
        for (int_t j = 0; j < NU; j++) nlp_u_mem[i][j] = ref[NX+j];      // resU(j, i)
    }
    for (int_t j = 0; j < NX; j++) nlp_x_mem[NN][j] = ref[j];  // resX(j, NN)

    int_t status;

    status = ocp_nlp_sqp(&nlp_in, &nlp_out, nlp_args, nlp_mem, nlp_work);
    printf("\n\nstatus = %i\n\n", status);

// for (int_t k =0; k < 3; k++) {
//     printf("x[%d] = \n", k);
//     d_print_mat(1, nx[k], nlp_out.x[k], 1);
//     printf("u[%d] = \n", k);
//     d_print_mat(1, nu[k], nlp_out.u[k], 1);
// }

#if 0
    real_t out_x[NX*(NN+1)];
    real_t out_u[NU*NN];
    for (int_t i = 0; i < NN; i++) {
        for (int_t j = 0; j < NX; j++) out_x[i*NX+j] = nlp_out.x[i][j];
        for (int_t j = 0; j < NU; j++) out_u[i*NU+j] = nlp_out.u[i][j];
    }
    for (int_t j = 0; j < NX; j++) out_x[NN*NX+j] = nlp_out.x[NN][j];
#endif

    d_free(W);
    d_free(WN);
    d_free(uref);
    d_free(x_end);
    d_free(u_end);

    int_free(idxb0);
    d_free(lb0);
    d_free(ub0);
    int_free(idxb1);
    for (jj = 0; jj < NN - 1; jj++) {
        d_free(lb1[jj]);
        d_free(ub1[jj]);
    }
    int_free(idxbN);
    d_free(lbN);
    d_free(ubN);

    // LS cost and path constraints
    for (int_t i = 0; i <= NN; i++) {
        // Least-squares cost
        free_ls_cost(ls_cost[i]);
        free_nlp_function(path_constraints[i]);
    }

    // Integrators
    for (jj = 0; jj < NN; jj++) {
        free(sim_in[jj].x);
        free(sim_in[jj].u);
        free(sim_in[jj].S_forw);
        free(sim_in[jj].S_adj);
        free(sim_in[jj].grad_K);
        free(sim_out[jj].xn);
        free(sim_out[jj].S_forw);
        free(sim_out[jj].grad);
    }

    // NLP memory and workspace
    ocp_nlp_sqp_destroy(nlp_mem, nlp_work);
    // NLP arguments
    free(nlp_args);
    // NLP output
    for (int_t i = 0; i <= NN; i++) {
        free(nlp_out.x[i]);
        free(nlp_out.u[i]);
        free(nlp_out.pi[i]);
        free(nlp_out.lam[i]);
    }
    free(nlp_out.x);
    free(nlp_out.u);
    free(nlp_out.lam);
    free(nlp_out.pi);
}
