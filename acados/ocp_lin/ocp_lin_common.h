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

#ifndef ACADOS_OCP_LIN_OCP_LIN_COMMON_H_
#define ACADOS_OCP_LIN_OCP_LIN_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/sim/sim_common.h"
#include "acados/utils/external_function.h"
#include "acados/utils/types.h"



typedef struct {
    int *nx;
    int *nu;
    int *np;
    int *nb;  // nbx + nbu
    int *nbx;
    int *nbu;
    int *ng;  // number of general linear constraints
    int *nh;  // number of path constraints - ONLY difference with ocp_qp_dims atm
    int *ns;  // number of soft constraints
    int *num_stages;
    int N;
} ocp_lin_dims;



typedef struct {
    real_t **x;
    real_t **u;
    real_t **p;
    real_t **pi;
    real_t **lam;
} ocp_lin_in;



typedef struct {
    real_t **hess_l;       // TODO(nielsvd): Hessians of stage-wise
                           // Lagrangians, document precise definition.
    real_t **grad_f;       // Gradients of stage-wise cost terms.
    real_t **jac_xp;       // Jacobians of stage-wise integration operator.
    real_t **jac_h;        // Jacobians of stage-wise path constraints.
    real_t **grad_pi_xp;   // Adjoint derivative of system dynamics
    real_t **grad_lam_h;   // Adjoint derivative of path constraints
    real_t **xp;           // Evaluation of stage-wise integration operator.
    real_t **h;            // Evaluation of stage-wise path constraints.
} ocp_lin_out;



typedef struct {
    int (*fun)(ocp_lin_in *qp_in, ocp_lin_out *qp_out, void *args, void *mem, void *work);
    int (*calculate_args_size)(ocp_lin_dims *dims, void *submodules);
    void *(*assign_args)(ocp_lin_dims *dims, void **submodules, void *raw_memory);
    void *(*copy_args)(ocp_lin_dims *dims, void *raw_memory, void *source);
    void (*initialize_default_args)(ocp_lin_dims *dims, void *args);
    int (*calculate_memory_size)(ocp_lin_dims *dims, void *args);
    void *(*assign_memory)(ocp_lin_dims *dims, void *args, void *raw_memory);
    int (*calculate_workspace_size)(ocp_lin_dims *dims, void *args);
    void *submodules;
} ocp_lin_method_fcn_ptrs;



//
int ocp_lin_dims_calculate_size(int N);
//
ocp_lin_dims *assign_ocp_lin_dims(int N, void *raw_memory);
//
int ocp_lin_in_calculate_size(ocp_lin_dims *dims);
//
ocp_lin_in *assign_ocp_lin_in(ocp_lin_dims *dims, void *raw_memory);
//
int ocp_lin_out_calculate_size(ocp_lin_dims *dims);
//
ocp_lin_out *assign_ocp_lin_out(ocp_lin_dims *dims, void *raw_memory);
//



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_LIN_OCP_LIN_COMMON_H_
