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


#ifndef ACADOS_OCP_LIN_OCP_LIN_GN_H_
#define ACADOS_OCP_LIN_OCP_LIN_GN_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_lin/ocp_lin_common.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/external_function.h"



typedef struct {
    external_function_fcn_ptrs **ls_res;
    external_function_fcn_ptrs **h;
    sim_solver_fcn_ptrs **xp;
} ocp_lin_gn_submodules;



typedef struct {
    ocp_lin_gn_submodules submodules;
    
    void **ls_res_args;
    void **h_args;
    void **xp_args;
} ocp_lin_gn_args;



typedef struct {
    void **ls_res_mem;
    void **h_mem;
    void **xp_mem;
} ocp_lin_gn_memory;



typedef struct {
    double **S_adj_in;
    double **S_forw_in;
    double **F;
    double **DF;
    double **DFT;
    double **H;
    double **DH;

    void **ls_res_work;
    void **h_work;
    void **xp_work;
} ocp_lin_gn_workspace;



//
int ocp_lin_gn_calculate_args_size(ocp_lin_dims *dims, void *submodules_);
//
void *ocp_lin_gn_assign_args(ocp_lin_dims *dims, void **submodules_, void *raw_memory);
//
void *ocp_lin_gn_copy_args(ocp_lin_dims *dims, void *raw_memory, void *source_);
//
void ocp_lin_gn_initialize_default_args(ocp_lin_dims *dims, void *args_);
//
int ocp_lin_gn_calculate_memory_size(ocp_lin_dims *dims, void *args_);
//
void *ocp_lin_gn_assign_memory(ocp_lin_dims *dims, void *args_, void *raw_memory);
//
int ocp_lin_gn_calculate_workspace_size(ocp_lin_dims *dims, void *args_);
//
int ocp_lin_gn(ocp_lin_in *lin_in, ocp_lin_out *lin_out, void *args_, void *memory_, void *work_);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_LIN_OCP_LIN_GN_H_