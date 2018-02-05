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



typedef struct ocp_lin_gn_args_ {

} ocp_lin_gn_args;



typedef struct ocp_lin_gn_memory_ {

} ocp_lin_gn_memory;



//
int ocp_lin_gn_calculate_args_size(ocp_lin_dims *dims, void *submodules_);
//
void *ocp_lin_gn_assign_args(ocp_lin_dims *dims, void *submodules_, void *raw_memory);
//
void ocp_lin_gn_initialize_default_args(void *args_);
//
int ocp_lin_gn_calculate_memory_size(ocp_lin_dims *dims, void *args_);
//
void *ocp_lin_gn_assign_memory(ocp_lin_dims *dims, void *args_, void *raw_memory);
//
int ocp_lin_gn_calculate_workspace_size(ocp_lin_dims *dims, void *args_);
//
int ocp_lin_gn(ocp_lin_in *qp_in, ocp_lin_out *qp_out, void *args_, void *memory_, void *work_);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_LIN_OCP_LIN_GN_H_