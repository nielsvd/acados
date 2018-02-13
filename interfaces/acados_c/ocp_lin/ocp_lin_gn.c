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

#include "acados_c/ocp_lin/ocp_lin_gn.h"

#include <assert.h>
#include <string.h>

#include "acados_c/external_function.h"
#include "acados_c/sim.h"



int ocp_lin_gn_calculate_submodules_size(ocp_lin_method_config *config, ocp_lin_dims *dims)
{
    int size = sizeof(ocp_lin_gn_submodules);

    int N = dims->N;

    size += 2*(N+1)*sizeof(external_function_fcn_ptrs);

    for (int i=0;i<=N;i++)
    {
        extern external_function_dims ocp_lin_gn_ls_res_dims;
        size += calculate_external_function_fcn_ptrs_size(&config->extfun_ls_res, &ocp_lin_gn_ls_res_dims);

        extern external_function_dims ocp_lin_gn_h_dims;
        size += calculate_external_function_fcn_ptrs_size(&config->extfun_h, &ocp_lin_gn_h_dims);        
    }

    size += N*sizeof(sim_solver_fcn_ptrs);

    for (int i=0;i<N;i++)
    {
        sim_dims sdi;
        sdi.num_stages = dims->num_stages[i];
        sdi.nx = dims->nx[i];
        sdi.nu = dims->nu[i];
        sdi.np = dims->np[i];

        size += calculate_sim_solver_fcn_ptrs_size(&config->simsol, &sdi);
    }

    return size;
}



void *ocp_lin_gn_assign_submodules(ocp_lin_method_config *config, ocp_lin_dims *dims, void *raw_memory)
{
    ocp_lin_gn_submodules *submodules;

    int N = dims->N;

    char *c_ptr = (char *) raw_memory;

    submodules = (ocp_lin_gn_submodules *) c_ptr;
    c_ptr += sizeof(ocp_lin_gn_submodules);

    submodules->ls_res = (external_function_fcn_ptrs **) c_ptr;
    c_ptr += (N+1)*sizeof(external_function_fcn_ptrs *);

    submodules->h = (external_function_fcn_ptrs **) c_ptr;
    c_ptr += (N+1)*sizeof(external_function_fcn_ptrs *);

    for (int i=0;i<=N;i++)
    {
        extern external_function_dims ocp_lin_gn_ls_res_dims;
        submodules->ls_res[i] = assign_external_function_fcn_ptrs(&config->extfun_ls_res, &ocp_lin_gn_ls_res_dims, c_ptr);
        c_ptr += calculate_external_function_fcn_ptrs_size(&config->extfun_ls_res, &ocp_lin_gn_ls_res_dims);

        extern external_function_dims ocp_lin_gn_h_dims;
        submodules->h[i] = assign_external_function_fcn_ptrs(&config->extfun_h, &ocp_lin_gn_h_dims, c_ptr);
        c_ptr += calculate_external_function_fcn_ptrs_size(&config->extfun_h, &ocp_lin_gn_h_dims);
    }

    submodules->xp = (sim_solver_fcn_ptrs **) c_ptr;
    c_ptr += N*sizeof(sim_solver_fcn_ptrs *);

    for (int i=0;i<N;i++)
    {
        sim_dims sdi;
        sdi.num_stages = dims->num_stages[i];
        sdi.nx = dims->nx[i];
        sdi.nu = dims->nu[i];
        sdi.np = dims->np[i];

        submodules->xp[i] = assign_sim_solver_fcn_ptrs(&config->simsol, &sdi, c_ptr);
        c_ptr += calculate_sim_solver_fcn_ptrs_size(&config->simsol, &sdi);
    }

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    assert((char*)raw_memory + ocp_lin_gn_calculate_submodules_size(config, dims) == c_ptr);

    return (void *)submodules;
}
