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

#include "acados/ocp_lin/ocp_lin_gn.h"

#include <stdlib.h>



#define LS_RES_NUMIN 3
#define LS_RES_NUMOUT 2
#define H_NUMIN 3
#define H_NUMOUT 2



int ocp_lin_gn_ls_res_input_dims[LS_RES_NUMIN] = {0};
int ocp_lin_gn_ls_res_output_dims[LS_RES_NUMOUT] = {0};
external_function_dims ocp_lin_gn_ls_res_dims = {
    LS_RES_NUMIN,
    LS_RES_NUMOUT,
    ocp_lin_gn_ls_res_input_dims,
    ocp_lin_gn_ls_res_output_dims
};



int ocp_lin_gn_h_input_dims[H_NUMIN] = {0};
int ocp_lin_gn_h_output_dims[H_NUMOUT] = {0};
external_function_dims ocp_lin_gn_h_dims = {
    H_NUMIN,
    H_NUMOUT,
    ocp_lin_gn_h_input_dims,
    ocp_lin_gn_h_output_dims
};



int ocp_lin_gn_calculate_args_size(ocp_lin_dims *dims, void *submodules_)
{
    ocp_lin_gn_submodules *submodules = (ocp_lin_gn_submodules *) submodules_;
    
    int N = dims->N;

    int size = 0;

    size += 2*(N+1)*sizeof(external_function_fcn_ptrs);

    for (int i=0;i<=N;i++)
    {
        if (submodules->ls_res[i] != NULL)
        {
            size += submodules->ls_res[i]->calculate_args_size(&ocp_lin_gn_ls_res_dims, submodules->ls_res[i]->submodules);
        }
    
        if (submodules->h != NULL) 
        {
            size += submodules->h[i]->calculate_args_size(&ocp_lin_gn_h_dims, submodules->h[i]->submodules);
        }
    }

    size += N*sizeof(sim_solver_fcn_ptrs);

    for (int i=0;i<N;i++)
    {
        sim_dims sdi;
        sdi.num_stages = 3;  // TODO(nielsvd): remove, should be an argument
        sdi.nx = dims->nx[i];
        sdi.nu = dims->nu[i];
        sdi.np = dims->np[i];
     
        size += submodules->xp[i]->calculate_args_size(&sdi, submodules->xp[i]->submodules);
    }

    return size;
}



void *ocp_lin_gn_assign_args(ocp_lin_dims *dims, void **submodules_, void *raw_memory)
{
    ocp_lin_gn_submodules *submodules = (ocp_lin_gn_submodules *) *submodules_;

    int N = dims->N;

    char *c_ptr = (char *) raw_memory;

    ocp_lin_gn_args *args = (ocp_lin_gn_args *) c_ptr;
    c_ptr += sizeof(ocp_lin_gn_args);

    args->submodules.ls_res = (external_function_fcn_ptrs **) c_ptr;
    c_ptr += (N+1) * sizeof(external_function_fcn_ptrs *);

    args->submodules.h = (external_function_fcn_ptrs **) c_ptr;
    c_ptr += (N + 1) * sizeof(external_function_fcn_ptrs *);

    for (int i=0;i<=N;i++)
    {
        if (submodules->ls_res[i] != NULL)
        {
            args->submodules.ls_res[i] = (external_function_fcn_ptrs *) c_ptr;
            c_ptr += sizeof(external_function_fcn_ptrs);

            void *ls_res_submodules = submodules->ls_res[i]->submodules;
            args->ls_res_args[i] = submodules->ls_res[i]->assign_args(&ocp_lin_gn_ls_res_dims, &(ls_res_submodules), c_ptr);
            c_ptr += submodules->ls_res[i]->calculate_args_size(&ocp_lin_gn_ls_res_dims, submodules->ls_res[i]->submodules);

            *(args->submodules.ls_res[i]) = *(submodules->ls_res[i]);
            args->submodules.ls_res[i]->submodules = ls_res_submodules;
        }
        else
        {
            args->submodules.ls_res[i] = NULL;
        }

        if (submodules->h[i] != NULL)
        {
            args->submodules.h[i] = (external_function_fcn_ptrs *) c_ptr;
            c_ptr += sizeof(external_function_fcn_ptrs);

            void *h_submodules = submodules->h[i]->submodules;
            args->h_args[i] = submodules->h[i]->assign_args(&ocp_lin_gn_h_dims, &(h_submodules), c_ptr);
            c_ptr += submodules->h[i]->calculate_args_size(&ocp_lin_gn_h_dims, submodules->h[i]->submodules);

            *(args->submodules.h[i]) = *(submodules->h[i]);
            args->submodules.h[i]->submodules = h_submodules;
        }
        else
        {
            args->submodules.h[i] = NULL;
        }
    }

    args->submodules.xp = (sim_solver_fcn_ptrs **) c_ptr;
    c_ptr += N * sizeof(sim_solver_fcn_ptrs *);

    for (int i=0;i<N;i++)
    {
        args->submodules.xp[i] = (sim_solver_fcn_ptrs *) c_ptr;
        c_ptr += sizeof(sim_solver_fcn_ptrs);

        sim_dims sdi;
        sdi.num_stages = 3;  // TODO(nielsvd): remove, should be an argument
        sdi.nx = dims->nx[i];
        sdi.nu = dims->nu[i];
        sdi.np = dims->np[i];

        void *xp_submodules = submodules->xp[i]->submodules;
        args->xp_args[i] = submodules->xp[i]->assign_args(&sdi, &xp_submodules, c_ptr);
        c_ptr += submodules->xp[i]->calculate_args_size(&sdi, submodules->xp[i]->submodules);

        *(args->submodules.xp[i]) = *(submodules->xp[i]);
        args->submodules.xp[i] = xp_submodules;
    }

    assert((char*)raw_memory + ocp_lin_gn_calculate_args_size(dims, *submodules_) >= c_ptr);

    // Update submodules pointer
    *submodules_ = (void *)&(args->submodules);

    return (void *)args;
}



void *ocp_lin_gn_copy_args(ocp_lin_dims *dims, void *raw_memory, void *source_)
{
    ocp_lin_gn_args *source = (ocp_lin_gn_args *) source_;
    ocp_lin_gn_args *dest;

    int N = dims->N;

    for (int i=0;i<=N;i++)
    {
        source->submodules.ls_res[i]->copy_args(&ocp_lin_gn_ls_res_dims, dest->ls_res_args[i], source->ls_res_args[i]);
        source->submodules.h[i]->copy_args(&ocp_lin_gn_h_dims, dest->h_args[i], source->h_args[i]);
    }

    for (int i=0;i<N;i++)
    {
        sim_dims sdi;
        sdi.num_stages = 3;  // TODO(nielsvd): remove, should be an argument
        sdi.nx = dims->nx[i];
        sdi.nu = dims->nu[i];
        sdi.np = dims->np[i];

        source->submodules.xp[i]->copy_args(&sdi, dest->xp_args[i], source->xp_args[i]);
    }

    return (void *)dest;
}



void ocp_lin_gn_initialize_default_args(ocp_lin_dims *dims, void *args_)
{
    ocp_lin_gn_args *args = (ocp_lin_gn_args *) args_;

    int N = dims->N;

    for (int i=0;i<=N;i++)
    {
        args->submodules.ls_res[i]->initialize_default_args(args->ls_res_args[i]);
        args->submodules.h[i]->initialize_default_args(args->h_args[i]);
    }

    for (int i=0;i<N;i++)
    {
        args->submodules.ls_res[i]->initialize_default_args(args->ls_res_args[i]);
    }
}



int ocp_lin_gn_calculate_memory_size(ocp_lin_dims *dims, void *args_)
{
    ocp_lin_gn_args *args = (ocp_lin_gn_args *) args_;

    int N = dims->N;

    int size = sizeof(ocp_lin_gn_memory);

    size += 2*(N+1)*sizeof(void *);

    for (int i=0;i<=N;i++)
    {
        size += args->submodules.ls_res[i]->calculate_memory_size(&ocp_lin_gn_ls_res_dims, args->ls_res_args[i]);
        size += args->submodules.h[i]->calculate_memory_size(&ocp_lin_gn_h_dims, args->h_args[i]);
    }

    size += N*sizeof(void *);

    for (int i=0;i<N;i++)
    {
        sim_dims sdi;
        sdi.num_stages = 3;  // TODO(nielsvd): remove, should be an argument
        sdi.nx = dims->nx[i];
        sdi.nu = dims->nu[i];
        sdi.np = dims->np[i];

        size += args->submodules.xp[i]->calculate_memory_size(&sdi, args->xp_args[i]);
    }

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



void *ocp_lin_gn_assign_memory(ocp_lin_dims *dims, void *args_, void *raw_memory)
{
    ocp_lin_gn_args *args = (ocp_lin_gn_args *) args_;

    ocp_lin_gn_memory *mem;

    int N = dims->N;

    char *c_ptr = (char *) raw_memory;

    mem = (ocp_lin_gn_memory *) c_ptr;
    c_ptr += sizeof(ocp_lin_gn_memory);

    mem->ls_res_mem = (void **) c_ptr;
    c_ptr += (N+1) * sizeof(void *);

    mem->h_mem = (void **) c_ptr;
    c_ptr += (N+1) * sizeof(void *);

    for (int i=0;i<=N;i++)
    {
        mem->ls_res_mem[i] = args->submodules.ls_res[i]->assign_memory(&ocp_lin_gn_ls_res_dims, args->ls_res_args[i], c_ptr);
        c_ptr += args->submodules.ls_res[i]->calculate_memory_size(&ocp_lin_gn_ls_res_dims, args->ls_res_args[i]);

        mem->h_mem[i] = args->submodules.h[i]->assign_memory(&ocp_lin_gn_h_dims, args->h_args[i], c_ptr);
        c_ptr += args->submodules.h[i]->calculate_memory_size(&ocp_lin_gn_h_dims, args->h_args);
    }

    mem->xp_mem = (void **) c_ptr;
    c_ptr += N * sizeof(void *);

    for (int i=0;i<N;i++)
    {
        sim_dims sdi;
        sdi.num_stages = 3;  // TODO(nielsvd): remove, should be an argument
        sdi.nx = dims->nx[i];
        sdi.nu = dims->nu[i];
        sdi.np = dims->np[i];

        mem->xp_mem[i] = args->submodules.xp[i]->assign_memory(&sdi, args->xp_args[i], c_ptr);
        c_ptr += args->submodules.xp[i]->calculate_memory_size(&sdi, args->xp_args[i]);
    }

    align_char_to(8, &c_ptr);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    assert((char*)raw_memory + ocp_lin_gn_calculate_memory_size(dims, args_) >= c_ptr);

    return (void *)mem;
}



int ocp_lin_gn_calculate_workspace_size(ocp_lin_dims *dims, void *args_)
{
    ocp_lin_gn_args *args = (ocp_lin_gn_args *) args_;

    int N = dims->N;

    int size = sizeof(ocp_lin_gn_workspace);

    size += 2*(N+1)*sizeof(void *);

    for (int i=0;i<=N;i++)
    {
        size += args->submodules.ls_res[i]->calculate_workspace_size(&ocp_lin_gn_ls_res_dims, args->ls_res_args[i]);
        size += args->submodules.h[i]->calculate_workspace_size(&ocp_lin_gn_h_dims, args->h_args[i]);
    }

    size += N*sizeof(void *);

    for (int i=0;i<N;i++)
    {
        sim_dims sdi;
        sdi.num_stages = 3;  // TODO(nielsvd): remove, should be an argument
        sdi.nx = dims->nx[i];
        sdi.nu = dims->nu[i];
        sdi.np = dims->np[i];

        size += args->submodules.xp[i]->calculate_workspace_size(&sdi, args->xp_args[i]);
    }

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



static void *cast_workspace(ocp_lin_dims *dims, void *args_, void *raw_memory)
{
    ocp_lin_gn_args *args = (ocp_lin_gn_args *) args_;

    ocp_lin_gn_workspace *work;

    int N = dims->N;

    char *c_ptr = (char *) raw_memory;

    work = (ocp_lin_gn_workspace *) c_ptr;
    c_ptr += sizeof(ocp_lin_gn_workspace);

    // S_adj_in
    assign_double_ptrs(N, &work->S_adj_in, &c_ptr);
    for (int i=0;i<N;i++)
    {
        assign_double(dims->nx[i+1], &work->S_adj_in[i], &c_ptr);
    }

    // S_forw_in
    assign_double_ptrs(N, &work->S_forw_in, &c_ptr);
    for (int i=0;i<N;i++)
    {
        int nxu = dims->nx[i] + dims->nu[i];
        assign_double(dims->nx[i+1]*nxu, &work->S_forw_in[i], &c_ptr);
    }

    // F
    assign_double_ptrs(N+1, &work->F, &c_ptr);
    for (int i=0;i<=N;i++)
    {
        assign_double(dims->ny[i], &work->F[i], &c_ptr);
    }

    // DF
    assign_double_ptrs(N+1, &work->DF, &c_ptr);
    for (int i=0;i<=N;i++)
    {
        int nxu = dims->nx[i] + dims->nu[i];
        assign_double(dims->ny[i]*nxu, &work->DF[i], &c_ptr);
    }

    // DFT
    assign_double_ptrs(N+1, &work->DFT, &c_ptr);
    for (int i=0;i<=N;i++)
    {
        int nxu = dims->nx[i] + dims->nu[i];
        assign_double(nxu*dims->ny[i], &work->DFT[i], &c_ptr);
    }

    // H
    assign_double_ptrs(N+1, &work->H, &c_ptr);
    for (int i=0;i<=N;i++)
    {
        assign_double(dims->nh[i], &work->H[i], &c_ptr);
    }

    // DH
    assign_double_ptrs(N+1, &work->DH, &c_ptr);
    for (int i=0;i<=N;i++)
    {
        int nxu = dims->nx[i] + dims->nu[i];
        assign_double(dims->nh[i]*nxu, &work->DFT[i], &c_ptr);
    }

    /*********************************
     * submodules
     *********************************/

    work->ls_res_work = (void **) c_ptr;
    c_ptr += (N+1) * sizeof(void *);

    work->h_work = (void **) c_ptr;
    c_ptr += (N+1) * sizeof(void *);

    for (int i=0;i<=N;i++)
    {
        work->ls_res_work[i] = (void *) c_ptr;
        c_ptr += args->submodules.ls_res[i]->calculate_workspace_size(&ocp_lin_gn_ls_res_dims, args->ls_res_args[i]);

        work->h_work[i] = (void *) c_ptr;
        c_ptr += args->submodules.h[i]->calculate_workspace_size(&ocp_lin_gn_h_dims, args->h_args);
    }

    work->xp_work = (void **) c_ptr;
    c_ptr += N * sizeof(void *);

    for (int i=0;i<N;i++)
    {
        sim_dims sdi;
        sdi.num_stages = 3;  // TODO(nielsvd): remove, should be an argument
        sdi.nx = dims->nx[i];
        sdi.nu = dims->nu[i];
        sdi.np = dims->np[i];

        work->xp_work[i] = (void *) c_ptr;
        c_ptr += args->submodules.xp[i]->calculate_workspace_size(&sdi, args->xp_args[i]);
    }

    align_char_to(8, &c_ptr);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    assert((char*)raw_memory + ocp_lin_gn_calculate_workspace_size(dims, args_) >= c_ptr);

    return (void *)work;
}



static void *cast_workspace(ocp_lin_dims *dims, void *args_, void *raw_memory)
{
    return NULL;
}



int ocp_lin_gn(ocp_lin_in *lin_in, ocp_lin_out *lin_out, void *args_, void *memory_, void *work_)
{
    int return_value = ACADOS_SUCCESS;

    ocp_lin_dims *dims = lin_in->dims;

    ocp_lin_gn_args *args = (ocp_lin_gn_args *) args_;
    ocp_lin_gn_memory *mem = (ocp_lin_gn_memory *) memory_;
    ocp_lin_gn_workspace *work = (ocp_lin_gn_workspace *) cast_workspace(dims, args_, work_);

    int N = dims->N;

    // Simulate
    for (int i=0;i<N;i++)
    {
        sim_in in;
        in.nx = dims->nx[i];
        in.nu = dims->nu[i];
        in.np = dims->np[i];
        in.x = lin_in->x[i];
        in.u = lin_in->u[i];
        in.p = lin_in->p[i];
        in.S_adj = work->S_adj_in[i];
        in.S_forw = work->S_forw_in[i];

        // Initialize forward seed
        assert(dims->nx[i] == dims->nx[i+1] && "No support for state transition maps yet.");
        memset(in.S_forw, 0, dims->nx[i+1]*(dims->nx[i]+dims->nu[i])*sizeof(double));
        for (int j=0;j<dims->nx[j];j++)
        {
            in.S_forw[j*dims->nx[i]+j] = 1;
        }

        // Initialize adjoint seed
        for (int j=0;j<dims->nx[i+1];j++)
        {
            in.S_adj[j] = -lin_in->pi[i][j];
        }

        sim_out out;
        sim_info info;
        out.xn = lin_out->xp[i];
        // out.grad = lin_out->grad_pi_xp[i];  // TODO(nielsvd): Find out what grad is...
        out.S_hess = NULL;
        out.S_forw = lin_out->jac_xp[i];
        out.S_adj = lin_out->grad_pi_xp[i];
        out.info = &info;

        // eval
        args->submodules.xp[i]->fun(&in, &out, args->xp_args[i], mem->xp_mem[i], work->xp_work[i]);
    }

    // Compute objective
    for (int i=0;i<=N;i++)
    {
        double *inputs [3];
        inputs[0] = lin_in->x[i];
        inputs[1] = lin_in->u[i];
        inputs[2] = lin_in->p[i];

        bool compute_outputs [3];
        compute_outputs[0] = true;
        compute_outputs[1] = true;
        compute_outputs[2] = false;

        double *outputs [3];
        outputs[0] = work->F[i];
        outputs[1] = work->DF[i];
        outputs[3] = NULL;

        external_function_in in;
        in.inputs = inputs;
        in.compute_output = compute_outputs;

        external_function_out out;
        out.outputs = outputs;

        // eval
        args->submodules.ls_res[i]->fun(&in, &out, args->ls_res_args[i], mem->ls_res_mem[i], work->ls_res_work[i]);
    }

    // Compute path constraints
    for (int i=0;i<=N;i++)
    {
        double *inputs [3];
        inputs[0] = lin_in->x[i];
        inputs[1] = lin_in->u[i];
        inputs[2] = lin_in->p[i];

        bool compute_outputs [3];
        compute_outputs[0] = true;
        compute_outputs[1] = true;
        compute_outputs[2] = false;

        double *outputs [3];
        outputs[0] = work->H[i];
        outputs[1] = work->DH[i];
        outputs[3] = NULL;

        external_function_in in;
        in.inputs = inputs;
        in.compute_output = compute_outputs;

        external_function_out out;
        out.outputs = outputs;

        // eval
        args->submodules.h[i]->fun(&in, &out, args->h_args[i], mem->h_mem[i], work->h_work[i]);
    }

    // Compute Hessian approximation and gradient of cost
    for (int i=0;i<=N;i++)
    {
        int ny = dims->ny[i];
        int nxu = dims->nx[i] + dims->nu[i];

        // Transpose of DF
        for (int_t j = 0; j < nxu; j++) {
            for (int_t k = 0; k < ny; k++)
                work->DFT[i][k * nxu + j] = work->DF[i][j * ny + k];
        }

        // Compute Gauss-Newton Hessian
        for (int_t j = 0; j < nxu*nxu; j++) lin_out->hess_l[i][j] = 0;
        dgemm_nn_3l(nxu, nxu, ny, work->DFT[i],
                    nxu, work->DF[i], ny, lin_out->hess_l[i], nxu);
        // Compute Gauss-Newton gradient
        for (int_t j = 0; j < nxu; j++) lin_out->grad_f[i][j] = 0;
        dgemv_n_3l(nxu, ny, work->DFT[i], nxu, work->F[i],
                   lin_out->grad_f[i]);
    }

    // Write to output

    return return_value;
}
