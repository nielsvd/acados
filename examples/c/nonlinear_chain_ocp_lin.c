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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// acados
#include <acados_c/ocp_lin.h>
#include <acados_c/options.h>
// Model
#include "chain_model/chain_model.h"
// for casting
#include <acados/ocp_lin/ocp_lin_gn.h>
#include <acados/utils/casadi_wrapper.h>
#include <acados/sim/sim_erk_integrator.h>


int main()
{
    int N = 15;
    double Tf = 3;

    int nmf = 3;  // Number of free masses

    int nx = 6 * nmf;
    int nu = 3;
    int ny = nx+nu;
    int nyN = nx;

    /**********************************
     * Specify dimensions
     **********************************/

    ocp_lin_dims *dims = create_ocp_lin_dims(N);
    for (int i=0; i <= N; i++)
    {
        dims->nx[i] = nx;
        dims->nbx[i] = nmf;
        dims->nbu[i] = nu;
        dims->nb[i] = dims->nbu[i]+dims->nbx[i];;
        dims->ng[i] = 0;
        dims->nh[i] = 0;
        dims->ns[i] = 0;
    }
    for (int i = 0; i < N; i++)
    {
        dims->num_stages[i] = 4;
        dims->nu[i] = nu;
        dims->np[i] = ny;
        dims->ny[i] = ny;
    }
    dims->nu[N] = 0;
    dims->np[N] = nyN;
    dims->ny[N] = nyN;

    /**********************************
     * Configure linearization method
     **********************************/

    ocp_lin_method_config config;
    config.lin_method = GAUSS_NEWTON;
    config.extfun.type = CASADI_WRAPPER;
    config.simsol.sim_solver = ERK;
    config.simsol.extfun.type = CASADI_WRAPPER;

    /**********************************
     * Generate function pointer tree
     **********************************/

    ocp_lin_method_fcn_ptrs *fcn_ptrs = create_ocp_lin_method_fcn_ptrs(&config, dims);

    /**********************************
     * Generate arguments struct
     **********************************/

    void *args = ocp_lin_create_args(fcn_ptrs, dims);

    /**********************************
     * Set arguments
     **********************************/

    printf("\nCompute GGN sensitivities of chain-problem using ERK:\n\n");

    ocp_lin_gn_args *gn_args = (ocp_lin_gn_args *) args;
    for (int i=0;i<N;i++)
    {
        // Least-squares cost
        ((casadi_wrapper_args *) gn_args->ls_res_args[i])->fun = &ls_cost_nm4;
        ((casadi_wrapper_args *) gn_args->ls_res_args[i])->dims = &ls_cost_nm4_work;
        ((casadi_wrapper_args *) gn_args->ls_res_args[i])->sparsity = &ls_cost_nm4_sparsity_out;
        // Inequality constraints
        ((casadi_wrapper_args *) gn_args->h_args[i])->fun = &pathcon_nm4;
        ((casadi_wrapper_args *) gn_args->h_args[i])->dims = &pathcon_nm4_work;
        ((casadi_wrapper_args *) gn_args->h_args[i])->sparsity = &pathcon_nm4_sparsity_out;
        // Simulator options
        sim_erk_integrator_args *erk_args = (sim_erk_integrator_args *) gn_args->xp_args[i];
        erk_args->num_steps = 4;
        erk_args->sens_forw = true;
        erk_args->sens_adj = false;
        erk_args->sens_hess = false;
        erk_args->num_forw_sens = dims->nx[i]+dims->nu[i];
        // Simulator: forward VDE
        ((casadi_wrapper_args *)erk_args->forward_vde_args)->fun = &vde_chain_nm4;
        ((casadi_wrapper_args *)erk_args->forward_vde_args)->dims = &vde_chain_nm4_work;
        ((casadi_wrapper_args *)erk_args->forward_vde_args)->sparsity = &vde_chain_nm4_sparsity_out;
        // Simulator: adjoint VDE
        ((casadi_wrapper_args *)erk_args->adjoint_vde_args)->fun = NULL;
        ((casadi_wrapper_args *)erk_args->adjoint_vde_args)->dims = NULL;
        ((casadi_wrapper_args *)erk_args->adjoint_vde_args)->sparsity = NULL;
        // Simulator: Hessian VDE
        ((casadi_wrapper_args *)erk_args->hess_vde_args)->fun = NULL;
        ((casadi_wrapper_args *)erk_args->hess_vde_args)->dims = NULL;
        ((casadi_wrapper_args *)erk_args->hess_vde_args)->sparsity = NULL;
        
    }
    // Least-squares cost
    ((casadi_wrapper_args *) gn_args->ls_res_args[N])->fun = &ls_costN_nm4;
    ((casadi_wrapper_args *) gn_args->ls_res_args[N])->dims = &ls_costN_nm4_work;
    ((casadi_wrapper_args *) gn_args->ls_res_args[N])->sparsity = &ls_costN_nm4_sparsity_out;
    // Inequality constraints
    ((casadi_wrapper_args *) gn_args->h_args[N])->fun = &pathconN_nm4;
    ((casadi_wrapper_args *) gn_args->h_args[N])->dims = &pathconN_nm4_work;
    ((casadi_wrapper_args *) gn_args->h_args[N])->sparsity = &pathconN_nm4_sparsity_out;


    /**********************************
     * Create linearization method
     **********************************/

    void *method = ocp_lin_create(fcn_ptrs, dims, args);

    free(fcn_ptrs);

    /**********************************
     * Create input-struct
     **********************************/

    ocp_lin_in *in = create_ocp_lin_in(dims);

    // set parameters
    double xref[nx];
    FILE *xn_file = fopen(XN_NM4_FILE, "r");
    for (int i = 0; i < nx; i++)
    if (!fscanf(xn_file, "%lf", &xref[i]))
        break;
    fclose(xn_file);

    for (int i=0;i<=N;i++)
    {
        for (int j=0;j<dims->nx[i];j++)
        {
            in->p[i][j] = xref[j];
            in->x[i][j] = xref[j];
        }
    }


    /**********************************
     * Create output struct
     **********************************/

    ocp_lin_out *out = create_ocp_lin_out(dims);

    /**********************************
     * Evaluate
     **********************************/

    ocp_lin_eval(method, in, out);

    /**********************************
     * Print output
     **********************************/

    for (int i=0;i<N;i++)
    {
        int nx = dims->nx[i];
        int nu = dims->nu[i];
        int nxu = nx + nu;

        printf("\nStage %d:\n=========\n\n",i);

        printf("\nHessian approximation:\n");
        for (int j=0;j<nxu;j++)
        {
            for (int k=0;k<nxu;k++)
            {
                printf("%8.5f ", out->hess_l[i][k*nxu+j]);
            }
            printf("\n");
        }

        printf("\nGradient of cost:\n");
        for (int j=0;j<nxu;j++)
        {
            printf("%8.5f ", out->grad_f[i][j]);
        }
        printf("\n");

        printf("\nSystem dynamics Jacobian:\n");
        for (int j=0;j<dims->nx[i+1];j++)
        {
            for (int k=0;k<nxu;k++)
            {
                printf("%8.5f ", out->jac_xp[i][k*dims->nx[i+1]+j]);
            }
            printf("\n");
        }

        printf("\nInequalities Jacobian:\n");
        for (int j=0;j<dims->nh[i];j++)
        {
            for (int k=0;k<nxu;k++)
            {
                printf("%8.5f ", out->jac_h[i][k*dims->nh[i]+j]);
            }
            printf("\n");
        }

        printf("\nForward simulation:\n");
        for (int j=0;j<dims->nx[i+1];j++)
        {
            printf("%8.5f ", out->xp[i][j]);
        }
        printf("\n");

        printf("\nEvaluation of inequalities:\n");
        for (int j=0;j<dims->nh[i];j++)
        {
            printf("%8.5f ", out->h[i][j]);
        }
        printf("\n");

    }

    return 0;
}