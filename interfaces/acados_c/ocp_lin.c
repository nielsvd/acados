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

#include "acados_c/ocp_lin.h"

// external
#include <assert.h>
#include <stdlib.h>
#include <string.h>
// acados_c
#include "acados_c/ocp_lin/ocp_lin_gn.h"



void ocp_lin_copy_dims(ocp_lin_dims *dest, ocp_lin_dims *src)
{
    dest->N = src->N;

    for (int ii = 0; ii < src->N + 1; ii++) {
        dest->nx[ii] = src->nx[ii];
        dest->nu[ii] = src->nu[ii];
        dest->np[ii] = src->np[ii];
        dest->ny[ii] = src->ny[ii];
        dest->nb[ii] = src->nb[ii];
        dest->nbx[ii] = src->nbx[ii];
        dest->nbu[ii] = src->nbu[ii];
        dest->ng[ii] = src->ng[ii];
        dest->nh[ii] = src->nh[ii];
        dest->ns[ii] = src->ns[ii];
        dest->nbu[ii] = src->nbu[ii];
        dest->nbx[ii] = src->nbx[ii];
    }

    for (int ii = 0; ii < src->N; ii++) {
        dest->num_stages[ii] = src->num_stages[ii];
    }
}



ocp_lin_dims *create_ocp_lin_dims(int N)
{
    int bytes = ocp_lin_dims_calculate_size(N);

    void *ptr = calloc(1, bytes);

    ocp_lin_dims *dims = assign_ocp_lin_dims(N, ptr);
    dims->N = N;

    return dims;
}



ocp_lin_in *create_ocp_lin_in(ocp_lin_dims *dims)
{
    int bytes = ocp_lin_in_calculate_size(dims);

    void *ptr = calloc(1, bytes);

    ocp_lin_in *in = assign_ocp_lin_in(dims, ptr);

    return in;
}



ocp_lin_out *create_ocp_lin_out(ocp_lin_dims *dims)
{
    int bytes = ocp_lin_out_calculate_size(dims);

    void *ptr = calloc(1, bytes);

    ocp_lin_out *out = assign_ocp_lin_out(dims, ptr);

    return out;
}



int ocp_lin_calculate_args_size(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims)
{
    return fcn_ptrs->calculate_args_size(dims, fcn_ptrs->submodules);
}



void *ocp_lin_assign_args(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims, void *raw_memory)
{
    void *submodules = fcn_ptrs->submodules;
    void *args = fcn_ptrs->assign_args(dims, &submodules, raw_memory);

    fcn_ptrs->initialize_default_args(dims, args);

    return args;
}



void *ocp_lin_create_args(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims)
{
    int bytes = ocp_lin_calculate_args_size(fcn_ptrs, dims);

    void *ptr = calloc(1, bytes);

    void *args = ocp_lin_assign_args(fcn_ptrs, dims, ptr);

    return args;
}



void *ocp_lin_copy_args(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims, void *raw_memory, void *source)
{
    return fcn_ptrs->copy_args(dims, raw_memory, source);
}



int ocp_lin_calculate_size(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims, void *args_)
{
    int bytes = 0;

    bytes += sizeof(ocp_lin_method);

    bytes += sizeof(ocp_lin_method_fcn_ptrs);

    bytes += ocp_lin_dims_calculate_size(dims->N);

    bytes += ocp_lin_calculate_args_size(fcn_ptrs, dims);

    bytes += fcn_ptrs->calculate_memory_size(dims, args_);

    bytes += fcn_ptrs->calculate_workspace_size(dims, args_);

    return bytes;
}

ocp_lin_method *ocp_lin_assign(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims, void *args_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_lin_method *method = (ocp_lin_method *) c_ptr;
    c_ptr += sizeof(ocp_lin_method);

    method->fcn_ptrs = (ocp_lin_method_fcn_ptrs *) c_ptr;
    c_ptr += sizeof(ocp_lin_method_fcn_ptrs);

    method->dims = assign_ocp_lin_dims(dims->N, c_ptr);
    c_ptr += ocp_lin_dims_calculate_size(dims->N);

    method->args = ocp_lin_copy_args(fcn_ptrs, dims, c_ptr, args_);
    c_ptr += ocp_lin_calculate_args_size(fcn_ptrs, dims);

    method->mem = fcn_ptrs->assign_memory(dims, args_, c_ptr);
    c_ptr += fcn_ptrs->calculate_memory_size(dims, args_);

    method->work = (void *) c_ptr;
    c_ptr += fcn_ptrs->calculate_workspace_size(dims, args_);

    assert((char*)raw_memory + ocp_lin_calculate_size(fcn_ptrs, dims, args_) == c_ptr);

    *method->fcn_ptrs = *fcn_ptrs;
    method->fcn_ptrs->submodules = NULL;

    ocp_lin_copy_dims(method->dims, dims);

    return method;
}



ocp_lin_method *ocp_lin_create(ocp_lin_method_fcn_ptrs *fcn_ptrs, ocp_lin_dims *dims, void *args_)
{
    int bytes = ocp_lin_calculate_size(fcn_ptrs, dims, args_);

    void *ptr = calloc(1, bytes);

    ocp_lin_method *method = ocp_lin_assign(fcn_ptrs, dims, args_, ptr);

    return method;
}



int ocp_lin_eval(ocp_lin_method *method, ocp_lin_in *qp_in, ocp_lin_out *qp_out)
{
    return method->fcn_ptrs->fun(qp_in, qp_out, method->args, method->mem, method->work);
}



int ocp_lin_calculate_submodules_size(ocp_lin_method_config *config, ocp_lin_dims *dims)
{
    ocp_lin_method_t method_name = config->lin_method;

    int size;

    switch (method_name)
    {
        case GAUSS_NEWTON:
            size = ocp_lin_gn_calculate_submodules_size(config, dims);
            break;
        default:
            size = 0;
    }

    return size;
}



void *ocp_lin_assign_submodules(ocp_lin_method_config *config, ocp_lin_dims *dims, void *raw_memory)
{
    ocp_lin_method_t method_name = config->lin_method;

    void *submodules;

    switch (method_name)
    {
        case GAUSS_NEWTON:
            submodules = ocp_lin_gn_assign_submodules(config, dims, raw_memory);
            break;
        default:
            submodules = NULL;
    }

    return submodules;
}



int calculate_ocp_lin_method_fcn_ptrs_size(ocp_lin_method_config *config, ocp_lin_dims *dims)
{
    int size = sizeof(ocp_lin_method_fcn_ptrs);

    size += ocp_lin_calculate_submodules_size(config, dims);

    return size;
}



void *assign_ocp_lin_method_fcn_ptrs(ocp_lin_method_config *config, ocp_lin_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *)raw_memory;

    ocp_lin_method_fcn_ptrs *fcn_ptrs = (ocp_lin_method_fcn_ptrs *)c_ptr;
    c_ptr += sizeof(ocp_lin_method_fcn_ptrs);

    set_ocp_lin_method_fcn_ptrs(config, fcn_ptrs);

    fcn_ptrs->submodules = ocp_lin_assign_submodules(config, dims, c_ptr);
    c_ptr += ocp_lin_calculate_submodules_size(config, dims);

    assert((char*)raw_memory + calculate_ocp_lin_method_fcn_ptrs_size(config, dims) == c_ptr);

    return (void *)fcn_ptrs;
}



void *create_ocp_lin_method_fcn_ptrs(ocp_lin_method_config *config, ocp_lin_dims *dims)
{
    int bytes = calculate_ocp_lin_method_fcn_ptrs_size(config, dims);

    void *ptr = malloc(bytes);

    ocp_lin_method_fcn_ptrs *fcn_ptrs = assign_ocp_lin_method_fcn_ptrs(config, dims, ptr);

    return fcn_ptrs;
}



int set_ocp_lin_method_fcn_ptrs(ocp_lin_method_config *config, ocp_lin_method_fcn_ptrs *fcn_ptrs)
{
    int return_value = ACADOS_SUCCESS;
    ocp_lin_method_t method_name = config->lin_method;

    switch (method_name)
    {
        case GAUSS_NEWTON:
            fcn_ptrs->fun = &ocp_lin_gn;
            fcn_ptrs->calculate_args_size = &ocp_lin_gn_calculate_args_size;
            fcn_ptrs->assign_args = &ocp_lin_gn_assign_args;
            fcn_ptrs->copy_args = &ocp_lin_gn_copy_args;
            fcn_ptrs->initialize_default_args = &ocp_lin_gn_initialize_default_args;
            fcn_ptrs->calculate_memory_size = &ocp_lin_gn_calculate_memory_size;
            fcn_ptrs->assign_memory = &ocp_lin_gn_assign_memory;
            fcn_ptrs->calculate_workspace_size = &ocp_lin_gn_calculate_workspace_size;
            break;
        default:
            return_value = ACADOS_FAILURE;
    }

    return return_value;
}
