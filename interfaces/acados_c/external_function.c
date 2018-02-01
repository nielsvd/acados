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

#include "acados_c/external_function.h"

//external
#include <stdlib.h>
#include <assert.h>
//acados
#include <acados/utils/mem.h>
//acados_c
#include "acados_c/utils/casadi_wrapper.h"



void external_function_copy_dims(external_function_dims *dest, external_function_dims *src)
{
    dest->num_inputs = src->num_inputs;

    dest->num_outputs = dest->num_outputs;

    for (int i=0; i<src->num_inputs; i++) {
        dest->input_dims[i] = src->input_dims[i];
    }

    for (int i=0; i<src->num_outputs; i++) {
        dest->output_dims[i] = src->output_dims[i];
    }
}



external_function_dims *create_external_function_dims(int num_inputs, int num_outputs)
{
    int bytes = external_function_dims_calculate_size(num_inputs, num_outputs);

    void *ptr = calloc(1, bytes);

    external_function_dims *dims = assign_external_function_dims(num_inputs, num_outputs, ptr);

    return dims;
}



external_function_in *create_external_function_in(external_function_dims *dims)
{
    int bytes = external_function_in_calculate_size(dims);

    void *ptr = acados_malloc(bytes, 1);

    external_function_in *in = assign_external_function_in(dims, ptr);

    return in;
}



external_function_out *create_external_function_out(external_function_dims *dims)
{
    int bytes = external_function_out_calculate_size(dims);

    void *ptr = malloc(bytes);

    external_function_out *out = assign_external_function_out(dims, ptr);

    return out;
}



int external_function_calculate_args_size(external_function_config *config, external_function_dims *dims)
{
    external_function_fcn_ptrs fcn_ptrs;

    set_external_function_fcn_ptrs(config, &fcn_ptrs);

    int size = fcn_ptrs.calculate_args_size(dims, NULL);

    return size;
}



void *external_function_assign_args(external_function_config *config, external_function_dims *dims, void *raw_memory)
{
    external_function_fcn_ptrs fcn_ptrs;

    set_external_function_fcn_ptrs(config, &fcn_ptrs);

    void *args = fcn_ptrs.assign_args(dims, NULL, raw_memory);

    fcn_ptrs.initialize_default_args(args);

    return args;
}



void *external_function_create_args(external_function_config *config, external_function_dims *dims)
{
    int bytes = external_function_calculate_args_size(config, dims);

    void *ptr = malloc(bytes);

    void *args = external_function_assign_args(config, dims, ptr);

    return args;
}



void *external_function_copy_args(external_function_config  *config, external_function_dims *dims, void *raw_memory, void *source)
{
    external_function_t function_type = config->type;

    void *args;

    switch (function_type)
    {
        case CASADI_WRAPPER:
            args = casadi_wrapper_copy_args(config, dims, raw_memory, source);
            break;
    }

    return args;
}



int external_function_calculate_size(external_function_config *config, external_function_dims *dims, void *args_)
{
    external_function_fcn_ptrs fcn_ptrs;

    set_external_function_fcn_ptrs(config, &fcn_ptrs);

    int bytes = 0;

    bytes += sizeof(external_function);

    bytes += sizeof(external_function_fcn_ptrs);

    bytes += external_function_dims_calculate_size(dims->num_inputs, dims->num_outputs);

    bytes += external_function_calculate_args_size(config, dims);

    bytes += fcn_ptrs.calculate_memory_size(dims, args_);

    bytes += fcn_ptrs.calculate_workspace_size(dims, args_);

    return bytes;
}



external_function *external_function_assign(external_function_config *config, external_function_dims *dims, void *args_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    external_function *ext_fun = (external_function *) c_ptr;
    c_ptr += sizeof(external_function);

    ext_fun->fcn_ptrs = (external_function_fcn_ptrs *) c_ptr;
    c_ptr += sizeof(external_function_fcn_ptrs);
    set_external_function_fcn_ptrs(config, ext_fun->fcn_ptrs);

    ext_fun->dims = assign_external_function_dims(dims->num_inputs, dims->num_outputs, c_ptr);
    c_ptr += external_function_dims_calculate_size(dims->num_inputs, dims->num_outputs);
    external_function_copy_dims(ext_fun->dims, dims);

    ext_fun->args = external_function_copy_args(config, dims, c_ptr, args_);
    c_ptr += external_function_calculate_args_size(config, dims);

    ext_fun->mem = ext_fun->fcn_ptrs->assign_memory(dims, args_, c_ptr);
    c_ptr += ext_fun->fcn_ptrs->calculate_memory_size(dims, args_);

    ext_fun->work = (void *) c_ptr;
    c_ptr += ext_fun->fcn_ptrs->calculate_workspace_size(dims, args_);

    assert((char*)raw_memory + external_function_calculate_size(config, dims, args_) == c_ptr);

    return ext_fun;
}



external_function *external_function_create(external_function_config *config, external_function_dims *dims, void *args_)
{
    int bytes = external_function_calculate_size(config, dims, args_);

    void *ptr = malloc(bytes);

    external_function *ext_fun = external_function_assign(config, dims, args_, ptr);

    return ext_fun;
}



int external_function_eval(external_function *ext_fun, external_function_in *ef_in, external_function_out *ef_out)
{
    return ext_fun->fcn_ptrs->fun(ef_in, ef_out, ext_fun->args, ext_fun->mem, ext_fun->work);
}



int set_external_function_fcn_ptrs(external_function_config *config, external_function_fcn_ptrs *fcn_ptrs)
{
    int return_value = ACADOS_SUCCESS;

    external_function_t function_type = config->type;

    switch (function_type)
    {
        case CASADI_WRAPPER:
            fcn_ptrs->fun = &casadi_wrapper;
            fcn_ptrs->calculate_args_size = &casadi_wrapper_calculate_args_size;
            fcn_ptrs->assign_args = &casadi_wrapper_assign_args;
            fcn_ptrs->initialize_default_args = &casadi_wrapper_initialize_default_args;
            fcn_ptrs->calculate_memory_size = &casadi_wrapper_calculate_memory_size;
            fcn_ptrs->assign_memory = &casadi_wrapper_assign_memory;
            fcn_ptrs->calculate_workspace_size = &casadi_wrapper_calculate_workspace_size;
            break;
        default:
            return_value = ACADOS_FAILURE;
    }

    return return_value;
}