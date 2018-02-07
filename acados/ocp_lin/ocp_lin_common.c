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

#include "acados/ocp_lin/ocp_lin_common.h"

#include <stdlib.h>



int ocp_lin_dims_calculate_size(int N)
{
    int size = sizeof(ocp_lin_dims);

    size += 9*(N+1)*sizeof(int);

    size += 8;  // initial align

    return size;
}



ocp_lin_dims *assign_ocp_lin_dims(int N, void *raw_memory)
{
    char *c_ptr = (char *)raw_memory;

    // initial align
    align_char_to(8, &c_ptr);

    // struct
    ocp_lin_dims *dims = (ocp_lin_dims *) c_ptr;
    c_ptr += sizeof(ocp_lin_dims);

    // nx
    assign_int(N + 1, &dims->nx, &c_ptr);
    // nu
    assign_int(N + 1, &dims->nu, &c_ptr);
    // np
    assign_int(N + 1, &dims->np, &c_ptr);
    // nb
    assign_int(N + 1, &dims->nb, &c_ptr);
    // nbx
    assign_int(N + 1, &dims->nbx, &c_ptr);
    // nbu
    assign_int(N + 1, &dims->nbu, &c_ptr);
    // ng
    assign_int(N + 1, &dims->ng, &c_ptr);
    // nh
    assign_int(N + 1, &dims->nh, &c_ptr);
    // ns
    assign_int(N + 1, &dims->ns, &c_ptr);

    // N
    dims->N = N;

    assert((char *) raw_memory + ocp_lin_dims_calculate_size(N) >= c_ptr);

    return dims;
}



int ocp_lin_in_calculate_size(ocp_lin_dims *dims)
{
    int size = sizeof(ocp_lin_in);
    int N = dims->N;

    // x
    size += (N+1) * sizeof(double *);

    // u

    // p

    // pi

    // lam

}



ocp_lin_in *assign_ocp_lin_in(ocp_lin_dims *dims, void *raw_memory)
{
    return NULL;
}



int ocp_lin_out_calculate_size(ocp_lin_dims *dims)
{
    return 0;
}



ocp_lin_out *assign_ocp_lin_out(ocp_lin_dims *dims, void *raw_memory)
{
    return NULL;
}
