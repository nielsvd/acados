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

#include "acados_c/ocp_qp/ocp_qp_qpdunes.h"



void *ocp_qp_qpdunes_copy_args(ocp_qp_solver_config *config, ocp_qp_dims *dims, void *raw_memory, void *source_)
{
    ocp_qp_qpdunes_args *source = (ocp_qp_qpdunes_args *) source_;
    ocp_qp_qpdunes_args *dest;

    dest = (ocp_qp_qpdunes_args *) ocp_qp_qpdunes_assign_args(dims, raw_memory);

    dest->options = source->options;
    dest->stageQpSolver = source->stageQpSolver;
    dest->warmstart = source->warmstart;
    dest->isLinearMPC = source->isLinearMPC;

    return (void *)dest;
}