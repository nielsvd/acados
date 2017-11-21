import matplotlib.pyplot as plt
from numpy import array, diag
from scipy.linalg import block_diag

from casadi import SX, Function, vertcat
from acados import ocp_nlp_function, ocp_nlp_ls_cost, ocp_nlp, ocp_nlp_solver
from models import chen_model

N = 10
ode_fun, nx, nu = chen_model()
nlp = ocp_nlp({'N': N, 'nx': nx, 'nu': nu, 'ng': N*[1] + [0]})

# ODE Model
step = 0.1
nlp.set_model(ode_fun, step)

# Cost function
Q = diag([1.0, 1.0])
R = 1e-2
x = SX.sym('x', nx)
u = SX.sym('u', nu)
u_N = SX.sym('u', 0)
f = ocp_nlp_function(Function('ls_cost', [x, u], [vertcat(x, u)]))
f_N = ocp_nlp_function(Function('ls_cost_N', [x, u_N], [x]))

stage_costs = []#N*[ocp_nlp_ls_cost(f)] + [ocp_nlp_ls_cost(f_N)]
for i in range(N):
    stage_costs += [ocp_nlp_ls_cost(f)]
    stage_costs[i].ls_cost_matrix = block_diag(Q, R)
stage_costs += [ocp_nlp_ls_cost(f_N)]
stage_costs[N].ls_cost_matrix = Q

nlp.set_cost(stage_costs)

# Constraints
g = ocp_nlp_function(Function('path_constraint', [x, u], [u]))
g_N = ocp_nlp_function(Function('path_constraintN', [x, u], [SX([])]))
nlp.set_path_constraints(N*[g] + [g_N])
for i in range(N):
    nlp.lg[i] = -0.5
    nlp.ug[i] = +0.5

solver = ocp_nlp_solver('sqp', nlp, {'integrator_steps': 2, 'qp_solver':'hpipm', 'sensitivity_method': 'gauss-newton'})

# Simulation
STATES = [array([0.1, 0.1])]
CONTROLS = []
for i in range(20):
    state_traj, control_traj = solver.evaluate(STATES[-1])
    STATES += [state_traj[1]]
    CONTROLS += [control_traj[0]]

plt.ion()
plt.plot([x[0] for x in STATES], [x[1] for x in STATES])
plt.axis('equal')
