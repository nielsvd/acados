# Author: Gianluca Frison



include ../../Makefile.rule



LIBS = -lacados_c -lacore -lhpipm -lqpOASES_e

ifeq ($(ACADOS_WITH_HPMPC), 1)
LIBS+=-lhpmpc
endif
ifeq ($(ACADOS_WITH_QPDUNES), 1)
CCC = $(CXX)
LIBS+= -lqpdunes -lqpOASES
else
CCC = $(CC)
endif
ifeq ($(ACADOS_WITH_QORE), 1)
LIBS+=-lqore_dense
endif
LIBS+=-lblasfeo -lm -lblas

# Comment this out to enable using gprof
# CFLAGS  += -pg
# LDFLAGS += -pg

EXAMPLES =
EXAMPLES += mass_spring_example
EXAMPLES += sim_crane
EXAMPLES += mass_spring_pcond_hpipm_split
EXAMPLES += mass_spring_fcond_hpipm_split
EXAMPLES += mass_spring_offline_fcond_qpoases_split
# EXAMPLES += nonlinear_chain_ocp_nlp
EXAMPLES += nonlinear_chain_ocp_lin

examples: $(EXAMPLES)



RUN_EXAMPLES =
RUN_EXAMPLES += run_mass_spring_example
RUN_EXAMPLES += run_sim_crane
RUN_EXAMPLES += run_mass_spring_pcond_hpipm_split
RUN_EXAMPLES += run_mass_spring_fcond_hpipm_split
RUN_EXAMPLES += run_mass_spring_offline_fcond_qpoases_split
# RUN_EXAMPLES += run_nonlinear_chain_ocp_nlp
RUN_EXAMPLES += run_nonlinear_chain_ocp_lin

run_examples: $(RUN_EXAMPLES)



CRANE_OBJS =
CRANE_OBJS += crane_model/jac_model.o
CRANE_OBJS += crane_model/ode_model.o
CRANE_OBJS += crane_model/vde_adj_model.o
CRANE_OBJS += crane_model/vde_forw_model.o
CRANE_OBJS += crane_model/vde_hess_model.o
CRANE_OBJS += crane_model/impl_ode.o
CRANE_OBJS += crane_model/impl_jac.o

sim_crane: $(CRANE_OBJS) sim_crane.o
	$(CCC) -o sim_crane.out sim_crane.o $(CRANE_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_crane build complete."
	@echo
	


mass_spring_example: mass_spring_example.o
	$(CCC) -o mass_spring_example.out mass_spring_example.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_example build complete."
	@echo



mass_spring_pcond_hpipm_split: mass_spring_pcond_hpipm_split.o
	$(CCC) -o mass_spring_pcond_hpipm_split.out mass_spring_pcond_hpipm_split.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_pcond_hpipm_split build complete."
	@echo



mass_spring_fcond_hpipm_split: mass_spring_fcond_hpipm_split.o
	$(CCC) -o mass_spring_fcond_hpipm_split.out mass_spring_fcond_hpipm_split.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_fcond_hpipm_split build complete."
	@echo


mass_spring_offline_fcond_qpoases_split: mass_spring_offline_fcond_qpoases_split.o
	$(CCC) -o mass_spring_offline_fcond_qpoases_split.out mass_spring_offline_fcond_qpoases_split.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_offline_fcond_qpoases_split build complete."
	@echo


CHAIN_OBJS =
# 2 masses
CHAIN_OBJS += chain_model/vde_chain_nm2.o
CHAIN_OBJS += chain_model/jac_chain_nm2.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm2.o
CHAIN_OBJS += chain_model/vde_hess_chain_nm2.o
CHAIN_OBJS += chain_model/ls_cost_nm2.o
CHAIN_OBJS += chain_model/ls_costN_nm2.o
CHAIN_OBJS += chain_model/pathcon_nm2.o
CHAIN_OBJS += chain_model/pathconN_nm2.o
# 3 masses
CHAIN_OBJS += chain_model/vde_chain_nm3.o
CHAIN_OBJS += chain_model/jac_chain_nm3.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm3.o
CHAIN_OBJS += chain_model/vde_hess_chain_nm3.o
CHAIN_OBJS += chain_model/ls_cost_nm3.o
CHAIN_OBJS += chain_model/ls_costN_nm3.o
CHAIN_OBJS += chain_model/pathcon_nm3.o
CHAIN_OBJS += chain_model/pathconN_nm3.o
# 4 masses
CHAIN_OBJS += chain_model/vde_chain_nm4.o
CHAIN_OBJS += chain_model/jac_chain_nm4.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm4.o
CHAIN_OBJS += chain_model/vde_hess_chain_nm4.o
CHAIN_OBJS += chain_model/ls_cost_nm4.o
CHAIN_OBJS += chain_model/ls_costN_nm4.o
CHAIN_OBJS += chain_model/pathcon_nm4.o
CHAIN_OBJS += chain_model/pathconN_nm4.o
# CHAIN_OBJS += nonlinear_chain_ocp_nlp.o

nonlinear_chain_ocp_nlp: $(CHAIN_OBJS) nonlinear_chain_ocp_nlp.o
	# $(CCC) -o nonlinear_chain_ocp_nlp.out nonlinear_chain_ocp_nlp.o $(CHAIN_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example nonlinear_chain_ocp_nlp build complete."
	@echo



nonlinear_chain_ocp_lin: $(CHAIN_OBJS) nonlinear_chain_ocp_lin.o
	$(CCC) -o nonlinear_chain_ocp_lin.out nonlinear_chain_ocp_lin.o $(CHAIN_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example nonlinear_chain_ocp_lin build complete."
	@echo



run_sim_crane:
	./sim_crane.out



run_mass_spring_example:
	./mass_spring_example.out



run_mass_spring_pcond_hpipm_split:
	./mass_spring_pcond_hpipm_split.out



run_mass_spring_fcond_hpipm_split:
	./mass_spring_fcond_hpipm_split.out



run_mass_spring_offline_fcond_qpoases_split:
	./mass_spring_offline_fcond_qpoases_split.out



run_nonlinear_chain_ocp_nlp:
	./nonlinear_chain_ocp_nlp.out

run_nonlinear_chain_ocp_lin:
	./nonlinear_chain_ocp_lin.out



clean:
	rm -f *.o
	rm -f *.out

deep_clean: clean
	rm -f chain_model/*.o

