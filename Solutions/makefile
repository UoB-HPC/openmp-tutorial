#
#  USAGE:
#     make          ... to build the program
#     make test     ... to run the default test case
#
include ../make.def

EXES=pi_spmd_final$(EXE) pi_loop$(EXE) vadd_par$(EXE) \
     jac_solv_dat_reg$(EXE) jac_solv_targ$(EXE)  \
     jac_solv_coal$(EXE) jac_solv_nobr$(EXE) \
     heat_target$(EXE) \
     pi_targ_v1$(EXE) pi_targ_v2$(EXE) pi_targ_v3$(EXE) pi_targ_v4$(EXE)
     # jac_solv_parfor$(EXE) jac_solv_par_for$(EXE)

JAC_PAR_FOR_OBJS  = jac_solv_par_for.$(OBJ) mm_utils.$(OBJ)

JAC_PARFOR_OBJS   = jac_solv_parfor.$(OBJ) mm_utils.$(OBJ)

JAC_COAL_OBJS  = jac_solv_coal.$(OBJ) mm_utils.$(OBJ)

JAC_NOBR_OBJS  = jac_solv_nobr.$(OBJ) mm_utils.$(OBJ)

JAC_DAT_REG_OBJS  = jac_solv_par_dat_reg.$(OBJ) mm_utils.$(OBJ)

JAC_DAT_TARG_OBJS = jac_solv_par_target.$(OBJ) mm_utils.$(OBJ)

all: $(EXES)

jac_solv_par_for$(EXE): $(JAC_PAR_FOR_OBJS) mm_utils.h
	$(CLINKER) $(OPTFLAGS) -o jac_solv_par_for$(EXE) $(JAC_PAR_FOR_OBJS) $(LIBS)

jac_solv_parfor$(EXE): $(JAC_PARFOR_OBJS) mm_utils.h
	$(CLINKER) $(OPTFLAGS) -o jac_solv_parfor$(EXE) $(JAC_PARFOR_OBJS) $(LIBS)

jac_solv_coal$(EXE): $(JAC_COAL_OBJS) mm_utils.h
	$(CLINKER) $(OPTFLAGS) -o jac_solv_coal$(EXE) $(JAC_COAL_OBJS) $(LIBS)

jac_solv_nobr$(EXE): $(JAC_NOBR_OBJS) mm_utils.h
	$(CLINKER) $(OPTFLAGS) -o jac_solv_nobr$(EXE) $(JAC_NOBR_OBJS) $(LIBS)

jac_solv_dat_reg$(EXE): $(JAC_DAT_REG_OBJS) mm_utils.h
	$(CLINKER) $(OPTFLAGS) -o jac_solv_dat_reg$(EXE) $(JAC_DAT_REG_OBJS) $(LIBS)

jac_solv_targ$(EXE): $(JAC_DAT_TARG_OBJS) mm_utils.h
	$(CLINKER) $(OPTFLAGS) -o jac_solv_targ$(EXE) $(JAC_DAT_TARG_OBJS) $(LIBS)

pi_spmd_final$(EXE): pi_spmd_final.$(OBJ)
	$(CLINKER) $(OPTFLAGS) -o pi_spmd_final$(EXE) pi_spmd_final.$(OBJ) $(LIBS)

pi_loop$(EXE): pi_loop.$(OBJ)
	$(CLINKER) $(OPTFLAGS) -o pi_loop$(EXE) pi_loop.$(OBJ) $(LIBS)

vadd_par$(EXE): vadd_par.$(OBJ)
	$(CLINKER) $(OPTFLAGS) -o vadd_par$(EXE) vadd_par.$(OBJ) $(LIBS)

heat_target$(EXE): heat_target.$(OBJ)
	$(CLINKER) $(OPTFLAGS) -o heat_target$(EXE) heat_target.$(OBJ) $(LIBS)

test: $(EXES)
	for i in $(EXES); do \
            $(PRE)$$i; \
        done

clean:
	$(RM) $(EXES) *.$(OBJ) *.ptx *.cub

jac_solv_par_dat_reg.$(OBJ): mm_utils.h
jac_solv_par_for.$(OBJ): mm_utils.h
jac_solv_par_target.$(OBJ): mm_utils.h
jac_solv_parfor.$(OBJ): mm_utils.h
mm_utils.$(OBJ): mm_utils.h

.SUFFIXES:
.SUFFIXES: .c .cpp  .$(OBJ)

.c.$(OBJ):
	$(CC) $(CFLAGS) -c $<

.cpp.$(OBJ):
	$(CC) $(CFLAGS) -c $<
