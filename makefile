#
#  USAGE:
#     make          ... to build the program
#     make test     ... to run the default test case
#

include make.def

EXES= pi$(EXE) jac_solv$(EXE) vadd$(EXE)

JAC_OBJS  = jac_solv.$(OBJ) mm_utils.$(OBJ)

all: $(EXES)

jac_solv$(EXE): $(JAC_OBJS) mm_utils.h
	$(CLINKER) $(CFLAGS) -o jac_solv$(EXE) $(JAC_OBJS) $(LIBS)

pi$(EXE): pi.$(OBJ)
	$(CLINKER) $(OPTFLAGS) -o pi pi.$(OBJ) $(LIBS)

vadd$(EXE): vadd.$(OBJ)
	$(CLINKER) $(OPTFLAGS) -o vadd vadd.$(OBJ) $(LIBS)


test: $(EXES)
	for i in $(EXES); do \
            $(PRE)$$i; \
        done

clean:
	$(RM) $(EXES) *.$(OBJ) *.ptx *.cub

jac_solv.$(OBJ) : mm_utils.h

.SUFFIXES:
.SUFFIXES: .c .cpp .$(OBJ)

.c.$(OBJ):
	$(CC) $(CFLAGS) -c $<

.cpp.$(OBJ):
	$(CC) $(CFLAGS) -c $<
