# for the NVIDIA nvhpc compiler (nvc)
# copy to make.def

CC          = clang
CLINKER     = $(CC)
OPTFLAGS    = -O3 -fopenmp --offload-arch=sm_75
LIBS        = 
PRE         = ./

CFLAGS    = $(OPTFLAGS)

OBJ=o
EXE=
RM=rm -f
