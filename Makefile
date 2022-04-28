# Makefile for various platforms
# Execute using Build csh-script only!
# Used together with Perl scripts in SRC/SCRIPT 
# (C) 2005 Marat Khairoutdinov
#------------------------------------------------------------------
# uncomment to disable timers:
#
#NOTIMERS=-DDISABLE_TIMERS
#-----------------------------------------------------------------

SAM = SAM_$(ADV_DIR)_$(SGS_DIR)_$(RAD_DIR)_$(MICRO_DIR)

# Determine platform 
PLATFORM := $(shell uname -s)


#----------------------------------------------------------------------
# Linux, Portland Group Compiler
#

# Yunyan:
#
ifeq ($(PLATFORM),Linux)
#
#  # Default compiler flags
  FFLAGS =  -Kieee -fastsse #-g -C -Ktrap=fp #
#
#  ifeq ($(HOSTNAME),hopper)
#
# # Determine platform
    PROCESSOR := $(shell uname -p)
#
# # DEFAULT MPI LOCATION AND COMPILER

    LIB_NETCDF = $(NETCDF_DIR)/lib
    INC_NETCDF = $(NETCDF_DIR)/include

    LIB_MPI = $(MPICH_DIR)/lib
    INC_MPI = $(MPICH_DIR)/include
    MPIF90 = ftn
#  endif

  FF77 = ${MPIF90} -c
  FF90 = ${MPIF90} -c
  CC = gcc -c -DLINUX -g

  LD = ${MPIF90} ${FFLAGS}
  FFLAGS += -I${INC_NETCDF}
  LDFLAGS += -L${LIB_NETCDF} -lnetcdf

#
## end Yunyan:


endif


#----------------------------------
#----------------------------------------------
# you dont need to edit below this line


#compute the search path
dirs := . $(shell cat Filepath)
VPATH    := $(foreach dir,$(dirs),$(wildcard $(dir))) 

.SUFFIXES:
.SUFFIXES: .f .f90 .c .o



all: $(SAM_DIR)/$(SAM)


SOURCES   := $(shell cat Srcfiles)

Depends: Srcfiles Filepath
	$(SAM_SRC)/SCRIPT/mkDepends Filepath Srcfiles > $@

Srcfiles: Filepath
	$(SAM_SRC)/SCRIPT/mkSrcfiles > $@

OBJS      := $(addsuffix .o, $(basename $(SOURCES))) 

$(SAM_DIR)/$(SAM): $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)


.f90.o:
	${FF90}  ${FFLAGS} $<
.f.o:
	${FF77}  ${FFLAGS} $<
.c.o:
	${CC}  ${CFLAGS} -I$(SAM_SRC)/TIMING $(NOTIMERS) $<



include Depends



clean: 
	rm ./OBJ/*


