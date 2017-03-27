################################################################################
#                               AtC_Elasticity
#
#     File        : Makefile
#     Description : This archive makefile using to generate the executable to
#                   sub-project Model
#     Author      : Fabiano F. Bargos
#     Date        : Sept/28/2011
#
################################################################################


# Defines the commands mainly.
GPP    = g++
AR     = ar
RM     = rm
CD     = cd
CP     = cp
RANLIB = ranlib


# Defines flags, defined symbols to preprocessor and the libraries that
#should be used by the compiler.
CFLAGS     = -w -ansi -pedantic -g -I../../lammps-14Aug11/src -I../../lammps-14Aug11/src/STUBS
DEFINES    = -D"__UNIX__" -D"__32BITS__" -D_REENTRANT
LIBS_FLAGS = -latcelastic -llmp_serial_debug
SRCS_MAINS = simple.cpp


# MPI library, can be src/STUBS dummy lib
# INC = path for mpi.h, MPI compiler settings
# PATH = path for MPI library
# LIB = name of MPI library

MPI_INC =   -I../../lammps-14Aug11/src/STUBS
MPI_PATH =	
MPI_LIB =	../../lammps-14Aug11/src/STUBS/libmpi.a

EXTRA_INC = $(MPI_INC)
EXTRA_LIB = $(MPI_LIB)

OBJS_MAINS = $(addsuffix .o, $(basename $(SRCS_MAINS)))


# Creates the variables to store source codes (*.cpp),
# binaries code (*.o), includes path and libraries path used.
SRCS     = $(wildcard *.cpp) $(wildcard */*.cpp)
SOURCES  = $(filter-out $(SRCS_MAINS), $(SRCS))
OBJECTS  = $(addsuffix .o, $(basename $(SOURCES)))
INCLUDES = -I$(PWD)/../ -I$(PWD)/../../lammps-14Aug11/src/ -I$(PWD)/../../lammps-14Aug11/src/STUBS      \

LIBS_DIR = -L$(PWD)/../ -L$(PWD)/../../lammps-14Aug11/src/ -L$(PWD)/../../lammps-14Aug11/src/STUBS

# Defines variable to store the path directory.
# Besides, stores the root directory path
DIR_ATC       = $(PWD)/../
DIR_LMP       = $(PWD)/../

# Creates variable to store the path of the static libraries
LIBRARY_ATC = $(PWD)/../libatcelastic.a
LIBRARY_LMP = $(PWD)/../../lammps-14Aug11/src/liblmp_serial_debug.a
LIBRARIES   = $(LIBRARY_ATC) $(LIBRARY_LMP)

# Defines the variable that represent the executable name of the sub-project
EXEC = simple

# Defines which names should be treated as rules and not as files.
.PHONY: all clean mrproper

# Creating libraries, objects and executable.
all: $(LIBRARIES) $(SOURCES) $(OBJECTS) \
$(EXEC)

# Create library if it is necessary.
$(LIBRARIES):
	$(CD) $(DIR_ATC)       && $(MAKE)

# Creating executable
$(EXEC): $(OBJECTS) simple.o
	@echo '> Creating $(EXEC)'
	$(GPP) $(LIBS_DIR)  -o $(EXEC)  \
simple.o $(OBJECTS) $(LIBS_FLAGS) $(EXTRA_LIB)

# Generate s binary codes (*.o) from the all sources codes (*.cpp)
%.o: %.cpp
	@echo '> Creating $@'
	@$(GPP) -c $< -o $@ $(INCLUDES) \
$(DEFINES) $(CFLAGS) $(EXTRA_INC)

# Deletes the all binaries code (*.o)
clean:
	@$(CD) $(DIR_ATC)       && $(MAKE) clean
	rm -rf *.o

# Deletes the library and binaries code.
mrproper:
	@$(CD) $(DIR_ATC)       && $(MAKE) clean
	rm -rf *.exe
