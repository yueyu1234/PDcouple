SHELL = /bin/sh

# ------ FILES ------

SRC = FE_Element.cpp \
FE_Mesh.cpp \
Material.cpp \
FE_Engine.cpp\
Matrix.cpp \
MatrixSolver.cpp \
Solver.cpp \


INC = Array2D.h \
Array.h \
CG.h \
CloneVector.h \
DenseMatrix.h \
DenseVector.h \
DiagonalMatrix.h \
FE_Element.h \
FE_Mesh.h \
Material.h \
FE_Engine.h\
MatrixDef.h \
Matrix.h \
MatrixLibrary.h \
MatrixSolver.h \
Quadrature.h \
Solver.h \
SparseMatrix.h \
SparseMatrix-inl.h \
SparseVector.h \
SparseVector-inl.h \
Utility.h \
Vector.h \
table_delaunay.h \

# ------ DEFINITIONS ------

LIB = libatcelastic.a
OBJ =   $(SRC:.cpp=.o)

# ------ SETTINGS ------

# include any MPI settings needed for the ATC library to build with
# the same MPI library that LAMMPS is built with

CC =	        g++
CCFLAGS =       -O -g -I../../src -I../../src/STUBS
ARCHIVE =	ar
ARCHFLAG =	-rc
DEPFLAGS =      -M
LINK =         	g++
LINKFLAGS =	-O
USRLIB =
SYSLIB =

# Defines which names should be treated as rules and not as files.
.PHONY: all clean mrproper

# ------ MAKE PROCEDURE ------

lib: 	$(OBJ)
	$(ARCHIVE) $(ARFLAGS) $(LIB) $(OBJ)

# ------ COMPILE RULES ------

%.o:%.cpp
	$(CC) $(CCFLAGS) -c $<
%.d:%.cpp
	$(CC) $(CCFLAGS) $(DEPFLAGS) $< > $@		

# ------ DEPENDENCIES ------

DEPENDS = $(OBJ:.o=.d)

# ------ CLEAN ------

clean:
	rm *.o $(LIB)
