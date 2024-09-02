LISROOT     = usr/local
SRC         = src
INC         = include
DIR_IVP     = ${INC}/IVPs
DIR_VSBF    = ${SRC}/VSBDF
OBJ_SOLVER         = IVP_Solver
SBDF_Solver = SBDF_Solver
COMP_TEST = Comparison_Test


IVP=IVP_ODE
IVP1=IVP_ODE_combustion
IVP2=IVP_ODE_advdiff1d
IVP3=IVP_ODE_cusp
IVP4=IVP_ODE_RD
IVP5=IVP_ODE_heat
IVP6=IVP_ODE_brusselator
IVP7=IVP_ODE_brusselator2D
IVP8=IVP_ODE_simpleadvdiff1d
ALL_IVPs=${DIR_IVP}/$(IVP1).h ${DIR_IVP}/$(IVP2).h ${DIR_IVP}/$(IVP3).h ${DIR_IVP}/$(IVP4).h ${DIR_IVP}/$(IVP5).h ${DIR_IVP}/$(IVP6).h ${DIR_IVP}/$(IVP7).h ${DIR_IVP}/${IVP8}.h ${DIR_IVP}/$(IVP).h

CCLINKLIBS = -L$(LISROOT)/lib  -llis -lopenblas 

# C++ Compiler:
CC=g++

# General compiler flags:
C_FLAGS= -I. -I$(LISROOT)/include/   -Wall  -O3 -m64 -ffast-math -fomit-frame-pointer -fopenmp


all: ${SRC}/$(OBJ_SOLVER).cpp ${SRC}/$(COMP_TEST).cpp   ${SBDF_Solver}.o ${ALL_IVPs}   
	$(CC)  -I${DIR_IVP} -I${INC} ${SRC}/$(OBJ_SOLVER).cpp ${SBDF_Solver}.o $(C_FLAGS) -o $(OBJ_SOLVER) $(CCLINKLIBS) 
	$(CC)  -I${DIR_IVP} -I${INC} ${SRC}/$(COMP_TEST).cpp ${SBDF_Solver}.o $(C_FLAGS) -o $(COMP_TEST) $(CCLINKLIBS)    

.PHONY: clean

clean:
	rm  $(OBJ_Solver) ${SBDF_Solver}.o 
cleanout:
	rm  *.txt 

	
	
${SBDF_Solver}.o: ${DIR_VSBF}/${SBDF_Solver}.cpp   ${INC}/${SBDF_Solver}.h
	$(CC)  -I${DIR_IVP} -I${INC} ${DIR_VSBF}/${SBDF_Solver}.cpp  $(C_FLAGS) -c     
 

