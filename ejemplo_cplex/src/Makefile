CPLEX_LIB       = /opt/ibm/ILOG/CPLEX_Studio124/cplex/lib/x86-64_sles10_4.1/static_pic/libcplex.a /opt/ibm/ILOG/CPLEX_Studio124/cplex/lib/x86-64_sles10_4.1/static_pic/libilocplex.a
LP_LIBS         = $(CPLEX_LIB)
INC_DIR		= -I/opt/ibm/ILOG/CPLEX_Studio124/cplex/include/ -I/opt/ibm/ILOG/CPLEX_Studio124/concert/include/

CCC = g++

CCFLAGS = -m64 -O2 -fPIC -fexceptions -DNDEBUG -DIL_STD

LDFLAGS = $(LP_LIBS) -lc -lm

SRCPRI = dieta.cpp 

SOURCES = $(SRCPRI)

OBJECTS = $(SOURCES:.cpp=.o)

all: $(OBJECTS) 
	$(CCC) -o dieta $(OBJECTS) $(CCFLAGS) $(LDFLAGS) -lpthread -lm -lnsl

#
# produce the .o files from the .cc and .c files
#

%.o: %.cpp
	$(CCC) $(CCFLAGS) $(INC_DIR) -c $<
