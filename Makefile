CXX=g++
#CXX=icpc
CC=gcc

# matlab path 
MATPATH=/opt/Matlab/R2011a# for local server
#MATPATH=/opt2/Matlab/R2011a# for servers of college

# Matlab link path
MATLINK=-Wl,-rpath=$(MATPATH)/bin/glnxa64

# Matlab link option
MATLIB=$(MATLINK) -L$(MATPATH)/bin/glnxa64 -lmx -leng

# Matlab link path
MATINC=-I$(MATPATH)/extern/include

SRC=src
CPPFLAGS=-Wall -g -DDEBUG=3 -DWITH_DENSITY #-DMATLAB_SIMULATION $(MATINC)
#CPPFLAGS=-Wall -g -DDEBUG=3 -DWETH_DENSITY -DMATLAB_SIMULATION $(MATINC)
LIB= -lm #$(MATLIB)


TEST=testCPML sine testMain
TEST_SRC_DIR=./test/
OBJS=cpml.o hpw3d.o fdtd.o InonizationFormula.o #datastruct.o
TEST_OBJ=sine.o testMain.o testcpml.o
projects=$(TEST) origProgram testCPML hpw3d orig cmain emain dmain tcpml sine#3DFormulaTransforming.pdf
.PHONY:all clean test

all:hpw3d 

hpw3d:$(OBJS)
	$(CXX) -o $@ $(OBJS) $(LIB)


# ==================================================
# test
# ==================================================

test: $(TEST)
sine:sine.o fdtd.o InonizationFormula.o
	$(CXX) -o $@ $? $(CPPFLAGS) $(LIB)
testMain:testMain.o
	$(CXX) -o $@ testMain.o $(CPPFLAGS) $(LIB)
testCPML:testcpml.o cpml.o
	$(CXX) $(CPPFLAGS) -o $@  cpml.o testcpml.o $(LIB)
$(TEST_OBJ):
	cd $(TEST_SRC_DIR) && make $@
$(OBJS):
	cd $(SRC) && make $@

# ==========================================
# 3DFormulaTransforming.pdf
# ==========================================
3DFormulaTransforming.pdf:3DFormulaTransforming.tex
	texi2pdf 3DFormulaTransforming.tex
clean:
	-rm -f *.o $(projects)
	sh clean.sh
	cd $(TEST_SRC_DIR) && make clean
	cd $(SRC) && make clean
