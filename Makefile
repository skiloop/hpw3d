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
#CPPFLAGS=-Wall -g -DDEBUG=3 -DWITH_DENSITY -DMATLAB_SIMULATION $(MATINC)
#LIB= $(MATLIB)

OBJS=cpml.o test.o  fdtd.o InonizationFormula.o #datastruct.o

projects=origProgram testCPML hpw3d orig cmain emain dmain tcpml #3DFormulaTransforming.pdf
.PHONY:all clean

all:hpw3d

hpw3d:$(OBJS)
	$(CXX) -o $@ $(OBJS) $(LIB)
origProgram:testMain.o
	$(CXX) -o $@ testMain.o $(CPPFLAGS) $(LIB)
cmain:cmain.o datastruct.o
	$(CXX) -o $@ cmain.o $(CPPFLAGS) datastruct.o $(LIB)
emain:emain.o datastruct.o
	$(CXX) -o $@ emain.o $(CPPFLAGS) datastruct.o $(LIB)
orig:orig.o
	$(CXX) -o $@ orig.o $(CPPFLAGS) $(LIB)
testCPML:testcpml.o cpml.o datastruct.o
	$(CXX) $(CPPFLAGS) -o $@  datastruct.o cpml.o testcpml.o $(LIB)
testcpml.o:$(SRC)/testcpml.cpp
	$(CXX) $(CPPFLAGS) -c $< 
testMain.o:$(SRC)/testMain.cpp
	$(CXX) $(CPPFLAGS) -c $< 
orig.o:$(SRC)/orig.cpp
	$(CXX) $(CPPFLAGS) -c $< 
test.o:$(SRC)/test.cpp fdtd.o 
	$(CXX) $(CPPFLAGS) -c $< 
cmain.o:$(SRC)/cmain.cpp
	$(CXX) $(CPPFLAGS) -c $< 
emain.o:$(SRC)/emain.cpp
	$(CXX) $(CPPFLAGS) -c $< 
fdtd.o:$(SRC)/fdtd.cpp $(SRC)/fdtd.h cpml.o $(SRC)/datastruct.h
	$(CXX) $(CPPFLAGS) -c $< 
cpml.o:$(SRC)/cpml.cpp $(SRC)/cpml.h $(SRC)/datastruct.h
	$(CXX) $(CPPFLAGS) -c $< 
cpmld.o:$(SRC)/cpmld.cpp $(SRC)/cpmld.h $(SRC)/datastruct.h
	$(CXX) $(CPPFLAGS) -c $< 
dmain:dmain.o datastruct.o cpmld.o
	$(CXX) -o $@ dmain.o $(CPPFLAGS) datastruct.o cpmld.o $(LIB)
dmain.o:$(SRC)/dmain.cpp
	$(CXX) $(CPPFLAGS) -c $< 
tcpml:tcpml.o datastruct.o cpml.o
	$(CXX) -o $@ tcpml.o $(CPPFLAGS) datastruct.o cpml.o $(LIB)
tcpml.o:$(SRC)/tcpml.cpp
	$(CXX) $(CPPFLAGS) -c $< 
datastruct.o:$(SRC)/datastruct.cpp $(SRC)/datastruct.h
	$(CXX) $(CPPFLAGS) -c $< 
InonizationFormula.o:$(SRC)/InonizationFormula.cpp $(SRC)/InonizationFormula.h
	$(CXX) $(CPPFLAGS) -c $< 
.cpp.o:
	$(CC) $(CPPFLAGS) -c $(SRC)/$< 
.c.o:$(SRC)/%.c $(SRC)/%.h
	$(CC) $(CPPFLAGS) -c $< 
# ==========================================
# 3DFormulaTransforming.pdf
# ==========================================
3DFormulaTransforming.pdf:3DFormulaTransforming.tex
	texi2pdf 3DFormulaTransforming.tex
clean:
	-rm -f *.o $(projects)
	sh clean.sh
