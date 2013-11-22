
include makefile.in

SRC=src
#TEST=testCPML sine testMain gaussian cpmlfdtd3d
#TEST_SRC_DIR=./test/
OBJS=cpml.o hpw3d.o fdtd.o InonizationFormula.o inputChecker.o source.o Point.o
#CPPOBJECT=sine.o testMain.o testcpml.o gaussian.o openmp.o
#TEST_OBJ=$(CPPOBJECTS) cpmlfdtd3d.o
VPATH = $(SRC):$(TEST_SRC_DIR)
projects=$(TEST) hpw3d#3DFormulaTransforming.pdf
.PHONY:all clean test

all:hpw3d 3DFormulaTransforming.pdf

hpw3d:$(OBJS)
	$(CXX) -o $@ $(OBJS) $(LIB)


# ==================================================
# test
# ==================================================

test: $(TEST)
gaussian:gaussian.o fdtd.o InonizationFormula.o source.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIB)
sine:sine.o fdtd.o InonizationFormula.o source.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIB)
testMain:testMain.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIB)
cpmlfdtd3d:cpmlfdtd3d.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIB)
testCPML:testcpml.o cpml.o
	$(CXX) $(CXXFLAGS) -o $@  $^ $(LIB)
$(CPPOBJECT): %.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $<	
cpmlfdtd3d.o:cpmlfdtd3d.c
	$(CC) -c $< $(CFLAGS)
openmp:openmp.o
	$(CXX) -o $@ openmp.o $(CXXFLAGS) $(LIB)
datastruct.h:datastruct.cpp	common.h
fdtd.o:cpml.o Point.o
cpml.o:cpml.cpp datastruct.h datastruct.cpp common.h cpml.h 
	$(CXX) -c $< $(CXXFLAGS)

# ==========================================
# 3DFormulaTransforming.pdf
# ==========================================
3DFormulaTransforming.pdf:3DFormulaTransforming.tex
	texi2pdf 3DFormulaTransforming.tex
clean:
	-rm -f *.o *.d $(projects) *.aux *.log
	sh clean.sh
#	cd $(TEST_SRC_DIR) && make clean
#	cd $(SRC) && make clean
