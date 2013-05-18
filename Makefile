
include makefile.in

SRC=src
TEST=testCPML sine testMain gaussian cpmlfdtd3d
TEST_SRC_DIR=./test/
CPPOBJECT=hpw3d.o fdtd.o InonizationFormula.o inputChecker.o source.o
OBJS=cpml.o $(CPPOBJECT)
TEST_OBJ=sine.o testMain.o testcpml.o gaussian.o openmp.o
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
$(TEST_OBJ): %.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $<	
$(CPPOBJECT): %.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $<	
cpmlfdtd3d.o:cpmlfdtd3d.c
	$(CC) $(CFLAGS) -c $<	
openmp:openmp.o
	$(CXX) -o $@ openmp.o $(CXXFLAGS) $(LIB)
datastruct.h:datastruct.cpp	common.h
fdtd.o:cpml.o
cpml.o:cpml.h cpml.cpp datastruct.h datastruct.cpp common.h
	$(CXX) -c $< $(CXXFLAGS)

# ==========================================
# 3DFormulaTransforming.pdf
# ==========================================
3DFormulaTransforming.pdf:3DFormulaTransforming.tex
	texi2pdf 3DFormulaTransforming.tex
clean:
	-rm -f *.o *.d $(projects)
	sh clean.sh
#	cd $(TEST_SRC_DIR) && make clean
#	cd $(SRC) && make clean
