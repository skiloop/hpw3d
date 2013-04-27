
include makefile.in

SRC=src
TEST=testCPML sine testMain gaussian
TEST_SRC_DIR=./test/
OBJS=cpml.o hpw3d.o fdtd.o InonizationFormula.o #datastruct.o
TEST_OBJ=sine.o testMain.o testcpml.o gaussian.o openmp.o
projects=$(TEST) hpw3d#3DFormulaTransforming.pdf
.PHONY:all clean test

all:hpw3d 3DFormulaTransforming.pdf

hpw3d:$(OBJS)
	$(CXX) -o $@ $(OBJS) $(LIB)


# ==================================================
# test
# ==================================================

test: $(TEST)
gaussian:gaussian.o fdtd.o InonizationFormula.o
	$(CXX) -o $@ $? $(CXXFLAGS) $(LIB)
sine:sine.o fdtd.o InonizationFormula.o
	$(CXX) -o $@ $? $(CXXFLAGS) $(LIB)
testMain:testMain.o
	$(CXX) -o $@ testMain.o $(CXXFLAGS) $(LIB)
testCPML:testcpml.o cpml.o
	$(CXX) $(CXXFLAGS) -o $@  cpml.o testcpml.o $(LIB)
$(TEST_OBJ):
	cd $(TEST_SRC_DIR) && make $@
$(OBJS):
	cd $(SRC) && make $@
openmp:openmp.o
	$(CXX) -o $@ openmp.o $(CXXFLAGS) $(LIB)

# ==========================================
# 3DFormulaTransforming.pdf
# ==========================================
3DFormulaTransforming.pdf:3DFormulaTransforming.tex
	texi2pdf 3DFormulaTransforming.tex
clean:
	-rm -f *.o *.d $(projects)
	sh clean.sh
	cd $(TEST_SRC_DIR) && make clean
	cd $(SRC) && make clean
