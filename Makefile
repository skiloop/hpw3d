CXX=/opt/intel/bin/icpc
#CXX=g++
CPPFLAGS=-Wall -g -DDEBUG=3
SRC=src
OBJS=cpml.o test.o datastruct.o fdtd.o InonizationFormula.o

projects=origProgram testCPML hpw3d orig
.PHONY:all clean

all:$(projects)

hpw3d:$(OBJS)
	$(CXX) -o $@ $(OBJS)
origProgram:testMain.o
	$(CXX) -o $@ testMain.o $(CPPFLAGS)
orig:orig.o
	$(CXX) -o $@ orig.o $(CPPFLAGS)
testCPML:testcpml.o cpml.o datastruct.o
	$(CXX) $(CPPFLAGS) -o $@  datastruct.o cpml.o testcpml.o
testcpml.o:$(SRC)/testcpml.cpp
	$(CXX) $(CPPFLAGS) -c $< 
testMain.o:$(SRC)/testMain.cpp
	$(CXX) $(CPPFLAGS) -c $< 
orig.o:$(SRC)/orig.cpp
	$(CXX) $(CPPFLAGS) -c $< 
test.o:$(SRC)/test.cpp
	$(CXX) $(CPPFLAGS) -c $< 
fdtd.o:$(SRC)/fdtd.cpp $(SRC)/fdtd.h cpml.o datastruct.o
	$(CXX) $(CPPFLAGS) -c $< 
cpml.o:$(SRC)/cpml.cpp $(SRC)/cpml.h datastruct.o
	$(CXX) $(CPPFLAGS) -c $< 
%.o:$(SRC)/%.cpp $(SRC)/%.h
	$(CXX) $(CPPFLAGS) -c $< 
clean:
	-rm -f *.o $(projects)
