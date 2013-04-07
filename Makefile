#CXX=g++
CXX=icpc
CPPFLAGS=-Wall -g -DDEBUG=3 -DWITH_DENSITY
SRC=src
OBJS=cpml.o test.o datastruct.o fdtd.o InonizationFormula.o

projects=origProgram testCPML hpw3d orig cmain emain dmain tcpml
.PHONY:all clean

all:$(projects)

hpw3d:$(OBJS)
	$(CXX) -o $@ $(OBJS)
origProgram:testMain.o
	$(CXX) -o $@ testMain.o $(CPPFLAGS)
cmain:cmain.o datastruct.o
	$(CXX) -o $@ cmain.o $(CPPFLAGS) datastruct.o
emain:emain.o datastruct.o
	$(CXX) -o $@ emain.o $(CPPFLAGS) datastruct.o
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
cmain.o:$(SRC)/cmain.cpp
	$(CXX) $(CPPFLAGS) -c $< 
emain.o:$(SRC)/emain.cpp
	$(CXX) $(CPPFLAGS) -c $< 
fdtd.o:$(SRC)/fdtd.cpp $(SRC)/fdtd.h cpml.o datastruct.o
	$(CXX) $(CPPFLAGS) -c $< 
cpml.o:$(SRC)/cpml.cpp $(SRC)/cpml.h datastruct.o
	$(CXX) $(CPPFLAGS) -c $< 
cpmld.o:$(SRC)/cpmld.cpp $(SRC)/cpmld.h datastruct.o
	$(CXX) $(CPPFLAGS) -c $< 
dmain:dmain.o datastruct.o cpmld.o
	$(CXX) -o $@ dmain.o $(CPPFLAGS) datastruct.o cpmld.o
dmain.o:$(SRC)/dmain.cpp
	$(CXX) $(CPPFLAGS) -c $< 
tcpml:tcpml.o datastruct.o cpml.o
	$(CXX) -o $@ tcpml.o $(CPPFLAGS) datastruct.o cpml.o
tcpml.o:$(SRC)/tcpml.cpp
	$(CXX) $(CPPFLAGS) -c $< 
%.o:$(SRC)/%.cpp $(SRC)/%.h
	$(CXX) $(CPPFLAGS) -c $< 
clean:
	-rm -f *.o $(projects)
	sh clean.sh
