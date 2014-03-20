
include makefile.in

EXCUTABLE:=hpw3d
PDFS=3DFormulaTransforming.pdf 1DFormulaTransforming.pdf
SRC_DIR:=src
TEST_SRC_DIR:=test
VPATH = $(SRC_DIR):$(TEST_SRC_DIR)
SOURCES=$(shell find $(SRC_DIR) -name "*.cpp")
OBJS:=$(patsubst $(SRC_DIR)/%.cpp,%.o,$(SOURCES))
DEPS:=$(patsubst $(SRC_DIR)/%.cpp,%.d,$(SOURCES))

#We don't need to clean up when we're making these targets
NODEPS:=clean tags svn
TEST:=sourceTest data3d
PROJECTS=$(TEST) $(EXCUTABLE) #3DFormulaTransforming.pdf
.PHONY:all clean test objs veryclean rebuild deps
	
all:$(PROJECTS) $(PDFS)

deps:$(DEPS)

objs:$(OBJS)

test:$(TEST)

#This is the rule for creating the dependency files
%.d:$(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -MM -MT $*.o $< -MF $@

%.o:%.cpp %.d
	$(CXX) $(CXXFLAGS) -o $@ -c $<

#Don't create dependencies when we're cleaning, for instance
ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
#Chances are, these files don't exist.  GMake will create them and
#clean up automatically afterwards
-include $(DEPS)
endif

$(EXCUTABLE):$(OBJS)
	$(CXX) -o $@ $(OBJS) $(LIB)
sourceTest:sourceTest.o SineWaveSource.o sourceType.o Point.o GaussianWaveSource.o CosineGaussianWave.o
	$(CXX) -o $@ $^ $(LIB)
data3d:data3d.o Point.o
	$(CXX) -o $@ $^ $(LIB)
# ==========================================
# 3DFormulaTransforming.pdf
# ==========================================
%.pdf:%.tex
	for f in $<; \
	   do \
	   $(PDFLATEX) $$f; \
	   done
clean:
	-rm -f $(DEPS) $(OBJS) $(PROJECTS) *.aux *.log *.o
veryclean:clean
	sh clean.sh

rebuild:veryclean all

