DEPDIR = .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td

CXX = g++ -std=c++11
include makefile.local
LIBS = -lmpfr -lgmp -lgsl -lgslcblas -lm -lglpk -lboost_iostreams -lboost_serialization


#CFLAGS += -Dno_output
#CFLAGS = -I . -I $(GMP_HOME) -g -O3
#CFLAGS += -Denablelog
#CXXFLAGS = -I . -I $(GMP_HOME) -g -O3
#LINK_FLAGS = -g -L $(GMP_LIB_HOME) -Wl,-rpath,$(GMP_LIB_HOME) 

#no profiling
#CFLAGS = -I . -g -O3
#CXXFLAGS = -I . -g -O3

#profiling
CFLAGS = -I . -g -O2
CXXFLAGS = -I . -g -O2

#CFLAGS += -pg
CFLAGS += -Denablelog
CXXFLAGS += -Denablelog
CFLAGS += -Ddo_checks
CXXFLAGS += -Ddo_checks
#CFLAGS += -Dno_output
#CXXFLAGS += -Dno_output
LINK_FLAGS = -g



# Uncomment for "perf --call-graph fp" profiling
# CFLAGS += -fno-omit-frame-pointer

# Uncomment for gprof profiling
#CFLAGS += -pg
#CXXFLAGS += -pg
#LINK_FLAGS += -pg


#if using local install for dependencies
#CFLAGS += -I $(DEP_HOME)/include
#LINK_FLAGS += -L $(DEP_LIB) -Wl,-rpath,$(DEP_LIB)


#if using double Interval
#CFLAGS += -DDInt
#CXXFLAGS += -DDInt
#OBJS = DoubleInterval.o
#if using mpfr Interval
OBJS = Interval.o

# PBJ: modelParser.o added.  Not clear why omitted.
# KL: removed modelParser again, it produces conflicts if using profiling (two main functions)
OBJS += Matrix.o Monomial.o Polynomial.o TaylorModel.o Continuous.o 
OBJS += Geometry.o Constraints.o Hybrid.o MyLogger.o ExtractedPicard.o 
OBJS += OutputWriter.o MyComponent.o  Exceptions.o ParsingUtils.o Utils.o 
OBJS += Transformer.o CompSolver.o PreconditionedTMV.o

# All sources
SRCS= $(OBJS:.o=.cpp)

# All dependency files
#DFILES= $(SRCS:.cpp=.d)
#DFILES= $(patsubst %.o,$(DEPDIR)/%.d,$(OBJS))
DFILES= $(patsubst %.o,$(DEPDIR)/%.d,$(OBJS) modelParser.o profiling.o)

all:
	make flowstar
	#make profiling

.PHONY: test profiling

test: 
	make testcompile
	./testcompile

testcompile: $(OBJS) unittesting.o lex.yy.o modelParser.tab.o 
	g++ $(LINK_FLAGS) -o $@ $^ $(LIBS)



.PHONY: fast
fast: 
	make fastcompile
	./fastcompile
  

fastcompile: $(OBJS) lex.yy.o ParsingUtils.o modelParser.tab.o unittestingfast.o
	g++ $(LINK_FLAGS) -o $@ $^ $(LIBS)
  
profiling: $(OBJS) lex.yy.o ParsingUtils.o modelParser.tab.o profiling.o
	g++ $(LINK_FLAGS) -o $@ $^ $(LIBS)

myapp: $(OBJS) lex.yy.o modelParser.tab.o modelParser.o
	g++ $(LINK_FLAGS) -o $@ $^ $(LIBS)
  
flowstar: $(OBJS) lex.yy.o modelParser.tab.o modelParser.o
	g++ $(LINK_FLAGS) -o $@ $^ $(LIBS)

%.o: %.cc %.d
	$(CXX) -c $(CFLAGS) -o $@ $<
%.o: %.cpp %.h $(DEPDIR)/%.d
#	@echo "printing" #why doesn't this print anything
	$(CXX) -c $(CFLAGS) -o $@ $<
unittestingfast.o: unittestingfast.cpp unittesting.h
	$(CXX) -c $(CFLAGS) -o $@ $<

#%.o: %.cpp $(DEPDIR)/%.d
#	$(CXX) -c $(DEPFLAGS) $(CFLAGS) -o $@ $<
#  $(POSTCOMPILE)
%.o: %.c
	$(CXX) -c $(CFLAGS) -o $@ $<

modelParser.tab.c: modelParser.y
	bison -d -v modelParser.y

lex.yy.c: modelLexer.l modelParser.tab.c
	flex modelLexer.l

$(DEPDIR)/%.d: %.cpp
	@echo Making dependency file $@
	@set -e; rm -f $@; \
	 $(CXX) -MM $(CPPFLAGS) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	 rm -f $@.$$$$

$(DEPDIR)/%.d: %.c
	@echo Making dependency file $@
	@set -e; rm -f $@; \
	 $(CXX) -MM $(CPPFLAGS) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	 rm -f $@.$$$$

include $(DFILES)
include $(DEPDIR)/modelParser.tab.d $(DEPDIR)/lex.yy.d

clean: 
	rm -f flowstar *.o *.d *~ modelParser.tab.c modelParser.tab.h \
	modelParser.output lex.yy.c myapp fastcompile testcompile .d/*
