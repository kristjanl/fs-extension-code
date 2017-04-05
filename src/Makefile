DEPDIR = .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td

CXX = g++
include makefile.local
LIBS = -lmpfr -lgmp -lgsl -lgslcblas -lm -lglpk
CFLAGS = -I . -I $(GMP_HOME) -g -O3
CXXFLAGS = -I . -I $(GMP_HOME) -g -O3
LINK_FLAGS = -g -L $(GMP_LIB_HOME) -Wl,-rpath,$(GMP_LIB_HOME)
OBJS = Interval.o Matrix.o Monomial.o Polynomial.o TaylorModel.o Continuous.o Geometry.o Constraints.o Hybrid.o MyLogger.o SimpleImplApp.o ExtractedPicard.o OutputWriter.o SimpleImpl.o SimpleComp.o SmallComp.o CompApp.o MyComponent.o  Exceptions.o ParsingUtils.o Utils.o

TESTOBJS = Interval.o Monomial.o Polynomial.o Matrix.o MyLogger.o TaylorModel.o MyComponent.o Exceptions.o

# All sources
SRCS= $(OBJS:.o=.cpp)

# All dependency files
DFILES= $(SRCS:.cpp=.d)

all: flowstar

.PHONY: test

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

myapp: $(OBJS) lex.yy.o modelParser.tab.o modelParser.o
	g++ $(LINK_FLAGS) -o $@ $^ $(LIBS)

flowstar: $(OBJS) lex.yy.o modelParser.tab.o modelParser.o
	g++ $(LINK_FLAGS) -o $@ $^ $(LIBS) -pg

%.o: %.cc %.d
	$(CXX) -c $(CFLAGS) -o $@ $<
%.o: %.cpp %.h %.d
	$(CXX) -c $(CFLAGS) -o $@ $< -pg
unittestingfast.o: unittestingfast.cpp unittesting.h
	$(CXX) -c $(CFLAGS) -o $@ $< -pg

#%.o: %.cpp $(DEPDIR)/%.d
#	$(CXX) -c $(DEPFLAGS) $(CFLAGS) -o $@ $<
#  $(POSTCOMPILE)
%.o: %.c
	$(CXX) -c $(CFLAGS) -o $@ $<

modelParser.tab.c: modelParser.y
	bison -d -v modelParser.y

lex.yy.c: modelLexer.l modelParser.tab.c
	flex modelLexer.l

%.d: %.cpp
	@echo Making dependency file $@
	@set -e; rm -f $@; \
	 $(CXX) -MM $(CPPFLAGS) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	 rm -f $@.$$$$

%.d: %.c
	@echo Making dependency file $@
	@set -e; rm -f $@; \
	 $(CXX) -MM $(CPPFLAGS) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	 rm -f $@.$$$$

 include $(DFILES)
 include modelParser.tab.d lex.yy.d

clean: 
	rm -f flowstar *.o *.d *~ modelParser.tab.c modelParser.tab.h \
	modelParser.output lex.yy.c myapp fastcompile testcompile
