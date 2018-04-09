CC := $(or $(CXX), g++)
SRCDIR := ./src
BUILDDIR := build
BINDIR := bin
TARGET := bin/blitzdg
TESTTARGET := bin/test

SRCEXT := cpp
SOURCES := $(wildcard $(SRCDIR)/*.cpp)
SOURCES += $(wildcard $(SRCDIR)/test/*.cpp)
ALLOBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS := $(patsubst build/test/tests.o,,$(ALLOBJECTS))
TESTOBJECTS := $(patsubst build/main.o,,$(ALLOBJECTS)) 

CFLAGS := -g -Wall -std=c++0x -fprofile-arcs -ftest-coverage 
LINKERFLAGS := -fprofile-arcs
LIB := -L lib -lblitz -lumfpack -lmetis -llapack -lblas -lgfortran -lamd -lcholmod -lcolamd -lsuitesparseconfig -lquadmath -lGKlib
INC := -I include -I /usr/include

$(TARGET): $(OBJECTS) $(TESTTARGET)
	@echo " Linking main binary..."
	@mkdir -p bin
	@echo " $(CC) $(LINKERFLAGS) $(OBJECTS) -o $(TARGET) $(LIB)"; $(CC) $(LINKERFLAGS) $(OBJECTS) -o $(TARGET) $(LIB)

$(TESTTARGET): $(TESTOBJECTS)
	@mkdir -p bin
	@echo " Linking tests..."
	@echo " $(CC) $(LINKERFLAGS) $(TESTOBJECTS) -o $(TESTTARGET) $(LIB)"; $(CC) $(LINKERFLAGS) $(TESTOBJECTS) -o $(TESTTARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@mkdir -p $(BUILDDIR)/test
	@echo " Building...";
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(BINDIR)"; $(RM) -r $(BUILDDIR) $(BINDIR)

test: bin/test
	@bin/test

get-deps:
	@sudo /bin/bash pull-deps.sh

docs:
	@doxygen doxygen.conf
	@cp -r doxygen/html/* docs/.

.PHONY: clean docs
