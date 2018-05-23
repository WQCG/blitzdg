CXX := $(or $(CXX), g++)
SRCDIR := src
BUILDDIR := build
BINDIR := bin
SRCEXT := cpp

TARGETS := $(patsubst src/,,$(patsubst src/%/,bin/%,$(sort $(dir $(wildcard src/**/)))))
BUILDDIRS := $(subst bin,build,$(TARGETS))

COMMONSOURCES := $(wildcard $(SRCDIR)/*.$(SRCEXT))
COMMONOBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(COMMONSOURCES:.$(SRCEXT)=.o))

SPECIFICSOURCES := $(wildcard $(SRCDIR)/**/*.$(SRCEXT))
SPECIFICOBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SPECIFICSOURCES:.$(SRCEXT)=.o))

ALLOBJECTS := $(COMMONOBJECTS)
ALLOBJECTS += $(SPECIFICOBJECTS)

CFLAGS := -g -Wall -std=c++0x -fprofile-arcs -ftest-coverage
LINKERFLAGS := -fprofile-arcs
INC := -I include -I /usr/include
LIB := -L lib -lblitz -lmetis -lumfpack -lcxsparse -llapack -lblas
EXPLICITLIBS := -lgfortran -lcholmod -lamd -lcolamd -lquadmath -lsuitesparseconfig

ifeq ($(OS), Windows_NT)
	LIB += $(EXPLICITLIBS)
endif

all: builddirs $(TARGETS)

$(TARGETS): $(ALLOBJECTS)
	@mkdir -p $(BINDIR)
	@echo " Linking $@..."
	@echo " $(CXX) $(LINKERFLAGS) $(COMMONOBJECTS) $(wildcard $(subst bin/,,build/$@/*.o)) -o $@ $(LIB)"; $(CXX) $(LINKERFLAGS) $(COMMONOBJECTS) $(wildcard $(subst bin/,,build/$@/*.o)) -o $@ $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT) builddirs
	@echo " $(CXX) $(CFLAGS) $(EXTRACFLAGS) $(INC) -c -o $@ $<"; $(CXX) $(CFLAGS) $(EXTRACFLAGS) $(INC) -c -o $@ $<

builddirs:
	@mkdir -p $(BUILDDIRS)

builddirs:

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(BINDIR)"; $(RM) -r $(BUILDDIR) $(BINDIR)

test: bin/test
	@bin/test

get-deps:
	@sudo /bin/bash pull-deps.sh

docs:
	@echo DOXYGEN VERSION:
	@doxygen -v
	@doxygen -u doxygen.conf
	@sed -ir "s,DOT_PATH.*=,DOT_PATH               = $$(which dot)," doxygen.conf
	@cat doxygen.conf | grep DOT_PATH
	@doxygen doxygen.conf
.PHONY: clean docs
