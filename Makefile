# Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
# See COPYING and LICENSE files at project root for more details.

CXX := $(or $(CXX), g++)
SRCDIR := src
BUILDDIR := build
BINDIR := bin
INCDIR := include

SRCEXT := cpp
DEPEXT := d
OBJECT := o

.DEFAULT_GOAL := all

TARGETS := $(patsubst src/,,$(patsubst src/%/,bin/%,$(sort $(dir $(wildcard src/**/)))))
BUILDDIRS := $(subst bin,build,$(TARGETS))
VALGRINDTARGETS :=$(patsubst bin/%,artifacts/%,$(TARGETS:%=%.log))

COMMONSOURCES := $(wildcard $(SRCDIR)/*.$(SRCEXT))
COMMONOBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(COMMONSOURCES:.$(SRCEXT)=.o))

SPECIFICSOURCES := $(wildcard $(SRCDIR)/**/*.$(SRCEXT))
SPECIFICOBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SPECIFICSOURCES:.$(SRCEXT)=.o))

ALLOBJECTS := $(COMMONOBJECTS)
ALLOBJECTS += $(SPECIFICOBJECTS)

ALLSOURCES := $(COMMONSOURCES)
ALLSOURCES += $(SPECIFICSOURCES)

DEPS := $(ALLOBJECTS:%.o=%.$(DEPEXT))
DEPPATH := $(BUILDDIR)/$(*F)

.SUFFIXES:
.DELETE_ON_ERROR:

# Pull in dependencies
ifneq ($(MAKECMDGOALS), clean)
-include $(DEPS)
endif

CFLAGS := -g -Wall -std=c++0x -fprofile-arcs -ftest-coverage -DBZ_DEBUG
LINKERFLAGS := -fprofile-arcs
INC := -I $(INCDIR) -I $(INCDIR)/igloo -I /opt/blitzpp/blitz-1.0.1
LIB := -L lib -lblitz -lmetis -lumfpack -lcxsparse -llapack -lblas
EXPLICITLIBS := -lgfortran -lcholmod -lamd -lcolamd -lquadmath -lsuitesparseconfig
VALGRIND := valgrind --error-exitcode=1 --leak-check=full --track-origins=yes

ifeq ($(OS), Windows_NT)
	LIB += $(EXPLICITLIBS)
endif

all: $(TARGETS)

$(TARGETS): $(ALLOBJECTS)
	@mkdir -p $(BINDIR)
	@echo " Linking $@..."
	@echo " $(CXX) $(LINKERFLAGS) $(COMMONOBJECTS) $(wildcard $(subst bin/,,build/$@/*.o)) -o $@ $(LIB)"; $(CXX) $(LINKERFLAGS) $(COMMONOBJECTS) $(wildcard $(subst bin/,,build/$@/*.o)) -o $@ $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIRS)
	@echo " $(CXX) $(CFLAGS) $(INC) -c -o $@ $<"; $(CXX) $(CFLAGS) $(INC) -c -o $@ $<

$(BUILDDIR)/%.$(DEPEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIRS)
	@echo "Building dependency $@..."
	$(CXX) $(CFLAGS) $(INC) -MM -MF $@ -MP -MT $(DEPPATH).$(OBJEXT) -MT $(DEPPATH).$(DEPEXT) $<
	@$(RM) [0-9]

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(BINDIR)" artifacts; $(RM) -r $(BUILDDIR) $(BINDIR) artifacts

test: $(BINDIR)/test
	@$(BINDIR)/test

artifacts/%.log: $(BINDIR)/%
	@mkdir -p artifacts
	@echo "$(VALGRIND) --log-file=$@ $<"; $(VALGRIND) --log-file="$@" $< > $@.out

valgrind: $(VALGRINDTARGETS)

valgrind-print: $(VALGRINDTARGETS)
	@cat artifacts/*

get-deps:
	@sudo /bin/bash pull-deps.sh

docs:
	@echo DOXYGEN VERSION:
	@doxygen -v
	@doxygen -u doxygen.conf
	@doxygen doxygen.conf

.PHONY: clean docs valgrind valgrind-print all
