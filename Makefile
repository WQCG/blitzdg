CC := g++
SRCDIR := src
BUILDDIR := build
TARGET := bin/blitzdg
TESTTARGET := bin/test

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
BINOBJECTS := $(shell echo $(OBJECTS) | tr -s " " "\012" | sed -r 's/build\/test\/.*?\.o//')
TESTOBJECTS := $(shell echo $(OBJECTS) | sed 's/build\/main\.o//')

CFLAGS := -g -Wall -std=c++0x -fprofile-arcs -ftest-coverage 
LINKERFLAGS := -fprofile-arcs
LIB := -pthread -L lib -lblitz -lumfpack -lmetis -lsuperlu -lblas -llapack
INC := -I include -I /usr/include

$(TARGET): $(TESTTARGET)
	@echo " Linking main binary..."
	@echo " $(CC) $(LINKERFLAGS) $(BINOBJECTS) -o $(TARGET) $(LIB)"; $(CC) $(LINKERFLAGS) $(BINOBJECTS) -o $(TARGET) $(LIB)

$(TESTTARGET): $(OBJECTS)
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
	@echo " $(RM) -r $(BUILDDIR) $(TARGET) $(TESTTARGET)"; $(RM) -r $(BUILDDIR) $(TARGET) $(TESTTARGET)

test: bin/test
	@bin/test

get-deps:
	@sudo /bin/bash pull-deps.sh

docs:
	@doxygen doxygen.conf
	@cp -r doxygen/html/* docs/.

.PHONY: clean docs
