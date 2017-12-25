CC := g++
SRCDIR := src
BUILDDIR := build
TARGET := bin/blitzdg
TESTTARGET := bin/test

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
BINOBJECTS := $(shell echo $(OBJECTS) | sed 's/build\/test\/tests.o//')
TESTOBJECTS := $(shell echo $(OBJECTS) | sed 's/build\/main.o//')

CFLAGS := -g -Wall -std=c++0x #do we want this last one?
LIB := -pthread -L lib -lblitz -lumfpack -lmetis
INC := -I include -I /usr/include

$(TARGET): $(TESTTARGET)
	@echo " Linking main binary..."
	@echo " $(CC) $(BINOBJECTS) -o $(TARGET) $(LIB)"; $(CC) $(BINOBJECTS) -o $(TARGET) $(LIB)

$(TESTTARGET): $(OBJECTS)
	@mkdir -p $(BUILDDIR)/bin
	@echo " Linking tests..."
	@echo " $(CC) $(TESTOBJECTS) -o $(TESTTARGET) $(LIB)"; $(CC) $(TESTOBJECTS) -o $(TESTTARGET) $(LIB)

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

.PHONY: clean