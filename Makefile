CC := g++
SRCDIR := src
TESTDIR := test
BUILDDIR := build
TARGET := bin/blitzdg

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS := -g -Wall -std=c++0x #do we want this last one?
LIB := -pthread -L lib -lblitz -lumfpack
INC := -I include

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " Building...";
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

tests:
	@echo " $(CC) $(CFLAGS) ${TESTDIR}/bdd_ex.cpp $(INC) $(LIB) -o bin/test"; $(CC) $(CFLAGS) ${TESTDIR}/bdd_ex.cpp $(INC) $(LIB) -o bin/test

.PHONY: clean