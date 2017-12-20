CC := g++
SRCDIR := src
TESTDIR := test
BUILDDIR := build
TARGET := bin/blitzdg

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS := -g -Wall -std=c++0x #do we want this last one?
LIB := -pthread -L lib -lblitz -lumfpack -lmetis -lgmsh
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

tests: $(shell echo $(OBJECTS) | sed 's/build\/main.o//' )
	@echo " $(CC) $(CFLAGS) ${TESTDIR}/tests.cpp $(INC) $(LIB) -c -o build/test.o"; $(CC) $(CFLAGS) ${TESTDIR}/tests.cpp $(INC) $(LIB) -c -o build/test.o
	#@echo " $(CC) build/test.o build/LUSolver.o -o bin/test $(LIB)"; $(CC) build/test.o build/LUSolver.o -o bin/test $(LIB)
	@echo " $(CC) build/test.o $^ -o bin/test $(LIB)"; $(CC) build/test.o $^ -o bin/test $(LIB)
.PHONY: clean