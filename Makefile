# Compiler and compilation flags
CXX := g++
OPTIMIZE := -O2
DEBUG := -g
BOOM := ./
CXXFLAGS := $(OPTIMIZE) -fpermissive -w -I$(BOOM) -std=c++11
LDFLAGS := -lgsl -lm -lgslcblas

# Directories
SRCDIR := src
OBJDIR := obj

# Targets
EXE := QuickBEAST

# Automatically find all source files in the source directory
SRC := $(wildcard $(SRCDIR)/*.C)
OBJ := $(SRC:$(SRCDIR)/%.C=$(OBJDIR)/%.o)

# Rule for the final executable
$(EXE): $(OBJ)
	$(CXX) $^ $(CXXFLAGS) -o $@ $(LDFLAGS) 

# Rule to compile each source file into an object file
$(OBJDIR)/%.o: $(SRCDIR)/%.C
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean

clean:
	rm -rf $(OBJDIR) $(EXE)