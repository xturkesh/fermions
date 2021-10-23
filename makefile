# Makefile; Author X. Turkeshi
FC    = ifort
MKL   = /home/xhek.turkeshi/intel/oneapi/mkl/2021.4.0/lib/intel64
FLAGS = -L$(MKL) -lmkl_rt -O3 -module $(MODDIR)

SRC = helper free_fermions
LIB = expokit    

SRCDIR = src
LIBDIR = lib

OBJDIR := obj
MODDIR := mod

TARGET = main
EXEC   = driver

DEBUG = -fcheck=all -Wall -g -fbacktrace

OBJSS := $(addprefix $(OBJDIR)/,$(addsuffix .o  , $(SRC)))
OBJLS := $(addprefix $(OBJDIR)/,$(addsuffix .o, $(LIB)))

.PHONY: all
all: $(OBJLS) $(OBJSS) $(EXEC)
	@echo "Prod. Done"

.PHONY: debug
debug: FLAGS += $(DEBUG)
debug: $(OBJLS) $(OBJSS) $(EXEC)
	@echo "Dev. Done"

.PHONY: clean
clean:
	@rm -f $(MODDIR)/*.mod $(OBJDIR)/*.o $(EXEC) slurm*.out
	@rm -rf $(RESDIR)
	@rmdir $(MODDIR) $(OBJDIR) 

$(OBJDIR)/%.o: $(SRCDIR)/%.f90 
	$(FC) $(FLAGS) -o $@ -c $<  

$(OBJDIR)/%.o: $(LIBDIR)/%.*
	$(FC) $(FLAGS) -o $@ -c $<  

$(EXEC): $(OBJSS) $(OBJLS) $(SRCDIR)/$(TARGET).f90
	$(FC) $(FLAGS) -o $@ $^ 

$(OBJSS): | $(OBJDIR) $(MODDIR)

$(OBJLS): | $(OBJDIR) $(MODDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(MODDIR):
	@mkdir -p $(MODDIR)