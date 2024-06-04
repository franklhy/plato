# Build plato

# Check SWIG version
MIN_SWIG_VERSION = 4.0
SWIG_VERSION := $(shell echo `swig -version | grep SWIG | awk '{print $$3}'`)
IS_SWIG_ABOVE_MIN := $(shell expr "$(SWIG_VERSION)" ">=" "$(MIN_SWIG_VERSION)")
ifeq ($(IS_SWIG_ABOVE_MIN),0)
    $(error Minimum version requirement of SWIG ($(MIN_SWIG_VERSION)) is not met)
endif

UNAME = $(shell uname)

# OS-depedent compiler, since MacOS g++ compiler does not support openMP
ifeq ($(UNAME), Darwin)
    CC = clang
else
    CC = g++
endif

SRCDIR = $(shell pwd)/code/src
HELPERDIR = $(shell pwd)/code/helper
OUTDIR = $(shell pwd)/plato
CCFLAG = -O3 -fPIC -std=c++11 -lstdc++ -lm
OMP = -fopenmp

SRC_H = \
	$(SRCDIR)/Container.h \
	$(SRCDIR)/ReadDump.h \
	$(SRCDIR)/WriteDump.h \
	$(SRCDIR)/ReadData.h \
	$(SRCDIR)/WriteData.h \
	$(SRCDIR)/Neigh.h \
	$(SRCDIR)/Topo.h \
	$(SRCDIR)/forswig.h
SRC_CPP = \
	$(SRCDIR)/Container.cpp \
	$(SRCDIR)/ReadDump.cpp \
	$(SRCDIR)/WriteDump.cpp \
	$(SRCDIR)/ReadData.cpp \
	$(SRCDIR)/WriteData.cpp \
	$(SRCDIR)/Neigh.cpp \
	$(SRCDIR)/Topo.cpp

SRC_PY = $(SRCDIR)/dump.py $(SRCDIR)/data.py $(SRCDIR)/neigh.py $(SRCDIR)/topo.py $(SRCDIR)/hbond.py
SRC_I = $(SRCDIR)/plato.i
HELPER_I = $(HELPERDIR)/numpy.i

PREQUEST = $(SRC_H) $(SRC_CPP) $(SRC_PY) $(SRC_I) $(HELPER_I)


PYTHON_INCLUDES = $(shell python3-config --includes)
PYTHON_PREFIX = $(shell python3-config --prefix)
PYTHON_LIBS = $(shell python3-config --libs)
NUMPY_INCLUDES = $(shell python3 -c "import numpy; print(numpy.get_include())")

# OS-depedent shared library flag, see https://github.com/dfroger/swig-python-skel/commit/28539f1273
ifeq ($(UNAME), Darwin)
    SHLIB_FLAGS = -shared -undefined dynamic_lookup
else
    SHLIB_FLAGS = -shared
endif

.PHONY: clean

all: $(OUTDIR)/plato.py $(OUTDIR)/_plato.so

$(OUTDIR)/plato.py: $(PREQUEST)
	swig -c++ -python -outdir $(OUTDIR) -o $(OUTDIR)/plato_wrap.cxx -I$(SRCDIR) -I$(HELPERDIR) $(SRC_I)

$(OUTDIR)/_plato.so: $(OUTDIR)/Container.o $(OUTDIR)/ReadDump.o $(OUTDIR)/WriteDump.o $(OUTDIR)/ReadData.o $(OUTDIR)/WriteData.o $(OUTDIR)/Neigh.o $(OUTDIR)/Topo.o $(OUTDIR)/plato_wrap.o
	$(CC) $(SHLIB_FLAGS) $(OMP) $(OUTDIR)/Container.o $(OUTDIR)/ReadDump.o $(OUTDIR)/WriteDump.o $(OUTDIR)/ReadData.o $(OUTDIR)/WriteData.o $(OUTDIR)/Neigh.o $(OUTDIR)/Topo.o $(OUTDIR)/plato_wrap.o -o $(OUTDIR)/_plato.so

$(OUTDIR)/Container.o: $(SRCDIR)/Container.cpp
	$(CC) $(CCFLAG) -c $(SRCDIR)/Container.cpp -o $(OUTDIR)/Container.o

$(OUTDIR)/ReadDump.o: $(SRCDIR)/ReadDump.cpp
	$(CC) $(CCFLAG) -c $(SRCDIR)/ReadDump.cpp -o $(OUTDIR)/ReadDump.o

$(OUTDIR)/WriteDump.o: $(SRCDIR)/WriteDump.cpp
	$(CC) $(CCFLAG) -c $(SRCDIR)/WriteDump.cpp -o $(OUTDIR)/WriteDump.o

$(OUTDIR)/ReadData.o: $(SRCDIR)/ReadData.cpp
	$(CC) $(CCFLAG) -c $(SRCDIR)/ReadData.cpp -o $(OUTDIR)/ReadData.o

$(OUTDIR)/WriteData.o: $(SRCDIR)/WriteData.cpp
	$(CC) $(CCFLAG) -c $(SRCDIR)/WriteData.cpp -o $(OUTDIR)/WriteData.o

$(OUTDIR)/Neigh.o: $(SRCDIR)/Neigh.cpp
	$(CC) $(CCFLAG) $(OMP) -c $(SRCDIR)/Neigh.cpp -o $(OUTDIR)/Neigh.o

$(OUTDIR)/Topo.o: $(SRCDIR)/Topo.cpp
	$(CC) $(CCFLAG) -c $(SRCDIR)/Topo.cpp -o $(OUTDIR)/Topo.o

$(OUTDIR)/plato_wrap.o: $(OUTDIR)/plato_wrap.cxx
	$(CC) $(CCFLAG) -c $(OUTDIR)/plato_wrap.cxx -o $(OUTDIR)/plato_wrap.o $(PYTHON_INCLUDES) -L$(PYTHON_PREFIX)/lib $(PYTHON_LIB) -I$(NUMPY_INCLUDES) -I$(SRCDIR)

$(OUTDIR)/plato_wrap.cxx: $(PREQUEST)
	swig -c++ -python -outdir $(OUTDIR) -o $(OUTDIR)/plato_wrap.cxx -I$(SRCDIR) -I$(HELPERDIR) $(SRC_I)

clean:
	rm -rf $(OUTDIR)/*
