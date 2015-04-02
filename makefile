#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

# g++ -L /usr/lib -l uhd -o e100test test_routines.cpp
# g++ -g -L /usr/lib -l uhd -o rxtest  receiver_test.cpp uhd_utilities.cpp
# g++ -g -L /usr/lib -l uhd -o serial_port_test serial_port_test.cpp
# g++ -pthread -o thread_test thread_test.cpp

# g : Indicates debug mode
# c : Indicates compilation only

IDIR =.
CC =g++
CXXFLAGS = -O2 -std=gnu++0x -I$(IDIR)
LINKFLAGS =

OBJDIR = obj
LIBDIR = /usr/lib
INCLUDEDIR = ../../src_libraries/dsp

# Location of the DSP source files
DSP_DIR = ../../src_libraries/dsp

LIBS= 


################ TEST OF DSP ROUTINES

_DEPS = files.h generators.h filters.h upsampling_filters.h correlators.h dsp_complex.h
DEPS = $(patsubst %, $(INCLUDEDIR)/%, $(_DEPS))

_OBJ = SrcDsp.o 
OBJ = $(patsubst %, $(OBJDIR)/%, $(_OBJ))


_DSP_OBJ = dsp_complex.o demodulators.o
DSP_OBJ = $(patsubst %, $(DSP_DIR)/%, $(_DSP_OBJ))

$(OBJDIR)/%.o : %.cpp $(DEPS) 
	$(CC) -c -o $@ $< $(CXXFLAGS)

SrcDsp:$(OBJ) $(DSP_OBJ)	
	$(CC)  -L $(LIBDIR)  -o $@ $^  $(LINKFLAGS) $(LIBS)


.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o *~
