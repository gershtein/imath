include Makefile.inc

default: test

test: test.o imath.o
	$(LD) -o $@ $< imath.o $(LDFLAGS) $(LIBS) 

test.o:  test.cc

imath.o: imath.cc

all: test

