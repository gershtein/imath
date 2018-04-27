include Makefile.inc

default: test

test: test.o imath.o
	$(LD) -o $@ $< imath.o $(LDFLAGS) $(LIBS) 

test1: test1.o imath.o
	$(LD) -o $@ $< imath.o $(LDFLAGS) $(LIBS) 

test1.o:  test1.cc

test.o:  test.cc

imath.o: imath.cc

all: test

