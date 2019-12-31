SHELL=/bin/sh

.SUFFIXES:
.SUFFIXES: .o .C

#OUTPUT_OPTION=

LDFLAGS=-L../lib -lrt
LDLIBS=-lsundials_cvode -lm
CC=g++
CPPFLAGS=-I../include -I../cvode/include/
CFLAGS=-O3

COMPILE.C=$(CC) $(CFLAGS) $(CPPFLAGS) -c
LINK.C=$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS)

RM=rm -f
MV=mv
MD=makedepend

.C.o:
	$(COMPILE.C) $(OUTPUT_OPTION) $<

OBJS = main.o simulation.o parameterSet.o
TARGET=../bin/simu_mRNA

$(TARGET): $(OBJS)
	$(LINK.C) -o $@ $(OBJS) $(LDLIBS)

main.o: main.C
	$(CC) $(CFLAGS) $(CPPFLAGS) -c main.C

parameterSet.o: parameterSet.C
	$(CC) $(CFLAGS) $(CPPFLAGS) -c parameterSet.C

simulation.o: simulation.C
	$(CC) $(CFLAGS) $(CPPFLAGS) -c simulation.C

clean:
	$(RM) $(OBJS)

depend:
	$(MD) -- $(CPPFLAGS) -- $(CSRCS)

# DO NOT DELETE

