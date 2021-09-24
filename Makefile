CC = g++
CCFLAGS = -g --std=c++14
EXE = main
SRCDIR = FE_poroelasticity
TARGS = Node.o ElementQ4.o ElementQ4Cohesive.o Problem.o main.o
# Show the dependencies
vpath %.cc ${SRCDIR}
all: ${EXE}

%.o: %.cc
	${CC} ${CCFLAGS} -c $?

${EXE}: ${TARGS}
	${CC} ${CCFLAGS} ${TARGS} -o $@

run:
	./${EXE}
clean:
	rm -rf ${TARGS} ${EXE} *.txt
	
cleantext:
	rm -rf *.txt
.PHONY:
	clean cleantext run

