CXX = g++
CXXFLAGS = -Wall -O3
OBJS = prpack_preprocessed_graph.o prpack_solver.o test_driver.o
PROG = test_driver

all: ${PROG}
	
${PROG}: ${OBJS}
	${CXX} ${CXXFLAGS} -o $@ ${OBJS}

clean:
	rm *.o ${PROG} -f
