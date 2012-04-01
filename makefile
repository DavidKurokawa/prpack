CXX = g++
CXXFLAGS = -Wall -O3 -std=c++0x
OBJS = prpack_adjacency_list.o prpack_preprocessed_gs_graph.o prpack_preprocessed_scc_graph.o prpack_solver.o test_driver.o
PROG = test_driver

all: ${PROG}
	
${PROG}: ${OBJS}
	${CXX} ${CXXFLAGS} -o $@ ${OBJS}

clean:
	rm *.o ${PROG} -f
