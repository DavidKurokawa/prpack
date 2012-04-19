CXX = g++
CXXFLAGS = -Wall -O3 -fopenmp
OBJS = prpack_utils.o prpack_base_graph.o prpack_preprocessed_gs_graph.o prpack_preprocessed_schur_graph.o prpack_preprocessed_scc_graph.o prpack_solver.o prpack_result.o prpack_driver.o
PROG = prpack_driver

all: ${PROG}
	
${PROG}: ${OBJS}
	${CXX} ${CXXFLAGS} -o $@ ${OBJS}

clean:
	rm *.o ${PROG} -f
