CXX = g++
CXXFLAGS += -std=c++14 -pedantic -Wall -Wextra -O3

bin/%: src/%.cpp
	${CXX} ${CXXFLAGS} -o '$@' $<

# vim: ft=make
