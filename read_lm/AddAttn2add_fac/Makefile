bin = bin/AddAttn2add_fac

CXX = g++
CXXFLAGS = -std=c++14

$(bin) : src/AddAttn2add_fac.o
	$(CXX) $(CXXFLAGS) -o $(bin) src/AddAttn2add_fac.o

AddAttn2add_fac.o : src/AddAttn2add_fac.cpp

.PHONY : clean

clean :
	-rm -v $(bin) src/AddAttn2add_fac.o
