bin = ../combine_listmode
bin_local = ../combine_listmode_local

CXX = g++
CXXFLAGS = -std=c++14

all : $(bin) $(bin_local)

$(bin) : main.o
	$(CXX) $(CXXFLAGS) -o $(bin) main.o

$(bin_local) : main_local.o
	$(CXX) $(CXXFLAGS) -o $(bin_local) main_local.o

main.o : main.cpp

main_local.o : main_local.cpp

.PHONY : clean

clean :
	-rm -v $(bin) main.o $(bin_local) main_local.o
