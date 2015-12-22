
CXXFLAGS = -W -Wall -O3 --std=c++11

.PHONY : clean default

default : TestMyFFT

TestMyFFT : TestMyFFT.cc MyFFT.h
	$(CXX) $(CXXFLAGS) $< -o $@

clean :
	$(RM) *~ TestMyFFT
