
CXXFLAGS = -W -Wall -O3 --std=c++11 -ffast-math

.PHONY : clean default

default : TestMyFFT

TestMyFFT : TestMyFFT.cc MyFFT.h
	$(CXX) $(CXXFLAGS) $< -o $@

clean :
	$(RM) *~ TestMyFFT
