
CXXFLAGS = -W -Wall -O3 --std=c++11 -ffast-math

.PHONY : clean default

default : TestMyFFT

TestMyFFT : TestMyFFT.cc MyFFT.h MyFFT_CosineLUT.h MyFFT_LowLevel.h MyFFT_BitReverse.h
	$(CXX) $(CXXFLAGS) $< -o $@

clean :
	$(RM) *~ TestMyFFT
