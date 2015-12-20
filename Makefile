
CXXFLAGS = -W -Wall -O3 --std=c++11

.PHONY : clean default

default : test_fft test_fft_performance

test_fft : test_fft.cc

test_fft_performance : test_fft_performance.cc

clean :
	$(RM) test_fft test_fft_performance
