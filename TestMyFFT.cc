
//////////////////
// TestMyFFT.cc //
//////////////////

#include <cassert>
#include <complex>
#include <iostream>
#include <random>

#include <sys/time.h>

#include "MyFFT.h"
#include "MyFFT_CosineLUT.h"

double gettime()
{
    struct timeval tv;
    int rc = gettimeofday(&tv, nullptr);
    assert(rc == 0);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

template <typename real_type>
void test_correctness()
{
    // data = Prime[10 + Range[16]] + I * Prime[20 + Range[16]];
    // Fourier[data, FourierParameters -> {1, -1}]

    typedef std::complex<real_type> complex_type;

    complex_type x[16] = {
        { 31,  73},
        { 37,  79},
        { 41,  83},
        { 43,  89},
        { 47,  97},
        { 53, 101},
        { 59, 103},
        { 61, 107},
        { 67, 109},
        { 71, 113},
        { 73, 127},
        { 79, 131},
        { 83, 137},
        { 89, 139},
        { 97, 149},
        {101, 151}
    };

    MyFFT<real_type, 16> fft;

    fft(x);

    for (int i = 0; i < 16; ++i)
    {
        std::cout << i << " " << x[i] << std::endl;
    }
}

template <typename real_type, unsigned n>
void test_performance(unsigned num_repeats)
{
    assert(num_repeats % 2 == 1);
    typedef std::complex<real_type> complex_type;

    MyFFT<real_type, n> fft;

    complex_type * x = new complex_type[num_repeats * n];

    std::random_device rd;
    std::mt19937 gen(rd());
 
    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
    std::normal_distribution<real_type> d(0, 1);

    double durations[num_repeats];
    for (unsigned i = 0; i < num_repeats * n; ++i)
    {
        x[i] = { d(gen), d(gen) };
    }

    for (unsigned i = 0; i < num_repeats; ++i)
    {
        double t1 = gettime();
        fft(&x[i * n]);
        double t2 = gettime();
        durations[i] = (t2 - t1) *1e6; // convert to us
    }

    std::sort(durations, durations + num_repeats);
    std::cout << "[test_performance] min " << durations[0] << " median " << durations[(num_repeats - 1) / 2] << " max " << durations[num_repeats - 1] << std::endl;

    delete [] x;
}

int main()
{
    test_correctness<float>();
    //test_performance<float, 1024>(1001);
    return 0;
}
