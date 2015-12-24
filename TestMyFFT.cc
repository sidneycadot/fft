
//////////////////
// TestMyFFT.cc //
//////////////////

#include <cassert>
#include <complex>
#include <iostream>
#include <iomanip>
#include <random>

#include <sys/time.h>

#include "MyFFT.h"
#include "MyFFT_CosineLUT.h"
#include "MyFFT_BitReverse.h"
#include "MyFFT_LowLevel.h"

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

    complex_type x_ref[16];
    std::copy(std::begin(x), std::end(x), x_ref);

    MyFFT<real_type, 16> fft_ref;
    MyFFT_BitReverse<real_type, 16> fft;

    fft_ref(x_ref);
    fft(x);

    std::cout << std::fixed;
    for (int i = 0; i < 16; ++i)
    {
        std::cout << std::setw(8) << i
                  << std::setw(16) << std::real(x[i]) << std::setw(16) << std::imag(x[i]);

        if (std::abs(x[i] - x_ref[i]) > 1e-10)
        {
            std::cout << " incorrect; should be" << std::setw(16) << std::real(x_ref[i]) << std::setw(16) << std::imag(x_ref[i]);
        }

        std::cout << std::endl;
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
    test_correctness<double>();
    //test_performance<float, 1024>(1001);
    return 0;
}


//       0       152.000000      324.000000
//       1       -20.000000       -4.000000
//       2        -8.000000      -12.000000
//       3         0.000000      -16.000000
