
#include <cassert>
#include <sys/time.h>

#include <iostream>
#include <complex>
#include <algorithm>
#include "reference_fft.h"

double gettime(void)
{
    struct timeval tv;
    int rc = gettimeofday(&tv, NULL);
    assert(rc == 0);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

#if 0
template <typename real_type>
class FasterFFT
{
    typedef std::complex<real_type> complex_type;

    public:

        FasterFFT(unsigned n) // constructor
        {
            nn = n / 2;
            w = new complex_type[nn];

            for (unsigned i = 0; i < nn; ++i)
            {
                const real_type turn = -2.0 * M_PI * i / n;
                w[i] = std::polar(static_cast<real_type>(1), turn);
            }
        }

        ~FasterFFT()
        {
            delete [] w;
        }

        void apply(const FourierTransformDirection direction, complex_type * z, unsigned n, unsigned stride)
        {
            if (n == 1)
            {
                return;
            }

            assert (n % 2 == 0); // verify that n is even.

            // Do sub-FFTs on even / odd entries.

            this->apply(direction, z         , n / 2, 2 * stride);
            this->apply(direction, z + stride, n / 2, 2 * stride);

            // Now, we need to combine the values of the sub-FFTs.

            // First, allocate temp[] working space.

            complex_type temp[n];

            // Calculate temp[] values as sum/difference of sub-FFT values.

            for (unsigned i = 0; i < n / 2; ++i)
            {
                assert(i * stride < nn);

                const complex_type term  = w[i * stride] * z[stride * (2 * i + 1)];

                temp[2 * i + 0] = z[stride * 2 * i] + term;
                temp[2 * i + 1] = z[stride * 2 * i] - term;

                // Note that we can undo the scaling by dividing the 'temp' values by 2, here,
                // in case of a forward FFT.
            }

            // Rearrange the temp[] entries into the z[] array.
            // This is the 'butterfly' step.

            for (unsigned i = 0; i < n; ++i)
            {
                const unsigned j = (i / 2) + (i % 2) * (n / 2);

                assert(j < n);

                z[j * stride] = temp[i];
            }

            // All done.
        }

    private:

        complex_type * w;
        unsigned nn;
};
#endif

template <typename real_type>
class FasterFFT2
{
    typedef std::complex<real_type> complex_type;

    public:

        FasterFFT(unsigned n) // constructor
        {
            nn = n / 2;
            w = new complex_type[nn];

            for (unsigned i = 0; i < nn; ++i)
            {
                const real_type turn = -2.0 * M_PI * i / n;
                w[i] = std::polar(static_cast<real_type>(1), turn);
            }
        }

        ~FasterFFT()
        {
            delete [] w;
        }

        void apply(complex_type * z, unsigned n, unsigned stride)
        {
            if (n == 1)
            {
                return;
            }

            assert (n % 2 == 0); // verify that n is even.

            // Do sub-FFTs on even / odd entries.

            this->apply(direction, z         , n / 2, 2 * stride);
            this->apply(direction, z + stride, n / 2, 2 * stride);

            // Now, we need to combine the values of the sub-FFTs.

            // First, allocate temp[] working space.

            complex_type temp[n];

            // Calculate temp[] values as sum/difference of sub-FFT values.

            for (unsigned i = 0; i < n / 2; ++i)
            {
                assert(i * stride < nn);

                const complex_type term  = w[i * stride] * z[stride * (2 * i + 1)];

                temp[2 * i + 0] = z[stride * 2 * i] + term;
                temp[2 * i + 1] = z[stride * 2 * i] - term;

                // Note that we can undo the scaling by dividing the 'temp' values by 2, here,
                // in case of a forward FFT.
            }

            // Rearrange the temp[] entries into the z[] array.
            // This is the 'butterfly' step.

            for (unsigned i = 0; i < n; ++i)
            {
                const unsigned j = (i / 2) + (i % 2) * (n / 2);

                assert(j < n);

                z[j * stride] = temp[i];
            }

            // All done.
        }

    private:

        complex_type * w;
        unsigned nn;
};

template <typename real_type>
class FasterFFT3
{
    typedef std::complex<real_type> complex_type;

    public:

        FasterFFT2(unsigned n) // constructor
        {
            nn = n / 2;
            w = new complex_type[nn];

            for (unsigned i = 0; i < nn; ++i)
            {
                const real_type turn = -2.0 * M_PI * i / n;
                w[i] = std::polar(static_cast<real_type>(1), turn);
            }
        }

        ~FasterFFT2()
        {
            delete [] w;
        }

        void apply(complex_type * z, unsigned n)
        {
            complex_type temp[n];

            unsigned num_ffts = n;
            unsigned fft_size = 1;


            while (fft_size < n)
            {
                num_ffts /= 2;
                fft_size *= 2;

                std::cout << "FFT loop: num_ffts = " << num_ffts << ", fft_size = " << fft_size << std::endl;

                //unsigned substride = fft_size;

                // We perform 'num_ffts' sub-ffts, where we may assume that
                // the smaller-size FFTs have been completed (including reshuffling).

                for (unsigned fi = 0; fi < num_ffts; ++fi)
                {
                    std::cout << "    fft #" << fi << std::endl;

                    for (unsigned i = 0; i < fft_size / 2; ++i)
                    {
                        unsigned S1 = fi + 2 * i;
                        unsigned S2 = 0;
                        unsigned D1 = 0;
                        unsigned D2 = 0;
                        printf("        (%u, %u) -> (%u, %u)\n", S1, S2, D1, D2);
                    }

                    for (unsigned i = 0; i < fft_size / 2; ++i)
                    {
#if 0
                        assert(i * fft_size < nn);

                        const complex_type term  = w[i * fft_size] * z[(2 * i + 1)];

                        assert(i < n);
                        assert((i + fft_size / 2) < n);

                        temp[i               ] = z[2 * i + fi] + term;
                        temp[i + fft_size / 2] = z[2 * i + fi] - term;
#endif
                    }

                }

                for (unsigned i = 0; i < n; ++i)
                {
                    z[i] = temp[i];
                }
            }
        }

    private:

        complex_type * w;
        unsigned nn;
};

template <typename real_type>
void test_pow2_fft(unsigned n, unsigned rep)
{
    assert(rep % 2 == 1);
    typedef std::complex<real_type> complex_type;

    complex_type * x = new complex_type[rep * n];

    std::random_device rd;
    std::mt19937 gen(rd());
 
    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
    std::normal_distribution<real_type> d(0, 1);

    double durations[rep];
    for (unsigned i = 0; i < rep * n; ++i)
    {
        x[i] = { d(gen), d(gen) };
    }

    for (unsigned i = 0; i < rep; ++i)
    {
        double t1 = gettime();
        pow2_fft(Forward, &x[i * n], n, 1);
        double t2 = gettime();
        durations[i] = (t2 - t1) *1e6; // convert to us
    }

    std::sort(durations, durations + rep);
    std::cout << "[test_pow2_fft] min " << durations[0] << " median " << durations[(rep - 1) / 2] << " max " << durations[rep - 1] << std::endl;

    delete [] x;
}

template <typename real_type>
void test_faster_fft2(unsigned n, unsigned rep)
{
    assert(rep % 2 == 1);
    typedef std::complex<real_type> complex_type;

    FasterFFT2<real_type> fft(n);

    complex_type * x = new complex_type[rep * n];

    std::random_device rd;
    std::mt19937 gen(rd());
 
    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
    std::normal_distribution<real_type> d(0, 1);

    double durations[rep];
    for (unsigned i = 0; i < rep * n; ++i)
    {
        x[i] = { d(gen), d(gen) };
    }

    for (unsigned i = 0; i < rep; ++i)
    {
        double t1 = gettime();
        fft.apply(&x[i * n], n);
        double t2 = gettime();
        durations[i] = (t2 - t1) *1e6; // convert to us
    }

    std::sort(durations, durations + rep);
    std::cout << "[test_faster_fft] min " << durations[0] << " median " << durations[(rep - 1) / 2] << " max " << durations[rep - 1] << std::endl;

    delete [] x;
}

template <typename real_type>
void test_correctness_faster2_fft()
{
    // Fourier[Prime[Range[5, 12]] + I * Prime[Range[13, 20]], FourierParameters -> {1, -1}]

    typedef std::complex<real_type> complex_type;

    complex_type x[8] = {
        {11, 41},
        {13, 43},
        {17, 47},
        {19, 53},
        {23, 59},
        {29, 61},
        {31, 67},
        {37, 71}
    };

    FasterFFT2<real_type> fft(8);

    fft.apply(x, 8);

    for (int i = 0; i < 8; ++i)
    {
        std::cout << i << " " << x[i] << std::endl;
    }
}

int main()
{
    test_correctness_faster2_fft<float>();
    //test_pow2_fft<float>(1024, 1001);
    //test_faster_fft2<float>(1024, 1001);
    return 0;
}
