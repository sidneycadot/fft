
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

                temp[i        ] = z[stride * 2 * i] + term;
                temp[i + n / 2] = z[stride * 2 * i] - term;

                // Note that we can undo the scaling by dividing the 'temp' values by 2, here,
                // in case of a forward FFT.
            }

            // Rearrange the temp[] entries into the z[] array.
            // This is the 'butterfly' step.

            for (unsigned i = 0; i < n; ++i)
            {
                z[i * stride] = temp[i];
            }

            // All done.
        }

    private:

        complex_type * w;
        unsigned nn;
};

template <typename real_type>
class FasterFFT2
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

        void apply_nonrec(complex_type * z, unsigned n)
        {
            unsigned count = 1;
            unsigned stride = n;

            while (count <= n)
            {
                count  *= 2;
                stride /= 2;

                // Walk over the FFTs at this level
                for (unsigned offset = 0; offset < stride; ++offset)
                {
                    // Combine the result of the two sub-FFTs.

                    // First, allocate temp[] working space.

                    complex_type temp[count];

                    // Calculate temp[] values as sum/difference of sub-FFT values.

                    for (unsigned i = 0; i < count / 2; ++i)
                    {
                        assert(i * stride < nn);

                        const complex_type term  = w[i * stride] * z[offset + stride * (2 * i + 1)];

                        temp[i            ] = z[offset + stride * 2 * i] + term;
                        temp[i + count / 2] = z[offset + stride * 2 * i] - term;

                        // Note that we can undo the scaling by dividing the 'temp' values by 2, here,
                        // in case of a forward FFT.
                    }

                    // Rearrange the temp[] entries into the z[] array.
                    // This is the 'butterfly' step.

                    for (unsigned i = 0; i < count; ++i)
                    {
                        z[offset + i * stride] = temp[i];
                    }

                    // All done.
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

    FasterFFT2<real_type> fft(16);

    fft.apply_nonrec(x, 16);

    for (int i = 0; i < 16; ++i)
    {
        std::cout << i << " " << x[i] << std::endl;
    }
}

int main()
{
    test_correctness_faster2_fft<double>();
    //test_pow2_f ft<float>(1024, 1001);
    //test_faster_fft2<float>(1024, 1001);
    return 0;
}

// =============== 16
//
// 0 2 8
// 1 2 8
// 2 2 8
// 3 2 8
// 4 2 8
// 5 2 8
// 6 2 8
// 7 2 8
//
// 0 4 4
// 1 4 4
// 2 4 4
// 3 4 4
//
// 0 8 2
// 1 8 2
//
// 0 16 1
