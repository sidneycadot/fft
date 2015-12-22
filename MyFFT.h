
/////////////
// MyFFT.h //
/////////////

#ifndef MyFFT_h
#define MyFFT_h

#include <complex>
#include <algorithm>
#include <vector>

template <typename real_type, unsigned n> // n must be power of two.
class MyFFT
{
    typedef std::complex<real_type> complex_type;

    public:

        MyFFT() : w(n / 2) // constructor
        {
            // Pre-calculate array of twiddle factors.
            // Only values in range [0, pi) are needed.

            for (unsigned i = 0; i < w.size(); ++i)
            {
                const real_type turn = -2.0 * M_PI * i / n;
                w[i] = std::polar(real_type(1), turn);
            }
        }

        void operator()(complex_type * z) const // perform n-point FFT
        {
            complex_type znext[n]; // storage for FFT result.

            // The FFT of a single number is the number itself.
            // First, we calculate 2-element FFTs by combining two 1-element FFTs.
            // Next, we calculate 4-element FFTs by combining two 2-element FFTs.
            // Next, we calculate 8-element FFTs by combining two 4-element FFTs.
            // ...
            // Repeat until we calculate the n-element FFTs.

            unsigned half_size = 1;
            unsigned fft_count = n;

            while (fft_count != 1)
            {
                fft_count /= 2;

                // Perform the 'fft_count' FFTs at this level.

                for (unsigned i = 0; i < fft_count; ++i)
                {
                    // Combine the result of the two sub-FFTs.

                    for (unsigned k = 0; k < half_size; ++k)
                    {
                        const complex_type even = z[i + fft_count * (2 * k + 0)];
                        const complex_type odd  = z[i + fft_count * (2 * k + 1)];

                        const complex_type term = w[fft_count * k] * odd;

                        znext[i + fft_count * (k            )] = even + term;
                        znext[i + fft_count * (k + half_size)] = even - term;
                    }
                }

                std::copy(znext, znext + n, z);

                half_size *= 2;
            }
        }

    private:

        std::vector<complex_type> w; // twiddle factors
};

#endif // MyFFT_h
