
/////////////
// MyFFT.h //
/////////////

#ifndef MyFFT_h
#define MyFFT_h

#include <complex>
#include <algorithm>

template <typename real_type, unsigned n>
class MyFFT
{
    static_assert((n & (n - 1)) == 0, "MyFFT<real_type, n> will only work if n is zero or a power of two.");

    typedef std::complex<real_type> complex_type;

    public:

        MyFFT() // Constructor.
        {
            // Pre-calculate the twiddle factors.
            // Only values in range [0, pi) are needed.

            for (unsigned i = 0; i < n / 2; ++i)
            {
                const real_type turn = -2.0 * M_PI * i / n;
                w[i] = std::polar(real_type(1), turn);
            }
        }

        void operator()(complex_type * z) const // Perform an n-point FFT.
        {
            // The FFT of a single complex number is the number itself.
            //
            // First, we calculate (n/2) 2-element FFTs by combining two 1-element FFTs.
            // Next,  we calculate (n/4) 4-element FFTs by combining two 2-element FFTs.
            // Next,  we calculate (n/8) 8-element FFTs by combining two 4-element FFTs.
            // ...
            // Repeat until we calculate the single n-element FFT.

            complex_type z_next[n]; // Storage for FFT results.

            unsigned half_size = 1;
            unsigned fft_count = n / 2;

            while (fft_count != 0)
            {
                // Perform the 'fft_count' FFTs at this level.

                for (unsigned i = 0; i < fft_count; ++i)
                {
                    // Combine the result of the two sub-FFTs.

                    for (unsigned k = 0; k < half_size; ++k)
                    {
                        const complex_type even = z[i + fft_count * (2 * k + 0)];
                        const complex_type odd  = z[i + fft_count * (2 * k + 1)];

                        const complex_type term = w[fft_count * k] * odd;

                        z_next[i + fft_count * (k            )] = (even + term);
                        z_next[i + fft_count * (k + half_size)] = (even - term);
                    }
                }

                std::copy(z_next, z_next + n, z); // Copy z_next to z.

                fft_count /= 2;
                half_size *= 2;
            }
        }

    private:

        complex_type w[n / 2]; // Twiddle factors.
};

#endif // MyFFT_h
