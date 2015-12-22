
///////////////////////
// MyFFT_CosineLUT.h //
///////////////////////

#ifndef MyFFT_CosineLUT_h
#define MyFFT_CosineLUT_h

#include <complex>
#include <algorithm>
#include <cmath>

template <typename real_type, unsigned n>
class MyFFT_CosineLUT
{
    static_assert((n & (n - 1)) == 0, "MyFFT_CosineLUT<real_type, n> will only work if n is zero or a power of two.");

    typedef std::complex<real_type> complex_type;

    public:

        MyFFT_CosineLUT() // Constructor.
        {
            // Pre-calculate cosine table for twiddle factors.
            // Only values in range [0, pi/2) are needed.

            for (unsigned i = 0; i < n / 4; ++i)
            {
                const real_type turn = 2.0 * M_PI * i / n;
                cosine_table[i] = cos(turn);
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

                        const complex_type term = twiddle(fft_count * k) * odd;

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

        static const unsigned n_div_2 = n / 2;
        static const unsigned n_div_4 = n / 4;

        real_type cosine_table[n_div_4]; // Cosine values.

        complex_type twiddle(unsigned i) const
        {
            int sin_index = i % n_div_2 - n_div_4;
            if (sin_index < 0) // abs
            {
                sin_index = -sin_index;
            }
            real_type sinval = (sin_index == n_div_4) ? 0.0 : cosine_table[sin_index];
            if ((i & n_div_2) == 0)
            {
                // Correct sign of sine.
                // Note that we invert the sign of the sign (for the forward FFT).
                sinval = -sinval;
            }

            int cos_index = n_div_4 - sin_index;
            real_type cosval = (cos_index == n_div_4) ? 0.0 : cosine_table[cos_index];
            if ((i + n_div_4) & n_div_2)
            {
                // Correct sign of cosine.
                cosval = -cosval;
            }

            return complex_type(cosval, sinval);
        }
};

#endif // MyFFT_CosineLUT_h
