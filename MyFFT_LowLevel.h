
//////////////////////
// MyFFT_LowLevel.h //
//////////////////////

#ifndef MyFFT_LowLevel_h
#define MyFFT_LowLevel_h

#include <complex>
#include <cstring>

template <typename real_type, unsigned n>
class MyFFT_LowLevel
{
    static_assert((n & (n - 1)) == 0, "MyFFT_LowLevel<real_type, n> will only work if n is zero or a power of two.");

    typedef std::complex<real_type> complex_type;

    public:

        MyFFT_LowLevel() // Constructor.
        {
            // Pre-calculate cosine table for twiddle factors.
            // Only values in range [0, pi/2) are needed.

            for (unsigned i = 0; i < n_div_4; ++i)
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

            real_type z_next[2 * n]; // Storage for FFT results.

            unsigned half_size = 1;
            unsigned fft_count = n_div_2;

            while (fft_count != 0)
            {
                // Perform the 'fft_count' FFTs at this level.

                for (unsigned i = 0; i < fft_count; ++i)
                {
                    // Combine the result of the two sub-FFTs.

                    real_type * z_next_lo = (z_next   ) + 2 * i;
                    real_type * z_next_hi = (z_next_lo) + 2 * n_div_2;

                    real_type * zp = (real_type *)z + 2 * i;

                    for (unsigned k = 0; k < half_size; ++k)
                    {
                        const real_type even_re = zp[0];
                        const real_type even_im = zp[1];
                        zp += 2 * fft_count;

                        const real_type odd_re = zp[0];
                        const real_type odd_im = zp[1];
                        zp += 2 * fft_count;

                        real_type twiddle_re;
                        real_type twiddle_im;
                        twiddle_func(k * fft_count, &twiddle_re, &twiddle_im);

                        const real_type term_re = twiddle_re * odd_re - twiddle_im * odd_im;
                        const real_type term_im = twiddle_re * odd_im + twiddle_im * odd_re;

                        z_next_lo[0] = (even_re + term_re);
                        z_next_hi[0] = (even_re - term_re);
                        z_next_lo[1] = (even_im + term_im);
                        z_next_hi[1] = (even_im - term_im);
                        z_next_lo += (2 * fft_count);
                        z_next_hi += (2 * fft_count);
                    }
                }

                memcpy(z, z_next, n * 2 * sizeof(real_type)); // Copy z_next to z.

                fft_count /= 2;
                half_size *= 2;
            }
        }

    private:

        static const unsigned n_div_2 = n / 2;
        static const unsigned n_div_4 = n / 4;

        real_type cosine_table[n_div_4]; // Cosine values.

        void twiddle_func(unsigned i, real_type * re, real_type * im) const
        {
            assert(i < n_div_2);

            const unsigned sin_index = (i <  n_div_4) ? (n_div_4 - i) : (i - n_div_4);
            const unsigned cos_index = n_div_4 - sin_index;

            // Invert the sign of the sign (for the forward FFT).
            real_type sinval = (sin_index == n_div_4) ? 0.0 : -cosine_table[sin_index];

            real_type cosval = (cos_index == n_div_4) ? 0.0 : +cosine_table[cos_index];
            if ((i + n_div_4) & n_div_2)
            {
                // Correct sign of cosine.
                cosval = -cosval;
            }

            *re = cosval;
            *im = sinval;
        }
};

#endif // MyFFT_LowLevel_h
