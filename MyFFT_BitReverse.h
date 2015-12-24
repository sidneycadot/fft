
////////////////////////
// MyFFT_BitReverse.h //
////////////////////////

#ifndef MyFFT_BitReverse_h
#define MyFFT_BitReverse_h

#include <complex>
#include <algorithm>

template <typename real_type, unsigned n>
class MyFFT_BitReverse
{
    static_assert((n & (n - 1)) == 0, "MyFFT_BitReverse<real_type, n> will only work if n is zero or a power of two.");

    typedef std::complex<real_type> complex_type;

    public:

        MyFFT_BitReverse() // Constructor.
        {
            // Pre-calculate the twiddle factors.
            // Only values in range [0, pi) are needed.

            for (unsigned i = 0; i < n_div_2; ++i)
            {
                const real_type turn = -2.0 * M_PI * i / n;
                w[i] = std::polar(real_type(1), turn);
            }
        }

        void operator()(complex_type * z) const // Perform an n-point FFT.
        {
            // Perform a bit-reversal permutation.
            // This makes it possible to do an in-place FFT because the pairs
            // of elements being combined in the inner loop below end up in the same
            // positions (leaving all other elements undisturbed).

            for (unsigned i = 0; i < n; ++i)
            {
                const unsigned j = bit_reverse(i);
                if (i < j)
                {
                    std::swap(z[i], z[j]);
                }
            }
 
            // The FFT of a single complex number is the number itself.
            //
            // First, we calculate (n/2) 2-element FFTs by combining two 1-element FFTs.
            // Next,  we calculate (n/4) 4-element FFTs by combining two 2-element FFTs.
            // Next,  we calculate (n/8) 8-element FFTs by combining two 4-element FFTs.
            // ...
            // Repeat until we calculate the single n-element FFT.

            unsigned fft_count = n_div_2;
            unsigned half_size = 1;

            while (fft_count != 0)
            {
                // Perform the 'fft_count' FFTs at this level.

                unsigned idx0 = 0;
                unsigned idx1 = half_size;

                for (unsigned i = 0; i != fft_count; ++i)
                {
                    // Combine the result of the two sub-FFTs.
 
                    for (unsigned k = 0; k != half_size; ++k)
                    {
                        const complex_type even = z[idx0];
                        const complex_type odd  = z[idx1];

                        const complex_type term = w[fft_count * k] * odd;

                        z[idx0] = (even + term);
                        z[idx1] = (even - term);

                        ++idx0;
                        ++idx1;
                    }

                    idx0 += half_size;
                    idx1 += half_size;
                }

                fft_count /= 2;
                half_size *= 2;
            }
        }

    private:

        static unsigned bit_reverse(unsigned i)
        {
            unsigned r = 0;
            for (unsigned q = n; q != 1; q /= 2)
            {
                r *= 2;
                r += i % 2;
                i /= 2;
            }
            return r;
        }

        static const unsigned n_div_2 = n / 2;

        complex_type w[n_div_2]; // Twiddle factors.
};

#endif // MyFFT_BitReverse_h
