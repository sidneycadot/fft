
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

            for (unsigned i = 0; i < n_div_2; ++i)
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

            {
                complex_type zc[n];      // Storage for FFT results.
                std::copy(z, z + n, zc); // Copy z to zc.
                for (unsigned i = 0; i < n; ++i)
                {
                    z[i] = zc[bit_reverse(i)];
                    std::cout << "@@@ " << i << " " << bit_reverse(i) << std::endl;
                }
            }

            unsigned half_size = 1;
            unsigned fft_count = n_div_2;

            while (fft_count != 0)
            {
                // Perform the 'fft_count' FFTs at this level.

                for (unsigned i = 0; i < fft_count; ++i)
                {
                    // Combine the result of the two sub-FFTs.

                    for (unsigned k = 0; k < half_size; ++k)
                    {
                        //const unsigned ie = (i * half_size + k) * 2 + 0;
                        //const unsigned io = (i * half_size + k) * 2 + 1;
                        //std::cout << "@@@ " << fft_count << " " << i << " " << k << " " << ie << " " << io << std::endl;

                        const complex_type even = z[(i * half_size + k) * 2 + 0];
                        const complex_type odd  = z[(i * half_size + k) * 2 + half_size];

                        const complex_type term = w[fft_count * k] * odd;

                        z[(i * half_size + k) * 2 + 0] = (even + term);
                        z[(i * half_size + k) * 2 + 1] = (even - term);
                    }
                }

                fft_count /= 2;`
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

#endif // MyFFT_h
