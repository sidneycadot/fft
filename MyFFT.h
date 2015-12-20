
/////////////
// MyFFT.h //
/////////////

#ifndef MyFFT_h
#define MyFFT_h

#include <complex>
#include <algorithm>

template <typename real_type>
class MyFFT
{
    typedef std::complex<real_type> complex_type;

    public:

        MyFFT(const unsigned & n) : n(n) // constructor
        {
            // Pre-calculate array of twiddle factors.
            // only values in range [0, pi) are needed.

            const unsigned nn = n / 2;

            w = new complex_type[nn];

            for (unsigned i = 0; i < nn; ++i)
            {
                const real_type turn = -2.0 * M_PI * i / n;
                w[i] = std::polar(static_cast<real_type>(1), turn);
            }
        }

        ~MyFFT()
        {
            delete [] w;
        }

        void operator()(complex_type * z) const
        {
            complex_type temp[n];

            unsigned half_size = 1;
            unsigned fft_count = n;

            while (fft_count != 1)
            {
                fft_count /= 2;

                // Walk over the FFTs at this level.

                for (unsigned i = 0; i < fft_count; ++i)
                {
                    // Combine the result of the two sub-FFTs.

                    for (unsigned k = 0; k < half_size; ++k)
                    {
                        const complex_type even = z[i + fft_count * (2 * k + 0)];
                        const complex_type odd  = z[i + fft_count * (2 * k + 1)];

                        const complex_type term  = w[fft_count * k] * odd;

                        temp[i + fft_count * (k            )] = even + term;
                        temp[i + fft_count * (k + half_size)] = even - term;
                    }
                }

                std::copy(temp, temp + n, z);

                half_size *= 2;
            }
        }

    private:

        const unsigned n;
        complex_type * w; // twiddle factors
};

#endif // MyFFT_h
