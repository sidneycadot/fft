
/////////////////////
// reference_fft.h //
/////////////////////

#ifndef reference_fft_h
#define reference_fft_h

#include <cassert>
#include <complex>

// Just about the simplest FFT imaginable.
// Handles complex 1D in-place transform of power-of-two points.

enum FourierTransformDirection { Forward, Inverse };

template <typename real_type>
void pow2_fft(const FourierTransformDirection direction, std::complex<real_type> z[], const unsigned n, const unsigned stride)
{
    typedef std::complex<real_type> complex_type;

    if (n == 1)
    {
        return;
    }

    assert (n % 2 == 0); // verify that n is even.

    // Do sub-FFTs on even / odd entries.

    pow2_fft(direction, z         , n / 2, 2 * stride);
    pow2_fft(direction, z + stride, n / 2, 2 * stride);

    // Now, we need to combine the values of the sub-FFTs.

    // First, allocate temp[] working space.

    complex_type temp[n];

    // Calculate temp[] values as sum/difference of sub-FFT values.

    for (unsigned i = 0; i < n / 2; ++i)
    {
        const real_type turn = ((direction == Forward) ? (-2.0 * M_PI) : (+2.0 * M_PI)) * i / n;

        const complex_type coeff = std::polar(static_cast<real_type>(1), turn);

        temp[2 * i + 0] = (z[stride * 2 * i] + coeff * z[stride * (2 * i + 1)]);
        temp[2 * i + 1] = (z[stride * 2 * i] - coeff * z[stride * (2 * i + 1)]);

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

#endif // reference_fft_h
