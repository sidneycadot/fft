
#include <complex>
#include <iostream>
#include "reference_fft.h"

typedef double real_type;
typedef std::complex<real_type> complex_type;

int main()
{
    // Fourier[Prime[Range[5, 12]] + I * Prime[Range[13, 20]], FourierParameters -> {1, -1}]

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

    pow2_fft<real_type>(Forward, x, 8, 1);

    for (int i = 0; i < 8; ++i)
    {
        std::cout << i << " " << x[i] << std::endl;
    }

    return 0;
}

