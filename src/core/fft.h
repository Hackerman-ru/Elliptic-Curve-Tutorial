#ifndef ECG_FFT_H
#define ECG_FFT_H

#include <complex>
#include <vector>

#pragma warning(disable : 4244)   // for size_t -> complex conversions
#pragma warning(disable : 4530)   // for resize of complex values
#pragma warning(disable : 4723)   // for division by size

namespace ECG {
    using complex = std::complex<long double>;

    static complex fast_pow(const complex& value, size_t power) {
        if (power == 0) {
            return 1;
        }

        if (power == 1) {
            return value;
        }

        if (power & 1) {
            return value * fast_pow(value, power - 1);
        } else {
            complex temp = fast_pow(value, power << 1);
            return temp * temp;
        }
    }

    using namespace std::complex_literals;

    static constexpr size_t bit_size(size_t value) {
        size_t result = 0;

        do {
            ++result;
            value >>= 1;
        } while (value);

        return result;
    }

    template<size_t size>
    static constexpr void fill_reverse(std::array<size_t, size>& reverse) {
        static constexpr size_t shift = bit_size(size);
        reverse[0] = 0;

        for (size_t i = 1; i < size; ++i) {
            reverse[i] = (reverse[i >> 1] | ((i & 1) << shift)) >> 1;
        }
    }

    template<size_t size>
    void evaluate(std::array<complex, size>& coeffs) {
        if constexpr (size <= 1) {
            return;
        }

        static const complex wn =
            exp(static_cast<complex>(2 * 3.14) * static_cast<complex>(1.0i) / static_cast<complex>(size));

        static constexpr size_t half_size = size >> 1;

        std::array<complex, half_size> even_coeffs;
        std::array<complex, half_size> odd_coeffs;

        for (size_t i = 0; i < size; ++i) {
            if (i & 1) {
                odd_coeffs[i >> 1] = coeffs[i];
            } else {
                even_coeffs[i >> 1] = coeffs[i];
            }
        }

        evaluate<half_size>(even_coeffs);
        evaluate<half_size>(odd_coeffs);

        complex w(1);

        for (size_t i = 0; i < half_size; ++i) {
            coeffs[i] = even_coeffs[i] + w * odd_coeffs[i];
            coeffs[i + half_size] = even_coeffs[i] - w * odd_coeffs[i];
            w *= wn;
        }
    }

    template<size_t size>
    void interpolate(std::array<complex, size>& values) {
        if constexpr (size <= 1) {
            return;
        }

        static const double angle = -2.0 * acos(-1) / static_cast<double>(size);
        static complex wn(cos(angle), sin(angle));
        wn /= size;

        static constexpr size_t half_size = size >> 1;

        std::array<complex, half_size> even_values;
        std::array<complex, half_size> odd_values;

        for (size_t i = 0; i < size; ++i) {
            if (i & 1) {
                odd_values[i >> 1] = values[i];
            } else {
                even_values[i >> 1] = values[i];
            }
        }

        interpolate<half_size>(even_values);
        interpolate<half_size>(odd_values);

        complex w(1);

        for (size_t i = 0; i < half_size; ++i) {
            values[i] = even_values[i] + w * odd_values[i];
            values[i + half_size] = even_values[i] - w * odd_values[i];
            w *= wn;
        }
    }

    template<size_t size>
    std::array<complex, size> fft(std::array<uint32_t, size>&& coeffs) {
        std::array<complex, size> complex_coeffs;

        for (size_t i = 0; i < size; ++i) {
            complex_coeffs[i] = static_cast<complex>(coeffs[i]);
        }

        evaluate(complex_coeffs);
        return complex_coeffs;
    }

    template<size_t size>
    std::array<uint32_t, size> ifft(std::array<complex, size>&& values) {
        interpolate<size>(values);
        std::array<uint32_t, size> coeffs;

        for (size_t i = 0; i < size; ++i) {
            coeffs[i] = std::abs(values[i]);
        }

        return coeffs;
    }

    template<size_t size>
    std::array<uint32_t, size> multiply(std::array<uint32_t, size> lhs, std::array<uint32_t, size> rhs) {
        std::array<complex, size> lhs_values = fft(std::move(lhs));
        std::array<complex, size> rhs_values = fft(std::move(rhs));

        std::array<complex, size> result_values;

        for (size_t i = 0; i < size; ++i) {
            result_values[i] = lhs_values[i] * rhs_values[i];
        }

        std::array<uint32_t, size> result = ifft(std::move(result_values));
        return result;
    }
}   // namespace ECG

#endif
