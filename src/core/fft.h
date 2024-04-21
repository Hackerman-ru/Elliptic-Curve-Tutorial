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
            value >>= 1;
            ++result;
        } while (value);

        return result - 1;
    }

    template<size_t size>
    static constexpr std::array<size_t, size> get_reverse() {
        constexpr size_t c_shift = bit_size(size);
        std::array<size_t, size> reverse;
        reverse[0] = 0;

        for (size_t i = 1; i < size; ++i) {
            reverse[i] = (reverse[i >> 1] | ((i & 1) << c_shift)) >> 1;
        }

        return reverse;
    }

    template<size_t size>
    static constexpr void fft(std::array<complex, size>& a) {
        if constexpr (size <= 1) {
            return;
        }

        static constexpr std::array<size_t, size> c_reverse = get_reverse<size>();

        for (size_t i = 0; i < size; ++i) {
            if (i < c_reverse[i]) {
                std::swap(a[i], a[c_reverse[i]]);
            }
        }

        static std::vector<complex> temp(2, 1);
        static std::vector<complex> root(2, 1);

        for (static size_t k = 2; k < size; k <<= 1) {
            temp.resize(size);
            root.resize(size);
            auto x = std::polar(1.0L, std::acos(-1.0L) / static_cast<double>(k));

            for (size_t i = k; i < (k << 1); ++i) {
                root[i] = temp[i] = temp[i >> 1] * (i & 1 ? x : 1);
            }
        }

        for (size_t k = 1; k < size; k <<= 1) {
            for (size_t i = 0; i < size; i += (k << 1)) {
                for (size_t j = 0; j < k; ++j) {
                    complex z = root[j + k] * a[i + j + k];
                    a[i + j + k] = a[i + j] - z;
                    a[i + j] = a[i + j] + z;
                }
            }
        }
    }

    template<size_t size>
    constexpr std::array<uint32_t, size>
        multiply(std::array<uint32_t, size> lhs, std::array<uint32_t, size> rhs) {
        std::array<complex, size> in;

        for (size_t i = 0; i < size; ++i) {
            in[i] = complex(lhs[i], rhs[i]);
        }

        fft<size>(in);

        for (complex& value : in) {
            value *= value;
        }

        std::array<complex, size> out;

        for (size_t i = 0; i < size; ++i) {
            out[i] = in[(size - i) & (size - 1)] - std::conj(in[i]);
        }

        fft<size>(out);
        static constexpr size_t shift = bit_size(size) + 2;
        std::array<uint32_t, size> result;
        uint64_t carry = 0;

        for (size_t i = 0; i < size; ++i) {
            uint64_t value = (static_cast<uint64_t>(std::imag(out[i])) >> shift) + carry;
            carry = value >> 32;
            result[i] = static_cast<uint32_t>(value);
        }

        return result;
    }
}   // namespace ECG

#endif
