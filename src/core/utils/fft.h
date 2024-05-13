#ifndef ECG_FFT_H
#define ECG_FFT_H

#include <cassert>
#include <complex>
#include <numbers>
#include <vector>

namespace elliptic_curve_guide {
    namespace algorithm {
        namespace fast_fourier_transform {
            using float_point_t = double;
            using complex = std::complex<float_point_t>;

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
            static const std::array<complex, size> get_roots() {
                constexpr size_t c_log_size = bit_size(size);
                std::array<complex, size> roots;
                roots[1] = complex(1, 0);

                for (size_t k = 1; k < c_log_size; ++k) {
                    float_point_t angle =
                        2 * std::numbers::pi_v<float_point_t> / static_cast<float_point_t>((1ULL << (k + 1)));
                    complex z(std::cos(angle), std::sin(angle));

                    for (size_t i = 1ULL << (k - 1); i < (1ULL << k); ++i) {
                        roots[i << 1] = roots[i];
                        roots[(i << 1) + 1] = roots[i] * z;
                    }
                }

                return roots;
            }

            template<size_t size>
            static constexpr size_t real_size(const std::array<uint32_t, size>& arr) {
                size_t result = size;

                while (arr[result - 1] == 0) {
                    if (result == 1) {
                        return 0;
                    }

                    --result;
                }

                return result - 1;
            }

            template<size_t size>
            static constexpr std::array<complex, size> fft(const std::array<complex, size>& arr) {
                if constexpr (size <= 1) {
                    return arr;
                } else {
                    static constexpr std::array<size_t, size> reverse = get_reverse<size>();
                    static const std::array<complex, size> roots = get_roots<size>();
                    std::array<complex, size> y;

                    for (size_t i = 0; i < size; ++i) {
                        y[i] = arr[reverse[i]];
                    }

                    for (size_t k = 1; k < size; k <<= 1) {
                        const size_t len = k << 1;
                        for (size_t i = 0; i < size; i += len) {
                            for (size_t j = 0; j < k; ++j) {
                                complex z = roots[j + k] * y[i + j + k];
                                y[i + j + k] = y[i + j] - z;
                                y[i + j] = y[i + j] + z;
                            }
                        }
                    }

                    return y;
                }
            }

            template<size_t size>
            constexpr std::array<uint32_t, size>
                multiply(std::array<uint32_t, size> lhs, std::array<uint32_t, size> rhs) {
                static constexpr size_t double_size = size << 1;
                std::array<complex, double_size> in_lhs;
                std::array<complex, double_size> in_rhs;

                for (size_t i = 0; i < size; ++i) {
                    in_lhs[i] = complex(lhs[i]);
                    in_rhs[i] = complex(rhs[i]);
                }

                std::array<complex, double_size> y_lhs = fft<double_size>(in_lhs);
                std::array<complex, double_size> y_rhs = fft<double_size>(in_rhs);

                std::array<complex, double_size> y;

                for (size_t i = 0; i < double_size; ++i) {
                    y[i] = y_lhs[i] * y_rhs[i] / static_cast<complex>(double_size);
                }

                std::reverse(y.begin(), y.end());

                static constexpr uint64_t c_32_bit_mask = UINT32_MAX;
                std::array<complex, double_size> out = fft<double_size>(y);
                std::array<uint32_t, size> result = {};
                uint64_t carry = 0;

                const size_t result_size = std::min(size, real_size<size>(lhs) + real_size<size>(rhs) + 2);

                for (size_t i = 0; i < result_size; ++i) {
                    float_point_t val = std::abs(out[i]);
                    assert(val <= static_cast<float_point_t>(UINT64_MAX));
                    uint64_t value = static_cast<uint64_t>(std::round(val)) + carry;
                    carry = value >> 32;
                    result[i] = static_cast<uint32_t>(value & c_32_bit_mask);
                }

                return result;
            }

            /*using complex = std::complex<long double>;

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
        constexpr std::array<uint32_t, size> multiply(std::array<uint32_t, size> lhs,
                                                      std::array<uint32_t, size> rhs) {
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
        }*/
        }   // namespace fast_fourier_transform
    }       // namespace algorithm
}   // namespace elliptic_curve_guide

#endif
