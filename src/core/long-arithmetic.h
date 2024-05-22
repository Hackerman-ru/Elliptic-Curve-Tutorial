#ifndef ECG_LONG_ARITHMETIC_H
#define ECG_LONG_ARITHMETIC_H

#include "utils/concepts.h"
#include "utils/fft.h"
#include "utils/string-parser.h"

#include <array>
#include <cassert>

namespace elliptic_curve_guide {
    template<size_t c_bits>
    class uint_t {
        using digit_t = uint32_t;
        using double_digit_t = uint64_t;

        static constexpr size_t c_bits_in_byte = 8;
        static constexpr size_t c_digit_size = sizeof(digit_t) * c_bits_in_byte;
        static constexpr size_t c_digit_number = c_bits / c_digit_size;
        static constexpr size_t c_double_digit_size = sizeof(double_digit_t) * c_bits_in_byte;
        static constexpr size_t c_double_digit_number = c_bits / c_double_digit_size;

        template<size_t V>
        friend class uint_t;

        using digits = std::array<digit_t, c_digit_number>;
        digits m_digits = {};

    public:
        constexpr uint_t() = default;

        template<typename T>
        constexpr uint_t(const T& value) : m_digits(split_into_digits<T>(value)) {}

        constexpr uint_t(const char* str) : m_digits(algorithm::parse_into_uint<uint_t>(str).m_digits) {};

        constexpr uint_t& operator=(const uint_t& value) = default;

        template<typename T>
        constexpr uint_t& operator=(const T& value) {
            m_digits = split_into_digits<T>(value);
            return *this;
        }

        constexpr uint_t& operator=(const char* str) {
            return *this = parse_into_uint<uint_t>(str);
        }

        friend constexpr std::strong_ordering operator<=>(const uint_t& lhs, const uint_t& rhs) {
            for (size_t i = c_digit_number; i > 0; --i) {
                if (lhs[i - 1] != rhs[i - 1]) {
                    return lhs[i - 1] <=> rhs[i - 1];
                }
            }

            return std::strong_ordering::equal;
        }

        friend constexpr bool operator==(const uint_t& lhs, const uint_t& rhs) {
            return lhs.m_digits == rhs.m_digits;
        }

        // operator+
        friend constexpr uint_t operator+(const uint_t& lhs, const uint_t& rhs) {
            uint_t result = lhs;
            result += rhs;
            return result;
        }

        friend constexpr uint_t operator+(uint_t&& lhs, const uint_t& rhs) {
            lhs += rhs;
            return lhs;
        }

        friend constexpr uint_t operator+(const uint_t& lhs, uint_t&& rhs) {
            rhs += lhs;
            return rhs;
        }

        friend constexpr uint_t operator+(uint_t&& lhs, uint_t&& rhs) {
            lhs += rhs;
            return lhs;
        }

        // operator-
        friend constexpr uint_t operator-(const uint_t& lhs, const uint_t& rhs) {
            uint_t result = lhs;
            result -= rhs;
            return result;
        }

        friend constexpr uint_t operator-(uint_t&& lhs, const uint_t& rhs) {
            lhs -= rhs;
            return lhs;
        }

        friend constexpr uint_t operator-(const uint_t& lhs, uint_t&& rhs) {
            rhs -= lhs;
            rhs.negative();
            return rhs;
        }

        friend constexpr uint_t operator-(uint_t&& lhs, uint_t&& rhs) {
            lhs -= rhs;
            return lhs;
        }

        // operator*
        friend constexpr uint_t operator*(const uint_t& lhs, const uint_t& rhs) {
            uint_t result;

            for (size_t i = 0; i < c_digit_number; ++i) {
                uint64_t u = 0;

                for (size_t j = 0; i + j < c_digit_number; ++j) {
                    u = static_cast<uint64_t>(result[i + j])
                      + static_cast<uint64_t>(lhs[i]) * static_cast<uint64_t>(rhs[j]) + (u >> c_digit_size);
                    result[i + j] = static_cast<uint32_t>(u);
                }
            }

            return result;
            //return algorithm::fast_fourier_transform::multiply<c_digit_number>(lhs.m_digits, rhs.m_digits);
        }

        // operator/
        friend constexpr uint_t operator/(const uint_t& lhs, const uint_t& rhs) {
            uint_t result = divide(lhs, rhs);
            return result;
        }

        // operator%
        friend constexpr uint_t operator%(const uint_t& lhs, const uint_t& rhs) {
            uint_t remainder;
            divide(lhs, rhs, &remainder);
            return remainder;
        }

        // operator>>
        friend constexpr uint_t operator>>(const uint_t& lhs, const size_t& rhs) {
            uint_t result = lhs;
            return result >>= rhs;
        }

        friend constexpr uint_t operator>>(uint_t&& lhs, const size_t& rhs) {
            return lhs >>= rhs;
        }

        // operator<<
        friend constexpr uint_t operator<<(const uint_t& lhs, const size_t& rhs) {
            uint_t result = lhs;
            return result <<= rhs;
        }

        friend constexpr uint_t operator<<(uint_t&& lhs, const size_t& rhs) {
            return lhs <<= rhs;
        }

        // operator^
        friend constexpr uint_t operator^(const uint_t& lhs, const uint_t& rhs) {
            uint_t result = lhs;
            return result ^= rhs;
        }

        friend constexpr uint_t operator^(uint_t&& lhs, const uint_t& rhs) {
            return lhs ^= rhs;
        }

        friend constexpr uint_t operator^(const uint_t& lhs, uint_t&& rhs) {
            return rhs ^= lhs;
        }

        friend constexpr uint_t operator^(uint_t&& lhs, uint_t&& rhs) {
            return lhs ^= rhs;
        }

        // operator|
        friend constexpr uint_t operator|(const uint_t& lhs, const uint_t& rhs) {
            uint_t result = lhs;
            result |= rhs;
            return result;
        }

        friend constexpr uint_t operator|(uint_t&& lhs, const uint_t& rhs) {
            lhs |= rhs;
            return lhs;
        }

        friend constexpr uint_t operator|(const uint_t& lhs, uint_t&& rhs) {
            rhs |= lhs;
            return rhs;
        }

        friend constexpr uint_t operator|(uint_t&& lhs, uint_t&& rhs) {
            lhs |= rhs;
            return lhs;
        }

        // operator&
        friend constexpr uint_t operator&(const uint_t& lhs, const uint_t& rhs) {
            uint_t result = lhs;
            result &= rhs;
            return result;
        }

        friend constexpr uint_t operator&(uint_t&& lhs, const uint_t& rhs) {
            lhs &= rhs;
            return lhs;
        }

        friend constexpr uint_t operator&(const uint_t& lhs, uint_t&& rhs) {
            rhs &= lhs;
            return rhs;
        }

        friend constexpr uint_t operator&(uint_t&& lhs, uint_t&& rhs) {
            lhs &= rhs;
            return lhs;
        }

        constexpr uint_t operator-() const {
            uint_t result = *this;
            result.negative();
            return result;
        }

        constexpr uint_t& operator+=(const uint_t& other) {
            digit_t carry = 0;

            for (size_t i = 0; i < c_digit_number; ++i) {
                digit_t sum = carry + other[i];
                m_digits[i] += sum;
                carry = (m_digits[i] < sum) || (sum < carry);
            }

            return *this;
        }

        constexpr uint_t& operator-=(const uint_t& other) {
            digit_t remainder = 0;

            for (size_t i = 0; i < c_digit_number; ++i) {
                digit_t prev = m_digits[i];
                digit_t sum = other[i] + remainder;
                m_digits[i] -= sum;
                remainder = (m_digits[i] > prev) || (sum < remainder);
            }

            return *this;
        }

        constexpr uint_t& operator*=(const uint_t& other) {
            return *this = *this * other;
        }

        constexpr uint_t& operator/=(const uint_t& other) {
            return *this = *this / other;
        }

        constexpr uint_t& operator%=(const uint_t& other) {
            return *this = *this % other;
        }

        constexpr uint_t& operator>>=(size_t shift_size) {
            size_t digit_shift = shift_size >> 5;

            if (digit_shift > 0) {
                for (size_t i = 0; i < c_digit_number; ++i) {
                    if (i + digit_shift < c_digit_number) {
                        m_digits[i] = m_digits[i + digit_shift];
                    } else {
                        m_digits[i] = 0;
                    }
                }
            }

            shift_size %= c_digit_size;

            if (shift_size == 0) {
                return *this;
            }

            for (size_t i = 0; i + digit_shift < c_digit_number; ++i) {
                m_digits[i] >>= shift_size;

                if (i + 1 < c_digit_number) {
                    m_digits[i] |= m_digits[i + 1] << (c_digit_size - shift_size);
                }
            }

            if constexpr (c_bits % c_digit_size != 0) {
                if constexpr (c_digit_number > 1) {
                    m_digits[c_digit_number - 2] |= m_digits[c_digit_number - 1]
                                                 << (c_digit_size - shift_size);
                }
                m_digits[c_digit_number - 1] >>= shift_size;
            }

            return *this;
        }

        constexpr uint_t& operator<<=(size_t shift_size) {
            size_t digit_shift = shift_size >> 5;

            if (digit_shift > 0) {
                for (size_t i = c_digit_number; i > 0; --i) {
                    if (i > digit_shift) {
                        m_digits[i - 1] = m_digits[i - digit_shift - 1];
                    } else {
                        m_digits[i - 1] = 0;
                    }
                }
            }

            shift_size %= c_digit_size;

            if (shift_size == 0) {
                return *this;
            }

            if constexpr (c_bits % c_digit_size != 0) {
                m_digits[c_digit_number - 1] <<= shift_size;

                if constexpr (c_digit_number > 1) {
                    m_digits[c_digit_number - 1] |=
                        m_digits[c_digit_number - 2] >> (c_digit_size - shift_size);
                }
            }

            for (size_t i = c_digit_number; i > digit_shift; --i) {
                m_digits[i - 1] <<= shift_size;

                if (i - 1 > 0) {
                    m_digits[i - 1] |= m_digits[i - 2] >> (c_digit_size - shift_size);
                }
            }

            return *this;
        }

        constexpr uint_t& operator^=(const uint_t& other) {
            for (size_t i = 0; i < c_digit_number; ++i) {
                m_digits[i] ^= other[i];
            }

            return *this;
        }

        constexpr uint_t& operator|=(const uint_t& other) {
            for (size_t i = 0; i < c_digit_number; ++i) {
                m_digits[i] |= other[i];
            }

            return *this;
        }

        constexpr uint_t& operator&=(const uint_t& other) {
            for (size_t i = 0; i < c_digit_number; ++i) {
                m_digits[i] &= other[i];
            }

            return *this;
        }

        [[nodiscard("Optimize unary operator usage")]]
        constexpr uint_t
            operator++(int) {
            uint_t result = *this;
            increment();
            return result;
        }

        constexpr uint_t& operator++() {
            increment();
            return *this;
        }

        [[nodiscard("Optimize unary operator usage")]]
        constexpr uint_t
            operator--(int) {
            uint_t result = *this;
            decrement();
            return result;
        }

        constexpr uint_t& operator--() {
            decrement();
            return *this;
        }

        template<typename T>
        constexpr T convert_to() const;

        template<typename T>
        requires concepts::is_convertible_to<T, digit_t>
        constexpr T convert_to() const {
            size_t shift_size = sizeof(T) * c_bits_in_byte;
            size_t digits_number = shift_size / c_digit_size;

            if (digits_number == 0) {
                return static_cast<T>(m_digits[0]);
            }

            T result = 0;

            for (size_t i = 0; i < c_digit_number && i < digits_number; ++i) {
                result |= static_cast<T>(m_digits[i]) << (i * c_digit_size);
            }

            return result;
        }

        template<>
        constexpr std::string convert_to() const {
            std::string result;
            uint_t clone_of_this = *this;

            do {
                uint_t remainder;
                clone_of_this = divide(clone_of_this, 10, &remainder);
                result.push_back(remainder.m_digits[0] + '0');
            } while (clone_of_this > 0);

            std::reverse(result.begin(), result.end());
            return result;
        }

    private:
        static constexpr size_t size() {
            return c_digit_number;
        }

        template<typename T>
        requires std::numeric_limits<T>::is_integer && concepts::is_upcastable_to<T, digit_t>
        static constexpr digits split_into_digits(T value) {
            return {static_cast<digit_t>(value)};
        }

        template<typename T>
        requires std::numeric_limits<T>::is_integer && concepts::is_downcastable_to<T, digit_t>
        static constexpr digits split_into_digits(T value) {
            digits result = {};

            for (size_t i = 0; i < c_digit_number; ++i) {
                result[i] = static_cast<digit_t>(value);
                value >>= c_digit_size;

                if (value == 0) {
                    break;
                }
            }

            return result;
        }

        template<typename T>
        requires concepts::is_convertible_container<T, digit_t> || requires(T x) {
            { elliptic_curve_guide::uint_t {x} } -> std::same_as<T>;
        }
        static constexpr digits split_into_digits(const T& other) {
            const size_t min_size = std::min(size(), other.size());
            digits result = {};

            for (size_t i = 0; i < min_size; i++) {
                result[i] = static_cast<digit_t>(other[i]);
            }

            return result;
        }

        static constexpr uint_t divide(const uint_t& lhs, const uint_t& rhs, uint_t* remainder = nullptr) {
            if (lhs < rhs) {
                if (remainder != nullptr) {
                    *remainder = lhs;
                }

                return uint_t(0);
            }

            if (rhs.actual_size() == 1) {
                return divide(lhs, rhs[0], remainder);
            }

            return d_divide(lhs, rhs, remainder);
        }

        static constexpr uint_t divide(const uint_t& lhs, const digit_t& rhs, uint_t* remainder = nullptr) {
            uint_t result;
            double_digit_t part = 0;

            for (size_t i = c_digit_number; i > 0; --i) {
                part = (part << (c_digit_size)) + static_cast<double_digit_t>(lhs[i - 1]);

                if (part < rhs) {
                    continue;
                }

                result[i - 1] = static_cast<digit_t>(part / rhs);
                part %= rhs;
            }

            if (remainder != nullptr) {
                *remainder = uint_t(static_cast<digit_t>(part));
            }

            return result;
        }

        static constexpr double_digit_t c_half_digit = static_cast<double_digit_t>(1) << (c_digit_size - 1);
        static constexpr double_digit_t c_digit = static_cast<double_digit_t>(1) << c_digit_size;

        static constexpr uint_t d_divide(const uint_t& lhs, const uint_t& rhs, uint_t* remainder = nullptr) {
            size_t dividend_size = lhs.actual_size();
            size_t divisor_size = rhs.actual_size();

            uint_t<c_bits + c_digit_size> dividend(lhs);
            uint_t divisor(rhs);
            uint_t quotient;

            size_t shift_size = 0;
            digit_t divisor_head = divisor[divisor_size - 1];

            while (divisor_head < c_half_digit) {
                ++shift_size;
                divisor_head <<= 1;
            }

            dividend <<= shift_size;
            divisor <<= shift_size;

            double_digit_t divisor_ = divisor[divisor_size - 1];

            for (size_t i = dividend_size - divisor_size + 1; i > 0; --i) {
                double_digit_t part =
                    (static_cast<double_digit_t>(dividend[i + divisor_size - 1]) << c_digit_size)
                    + static_cast<double_digit_t>(dividend[i + divisor_size - 2]);
                double_digit_t quotient_temp = part / divisor_;
                part %= divisor_;

                if (quotient_temp == c_digit) {
                    --quotient_temp;
                    part += divisor_;
                }

                while (part < c_digit
                       && (quotient_temp * divisor[divisor_size - 2]
                           > (part << c_digit_size) + dividend[i + divisor_size - 3])) {
                    --quotient_temp;
                    part += divisor_;
                }

                int64_t carry = 0;
                int64_t widedigit = 0;

                for (size_t j = 0; j < divisor_size; ++j) {
                    double_digit_t product =
                        static_cast<digit_t>(quotient_temp) * static_cast<double_digit_t>(divisor[j]);
                    widedigit = (static_cast<int64_t>(dividend[i + j - 1]) + carry) - (product & UINT32_MAX);
                    dividend[i + j - 1] = static_cast<digit_t>(widedigit);
                    carry = (widedigit >> c_digit_size) - static_cast<int64_t>(product >> c_digit_size);
                }

                widedigit = static_cast<int64_t>(dividend[i + divisor_size - 1]) + carry;
                dividend[i + divisor_size - 1] = static_cast<digit_t>(widedigit);

                quotient[i - 1] = static_cast<digit_t>(quotient_temp);

                if (widedigit < 0) {
                    --quotient[i - 1];
                    widedigit = 0;

                    for (size_t j = 0; j < divisor_size; ++j) {
                        widedigit += static_cast<double_digit_t>(dividend[i + j - 1]) + divisor[j];
                        dividend[i + j - 1] = static_cast<digit_t>(widedigit);
                        widedigit >>= 32;
                    }
                }
            }

            if (remainder != nullptr) {
                *remainder = uint_t(0);

                for (size_t i = 0; i < divisor_size - 1; ++i) {
                    (*remainder)[i] =
                        (dividend[i] >> shift_size)
                        | (static_cast<double_digit_t>(dividend[i + 1]) << (c_digit_size - shift_size));
                }

                (*remainder)[divisor_size - 1] = dividend[divisor_size - 1] >> shift_size;
            }

            return quotient;
        }

        constexpr void negative() {
            for (size_t i = 0; i < c_digit_number; ++i) {
                m_digits[i] = ~(m_digits[i]);
            }

            ++*this;
        }

        constexpr void increment() {
            for (size_t i = 0; i < c_digit_number; ++i) {
                m_digits[i] += 1;

                if (m_digits[i] != 0) {
                    break;
                }
            }
        }

        constexpr void decrement() {
            for (size_t i = 0; i < c_digit_number; ++i) {
                digit_t temp = m_digits[i];
                m_digits[i] -= 1;

                if (temp >= m_digits[i]) {
                    break;
                }
            }
        }

        constexpr uint_t operator*(digit_t other) const {
            digit_t remainder = 0;
            uint_t result;

            for (size_t i = 0; i < c_digit_number; ++i) {
                double_digit_t prod = m_digits[i] * static_cast<double_digit_t>(other);
                result[i] = static_cast<digit_t>(prod) + remainder;
                remainder = static_cast<digit_t>(result[i] < remainder) + static_cast<digit_t>(prod >> 32);
            }

            return result;
        }

        constexpr const digit_t& operator[](size_t pos) const {
            return m_digits[pos];
        }

        constexpr digit_t& operator[](size_t pos) {
            return m_digits[pos];
        }

        constexpr size_t actual_size() const {
            size_t result = c_digit_number;

            while (result > 0 && m_digits[result - 1] == 0) {
                --result;
            }

            return result;
        }
    };
}   // namespace elliptic_curve_guide

#endif
