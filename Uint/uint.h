#ifndef ECG_LONGINT_H
#define ECG_LONGINT_H

#ifndef CHAR_BIT
    #define CHAR_BIT 8
#endif

#include "../util.h"

#include <array>
#include <cassert>
#include <complex>
#include <functional>
#include <string>

namespace ECG {
    template<size_t bits>
    class uint_t {
        using bucket = uint32_t;
        using double_bucket = uint64_t;
        using shift = size_t;

    public:
        constexpr uint_t() = default;

        template<is_convertible<bucket> T>
        constexpr uint_t(T value) : m_buckets({static_cast<bucket>(value)}) {};

        template<>
        constexpr uint_t(bucket value) : m_buckets({value}) {};

        template<>
        constexpr uint_t(size_t value) {
            for (size_t i = 0; i < c_BUCKET_NUMBER; ++i) {
                m_buckets[i] = static_cast<bucket>(value);
                value >>= c_BUCKET_SIZE;

                if (value == 0) {
                    break;
                }
            }
        }

        template<typename T>
        requires requires(T t, size_t i) {
            { t[i] } -> is_convertible<bucket>;
            { t.size() } -> std::same_as<size_t>;
        }
        uint_t(const T& other) {
            const size_t min_size = std::min(size(), other.size());

            for (size_t i = 0; i < min_size; i++) {
                m_buckets[i] = static_cast<bucket>(other[i]);
            }
        }

        uint_t(const char* str) {
            if (str == nullptr) {
                return;
            }

            while (*str != '\0') {
                uint_t temp = (*this <<= 1);
                *this <<= 2;
                *this += temp;
                *this += uint_t(*str - '0');
                ++str;
            }
        }

        uint_t(const std::string& str) {
            *this = uint_t(str, find_type(str));
        }

        uint_t(const std::string& str, StringType str_type) {
            switch (str_type) {
            case StringType::BINARY :
                for (size_t i = 0; i < str.size() && i < bits; ++i) {
                    size_t bucket_pos = i / c_BUCKET_SIZE;

                    m_buckets[bucket_pos] |= bucket(str[str.size() - 1 - i] - '0') << (i % c_BUCKET_SIZE);
                }
                break;
            case StringType::DECIMAL :
                for (char c : str) {
                    uint_t temp = (*this <<= 1);
                    *this <<= 2;
                    *this += temp;
                    *this += uint_t(c - '0');
                }
                break;
            case StringType::HEXADECIMAL :
                for (char c : str) {
                    *this <<= 4;

                    if (isdigit(c)) {
                        *this += uint_t(c - '0');
                    } else {
                        *this += uint_t((c - 'a') + 10);
                    }
                }
                break;
            }
        }

        uint_t(const std::string& str, std::function<bucket(char)> map, shift shift_size) {
            for (char c : str) {
                *this <<= shift_size;
                *this += uint_t(map(c));
            }
        }

        // operator<=>
        auto operator<=>(const uint_t& other) const = default;

        // operator==
        friend bool operator==(const uint_t& lhs, const uint_t& rhs) {
            return lhs.m_buckets == rhs.m_buckets;
        }

        // operator+
        friend uint_t operator+(const uint_t& lhs, const uint_t& rhs) {
            uint_t result = lhs;
            return result += rhs;
        }

        friend uint_t operator+(uint_t&& lhs, const uint_t& rhs) {
            return lhs += rhs;
        }

        friend uint_t operator+(const uint_t& lhs, uint_t&& rhs) {
            return rhs += lhs;
        }

        friend uint_t operator+(uint_t&& lhs, uint_t&& rhs) {
            return lhs += rhs;
        }

        // operator-
        friend uint_t operator-(const uint_t& lhs, const uint_t& rhs) {
            uint_t result = lhs;
            return result -= rhs;
        }

        friend uint_t operator-(uint_t&& lhs, const uint_t& rhs) {
            return lhs -= rhs;
        }

        friend uint_t operator-(const uint_t& lhs, uint_t&& rhs) {
            return -(rhs -= lhs);   // or result version, if operator-() will be expensive
        }

        friend uint_t operator-(uint_t&& lhs, uint_t&& rhs) {
            return lhs -= rhs;
        }

        // operator*
        friend uint_t operator*(const uint_t& lhs, const uint_t& rhs) {   // FFT will change this
            uint_t<bits> result;

            for (size_t i = 0; i < c_BUCKET_NUMBER; ++i) {
                result += (lhs * rhs[i]) << (c_BUCKET_SIZE * i);
            }

            return result;
        }

        // operator/
        friend uint_t operator/(const uint_t& lhs, const uint_t& rhs) {   // FFT will change this
            return divide(lhs, rhs);
        }

        // operator%
        friend uint_t operator%(const uint_t& lhs, const uint_t& rhs) {   // FFT will change this
            uint_t remainder;
            divide(lhs, rhs, &remainder);
            return remainder;
        }

        // operator>>
        friend uint_t operator>>(const uint_t& lhs, const size_t& rhs) {
            uint_t result = lhs;
            return result >>= rhs;
        }

        friend uint_t operator>>(uint_t&& lhs, const size_t& rhs) {
            return lhs >>= rhs;
        }

        // operator<<
        friend uint_t operator<<(const uint_t& lhs, const size_t& rhs) {
            uint_t result = lhs;
            return result <<= rhs;
        }

        friend uint_t operator<<(uint_t&& lhs, const size_t& rhs) {
            return lhs <<= rhs;
        }

        // operator^
        friend uint_t operator^(const uint_t& lhs, const uint_t& rhs) {
            uint_t result = lhs;
            return result ^= rhs;
        }

        friend uint_t operator^(uint_t&& lhs, const uint_t& rhs) {
            return lhs ^= rhs;
        }

        friend uint_t operator^(const uint_t& lhs, uint_t&& rhs) {
            return rhs ^= lhs;
        }

        friend uint_t operator^(uint_t&& lhs, uint_t&& rhs) {
            return lhs ^= rhs;
        }

        // operator|
        friend uint_t operator|(const uint_t& lhs, const uint_t& rhs) {
            uint_t result = lhs;
            return result |= rhs;
        }

        friend uint_t operator|(uint_t&& lhs, const uint_t& rhs) {
            return lhs |= rhs;
        }

        friend uint_t operator|(const uint_t& lhs, uint_t&& rhs) {
            return rhs |= lhs;
        }

        friend uint_t operator|(uint_t&& lhs, uint_t&& rhs) {
            return lhs |= rhs;
        }

        // operator&
        friend uint_t operator&(const uint_t& lhs, const uint_t& rhs) {
            uint_t result = lhs;
            return result &= rhs;
        }

        friend uint_t operator&(uint_t&& lhs, const uint_t& rhs) {
            return lhs &= rhs;
        }

        friend uint_t operator&(const uint_t& lhs, uint_t&& rhs) {
            return rhs &= lhs;
        }

        friend uint_t operator&(uint_t&& lhs, uint_t&& rhs) {
            return lhs &= rhs;
        }

        uint_t operator-() const {
            uint_t result;

            for (size_t i = 0; i < c_BUCKET_NUMBER; ++i) {
                result[i] = ~(m_buckets[i]);
            }

            result++;

            return result;
        }

        uint_t& operator+=(const uint_t& other) {
            bucket carry = 0;

            for (size_t i = 0; i < c_BUCKET_NUMBER; ++i) {
                bucket sum = carry + other[i];
                m_buckets[i] += sum;
                carry = (m_buckets[i] < sum) || (sum < carry);
            }

            return *this;
        }

        uint_t& operator-=(const uint_t& other) {
            bucket remainder = 0;

            for (size_t i = 0; i < c_BUCKET_NUMBER; ++i) {
                bucket temp = m_buckets[i];
                m_buckets[i] -= other[i] + remainder;
                remainder = (temp < m_buckets[i]);
            }

            return *this;
        }

        uint_t& operator*=(const uint_t& other) {
            return (*this = *this * other);
        }

        uint_t& operator/=(const uint_t& other) {
            return (*this = *this / other);
        }

        uint_t& operator%=(const uint_t& other) {
            return (*this = *this % other);
        }

        uint_t& operator>>=(shift shift_size) {
            shift bucket_shift = shift_size / c_BUCKET_SIZE;

            for (size_t i = 0; i + bucket_shift < c_BUCKET_NUMBER; ++i) {
                m_buckets[i] = m_buckets[i + bucket_shift];
            }

            shift_size %= c_BUCKET_SIZE;

            if (shift_size == 0) {
                return *this;
            }

            for (size_t i = 0; i + bucket_shift < c_BUCKET_NUMBER; ++i) {
                m_buckets[i] >>= shift_size;

                if (i + 1 < c_BUCKET_NUMBER) {
                    m_buckets[i] |= m_buckets[i + 1] << (c_BUCKET_SIZE - shift_size);
                }
            }

            return *this;
        }

        uint_t& operator<<=(shift shift_size) {
            shift bucket_shift = shift_size / c_BUCKET_SIZE;

            for (size_t i = c_BUCKET_NUMBER; i > 0; --i) {
                if (i > bucket_shift) {
                    m_buckets[i - 1] = m_buckets[i - bucket_shift - 1];
                } else {
                    m_buckets[i - 1] = 0;
                }
            }

            shift_size %= c_BUCKET_SIZE;

            if (shift_size == 0) {
                return *this;
            }

            for (size_t i = c_BUCKET_NUMBER; i > bucket_shift; --i) {
                m_buckets[i - 1] <<= shift_size;

                if (i - 1 > 0) {
                    m_buckets[i - 1] |= m_buckets[i - 2] >> (c_BUCKET_SIZE - shift_size);
                }
            }

            return *this;
        }

        uint_t& operator^=(const uint_t& other) {
            for (size_t i = 0; i < c_BUCKET_NUMBER; ++i) {
                m_buckets[i] ^= other[i];
            }

            return *this;
        }

        uint_t& operator|=(const uint_t& other) {
            for (size_t i = 0; i < c_BUCKET_NUMBER; ++i) {
                m_buckets[i] |= other[i];
            }

            return *this;
        }

        uint_t& operator&=(const uint_t& other) {
            for (size_t i = 0; i < c_BUCKET_NUMBER; ++i) {
                m_buckets[i] &= other[i];
            }

            return *this;
        }

        uint_t& operator++(int) {
            return (*this += 1);
        }

        uint_t& operator++() {
            return (*this += 1);
        }

        uint_t& operator--() {
            return (*this -= 1);
        }

        const bucket& operator[](size_t pos) const {
            return m_buckets[pos];
        }

        bucket& operator[](size_t pos) {
            return m_buckets[pos];
        }

        template<is_convertible_reverse<bucket> T>
        T convert_to() const {
            shift shift_size = sizeof(T) * CHAR_BIT;
            size_t bucket_number = shift_size / c_BUCKET_SIZE;
            T result = 0;

            for (size_t i = 0; i < c_BUCKET_NUMBER && i < bucket_number; ++i) {
                result |= T(m_buckets[i]) << (i * c_BUCKET_SIZE);
            }

            shift_size %= c_BUCKET_SIZE;

            result |= ((T(m_buckets[bucket_number]) << (c_BUCKET_SIZE - shift_size))
                       >> (c_BUCKET_SIZE - shift_size))
                   << (bucket_number * c_BUCKET_SIZE);

            return result;
        }

        std::string into_string(StringType str_type = StringType::DECIMAL) const {
            std::string result;
            uint_t clone_of_this = *this;

            switch (str_type) {
            case StringType::BINARY :
                do {
                    result.push_back((clone_of_this[0] & static_cast<bucket>(StringType::BINARY)) + '0');
                    clone_of_this >>= 1;
                } while (clone_of_this > 0);
                break;
            case StringType::DECIMAL :
                do {
                    uint_t remainder;
                    clone_of_this =
                        divide(clone_of_this, static_cast<bucket>(StringType::DECIMAL), &remainder);
                    result.push_back(remainder[0] + '0');
                } while (clone_of_this > 0);
                break;
            case StringType::HEXADECIMAL :
                do {
                    bucket value = clone_of_this[0] & static_cast<bucket>(StringType::HEXADECIMAL);

                    if (value >= 10) {
                        value -= 10;
                        result.push_back(value + 'a');
                    } else {
                        result.push_back(value + '0');
                    }

                    clone_of_this >>= 4;
                } while (clone_of_this > 0);
                break;
            }

            std::reverse(result.begin(), result.end());
            return result;
        }

        std::string into_string(std::function<char(bucket)> map, shift shift_size) const {
            std::string result;
            uint_t clone = *this;

            do {
                result.push_back(
                    map((bucket(clone) << (c_BUCKET_SIZE - shift_size)) >> (c_BUCKET_SIZE - shift_size)));
                clone >>= shift_size;
            } while (clone > 0);

            std::reverse(result.begin(), result.end());
            return result;
        }

        static constexpr size_t size() {
            return c_BUCKET_NUMBER;
        }

    private:
        static uint_t divide(const uint_t& lhs, const uint_t& rhs, uint_t* remainder = nullptr) {
            size_t dividend_size = lhs.clz();
            size_t divisor_size = rhs.clz();

            // CASE 0:
            if (dividend_size < divisor_size) {
                if (remainder != nullptr) {
                    *remainder = lhs;
                }

                return uint_t(0);
            }

            // CASE 1:
            if (divisor_size == 1) {
                return divide(lhs, rhs[0], remainder);
            }

            // CASE 2:
            return d_divide(lhs, rhs, remainder);
        }

        static uint_t divide(const uint_t& lhs, const bucket& rhs, uint_t* remainder = nullptr) {
            uint_t result;
            double_bucket part = 0;

            for (size_t i = c_BUCKET_NUMBER; i > 0; --i) {
                part = (part << (c_BUCKET_SIZE)) + static_cast<double_bucket>(lhs[i - 1]);

                if (part < rhs) {
                    continue;
                }

                result[i - 1] = static_cast<bucket>(part / rhs);
                part %= rhs;
            }

            if (remainder != nullptr) {
                *remainder = uint_t(static_cast<bucket>(part));
            }

            return result;
        }

        static uint_t d_divide(const uint_t& lhs, const uint_t& rhs, uint_t* remainder = nullptr) {
            size_t dividend_size = lhs.clz();
            size_t divisor_size = rhs.clz();

            uint_t<bits + c_BUCKET_SIZE> dividend(lhs);
            uint_t divisor(rhs);
            uint_t quotient;

            shift shift_size = 0;
            bucket divisor_head = divisor[divisor_size - 1];
            static constexpr double_bucket c_HALF_BUCKET = static_cast<double_bucket>(1)
                                                        << (c_BUCKET_SIZE - 1);

            while (divisor_head < c_HALF_BUCKET) {
                ++shift_size;
                divisor_head <<= 1;
            }

            dividend <<= shift_size;
            divisor <<= shift_size;

            double_bucket divisor_ = divisor[divisor_size - 1];
            static constexpr double_bucket c_BUCKET = static_cast<double_bucket>(1) << c_BUCKET_SIZE;

            for (size_t i = dividend_size - divisor_size + 1; i > 0; --i) {
                double_bucket part =
                    (static_cast<double_bucket>(dividend[i + divisor_size - 1]) << c_BUCKET_SIZE)
                    + static_cast<double_bucket>(dividend[i + divisor_size - 2]);
                double_bucket quotient_t = part / divisor_;
                part %= divisor_;

                if (quotient_t == c_BUCKET) {
                    --quotient_t;
                    part += divisor_;
                }

                while (part < c_BUCKET
                       && (quotient_t * divisor[divisor_size - 2]
                           > (part << c_BUCKET_SIZE) + dividend[i + divisor_size - 3])) {
                    --quotient_t;
                    part += divisor_;
                }

                int64_t carry = 0;
                int64_t widedigit = 0;

                for (size_t j = 0; j < divisor_size; ++j) {
                    double_bucket product =
                        static_cast<bucket>(quotient_t) * static_cast<double_bucket>(divisor[j]);
                    widedigit =
                        (static_cast<int64_t>(dividend[i + j - 1]) + carry) - (product & 0xffffffffLL);
                    dividend[i + j - 1] = static_cast<bucket>(widedigit);
                    carry = (widedigit >> c_BUCKET_SIZE) - static_cast<int64_t>(product >> c_BUCKET_SIZE);
                }

                widedigit = static_cast<int64_t>(dividend[i + divisor_size - 1]) + carry;
                dividend[i + divisor_size - 1] = static_cast<bucket>(widedigit);

                quotient[i - 1] = static_cast<bucket>(quotient_t);

                if (widedigit < 0) {
                    --quotient[i - 1];
                    widedigit = 0;

                    for (size_t j = 0; j < divisor_size; ++j) {
                        widedigit += static_cast<double_bucket>(dividend[i + j - 1]) + divisor[j];
                        dividend[i + j - 1] = static_cast<bucket>(widedigit);
                        widedigit >>= 32;
                    }
                }
            }

            if (remainder != nullptr) {
                *remainder = uint_t(0);

                for (size_t i = 0; i < divisor_size - 1; ++i) {
                    (*remainder)[i] =
                        (dividend[i] >> shift_size)
                        | (static_cast<double_bucket>(dividend[i + 1]) << (c_BUCKET_SIZE - shift_size));
                }

                (*remainder)[divisor_size - 1] = dividend[divisor_size - 1] >> shift_size;
            }

            return quotient;
        }

        static StringType find_type(const std::string& str) {
            if (str.size() < 2) {
                return StringType::DECIMAL;
            }

            switch (str[1]) {
            case 'x' :
                return StringType::HEXADECIMAL;
            case 'b' :
                return StringType::BINARY;
            default :
                return StringType::DECIMAL;
            }
        }

        uint_t operator*(bucket other) const {   // FFT will delete this
            bucket remainder = 0;
            uint_t result;

            for (size_t i = 0; i < c_BUCKET_NUMBER; ++i) {
                double_bucket prod = m_buckets[i] * static_cast<double_bucket>(other);
                result[i] = static_cast<bucket>(prod) + remainder;
                remainder = static_cast<bucket>(result[i] < remainder) + static_cast<bucket>(prod >> 32);
            }

            return result;
        }

        size_t clz() const {
            size_t result = c_BUCKET_NUMBER;

            while (result > 0 && m_buckets[result - 1] == 0) {
                --result;
            }

            return result;
        }

        static constexpr size_t c_BUCKET_SIZE = sizeof(bucket) * CHAR_BIT;
        static constexpr size_t c_BUCKET_NUMBER = bits / c_BUCKET_SIZE;

        std::array<bucket, c_BUCKET_NUMBER> m_buckets = {};
    };
}   // namespace ECG

#endif
