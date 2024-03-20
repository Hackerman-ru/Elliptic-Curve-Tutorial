#ifndef ECG_LONGINT_H
#define ECG_LONGINT_H

#include "../util.h"

#include <array>
#include <bitset>
#include <cassert>
#include <complex>
#include <functional>
#include <string>

namespace ECG {
    template<size_t bits>
    class uint_t {
        using block = uint32_t;
        using double_block = uint64_t;   // TEMPORARY

        static constexpr size_t c_BITS_IN_BYTE = 8;
        static constexpr size_t c_BLOCK_SIZE = sizeof(block) * c_BITS_IN_BYTE;
        static constexpr size_t c_BLOCK_NUMBER = bits / c_BLOCK_SIZE;

        std::array<block, c_BLOCK_NUMBER> m_blocks = {};

    public:
        constexpr uint_t() = default;

        template<is_convertible_to<block> T>
        constexpr uint_t(T value) : m_blocks(split_into_blocks<T>(value)) {
            if constexpr (sizeof(T) * c_BITS_IN_BYTE <= c_BLOCK_SIZE) {
                m_blocks[0] = static_cast<block>(value);
            } else {
                for (size_t i = 0; i < c_BLOCK_NUMBER; ++i) {
                    m_blocks[i] = static_cast<block>(value);
                    value >>= c_BLOCK_SIZE;

                    if (value == 0) {
                        break;
                    }
                }
            }
        }

        template<ConvertibleContainer<block> T>
        uint_t(const T& other) {
            const size_t min_size = std::min(size(), other.size());

            for (size_t i = 0; i < min_size; i++) {
                m_blocks[i] = static_cast<block>(other[i]);
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
                    size_t bucket_pos = i / c_BLOCK_SIZE;

                    m_blocks[bucket_pos] |= block(str[str.size() - 1 - i] - '0') << (i % c_BLOCK_SIZE);
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

        uint_t(const std::string& str, std::function<block(char)> map, size_t shift_size) {
            for (char c : str) {
                *this <<= shift_size;
                *this += uint_t(map(c));
            }
        }

        // operator<=>
        auto operator<=>(const uint_t& other) const = default;

        // operator==
        friend bool operator==(const uint_t& lhs, const uint_t& rhs) {
            return lhs.m_blocks == rhs.m_blocks;
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

            for (size_t i = 0; i < c_BLOCK_NUMBER; ++i) {
                result += (lhs * rhs[i]) << (c_BLOCK_SIZE * i);
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

            for (size_t i = 0; i < c_BLOCK_NUMBER; ++i) {
                result[i] = ~(m_blocks[i]);
            }

            result++;

            return result;
        }

        uint_t& operator+=(const uint_t& other) {
            block carry = 0;

            for (size_t i = 0; i < c_BLOCK_NUMBER; ++i) {
                block sum = carry + other[i];
                m_blocks[i] += sum;
                carry = (m_blocks[i] < sum) || (sum < carry);
            }

            return *this;
        }

        uint_t& operator-=(const uint_t& other) {
            block remainder = 0;

            for (size_t i = 0; i < c_BLOCK_NUMBER; ++i) {
                block temp = m_blocks[i];
                m_blocks[i] -= other[i] + remainder;
                remainder = (temp < m_blocks[i]);
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

        uint_t& operator>>=(size_t shift_size) {
            static constexpr size_t BUCKET_SIZE = 64;
            static constexpr size_t BUCKET_NUMBER = bits / BUCKET_SIZE;

            size_t bucket_shift = shift_size >> 6;
            auto data = static_cast<uint64_t*>(static_cast<void*>(m_blocks.data()));

            if (bucket_shift > 0) {
                for (size_t i = 0; i < BUCKET_NUMBER; ++i) {
                    if (i + bucket_shift < BUCKET_NUMBER) {
                        data[i] = data[i + bucket_shift];
                    } else {
                        data[i] = 0;
                    }
                }
            }

            shift_size %= BUCKET_SIZE;

            if (shift_size == 0) {
                return *this;
            }

            for (size_t i = 0; i + bucket_shift < BUCKET_NUMBER; ++i) {
                data[i] >>= shift_size;

                if (i + 1 < BUCKET_NUMBER) {
                    data[i] |= data[i + 1] << (BUCKET_SIZE - shift_size);
                }
            }

            if constexpr (bits % BUCKET_SIZE != 0) {
                if constexpr (c_BLOCK_NUMBER > 1) {
                    m_blocks[c_BLOCK_NUMBER - 2] |= m_blocks[c_BLOCK_NUMBER - 1]
                                                 << (c_BLOCK_SIZE - shift_size);
                }
                m_blocks[c_BLOCK_NUMBER - 1] >>= shift_size;
            }

            return *this;
        }

        uint_t& operator<<=(size_t shift_size) {
            static constexpr size_t BUCKET_SIZE = 64;
            static constexpr size_t BUCKET_NUMBER = bits / BUCKET_SIZE;

            size_t bucket_shift = shift_size >> 6;
            auto data = static_cast<uint64_t*>(static_cast<void*>(m_blocks.data()));

            if (bucket_shift > 0) {
                for (size_t i = BUCKET_NUMBER; i > 0; --i) {
                    if (i > bucket_shift) {
                        data[i - 1] = data[i - bucket_shift - 1];
                    } else {
                        data[i - 1] = 0;
                    }
                }
            }

            shift_size %= BUCKET_SIZE;

            if (shift_size == 0) {
                return *this;
            }

            if constexpr (bits % BUCKET_SIZE != 0) {
                m_blocks[c_BLOCK_NUMBER - 1] <<= shift_size;

                if constexpr (c_BLOCK_NUMBER > 1) {
                    m_blocks[c_BLOCK_NUMBER - 1] |=
                        m_blocks[c_BLOCK_NUMBER - 2] >> (c_BLOCK_SIZE - shift_size);
                }
            }

            for (size_t i = BUCKET_NUMBER; i > bucket_shift; --i) {
                data[i - 1] <<= shift_size;

                if (i - 1 > 0) {
                    data[i - 1] |= data[i - 2] >> (BUCKET_SIZE - shift_size);
                }
            }

            return *this;
        }

        uint_t& operator^=(const uint_t& other) {
            for (size_t i = 0; i < c_BLOCK_NUMBER; ++i) {
                m_blocks[i] ^= other[i];
            }

            return *this;
        }

        uint_t& operator|=(const uint_t& other) {
            for (size_t i = 0; i < c_BLOCK_NUMBER; ++i) {
                m_blocks[i] |= other[i];
            }

            return *this;
        }

        uint_t& operator&=(const uint_t& other) {
            for (size_t i = 0; i < c_BLOCK_NUMBER; ++i) {
                m_blocks[i] &= other[i];
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

        const block& operator[](size_t pos) const {
            return m_blocks[pos];
        }

        block& operator[](size_t pos) {
            return m_blocks[pos];
        }

        template<is_convertible_from<block> T>
        T convert_to() const {
            size_t shift_size = sizeof(T) * c_BITS_IN_BYTE;
            size_t bucket_number = shift_size / c_BLOCK_SIZE;
            T result = 0;

            for (size_t i = 0; i < c_BLOCK_NUMBER && i < bucket_number; ++i) {
                result |= T(m_blocks[i]) << (i * c_BLOCK_SIZE);
            }

            return result;
        }

        std::string into_string(StringType str_type = StringType::DECIMAL) const {
            std::string result;
            uint_t clone_of_this = *this;

            switch (str_type) {
            case StringType::BINARY :
                do {
                    result.push_back(((clone_of_this & 1) != 0) + '0');
                    clone_of_this >>= 1;
                } while (clone_of_this > 0);
                break;
            case StringType::DECIMAL :
                do {
                    uint_t remainder;
                    clone_of_this =
                        divide(clone_of_this, static_cast<block>(StringType::DECIMAL), &remainder);
                    result.push_back(remainder[0] + '0');
                } while (clone_of_this > 0);
                break;
            case StringType::HEXADECIMAL :
                do {
                    block value = clone_of_this[0] & static_cast<block>(StringType::HEXADECIMAL);

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

        std::string into_string(std::function<char(block)> map, size_t shift_size) const {
            std::string result;
            uint_t clone = *this;

            do {
                result += map((clone[0] << (c_BLOCK_SIZE - shift_size)) >> (c_BLOCK_SIZE - shift_size));
                clone >>= shift_size;
            } while (clone > 0);

            std::reverse(result.begin(), result.end());
            return result;
        }

        static constexpr size_t size() {
            return c_BLOCK_NUMBER;
        }

    private:
        template<typename T>
        static std::array<block, c_BLOCK_NUMBER> split_into_blocks(T value);

        template<is_upcastable_to<block> T>
        static std::array<block, c_BLOCK_NUMBER> split_into_blocks(T value) {
            return std::array<block, c_BLOCK_NUMBER>(value);
        }

        template<is_downcastable_to<block> T>
        static std::array<block, c_BLOCK_NUMBER> split_into_blocks(T value) {
            std::array<block, c_BLOCK_NUMBER> result;

            for (size_t i = 0; i < c_BLOCK_NUMBER; ++i) {
                result[i] = static_cast<block>(value);
                value >>= c_BLOCK_SIZE;

                if (value == 0) {
                    break;
                }
            }

            return result;
        }

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

        static uint_t divide(const uint_t& lhs, const block& rhs, uint_t* remainder = nullptr) {
            uint_t result;
            double_block part = 0;

            for (size_t i = c_BLOCK_NUMBER; i > 0; --i) {
                part = (part << (c_BLOCK_SIZE)) + static_cast<double_block>(lhs[i - 1]);

                if (part < rhs) {
                    continue;
                }

                result[i - 1] = static_cast<block>(part / rhs);
                part %= rhs;
            }

            if (remainder != nullptr) {
                *remainder = uint_t(static_cast<block>(part));
            }

            return result;
        }

        static uint_t d_divide(const uint_t& lhs, const uint_t& rhs, uint_t* remainder = nullptr) {
            size_t dividend_size = lhs.clz();
            size_t divisor_size = rhs.clz();

            uint_t<bits + c_BLOCK_SIZE> dividend(lhs);
            uint_t divisor(rhs);
            uint_t quotient;

            size_t shift_size = 0;
            block divisor_head = divisor[divisor_size - 1];
            static constexpr double_block c_HALF_BUCKET = static_cast<double_block>(1) << (c_BLOCK_SIZE - 1);

            while (divisor_head < c_HALF_BUCKET) {
                ++shift_size;
                divisor_head <<= 1;
            }

            dividend <<= shift_size;
            divisor <<= shift_size;

            double_block divisor_ = divisor[divisor_size - 1];
            static constexpr double_block c_BUCKET = static_cast<double_block>(1) << c_BLOCK_SIZE;

            for (size_t i = dividend_size - divisor_size + 1; i > 0; --i) {
                double_block part =
                    (static_cast<double_block>(dividend[i + divisor_size - 1]) << c_BLOCK_SIZE)
                    + static_cast<double_block>(dividend[i + divisor_size - 2]);
                double_block quotient_t = part / divisor_;
                part %= divisor_;

                if (quotient_t == c_BUCKET) {
                    --quotient_t;
                    part += divisor_;
                }

                while (part < c_BUCKET
                       && (quotient_t * divisor[divisor_size - 2]
                           > (part << c_BLOCK_SIZE) + dividend[i + divisor_size - 3])) {
                    --quotient_t;
                    part += divisor_;
                }

                int64_t carry = 0;
                int64_t widedigit = 0;

                for (size_t j = 0; j < divisor_size; ++j) {
                    double_block product =
                        static_cast<block>(quotient_t) * static_cast<double_block>(divisor[j]);
                    widedigit =
                        (static_cast<int64_t>(dividend[i + j - 1]) + carry) - (product & 0xffffffffLL);
                    dividend[i + j - 1] = static_cast<block>(widedigit);
                    carry = (widedigit >> c_BLOCK_SIZE) - static_cast<int64_t>(product >> c_BLOCK_SIZE);
                }

                widedigit = static_cast<int64_t>(dividend[i + divisor_size - 1]) + carry;
                dividend[i + divisor_size - 1] = static_cast<block>(widedigit);

                quotient[i - 1] = static_cast<block>(quotient_t);

                if (widedigit < 0) {
                    --quotient[i - 1];
                    widedigit = 0;

                    for (size_t j = 0; j < divisor_size; ++j) {
                        widedigit += static_cast<double_block>(dividend[i + j - 1]) + divisor[j];
                        dividend[i + j - 1] = static_cast<block>(widedigit);
                        widedigit >>= 32;
                    }
                }
            }

            if (remainder != nullptr) {
                *remainder = uint_t(0);

                for (size_t i = 0; i < divisor_size - 1; ++i) {
                    (*remainder)[i] =
                        (dividend[i] >> shift_size)
                        | (static_cast<double_block>(dividend[i + 1]) << (c_BLOCK_SIZE - shift_size));
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

        uint_t operator*(block other) const {   // FFT will delete this
            block remainder = 0;
            uint_t result;

            for (size_t i = 0; i < c_BLOCK_NUMBER; ++i) {
                double_block prod = m_blocks[i] * static_cast<double_block>(other);
                result[i] = static_cast<block>(prod) + remainder;
                remainder = static_cast<block>(result[i] < remainder) + static_cast<block>(prod >> 32);
            }

            return result;
        }

        size_t clz() const {
            size_t result = c_BLOCK_NUMBER;

            while (result > 0 && m_blocks[result - 1] == 0) {
                --result;
            }

            return result;
        }
    };
}   // namespace ECG

#endif
