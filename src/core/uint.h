#ifndef ECG_LONGINT_H
#define ECG_LONGINT_H

#include "fft.h"
#include "string-parser.h"
#include "util.h"

#include <array>
#include <bitset>
#include <cassert>
#include <string>

namespace ECG {
    template<size_t bits>
    class uint_t {
        using block_t = uint32_t;
        using double_block_t = uint64_t;   // FFT will delete this

        static constexpr size_t c_BITS_IN_BYTE = 8;
        static constexpr size_t c_BLOCK_SIZE = sizeof(block_t) * c_BITS_IN_BYTE;
        static constexpr size_t c_BLOCK_NUMBER = bits / c_BLOCK_SIZE;

        using blocks = std::array<block_t, c_BLOCK_NUMBER>;
        blocks m_blocks = {};

        template<size_t V>
        friend class uint_t;

    public:
        constexpr uint_t() = default;

        template<typename T>
        constexpr uint_t(const T& value) : m_blocks(split_into_blocks<T>(value)) {}

        constexpr uint_t& operator=(const char* str) {
            *this = parse_into<uint_t>(str);
            return *this;
        }

        // operator<=>
        friend auto operator<=>(const uint_t& lhs, const uint_t& rhs) = default;

        // operator+
        friend uint_t operator+(const uint_t& lhs, const uint_t& rhs) {
            uint_t result = lhs;
            lhs += rhs;
            return result;
        }

        friend uint_t operator+(uint_t&& lhs, const uint_t& rhs) {
            lhs += rhs;
            return lhs;
        }

        friend uint_t operator+(const uint_t& lhs, uint_t&& rhs) {
            rhs += lhs;
            return rhs;
        }

        friend uint_t operator+(uint_t&& lhs, uint_t&& rhs) {
            lhs += rhs;
            return lhs;
        }

        // operator-
        friend uint_t operator-(const uint_t& lhs, const uint_t& rhs) {
            uint_t result = lhs;
            result -= rhs;
            return result;
        }

        friend uint_t operator-(uint_t&& lhs, const uint_t& rhs) {
            lhs -= rhs;
            return lhs;
        }

        friend uint_t operator-(const uint_t& lhs, uint_t&& rhs) {
            rhs -= lhs;
            rhs.neg();
            return rhs;
        }

        friend uint_t operator-(uint_t&& lhs, uint_t&& rhs) {
            lhs -= rhs;
            return lhs;
        }

        // operator*
        friend uint_t operator*(const uint_t& lhs, const uint_t& rhs) {
            return multiply<c_BLOCK_NUMBER>(lhs.m_blocks, rhs.m_blocks);
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
            uint_t result = *this;
            result.neg();
            return result;
        }

        uint_t& operator+=(const uint_t& other) {
            block_t carry = 0;

            for (size_t i = 0; i < c_BLOCK_NUMBER; ++i) {
                block_t sum = carry + other[i];
                m_blocks[i] += sum;
                carry = (m_blocks[i] < sum) || (sum < carry);
            }

            return *this;
        }

        uint_t& operator-=(const uint_t& other) {
            block_t remainder = 0;

            for (size_t i = 0; i < c_BLOCK_NUMBER; ++i) {
                block_t temp = m_blocks[i];
                m_blocks[i] -= other[i] + remainder;
                remainder = (temp < m_blocks[i]);
            }

            return *this;
        }

        uint_t& operator*=(const uint_t& other) {
            return *this = *this * other;
        }

        uint_t& operator/=(const uint_t& other) {
            return *this = *this / other;
        }

        uint_t& operator%=(const uint_t& other) {
            return *this = *this % other;
        }

        uint_t& operator>>=(size_t shift_size) {
            static constexpr size_t BUCKET_SIZE = 64;
            static constexpr size_t BUCKET_NUMBER = bits / BUCKET_SIZE;

            size_t bucket_shift = shift_size >> 6;
            auto data = reinterpret_cast<double_block_t*>(m_blocks.data());

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
            auto data = reinterpret_cast<double_block_t*>(m_blocks.data());

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

        uint_t operator++(int) {
            uint_t result = *this;
            *this += 1;
            return result;
        }

        uint_t& operator++() {
            return *this += 1;
        }

        uint_t operator--(int) {
            uint_t result = *this;
            *this -= 1;
            return result;
        }

        uint_t& operator--() {
            return *this -= 1;
        }

        template<typename T>
        T convert_to() const;

        template<typename T>
        requires is_convertible_to<T, block_t>
        T convert_to() const {
            size_t shift_size = sizeof(T) * c_BITS_IN_BYTE;
            size_t blocks_number = shift_size / c_BLOCK_SIZE;
            T result = 0;

            for (size_t i = 0; i < c_BLOCK_NUMBER && i < blocks_number; ++i) {
                result |= T(m_blocks[i]) << (i * c_BLOCK_SIZE);
            }

            return result;
        }

        template<>
        std::string convert_to() const {
            std::string result;
            uint_t clone_of_this = *this;

            do {
                uint_t remainder;
                clone_of_this = divide(clone_of_this, 10, &remainder);
                result.push_back(remainder.m_blocks[0] + '0');
            } while (clone_of_this > 0);

            std::reverse(result.begin(), result.end());
            return result;
        }

    private:
        static constexpr size_t size() {
            return c_BLOCK_NUMBER;
        }

        template<typename T>
        requires std::numeric_limits<T>::is_integer && is_upcastable_to<T, block_t>
        static constexpr blocks split_into_blocks(T value) {
            return {static_cast<block_t>(value)};
        }

        template<typename T>
        requires std::numeric_limits<T>::is_integer && is_downcastable_to<T, block_t>
        static constexpr blocks split_into_blocks(T value) {
            blocks result = {};

            for (size_t i = 0; i < c_BLOCK_NUMBER; ++i) {
                result[i] = static_cast<block_t>(value);
                value >>= c_BLOCK_SIZE;

                if (value == 0) {
                    break;
                }
            }

            return result;
        }

        template<typename T>
        requires is_convertible_container<T, block_t> || requires(T x) {
            { ECG::uint_t {x} } -> std::same_as<T>;
        }
        static constexpr blocks split_into_blocks(const T& other) {
            const size_t min_size = std::min(size(), other.size());
            blocks result = {};

            for (size_t i = 0; i < min_size; i++) {
                result[i] = static_cast<block_t>(other[i]);
            }

            return result;
        }

        static uint_t divide(const uint_t& lhs,
                             const uint_t& rhs,
                             uint_t* remainder = nullptr) {   // FFT will delete this
            size_t dividend_size = lhs.actual_size();
            size_t divisor_size = rhs.actual_size();

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

        static uint_t divide(const uint_t& lhs,
                             const block_t& rhs,
                             uint_t* remainder = nullptr) {   // FFT will delete this
            uint_t result;
            double_block_t part = 0;

            for (size_t i = c_BLOCK_NUMBER; i > 0; --i) {
                part = (part << (c_BLOCK_SIZE)) + static_cast<double_block_t>(lhs[i - 1]);

                if (part < rhs) {
                    continue;
                }

                result[i - 1] = static_cast<block_t>(part / rhs);
                part %= rhs;
            }

            if (remainder != nullptr) {
                *remainder = uint_t(static_cast<block_t>(part));
            }

            return result;
        }

        static uint_t d_divide(const uint_t& lhs,
                               const uint_t& rhs,
                               uint_t* remainder = nullptr) {   // FFT will delete this
            size_t dividend_size = lhs.actual_size();
            size_t divisor_size = rhs.actual_size();

            uint_t<bits + c_BLOCK_SIZE> dividend(lhs);
            uint_t divisor(rhs);
            uint_t quotient;

            size_t shift_size = 0;
            block_t divisor_head = divisor[divisor_size - 1];
            static constexpr double_block_t c_HALF_BUCKET = static_cast<double_block_t>(1)
                                                         << (c_BLOCK_SIZE - 1);

            while (divisor_head < c_HALF_BUCKET) {
                ++shift_size;
                divisor_head <<= 1;
            }

            dividend <<= shift_size;
            divisor <<= shift_size;

            double_block_t divisor_ = divisor[divisor_size - 1];
            static constexpr double_block_t c_BUCKET = static_cast<double_block_t>(1) << c_BLOCK_SIZE;

            for (size_t i = dividend_size - divisor_size + 1; i > 0; --i) {
                double_block_t part =
                    (static_cast<double_block_t>(dividend[i + divisor_size - 1]) << c_BLOCK_SIZE)
                    + static_cast<double_block_t>(dividend[i + divisor_size - 2]);
                double_block_t quotient_t = part / divisor_;
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
                    double_block_t product =
                        static_cast<block_t>(quotient_t) * static_cast<double_block_t>(divisor[j]);
                    widedigit =
                        (static_cast<int64_t>(dividend[i + j - 1]) + carry) - (product & 0xffffffffLL);
                    dividend[i + j - 1] = static_cast<block_t>(widedigit);
                    carry = (widedigit >> c_BLOCK_SIZE) - static_cast<int64_t>(product >> c_BLOCK_SIZE);
                }

                widedigit = static_cast<int64_t>(dividend[i + divisor_size - 1]) + carry;
                dividend[i + divisor_size - 1] = static_cast<block_t>(widedigit);

                quotient[i - 1] = static_cast<block_t>(quotient_t);

                if (widedigit < 0) {
                    --quotient[i - 1];
                    widedigit = 0;

                    for (size_t j = 0; j < divisor_size; ++j) {
                        widedigit += static_cast<double_block_t>(dividend[i + j - 1]) + divisor[j];
                        dividend[i + j - 1] = static_cast<block_t>(widedigit);
                        widedigit >>= 32;
                    }
                }
            }

            if (remainder != nullptr) {
                *remainder = uint_t(0);

                for (size_t i = 0; i < divisor_size - 1; ++i) {
                    (*remainder)[i] =
                        (dividend[i] >> shift_size)
                        | (static_cast<double_block_t>(dividend[i + 1]) << (c_BLOCK_SIZE - shift_size));
                }

                (*remainder)[divisor_size - 1] = dividend[divisor_size - 1] >> shift_size;
            }

            return quotient;
        }

        void neg() {
            for (size_t i = 0; i < c_BLOCK_NUMBER; ++i) {
                m_blocks[i] = ~(m_blocks[i]);
            }

            ++*this;
        }

        uint_t operator*(block_t other) const {   // FFT will delete this
            block_t remainder = 0;
            uint_t result;

            for (size_t i = 0; i < c_BLOCK_NUMBER; ++i) {
                double_block_t prod = m_blocks[i] * static_cast<double_block_t>(other);
                result[i] = static_cast<block_t>(prod) + remainder;
                remainder = static_cast<block_t>(result[i] < remainder) + static_cast<block_t>(prod >> 32);
            }

            return result;
        }

        const block_t& operator[](size_t pos) const {
            return m_blocks[pos];
        }

        block_t& operator[](size_t pos) {
            return m_blocks[pos];
        }

        size_t actual_size() const {   // FFT will delete this
            size_t result = c_BLOCK_NUMBER;

            while (result > 0 && m_blocks[result - 1] == 0) {
                --result;
            }

            return result;
        }
    };
}   // namespace ECG

#endif
