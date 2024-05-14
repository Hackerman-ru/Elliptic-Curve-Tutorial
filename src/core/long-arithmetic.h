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
        using block_t = uint32_t;
        using double_block_t = uint64_t;

        static constexpr size_t c_bits_in_byte = 8;
        static constexpr size_t c_block_size = sizeof(block_t) * c_bits_in_byte;
        static constexpr size_t c_block_number = c_bits / c_block_size;

        template<size_t V>
        friend class uint_t;

        using blocks = std::array<block_t, c_block_number>;
        blocks m_blocks = {};

    public:
        constexpr uint_t() = default;

        template<typename T>
        constexpr uint_t(const T& value) : m_blocks(split_into_blocks<T>(value)) {}

        constexpr uint_t(const char* str) : m_blocks(algorithm::parse_into<uint_t>(str).m_blocks) {};

        constexpr uint_t& operator=(const uint_t& value) = default;

        template<typename T>
        constexpr uint_t& operator=(const T& value) {
            m_blocks = split_into_blocks<T>(value);
            return *this;
        }

        constexpr uint_t& operator=(const char* str) {
            return *this = parse_into<uint_t>(str);
        }

        friend constexpr std::strong_ordering operator<=>(const uint_t& lhs, const uint_t& rhs) {
            for (size_t i = c_block_number; i > 0; --i) {
                if (lhs[i - 1] != rhs[i - 1]) {
                    return lhs[i - 1] <=> rhs[i - 1];
                }
            }

            return std::strong_ordering::equal;
        }

        friend constexpr bool operator==(const uint_t& lhs, const uint_t& rhs) {
            return lhs.m_blocks == rhs.m_blocks;
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

            for (size_t i = 0; i < c_block_number; ++i) {
                result += (lhs * rhs[i]) << (c_block_size * i);
            }

            return result;
            //return algorithm::fast_fourier_transform::multiply<c_block_number>(lhs.m_blocks, rhs.m_blocks);
        }

        // operator/
        friend constexpr uint_t operator/(const uint_t& lhs,
                                          const uint_t& rhs) {   // fast_fourier_transform will change this
            uint_t result = divide(lhs, rhs);
            uint_t less = result * rhs;
            uint_t greater = (result + 1) * rhs;
            if (less > lhs || greater <= lhs) {
                result = 0;
            }
            assert(result * rhs <= lhs && (result + 1) * rhs > lhs
                   && "uint_t::operator/ : remainder must be less than divisor");
            return result;
        }

        // operator%
        friend constexpr uint_t operator%(const uint_t& lhs,
                                          const uint_t& rhs) {   // fast_fourier_transform will change this
            uint_t remainder;
            divide(lhs, rhs, &remainder);
            assert(rhs > remainder && "uint_t::operator% : remainder must be less than divisor");
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
            block_t carry = 0;

            for (size_t i = 0; i < c_block_number; ++i) {
                block_t sum = carry + other[i];
                m_blocks[i] += sum;
                carry = (m_blocks[i] < sum) || (sum < carry);
            }

            return *this;
        }

        constexpr uint_t& operator-=(const uint_t& other) {
            block_t remainder = 0;

            for (size_t i = 0; i < c_block_number; ++i) {
                block_t prev = m_blocks[i];
                block_t sum = other[i] + remainder;
                m_blocks[i] -= sum;
                remainder = (m_blocks[i] > prev) || (sum < remainder);
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
            static constexpr size_t c_double_bucket_size = sizeof(double_block_t) * c_bits_in_byte;
            static constexpr size_t c_double_bucket_number = c_bits / c_double_bucket_size;

            size_t bucket_shift = shift_size >> 6;
            auto data = reinterpret_cast<double_block_t*>(m_blocks.data());

            if (bucket_shift > 0) {
                for (size_t i = 0; i < c_double_bucket_number; ++i) {
                    if (i + bucket_shift < c_double_bucket_number) {
                        data[i] = data[i + bucket_shift];
                    } else {
                        data[i] = 0;
                    }
                }
            }

            shift_size %= c_double_bucket_size;

            if (shift_size == 0) {
                return *this;
            }

            for (size_t i = 0; i + bucket_shift < c_double_bucket_number; ++i) {
                data[i] >>= shift_size;

                if (i + 1 < c_double_bucket_number) {
                    data[i] |= data[i + 1] << (c_double_bucket_size - shift_size);
                }
            }

            if constexpr (c_bits % c_double_bucket_size != 0) {
                if constexpr (c_block_number > 1) {
                    m_blocks[c_block_number - 2] |= m_blocks[c_block_number - 1]
                                                 << (c_block_size - shift_size);
                }
                m_blocks[c_block_number - 1] >>= shift_size;
            }

            return *this;
        }

        constexpr uint_t& operator<<=(size_t shift_size) {
            static constexpr size_t c_double_bucket_size = sizeof(double_block_t) * c_bits_in_byte;
            static constexpr size_t c_double_bucket_number = c_bits / c_double_bucket_size;

            size_t bucket_shift = shift_size >> 6;
            auto data = reinterpret_cast<double_block_t*>(m_blocks.data());

            if (bucket_shift > 0) {
                for (size_t i = c_double_bucket_number; i > 0; --i) {
                    if (i > bucket_shift) {
                        data[i - 1] = data[i - bucket_shift - 1];
                    } else {
                        data[i - 1] = 0;
                    }
                }
            }

            shift_size %= c_double_bucket_size;

            if (shift_size == 0) {
                return *this;
            }

            if constexpr (c_bits % c_double_bucket_size != 0) {
                m_blocks[c_block_number - 1] <<= shift_size;

                if constexpr (c_block_number > 1) {
                    m_blocks[c_block_number - 1] |=
                        m_blocks[c_block_number - 2] >> (c_block_size - shift_size);
                }
            }

            for (size_t i = c_double_bucket_number; i > bucket_shift; --i) {
                data[i - 1] <<= shift_size;

                if (i - 1 > 0) {
                    data[i - 1] |= data[i - 2] >> (c_double_bucket_size - shift_size);
                }
            }

            return *this;
        }

        constexpr uint_t& operator^=(const uint_t& other) {
            for (size_t i = 0; i < c_block_number; ++i) {
                m_blocks[i] ^= other[i];
            }

            return *this;
        }

        constexpr uint_t& operator|=(const uint_t& other) {
            for (size_t i = 0; i < c_block_number; ++i) {
                m_blocks[i] |= other[i];
            }

            return *this;
        }

        constexpr uint_t& operator&=(const uint_t& other) {
            for (size_t i = 0; i < c_block_number; ++i) {
                m_blocks[i] &= other[i];
            }

            return *this;
        }

        constexpr uint_t operator++(int) {
            uint_t result = *this;
            increment();
            return result;
        }

        constexpr uint_t& operator++() {
            increment();
            return *this;
        }

        constexpr uint_t operator--(int) {
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
        requires concepts::is_convertible_to<T, block_t>
        constexpr T convert_to() const {
            size_t shift_size = sizeof(T) * c_bits_in_byte;
            size_t blocks_number = shift_size / c_block_size;

            if (blocks_number == 0) {
                return static_cast<T>(m_blocks[0]);
            }

            T result = 0;

            for (size_t i = 0; i < c_block_number && i < blocks_number; ++i) {
                result |= static_cast<T>(m_blocks[i]) << (i * c_block_size);
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
                result.push_back(remainder.m_blocks[0] + '0');
            } while (clone_of_this > 0);

            std::reverse(result.begin(), result.end());
            return result;
        }

    private:
        static constexpr size_t size() {
            return c_block_number;
        }

        template<typename T>
        requires std::numeric_limits<T>::is_integer && concepts::is_upcastable_to<T, block_t>
        static constexpr blocks split_into_blocks(T value) {
            return {static_cast<block_t>(value)};
        }

        template<typename T>
        requires std::numeric_limits<T>::is_integer && concepts::is_downcastable_to<T, block_t>
        static constexpr blocks split_into_blocks(T value) {
            blocks result = {};

            for (size_t i = 0; i < c_block_number; ++i) {
                result[i] = static_cast<block_t>(value);
                value >>= c_block_size;

                if (value == 0) {
                    break;
                }
            }

            return result;
        }

        template<typename T>
        requires concepts::is_convertible_container<T, block_t> || requires(T x) {
            { elliptic_curve_guide::uint_t {x} } -> std::same_as<T>;
        }
        static constexpr blocks split_into_blocks(const T& other) {
            const size_t min_size = std::min(size(), other.size());
            blocks result = {};

            for (size_t i = 0; i < min_size; i++) {
                result[i] = static_cast<block_t>(other[i]);
            }

            return result;
        }

        static constexpr uint_t divide(const uint_t& lhs, const uint_t& rhs, uint_t* remainder = nullptr) {
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

        static constexpr uint_t divide(const uint_t& lhs, const block_t& rhs, uint_t* remainder = nullptr) {
            uint_t result;
            double_block_t part = 0;

            for (size_t i = c_block_number; i > 0; --i) {
                part = (part << (c_block_size)) + static_cast<double_block_t>(lhs[i - 1]);

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

        static constexpr uint_t d_divide(const uint_t& lhs, const uint_t& rhs, uint_t* remainder = nullptr) {
            size_t dividend_size = lhs.actual_size();
            size_t divisor_size = rhs.actual_size();

            uint_t<c_bits + c_block_size> dividend(lhs);
            uint_t divisor(rhs);
            uint_t quotient;

            size_t shift_size = 0;
            block_t divisor_head = divisor[divisor_size - 1];
            static constexpr double_block_t c_HalfBlock = static_cast<double_block_t>(1)
                                                       << (c_block_size - 1);

            while (divisor_head < c_HalfBlock) {
                ++shift_size;
                divisor_head <<= 1;
            }

            dividend <<= shift_size;
            divisor <<= shift_size;

            double_block_t divisor_ = divisor[divisor_size - 1];
            static constexpr double_block_t c_Block = static_cast<double_block_t>(1) << c_block_size;

            for (size_t i = dividend_size - divisor_size + 1; i > 0; --i) {
                double_block_t part =
                    (static_cast<double_block_t>(dividend[i + divisor_size - 1]) << c_block_size)
                    + static_cast<double_block_t>(dividend[i + divisor_size - 2]);
                double_block_t quotient_temp = part / divisor_;
                part %= divisor_;

                if (quotient_temp == c_Block) {
                    --quotient_temp;
                    part += divisor_;
                }

                while (part < c_Block
                       && (quotient_temp * divisor[divisor_size - 2]
                           > (part << c_block_size) + dividend[i + divisor_size - 3])) {
                    --quotient_temp;
                    part += divisor_;
                }

                int64_t carry = 0;
                int64_t widedigit = 0;

                for (size_t j = 0; j < divisor_size; ++j) {
                    double_block_t product =
                        static_cast<block_t>(quotient_temp) * static_cast<double_block_t>(divisor[j]);
                    widedigit = (static_cast<int64_t>(dividend[i + j - 1]) + carry) - (product & UINT32_MAX);
                    dividend[i + j - 1] = static_cast<block_t>(widedigit);
                    carry = (widedigit >> c_block_size) - static_cast<int64_t>(product >> c_block_size);
                }

                widedigit = static_cast<int64_t>(dividend[i + divisor_size - 1]) + carry;
                dividend[i + divisor_size - 1] = static_cast<block_t>(widedigit);

                quotient[i - 1] = static_cast<block_t>(quotient_temp);

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
                        | (static_cast<double_block_t>(dividend[i + 1]) << (c_block_size - shift_size));
                }

                (*remainder)[divisor_size - 1] = dividend[divisor_size - 1] >> shift_size;
            }

            return quotient;
        }

        void constexpr negative() {
            for (size_t i = 0; i < c_block_number; ++i) {
                m_blocks[i] = ~(m_blocks[i]);
            }

            ++*this;
        }

        void constexpr increment() {
            for (size_t i = 0; i < c_block_number; ++i) {
                m_blocks[i] += 1;

                if (m_blocks[i] != 0) {
                    break;
                }
            }
        }

        void constexpr decrement() {
            for (size_t i = 0; i < c_block_number; ++i) {
                block_t temp = m_blocks[i];
                m_blocks[i] -= 1;

                if (temp >= m_blocks[i]) {
                    break;
                }
            }
        }

        constexpr uint_t operator*(block_t other) const {
            block_t remainder = 0;
            uint_t result;

            for (size_t i = 0; i < c_block_number; ++i) {
                double_block_t prod = m_blocks[i] * static_cast<double_block_t>(other);
                result[i] = static_cast<block_t>(prod) + remainder;
                remainder = static_cast<block_t>(result[i] < remainder) + static_cast<block_t>(prod >> 32);
            }

            return result;
        }

        constexpr const block_t& operator[](size_t pos) const {
            return m_blocks[pos];
        }

        constexpr block_t& operator[](size_t pos) {
            return m_blocks[pos];
        }

        constexpr size_t actual_size() const {
            size_t result = c_block_number;

            while (result > 0 && m_blocks[result - 1] == 0) {
                --result;
            }

            return result;
        }
    };
}   // namespace elliptic_curve_guide

#endif
