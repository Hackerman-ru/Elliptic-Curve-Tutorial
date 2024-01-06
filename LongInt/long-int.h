#ifndef ECG_LONGINT_H
#define ECG_LONGINT_H

#ifndef CHAR_BIT
    #define CHAR_BIT 8
#endif

#include <array>
#include <cassert>
#include <complex>
#include <functional>
#include <string>

namespace ECG {
    template<size_t bits>
    class uint_t {
        using bucket_type = uint32_t;
        static_assert(std::is_unsigned_v<bucket_type>, "The bucket type must be unsigned int");
        static_assert(bits % sizeof(bucket_type) == 0,
                      "The bits must be an integer divisible by the size of the bucket type");
        static_assert(
            std::is_convertible_v<bucket_type, char>,
            "The bucket type must be convertible to char for string initialization and string representation");

    public:
        enum class StringType {
            BINARY = 0b1,
            DECIMAL = 10,
            HEXADECIMAL = 0xF,
        };

        constexpr uint_t() = default;
        constexpr explicit uint_t(bucket_type value);
        template<size_t other_bits>
        uint_t(const uint_t<other_bits>& other);
        uint_t(const std::string& str, StringType str_type);
        uint_t(const std::string& str, std::function<bucket_type(char)> map, size_t shift);

        auto operator<=>(const uint_t& other) const {
            return other.m_value <=> m_value;
        }

        uint_t operator+(const uint_t& other) const;
        uint_t operator-(const uint_t& other) const;
        uint_t operator*(const uint_t& other) const;
        uint_t operator/(const uint_t& other) const;
        uint_t operator%(const uint_t& other) const;
        uint_t operator/(bucket_type value) const;
        uint_t operator%(bucket_type value) const;
        uint_t operator>>(size_t shift) const;
        uint_t operator<<(size_t shift) const;
        uint_t operator^(const uint_t& other) const;
        uint_t operator|(const uint_t& other) const;
        uint_t operator&(const uint_t& other) const;

        uint_t operator-() const;

        uint_t& operator+=(const uint_t& other);
        uint_t& operator-=(const uint_t& other);
        uint_t& operator*=(const uint_t& other);
        uint_t& operator/(const uint_t& other);
        uint_t& operator%(const uint_t& other);
        uint_t& operator/=(bucket_type value);
        uint_t& operator%=(bucket_type value);
        uint_t& operator>>=(size_t shift);
        uint_t& operator<<=(size_t shift);
        uint_t& operator^=(const uint_t& other);
        uint_t& operator|=(const uint_t& other);
        uint_t& operator&=(const uint_t& other);

        uint_t& operator++(int);
        uint_t& operator++();
        uint_t& operator--();

        explicit operator bucket_type() const;

        // Experimental feature //
        template<typename T>
        requires std::is_nothrow_convertible_v<bucket_type, T>
        T convert_to() const {
            size_t shift = sizeof(T) * CHAR_BIT;
            size_t bucket_number = shift / c_BUCKET_SIZE;
            T result = {};

            for (size_t i = 0; i < m_value.size() && i < bucket_number; ++i) {
                result |= T(m_value[i]) << (i * c_BUCKET_SIZE);
            }

            shift %= c_BUCKET_SIZE;

            result |= ((T(m_value[bucket_number]) << (c_BUCKET_SIZE - shift)) >> (c_BUCKET_SIZE - shift))
                   << (bucket_number * c_BUCKET_SIZE);

            return result;
        }

        //

        std::string into_string(StringType str_type) const;
        std::string into_string(std::function<char(bucket_type)> map, size_t shift) const;

    private:
        static constexpr size_t c_BUCKET_SIZE = sizeof(bucket_type) * CHAR_BIT;
        std::array<bucket_type, bits / c_BUCKET_SIZE> m_value = {};
    };

    template<size_t bits>
    constexpr uint_t<bits>::uint_t(bucket_type value) : m_value({value}) {};

    template<size_t bits>
    template<size_t other_bits>
    inline uint_t<bits>::uint_t(const uint_t<other_bits>& other) {
        constexpr size_t min_buckets = min(bits, other_bits) / c_BUCKET_SIZE;
        *this = uint_t();

        for (size_t i = 0; i < min_buckets; i++) {
            m_value[i] = other.m_value[i];
        }
    }

    template<size_t bits>
    uint_t<bits>::uint_t(const std::string& str, StringType str_type) {
        *this = uint_t();

        switch (str_type) {
        case uint_t::StringType::BINARY :
            for (size_t i = 0; i < str.size() && i < c_BITS; ++i) {
                size_t bucket_pos = i / c_BUCKET_SIZE;

                m_value[bucket_pos] |= bucket_type(str[str.size() - 1 - i] - '0') << (i % c_BUCKET_SIZE);
            }
            break;
        case uint_t::StringType::DECIMAL :
            for (char c : str) {
                uint_t temp = (*this <<= 1);
                *this <<= 2;
                *this += temp;
                *this += uint_t(c - '0');
            }
            break;
        case uint_t::StringType::HEXADECIMAL :
            for (char c : str) {
                *this <<= 4;
                *this += uint_t(c - (isdigit(c) ? '0' : 'a'));
            }
            break;
        }
    }

    template<size_t bits>
    uint_t<bits>::uint_t(const std::string& str, std::function<bucket_type(char)> map, size_t shift) {
        *this = uint_t();

        for (char c : str) {
            if (isdigit(c)) {
                *this <<= shift;
                *this += uint_t(map(c));
            }
        }
    }

    template<size_t bits>
    uint_t<bits> uint_t<bits>::operator+(const uint_t& other) const {
        uint_t result;

        for (size_t i = 0; i < m_value.size(); ++i) {
            result.m_value[i] += m_value[i] + other.m_value[i];

            if (i + 1 < m_value.size() && result.m_value[i] < m_value[i]) {
                ++result.m_value[i + 1];
            }
        }

        return result;
    }

    template<size_t bits>
    uint_t<bits> uint_t<bits>::operator-(const uint_t& other) const {
        return *this + (-other);
    }

    namespace {
        using complex_double = std::complex<double>;
        static const double PI = acos(-1);

        void fft(std::vector<complex_double>& a, bool invert) {
            size_t n = a.size();

            if (n == 1) {
                return;
            }

            std::vector<complex_double> a0(n / 2);
            std::vector<complex_double> a1(n / 2);

            for (size_t i = 0; 2 * i < n; i++) {
                a0[i] = a[2 * i];
                a1[i] = a[2 * i + 1];
            }

            fft(a0, invert);
            fft(a1, invert);

            double ang = 2 * PI / n * (invert ? -1 : 1);
            complex_double w(1), wn(cos(ang), sin(ang));

            for (size_t i = 0; 2 * i < n; i++) {
                a[i] = a0[i] + w * a1[i];
                a[i + n / 2] = a0[i] - w * a1[i];

                if (invert) {
                    a[i] /= 2;
                    a[i + n / 2] /= 2;
                }

                w *= wn;
            }
        }

        template<typename T>
        std::vector<T> multiply(std::vector<T> const& a, std::vector<T> const& b) {
            std::vector<complex_double> fa(a.begin(), a.end());
            std::vector<complex_double> fb(b.begin(), b.end());
            size_t n = 1;

            while (n < a.size() + b.size()) {
                n <<= 1;
            }

            fa.resize(n);
            fb.resize(n);

            fft(fa, false);
            fft(fb, false);

            for (size_t i = 0; i < n; i++) {
                fa[i] *= fb[i];
            }

            fft(fa, true);

            std::vector<T> result(n);

            for (size_t i = 0; i < n; i++) {
                result[i] = T(round(fa[i].real()));
            }

            return result;
        }
    }   // namespace

    template<size_t bits>
    uint_t<bits> uint_t<bits>::operator*(const uint_t& other) const {
        std::vector<bucket_type> lhs(m_value.begin(), m_value.end());
        std::vector<bucket_type> rhs(other.m_value.begin(), other.m_value.end());
        std::vector<bucket_type> product = multiply(lhs, rhs);
        uint_t result;

        for (size_t i = 0; i < m_value.size(); ++i) {
            result.m_value[i] = product[i];
        }

        return result;

        // If FFT won't work, then long multiplication should be completed
        /*static constexpr size_t HALF_BUCKET_SIZE = c_BUCKET_SIZE >> 1;
    static_assert((HALF_BUCKET_SIZE & 0b1) == 0, "Integer type must have even number of bits");

    uint512_t result {};

    auto get_lower_part = [](T value) {
        return (value << HALF_BUCKET_SIZE) >> HALF_BUCKET_SIZE;
    };

    auto get_higher_part = [](T value) {
        return value >> (c_BUCKET_SIZE - HALF_BUCKET_SIZE);
    };

    T carry = 0;

    for (size_t i = 0; i < m_value.size(); ++i) {
        T lhs_high = get_higher_part(m_value[i]);
        T lhs_low = get_lower_part(m_value[i]);
        T rhs_high = get_higher_part(other.m_value[i]);
        T rhs_low = get_lower_part(other.m_value[i]);

        T low_low = lhs_low * rhs_low;
        T low_high = lhs_low * rhs_high;
        T high_low = lhs_high * rhs_low;
        T high_high = lhs_high * rhs_high;

        T low_carry = (((low_high << HALF_BUCKET_SIZE) >> HALF_BUCKET_SIZE)
                       + ((high_low << HALF_BUCKET_SIZE) >> HALF_BUCKET_SIZE) + (low_low >> HALF_BUCKET_SIZE))
                   >> HALF_BUCKET_SIZE;

        result.m_value[i] = carry;
        carry = high_high + (low_high >> HALF_BUCKET_SIZE) + (high_low >> HALF_BUCKET_SIZE) + low_carry;
        T add = m_value[i] * other.m_value[i];
        result.m_value[i] += add;

        if (result.m_value[i] < add) {
            carry++;
        }
        ...
    }*/
    }

    template<size_t bits>
    uint_t<bits> uint_t<bits>::operator/(const uint_t& other) const {
        assert(other != uint_t(0) && "Division by 0");

        size_t new_bits = bits + c_BUCKET_SIZE;
        size_t bucket_number = new_bits / c_BUCKET_SIZE;

        uint_t<new_bits> dividend(*this);
        uint_t<new_bits> divisor_t(other);
        uint_t<new_bits> quotient;

        size_t divisor_head = 0;

        for (size_t i = bucket_number; i > 0; --i) {
            if (divisor_t.m_value[i - 1] != 0) {
                divisor_head = i - 1;
                break;
            }
        }

        static constexpr bucket_type half_of_bucket = 1 << (c_BUCKET_SIZE >> 1);

        while (divisor_t.m_value[divisor_head] < half_of_bucket) {
            dividend <<= 1;
            divisor_t <<= 1;
        }

        dividend >>= divisor_head * c_BUCKET_SIZE;
        divisor_t >>= divisor_head * c_BUCKET_SIZE;

        bucket_type divisor = static_cast<bucket_type>(divisor_t);
    }

    template<size_t bits>
    uint_t<bits> uint_t<bits>::operator%(const uint_t& other) const {
        return uint_t();
    }

    template<size_t bits>
    uint_t<bits> uint_t<bits>::operator/(bucket_type value) const {
        assert(value != 0 && "Division by 0");
        uint_t result;

        if constexpr (std::is_same_v<bucket_type, uint32_t>) {
            uint64_t remainder = 0;
            uint32_t divisor = value;

            for (size_t i = m_value.size(); i > 0; --i) {
                remainder += m_value[i - 1];

                if (remainder < divisor) {
                    continue;
                }

                result.m_value[i - 1] = static_cast<bucket_type>(remainder / divisor);
                remainder %= divisor;
            }
        } else {
            // Unlikely, but should be done to unify the code
        }

        return result;
    }

    template<size_t bits>
    uint_t<bits> uint_t<bits>::operator%(bucket_type value) const {
        assert(value != 0 && "Division by 0");
        return *this - (*this / value) * uint_t(value);
    }

    template<size_t bits>
    uint_t<bits> uint_t<bits>::operator>>(size_t shift) const {
        size_t bucket_shift = shift / c_BUCKET_SIZE;
        uint_t result;

        for (size_t i = 0; i + bucket_shift < m_value.size(); ++i) {
            result.m_value[i] = m_value[i + bucket_shift];
        }

        shift %= c_BUCKET_SIZE;

        for (size_t i = 0; i + bucket_shift < m_value.size(); ++i) {
            result.m_value[i] >>= shift;

            if (i + 1 < m_value.size()) {
                result.m_value[i] |= m_value[i + 1] << (c_BUCKET_SIZE - shift);
            }
        }

        return result;
    }

    template<size_t bits>
    uint_t<bits> uint_t<bits>::operator<<(size_t shift) const {
        size_t bucket_shift = shift / c_BUCKET_SIZE;
        uint_t result;

        for (size_t i = 0; i + bucket_shift < m_value.size(); ++i) {
            result.m_value[i + bucket_shift] = m_value[i];
        }

        shift %= c_BUCKET_SIZE;

        for (size_t i = m_value.size(); i > bucket_shift; --i) {
            result.m_value[i - 1] <<= shift;

            if (i - 1 > 0) {
                result.m_value[i - 1] |= m_value[i - 2] >> (c_BUCKET_SIZE - shift);
            }
        }

        return result;
    }

    template<size_t bits>
    uint_t<bits> uint_t<bits>::operator^(const uint_t& other) const {
        uint_t result;

        for (size_t i = 0; i < m_value.size(); ++i) {
            result.m_value[i] = m_value[i] ^ other.m_value[i];
        }

        return result;
    }

    template<size_t bits>
    uint_t<bits> uint_t<bits>::operator|(const uint_t& other) const {
        uint_t result;

        for (size_t i = 0; i < m_value.size(); ++i) {
            result.m_value[i] = m_value[i] | other.m_value[i];
        }

        return result;
    }

    template<size_t bits>
    uint_t<bits> uint_t<bits>::operator&(const uint_t& other) const {
        uint_t result;

        for (size_t i = 0; i < m_value.size(); ++i) {
            result.m_value[i] = m_value[i] & other.m_value[i];
        }

        return result;
    }

    template<size_t bits>
    uint_t<bits> uint_t<bits>::operator-() const {
        uint_t result;

        for (size_t i = 0; i < m_value.size(); ++i) {
            result.m_value[i] = ~(m_value[i]);
        }

        result++;

        return result;
    }

    template<size_t bits>
    uint_t<bits>& uint_t<bits>::operator+=(const uint_t& other) {
        return (*this = *this + other);
    }

    template<size_t bits>
    uint_t<bits>& uint_t<bits>::operator-=(const uint_t& other) {
        return (*this = *this - other);
    }

    template<size_t bits>
    uint_t<bits>& uint_t<bits>::operator*=(const uint_t& other) {
        return (*this = *this * other);
    }

    template<size_t bits>
    uint_t<bits>& uint_t<bits>::operator/=(bucket_type value) {
        return (*this = *this / value);
    }

    template<size_t bits>
    uint_t<bits>& uint_t<bits>::operator%=(bucket_type value) {
        return (*this = *this % value);
    }

    template<size_t bits>
    uint_t<bits>& uint_t<bits>::operator>>=(size_t shift) {
        return (*this = *this >> shift);
    }

    template<size_t bits>
    uint_t<bits>& uint_t<bits>::operator<<=(size_t shift) {
        return (*this = *this << shift);
    }

    template<size_t bits>
    uint_t<bits>& uint_t<bits>::operator^=(const uint_t& other) {
        return (*this = *this ^ other);
    }

    template<size_t bits>
    uint_t<bits>& uint_t<bits>::operator|=(const uint_t& other) {
        return (*this = *this | other);
    }

    template<size_t bits>
    uint_t<bits>& uint_t<bits>::operator&=(const uint_t& other) {
        return (*this = *this & other);
    }

    template<size_t bits>
    uint_t<bits>& uint_t<bits>::operator++(int) {
        return (*this += uint_t(1));
    }

    template<size_t bits>
    uint_t<bits>& uint_t<bits>::operator++() {
        return (*this += uint_t(1));
    }

    template<size_t bits>
    uint_t<bits>& uint_t<bits>::operator--() {
        return (*this -= uint_t(1));
    }

    template<size_t bits>
    uint_t<bits>::operator bucket_type() const {
        return m_value[0];
    }

    template<size_t bits>
    std::string uint_t<bits>::into_string(StringType str_type) const {
        std::string result;
        uint_t clone_of_this = *this;

        switch (str_type) {
        case uint_t::StringType::BINARY :
            do {
                result.push_back((static_cast<bucket_type>(clone_of_this)
                                  & static_cast<bucket_type>(uint_t::StringType::BINARY))
                                 + '0');
                clone_of_this >>= 1;
            } while (clone_of_this > uint_t(0));
            break;
        case uint_t::StringType::DECIMAL :
            do {
                result.push_back((static_cast<bucket_type>(clone_of_this)
                                  % static_cast<bucket_type>(uint_t::StringType::DECIMAL))
                                 + '0');
                clone_of_this /= bucket_type(uint_t::StringType::DECIMAL);
            } while (clone_of_this > uint_t(0));
            break;
        case uint_t::StringType::HEXADECIMAL :
            do {
                result.push_back((static_cast<bucket_type>(clone_of_this)
                                  & static_cast<bucket_type>(uint_t::StringType::HEXADECIMAL))
                                 + '0');
                clone_of_this >>= 4;
            } while (clone_of_this > uint_t(0));
            break;
        }

        std::reverse(result.begin(), result.end());
        return result;
    }

    template<size_t bits>
    std::string uint_t<bits>::into_string(std::function<char(bucket_type)> map, size_t shift) const {
        std::string result;
        uint_t clone = *this;

        do {
            result.push_back(map((bucket_type(clone) << (c_BUCKET_SIZE - shift)) >> (c_BUCKET_SIZE - shift)));
            clone >>= shift;
        } while (clone > uint_t(0));

        std::reverse(result.begin(), result.end());
        return result;
    }
}   // namespace ECG

#endif
