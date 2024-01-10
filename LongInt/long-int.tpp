#ifndef ECG_LONG_INT_HPP
#define ECG_LONG_INT_HPP
#include "long-int.h"

using ECG::uint_t;

template<size_t bits>
inline constexpr uint_t<bits>::uint_t(const bucket_type& value) {
    m_buckets[0] = value;
}

template<size_t bits>
uint_t<bits>::uint_t(const std::string& str, StringType str_type) {
    switch (str_type) {
    case uint_t::StringType::BINARY :
        for (size_t i = 0; i < str.size() && i < bits; ++i) {
            size_t bucket_pos = i / c_BUCKET_SIZE;

            m_buckets[bucket_pos] |= bucket_type(str[str.size() - 1 - i] - '0') << (i % c_BUCKET_SIZE);
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

            if (isdigit(c)) {
                *this += uint_t(c - '0');
            } else {
                *this += uint_t((c - 'a') + 10);
            }
        }
        break;
    }
}

template<size_t bits>
uint_t<bits>::uint_t(const std::string& str, std::function<bucket_type(char)> map, size_t shift) {
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

    for (size_t i = 0; i < c_BUCKET_NUMBER; ++i) {
        result[i] += m_buckets[i] + other[i];

        if (i + 1 < c_BUCKET_NUMBER && result[i] < m_buckets[i]) {
            ++result[i + 1];
        }
    }

    return result;
}

template<size_t bits>
uint_t<bits> uint_t<bits>::operator-(const uint_t& other) const {
    uint_t result;
    bucket_type remainder = 0;

    for (size_t i = 0; i < c_BUCKET_NUMBER; ++i) {
        result[i] += m_buckets[i] - other[i] + remainder;
        remainder = static_cast<bucket_type>(result[i] > (m_buckets[i] + remainder));
    }

    return result;
}

template<size_t bits>
uint_t<bits> uint_t<bits>::operator*(const uint_t& other) const {
    uint_t result;
    uint32_t remainder = 0;

    for (size_t i = 0; i < c_BUCKET_NUMBER; ++i) {
        result += (*this * other[i]) << (c_BUCKET_SIZE * i);
    }

    return result;
}

template<size_t bits>
uint_t<bits> ECG::uint_t<bits>::operator*(bucket_type other) const {
    uint_t result;
    uint32_t remainder = 0;

    for (size_t i = 0; i < c_BUCKET_NUMBER; ++i) {
        uint64_t prod = m_buckets[i] * static_cast<uint64_t>(other);
        result[i] = static_cast<uint32_t>(prod) + remainder;
        remainder = static_cast<uint32_t>(result[i] < remainder) + static_cast<uint32_t>(prod >> 32);
    }

    return result;
}

template<size_t bits>
uint_t<bits> uint_t<bits>::operator/(const uint_t& other) const {
    return divide(*this, other);
}

template<size_t bits>
uint_t<bits> uint_t<bits>::operator%(const uint_t& other) const {
    uint_t result;
    divide(*this, other, &result);
    return result;
}

template<size_t bits>
uint_t<bits> uint_t<bits>::operator/(bucket_type other) const {
    return divide(*this, other);
}

template<size_t bits>
uint_t<bits> uint_t<bits>::operator%(bucket_type other) const {
    uint_t result;
    divide(*this, other, &result);
    return result;
}

template<size_t bits>
uint_t<bits> uint_t<bits>::operator>>(size_t shift) const {
    size_t bucket_shift = shift / c_BUCKET_SIZE;
    uint_t result;

    for (size_t i = 0; i + bucket_shift < c_BUCKET_NUMBER; ++i) {
        result[i] = m_buckets[i + bucket_shift];
    }

    shift %= c_BUCKET_SIZE;

    if (shift == 0) {
        return result;
    }

    for (size_t i = 0; i + bucket_shift < c_BUCKET_NUMBER; ++i) {
        result[i] >>= shift;

        if (i + 1 < c_BUCKET_NUMBER) {
            result[i] |= m_buckets[i + 1] << (c_BUCKET_SIZE - shift);
        }
    }

    return result;
}

template<size_t bits>
uint_t<bits> uint_t<bits>::operator<<(size_t shift) const {
    size_t bucket_shift = shift / c_BUCKET_SIZE;
    uint_t result;

    for (size_t i = 0; i + bucket_shift < c_BUCKET_NUMBER; ++i) {
        result[i + bucket_shift] = m_buckets[i];
    }

    shift %= c_BUCKET_SIZE;

    if (shift == 0) {
        return result;
    }

    for (size_t i = c_BUCKET_NUMBER; i > bucket_shift; --i) {
        result[i - 1] <<= shift;

        if (i - 1 > 0) {
            result[i - 1] |= m_buckets[i - 2] >> (c_BUCKET_SIZE - shift);
        }
    }

    return result;
}

template<size_t bits>
uint_t<bits> uint_t<bits>::operator^(const uint_t& other) const {
    uint_t result;

    for (size_t i = 0; i < c_BUCKET_NUMBER; ++i) {
        result[i] = m_buckets[i] ^ other[i];
    }

    return result;
}

template<size_t bits>
uint_t<bits> uint_t<bits>::operator|(const uint_t& other) const {
    uint_t result;

    for (size_t i = 0; i < c_BUCKET_NUMBER; ++i) {
        result[i] = m_buckets[i] | other[i];
    }

    return result;
}

template<size_t bits>
uint_t<bits> uint_t<bits>::operator&(const uint_t& other) const {
    uint_t result;

    for (size_t i = 0; i < c_BUCKET_NUMBER; ++i) {
        result[i] = m_buckets[i] & other[i];
    }

    return result;
}

template<size_t bits>
uint_t<bits> uint_t<bits>::operator-() const {
    uint_t result;

    for (size_t i = 0; i < c_BUCKET_NUMBER; ++i) {
        result[i] = ~(m_buckets[i]);
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
uint_t<bits>& ECG::uint_t<bits>::operator*=(bucket_type other) {
    return (*this = *this * other);
}

template<size_t bits>
inline uint_t<bits>& uint_t<bits>::operator/=(const uint_t& other) {
    return (*this = *this / other);
}

template<size_t bits>
inline uint_t<bits>& uint_t<bits>::operator%=(const uint_t& other) {
    return (*this = *this % other);
}

template<size_t bits>
inline uint_t<bits>& uint_t<bits>::operator/=(bucket_type other) {
    return (*this = *this / other);
}

template<size_t bits>
inline uint_t<bits>& uint_t<bits>::operator%=(bucket_type other) {
    return (*this = *this % other);
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
    return m_buckets[0];
}

template<size_t bits>
std::string uint_t<bits>::into_string(StringType str_type) const {
    static const uint_t ZERO = uint_t(0);
    std::string result;
    uint_t clone_of_this = *this;

    switch (str_type) {
    case uint_t::StringType::BINARY :
        do {
            result.push_back((static_cast<bucket_type>(clone_of_this)
                              & static_cast<bucket_type>(uint_t::StringType::BINARY))
                             + '0');
            clone_of_this >>= 1;
        } while (clone_of_this > ZERO);
        break;
    case uint_t::StringType::DECIMAL :
        do {
            result.push_back((static_cast<bucket_type>(clone_of_this)
                              % static_cast<bucket_type>(uint_t::StringType::DECIMAL))
                             + '0');
            clone_of_this /= bucket_type(uint_t::StringType::DECIMAL);
        } while (clone_of_this > ZERO);
        break;
    case uint_t::StringType::HEXADECIMAL :
        do {
            bucket_type value = static_cast<bucket_type>(clone_of_this)
                              & static_cast<bucket_type>(uint_t::StringType::HEXADECIMAL);

            if (value >= 10) {
                value -= 10;
                result.push_back(value + 'a');
            } else {
                result.push_back(value + '0');
            }

            clone_of_this >>= 4;
        } while (clone_of_this > ZERO);
        break;
    }

    std::reverse(result.begin(), result.end());
    return result;
}

template<size_t bits>
std::string uint_t<bits>::into_string(std::function<char(bucket_type)> map, size_t shift) const {
    static constexpr uint_t ZERO = uint_t(0);
    std::string result;
    uint_t clone = *this;

    do {
        result.push_back(map((bucket_type(clone) << (c_BUCKET_SIZE - shift)) >> (c_BUCKET_SIZE - shift)));
        clone >>= shift;
    } while (clone > ZERO);

    std::reverse(result.begin(), result.end());
    return result;
}

template<size_t bits>
inline uint_t<bits> uint_t<bits>::divide(const uint_t& lhs, const uint_t& rhs, uint_t* remainder) {
    static constexpr uint_t ZERO = uint_t(0);

    assert(rhs != ZERO && "Division by 0");

    size_t dividend_size = lhs.clz();
    size_t divisor_size = rhs.clz();

    // CASE 0:
    if (dividend_size < divisor_size) {
        if (remainder != nullptr) {
            *remainder = lhs;
        }

        return ZERO;
    }

    // CASE 1:
    if (divisor_size == 1) {
        return divide(lhs, rhs[0], remainder);
    }

    // CASE 2:
    return d_divide(lhs, rhs, remainder);
}

template<size_t bits>
inline uint_t<bits> uint_t<bits>::divide(const uint_t& lhs, const bucket_type& rhs, uint_t* remainder) {
    assert(rhs != 0 && "Division by 0");

    uint_t result;
    uint64_t part = 0;

    for (size_t i = c_BUCKET_NUMBER; i > 0; --i) {
        part = (part << (c_BUCKET_SIZE)) + static_cast<uint64_t>(lhs[i - 1]);

        if (part < rhs) {
            continue;
        }

        result[i - 1] = static_cast<uint32_t>(part / rhs);
        part %= rhs;
    }

    if (remainder != nullptr) {
        *remainder = uint_t(static_cast<bucket_type>(part));
    }

    return result;
}

// Algorithm D implementation, Art of Computer Programming, Knuth, Vol. 2
template<size_t bits>
inline uint_t<bits> uint_t<bits>::d_divide(const uint_t& lhs, const uint_t& rhs, uint_t* remainder) {
    size_t dividend_size = lhs.clz();
    size_t divisor_size = rhs.clz();

    uint_t<bits + c_BUCKET_SIZE> dividend(lhs);
    uint_t divisor(rhs);
    uint_t quotient;

    size_t shift = 0;
    uint32_t divisor_head = divisor[divisor_size - 1];
    static constexpr uint64_t c_HALF_BUCKET = static_cast<uint64_t>(1) << (c_BUCKET_SIZE - 1);

    while (divisor_head < c_HALF_BUCKET) {
        ++shift;
        divisor_head <<= 1;
    }

    dividend <<= shift;
    divisor <<= shift;

    uint64_t divisor_ = divisor_head;
    static constexpr uint64_t c_BUCKET = static_cast<uint64_t>(1) << c_BUCKET_SIZE;

    for (size_t i = dividend_size - divisor_size + 1; i > 0; --i) {
        uint64_t part = (static_cast<uint64_t>(dividend[i + divisor_size - 1]) << c_BUCKET_SIZE)
                      + static_cast<uint64_t>(dividend[i + divisor_size - 2]);
        uint64_t quotient_t = part / divisor_;
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
            uint64_t product = static_cast<uint32_t>(quotient_t) * static_cast<uint64_t>(divisor[j]);
            widedigit = (static_cast<int64_t>(dividend[i + j - 1]) + carry) - (product & 0xffffffffLL);
            dividend[i + j - 1] = static_cast<uint32_t>(widedigit);
            carry = (widedigit >> c_BUCKET_SIZE) - static_cast<int64_t>(product >> c_BUCKET_SIZE);
        }

        widedigit = static_cast<int64_t>(dividend[i + divisor_size - 1]) + carry;
        dividend[i + divisor_size - 1] = static_cast<uint32_t>(widedigit);

        quotient[i - 1] = static_cast<uint32_t>(quotient_t);

        if (widedigit < 0) {
            --quotient[i - 1];
            widedigit = 0;

            for (size_t j = 0; j < divisor_size; ++j) {
                widedigit += static_cast<uint64_t>(dividend[i + j - 1]) + divisor[j];
                dividend[i + j - 1] = static_cast<uint32_t>(widedigit);
                widedigit >>= 32;
            }
        }
    }

    if (remainder != nullptr) {
        *remainder = uint_t(0);

        for (size_t i = 0; i < divisor_size - 1; ++i) {
            (*remainder)[i] =
                (dividend[i] >> shift) | (static_cast<uint64_t>(dividend[i + 1]) << (c_BUCKET_SIZE - shift));
        }

        (*remainder)[divisor_size - 1] = dividend[divisor_size - 1] >> shift;
    }

    return quotient;
}

template<size_t bits>
constexpr inline size_t uint_t<bits>::size() const {
    return c_BUCKET_NUMBER;
}

template<size_t bits>
inline size_t uint_t<bits>::clz() const {
    size_t result = c_BUCKET_NUMBER;

    while (result > 0 && m_buckets[result - 1] == 0) {
        --result;
    }

    return result;
}

#endif
