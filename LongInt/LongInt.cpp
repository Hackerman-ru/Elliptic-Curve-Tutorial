#include "LongInt.h"

#include <cassert>
#include <complex>

using ECG::uint512_t;

constexpr uint512_t::uint512_t(T value) : m_value({value}) {};

uint512_t::uint512_t(const std::string& str, StringType str_type) {
    *this = uint512_t {};

    switch (str_type) {
    case uint512_t::StringType::BINARY :
        for (size_t i = 0; i < str.size() && i < c_BITS; ++i) {
            size_t bucket_pos = i / c_BUCKET_SIZE;

            m_value[bucket_pos] |= T(str[str.size() - 1 - i] - '0') << (i % c_BUCKET_SIZE);
        }
        break;
    case uint512_t::StringType::DECIMAL :
        for (char c : str) {
            uint512_t temp = (*this <<= 1);
            *this <<= 2;
            *this += temp;
            *this += uint512_t(c - '0');
        }
        break;
    case uint512_t::StringType::HEXADECIMAL :
        for (char c : str) {
            *this <<= 4;
            *this += uint512_t(c - (isdigit(c) ? '0' : 'a'));
        }
        break;
    }
}

uint512_t::uint512_t(const std::string& str, std::function<T(char)> map, size_t shift) {
    *this = uint512_t {};

    for (char c : str) {
        if (isdigit(c)) {
            *this <<= shift;
            *this += uint512_t(map(c));
        }
    }
}

uint512_t uint512_t::operator+(const uint512_t& other) const {
    uint512_t result {};

    for (size_t i = 0; i < m_value.size(); ++i) {
        result.m_value[i] += m_value[i] + other.m_value[i];

        if (i + 1 < m_value.size() && result.m_value[i] < m_value[i]) {
            ++result.m_value[i + 1];
        }
    }

    return result;
}

uint512_t uint512_t::operator-(const uint512_t& other) const {
    return *this + (-other);
}

namespace {
    using cd = std::complex<double>;
    static const double PI = acos(-1);

    void fft(std::vector<cd>& a, bool invert) {
        size_t n = a.size();

        if (n == 1) {
            return;
        }

        std::vector<cd> a0(n / 2), a1(n / 2);

        for (size_t i = 0; 2 * i < n; i++) {
            a0[i] = a[2 * i];
            a1[i] = a[2 * i + 1];
        }

        fft(a0, invert);
        fft(a1, invert);

        double ang = 2 * PI / n * (invert ? -1 : 1);
        cd w(1), wn(cos(ang), sin(ang));

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
        std::vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
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
            result[i] = round(fa[i].real());
        }

        return result;
    }
}   // namespace

uint512_t uint512_t::operator*(const uint512_t& other) const {
    std::vector<T> lhs(m_value.begin(), m_value.end());
    std::vector<T> rhs(other.m_value.begin(), other.m_value.end());
    std::vector<T> product = multiply(lhs, rhs);
    uint512_t result {};

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

uint512_t uint512_t::operator/(T value) const {
    assert(value != 0 && "Division by 0");
    uint512_t result {};

    if constexpr (std::is_same_v<T, uint32_t>) {
        uint64_t remainder = 0;
        uint32_t divisor = value;

        for (size_t i = m_value.size(); i > 0; --i) {
            remainder += m_value[i - 1];

            if (remainder < divisor) {
                continue;
            }

            result.m_value[i - 1] = remainder / divisor;
            remainder %= divisor;
        }
    } else {
        // Unlikely, but should be done to unify the code
    }

    return result;
}

uint512_t uint512_t::operator%(T value) const {
    assert(value != 0 && "Division by 0");
    return *this - (*this / value) * uint512_t(value);
}

uint512_t uint512_t::operator>>(size_t shift) const {
    size_t bucket_shift = shift / c_BUCKET_SIZE;
    uint512_t result {};

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

uint512_t uint512_t::operator<<(size_t shift) const {
    size_t bucket_shift = shift / c_BUCKET_SIZE;
    uint512_t result {};

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

uint512_t uint512_t::operator^(const uint512_t& other) const {
    uint512_t result {};

    for (size_t i = 0; i < m_value.size(); ++i) {
        result.m_value[i] = m_value[i] ^ other.m_value[i];
    }

    return result;
}

uint512_t uint512_t::operator|(const uint512_t& other) const {
    uint512_t result {};

    for (size_t i = 0; i < m_value.size(); ++i) {
        result.m_value[i] = m_value[i] | other.m_value[i];
    }

    return result;
}

uint512_t uint512_t::operator&(const uint512_t& other) const {
    uint512_t result {};

    for (size_t i = 0; i < m_value.size(); ++i) {
        result.m_value[i] = m_value[i] & other.m_value[i];
    }

    return result;
}

uint512_t ECG::uint512_t::operator-() const {
    uint512_t result {};

    for (size_t i = 0; i < m_value.size(); ++i) {
        result.m_value[i] = ~(m_value[i]);
    }

    result++;

    return result;
}

uint512_t& uint512_t::operator+=(const uint512_t& other) {
    return (*this = *this + other);
}

uint512_t& uint512_t::operator-=(const uint512_t& other) {
    return (*this = *this - other);
}

uint512_t& uint512_t::operator*=(const uint512_t& other) {
    return (*this = *this * other);
}

uint512_t& uint512_t::operator/=(T value) {
    return (*this = *this / value);
}

uint512_t& uint512_t::operator%=(T value) {
    return (*this = *this % value);
}

uint512_t& uint512_t::operator>>=(size_t shift) {
    return (*this = *this >> shift);
}

uint512_t& uint512_t::operator<<=(size_t shift) {
    return (*this = *this << shift);
}

uint512_t& uint512_t::operator^=(const uint512_t& other) {
    return (*this = *this ^ other);
}

uint512_t& uint512_t::operator|=(const uint512_t& other) {
    return (*this = *this | other);
}

uint512_t& uint512_t::operator&=(const uint512_t& other) {
    return (*this = *this & other);
}

uint512_t& uint512_t::operator++(int) {
    return (*this += uint512_t(1));
}

uint512_t& uint512_t::operator++() {
    return (*this += uint512_t(1));
}

uint512_t& uint512_t::operator--() {
    return (*this -= uint512_t(1));
}

uint512_t::operator T() const {
    return m_value[0];
}

std::string uint512_t::into_string(StringType str_type) const {
    std::string result;
    uint512_t clone = *this;

    switch (str_type) {
    case uint512_t::StringType::BINARY :
        do {
            result.push_back((T(clone) & T(uint512_t::StringType::BINARY)) + '0');
            clone >>= 1;
        } while (clone > uint512_t(0));
        break;
    case uint512_t::StringType::DECIMAL :
        do {
            result.push_back((T(clone) % T(uint512_t::StringType::DECIMAL)) + '0');
            clone /= T(uint512_t::StringType::DECIMAL);
        } while (clone > uint512_t(0));
        break;
    case uint512_t::StringType::HEXADECIMAL :
        do {
            result.push_back((T(clone) & T(uint512_t::StringType::HEXADECIMAL)) + '0');
            clone >>= 4;
        } while (clone > uint512_t(0));
        break;
    }

    std::reverse(result.begin(), result.end());
    return result;
}

std::string uint512_t::into_string(std::function<char(T)> map, size_t shift) const {
    std::string result;
    uint512_t clone = *this;

    do {
        result.push_back(map((T(clone) << (c_BUCKET_SIZE - shift)) >> (c_BUCKET_SIZE - shift)));
        clone >>= shift;
    } while (clone > uint512_t(0));

    std::reverse(result.begin(), result.end());
    return result;
}
