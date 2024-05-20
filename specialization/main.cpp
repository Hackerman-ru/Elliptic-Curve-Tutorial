#include <array>
#include <cassert>
#include <iostream>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <vector>

template<typename From, typename To>
concept is_convertible_to = requires(From f) {
    { static_cast<To>(f) } noexcept;
};

template<typename From, typename To>
concept is_upcastable_to = sizeof(From) <= sizeof(To) && is_convertible_to<From, To>;

template<typename From, typename To>
concept is_downcastable_to = sizeof(From) > sizeof(To) && is_convertible_to<From, To>;

class F_256 {
    static constexpr size_t c_bits = 512;
    static constexpr size_t c_bits_in_byte = 8;
    static constexpr size_t c_digit_size = 32;
    static constexpr size_t c_digit_number = 16;
    static constexpr size_t c_half_digit_number = 8;

    using digits = std::array<uint32_t, c_digit_number>;
    digits m_digits = {};

public:
    static constexpr F_256 inverse(const F_256& element) {
        F_256 result = element;
        result.inverse();
        return result;
    }

    static constexpr F_256 inverse(F_256&& element) {
        element.inverse();
        return element;
    }

    static constexpr F_256 pow(const F_256& element, const F_256& power) {
        if (!power.is_even()) {
            if (power == 1) {
                return element;
            }

            return element * pow(element, power - 1);
        }

        F_256 temp = pow(element, power >> 1);
        return temp * temp;
    }

    constexpr F_256() = default;

    template<typename T>
    constexpr F_256(const T& value) : m_digits(split_into_digits<T>(value)) {}

    constexpr F_256(digits&& value) : m_digits(std::move(value)) {}

    constexpr F_256(const digits& value) : m_digits(value) {}

    constexpr F_256(const char* str) : m_digits(parse_into_uint(str).m_digits) {}

    constexpr F_256& operator=(const F_256& value) = default;

    template<typename T>
    constexpr F_256& operator=(const T& value) {
        m_digits = split_into_digits<T>(value);
        return *this;
    }

    constexpr F_256& operator=(const char* str) {
        return *this = parse_into_uint(str);
    }

    friend constexpr std::strong_ordering operator<=>(const F_256& lhs, const F_256& rhs) {
        for (size_t i = c_digit_number; i > 0; --i) {
            if (lhs[i - 1] != rhs[i - 1]) {
                return lhs[i - 1] <=> rhs[i - 1];
            }
        }

        return std::strong_ordering::equal;
    }

    friend constexpr bool operator==(const F_256& lhs, const F_256& rhs) {
        return lhs.m_digits == rhs.m_digits;
    }

    // operator+
    friend constexpr F_256 operator+(const F_256& lhs, const F_256& rhs) {
        F_256 result = lhs;
        result += rhs;
        return result;
    }

    friend constexpr F_256 operator+(F_256&& lhs, const F_256& rhs) {
        lhs += rhs;
        return lhs;
    }

    friend constexpr F_256 operator+(const F_256& lhs, F_256&& rhs) {
        rhs += lhs;
        return rhs;
    }

    friend constexpr F_256 operator+(F_256&& lhs, F_256&& rhs) {
        lhs += rhs;
        return lhs;
    }

    // operator-
    friend constexpr F_256 operator-(const F_256& lhs, const F_256& rhs) {
        F_256 result = lhs;
        result -= rhs;
        return result;
    }

    friend constexpr F_256 operator-(F_256&& lhs, const F_256& rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend constexpr F_256 operator-(const F_256& lhs, F_256&& rhs) {
        F_256 result = lhs;
        result -= rhs;
        return result;
    }

    friend constexpr F_256 operator-(F_256&& lhs, F_256&& rhs) {
        lhs -= rhs;
        return lhs;
    }

    // operator*
    friend constexpr F_256 operator*(const F_256& lhs, const F_256& rhs) {
        F_256 result;

        for (size_t i = 0; i < c_half_digit_number; ++i) {
            uint64_t u = 0;

            for (size_t j = 0; j < c_half_digit_number; ++j) {
                u = static_cast<uint64_t>(result[i + j])
                  + static_cast<uint64_t>(lhs[i]) * static_cast<uint64_t>(rhs[j]) + (u >> c_digit_size);
                result[i + j] = static_cast<uint32_t>(u);
            }

            result[i + c_half_digit_number] = static_cast<uint32_t>(u >> c_digit_size);
        }

        result.reduce();
        return result;
    }

    // operator/
    friend constexpr F_256 operator/(const F_256& lhs, const F_256& rhs) {
        return lhs * inverse(rhs);
    }

    // operator>>
    friend constexpr F_256 operator>>(const F_256& lhs, const size_t& rhs) {
        F_256 result = lhs;
        return result >>= rhs;
    }

    friend constexpr F_256 operator>>(F_256&& lhs, const size_t& rhs) {
        return lhs >>= rhs;
    }

    // operator<<
    friend constexpr F_256 operator<<(const F_256& lhs, const size_t& rhs) {
        F_256 result = lhs;
        return result <<= rhs;
    }

    friend constexpr F_256 operator<<(F_256&& lhs, const size_t& rhs) {
        return lhs <<= rhs;
    }

    constexpr F_256 operator-() const {
        if(!is_invertible) {
            return *this;
        }
    
        blocks result = p_values;
        uint32_t remainder = 0;
    
        for (size_t i = 0; i < c_digit_number; ++i) {
            uint32_t prev = result[i];
            uint32_t sum = m_digits[i] + remainder;
            result[i] -= sum;
            remainder = (result[i] > prev) || (sum < remainder);
        }
    
        return F_256(result);
    }

    constexpr F_256& operator+=(const F_256& other) {
        uint32_t carry = 0;

        for (size_t i = 0; i < c_digit_number; ++i) {
            uint32_t sum = carry + other[i];
            m_digits[i] += sum;
            carry = (m_digits[i] < sum) || (sum < carry);
        }

        if (!is_valid()) {
            subtract_p_uncheck();
        }
        assert(is_valid());
        return *this;
    }

    constexpr F_256& operator-=(const F_256& other) {
        uint32_t remainder = 0;

        for (size_t i = 0; i < c_digit_number; ++i) {
            uint32_t prev = m_digits[i];
            uint32_t sum = other[i] + remainder;
            m_digits[i] -= sum;
            remainder = (m_digits[i] > prev) || (sum < remainder);
        }

        if (!is_valid()) {
            add_p_uncheck();
        }
        assert(is_valid());
        return *this;
    }

    constexpr F_256& operator*=(const F_256& other) {
        return *this = *this * other;
    }

    constexpr F_256& operator/=(const F_256& other) {
        return *this = *this / other;
    }

    constexpr F_256& operator>>=(size_t shift_size) {
        size_t digit_shift = shift_size >> 6;

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

        return *this;
    }

    constexpr F_256& operator<<=(size_t shift_size) {
        size_t digit_shift = shift_size >> 6;

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

        for (size_t i = c_digit_number; i > digit_shift; --i) {
            m_digits[i - 1] <<= shift_size;

            if (i - 1 > 0) {
                m_digits[i - 1] |= m_digits[i - 2] >> (c_digit_size - shift_size);
            }
        }

        while (!is_valid()) {
            subtract_p_uncheck();
        }

        return *this;
    }

    [[nodiscard("Optimize unary operator usage")]]
    constexpr F_256
        operator++(int) {
        F_256 result = *this;
        increment();
        assert(is_valid());
        return result;
    }

    constexpr F_256& operator++() {
        increment();
        assert(is_valid());
        return *this;
    }

    [[nodiscard("Optimize unary operator usage")]]
    constexpr F_256
        operator--(int) {
        F_256 result = *this;
        decrement();
        assert(is_valid());
        return result;
    }

    constexpr F_256& operator--() {
        decrement();
        assert(is_valid());
        return *this;
    }

    constexpr std::string convert_to_hex_string() const {
        std::string result;
        size_t pos = c_digit_number;

        while (pos > 0 && m_digits[pos - 1] == 0) {
            --pos;
        }

        if (pos == 0) {
            return "0x0";
        }

        result += "0x";

        auto mini_push = [&](const uint8_t& value) {
            if (value < 10) {
                result.push_back(value + '0');
            } else {
                result.push_back(value - 10 + 'A');
            }
        };

        auto push = [&](const uint32_t& value, bool first_time = false) {
            size_t shift = 32;
            uint8_t temp = 0;

            if (first_time) {
                while (shift > 0 && temp == 0) {
                    shift -= 4;
                    temp = (value >> shift) & 0xF;
                }

                do {
                    mini_push(temp);

                    if (shift == 0) {
                        break;
                    }

                    shift -= 4;
                    temp = (value >> shift) & 0xF;
                } while (shift >= 0);

                return;
            }

            while (shift > 0) {
                shift -= 4;
                temp = (value >> shift) & 0xF;
                mini_push(temp);
            }
        };

        push(m_digits[pos - 1], true);
        --pos;

        while (pos > 0) {
            push(m_digits[pos - 1]);
            --pos;
        }

        return result;
    }

    constexpr bool is_invertible() const {
        constexpr F_256 zero;
        return *this != zero;
    }

    constexpr void inverse() {
        constexpr F_256 zero;
        constexpr F_256 one = 1;

        F_256 u = *this;
        F_256 v(p_values);

        F_256 x_1 = one;
        F_256 x_2;

        while (u != one && v != one) {
            while (u.is_even()) {
                u >>= 1;

                if (x_1.is_even()) {
                    x_1 >>= 1;
                } else {
                    x_1.add_p_uncheck();
                    x_1 >>= 1;
                    assert(x_1.is_valid());
                }
            }

            while (v.is_even()) {
                v >>= 1;

                if (x_2.is_even()) {
                    x_2 >>= 1;
                } else {
                    x_2.add_p_uncheck();
                    x_2 >>= 1;
                    assert(x_2.is_valid());
                }
            }

            if (u >= v) {
                u -= v;
                x_1 -= x_2;
            } else {
                v -= u;
                x_2 -= x_1;
            }
        }

        if (u == 1) {
            *this = x_1;
        } else {
            *this = x_2;
        }

        assert(is_valid());
    }

    constexpr bool is_even() const {
        return (m_digits[0] & 0b1) == 0;
    }

    const uint32_t& first_digit() const {
        return m_digits[0];
    }

    constexpr const uint32_t& operator[](size_t pos) const {
        return m_digits[pos];
    }

    constexpr uint32_t& operator[](size_t pos) {
        return m_digits[pos];
    }

private:
    static constexpr F_256 parse_into_uint(const char* str) {
        F_256 value = 0;
        uint16_t radix = 10;

        if (str[0] == '0' && str[1] == 'x') {
            radix = 16;
            str += 2;
        } else if (str[0] == '0' && str[1] == 'b') {
            radix = 2;
            str += 2;
        } else if (str[0] == '0') {
            radix = 8;
            ++str;
        }

        while (*str != '\0') {
            switch (radix) {
            case 16 :
                value <<= 4;
                break;
            case 8 :
                value <<= 3;
                break;
            case 2 :
                value <<= 1;
                break;
            default :
                value *= radix;
                break;
            }

            uint32_t symbol_value = 0;

            if (*str >= '0' && *str <= '9') {
                symbol_value = static_cast<uint32_t>(*str - '0');
            } else if (*str >= 'a' && *str <= 'f') {
                symbol_value = static_cast<uint32_t>(*str - 'a') + 10;
            } else if (*str >= 'A' && *str <= 'F') {
                symbol_value = static_cast<uint32_t>(*str - 'A') + 10;
            }

            value += symbol_value;
            ++str;
        }

        return value;
    }

    template<typename T>
    requires std::numeric_limits<T>::is_integer && is_upcastable_to<T, uint32_t>
    static constexpr digits split_into_digits(T value) {
        return {static_cast<uint32_t>(value)};
    }

    template<typename T>
    requires std::numeric_limits<T>::is_integer && is_downcastable_to<T, uint32_t>
    static constexpr digits split_into_digits(T value) {
        digits result = {};

        for (size_t i = 0; i < c_digit_number; ++i) {
            result[i] = static_cast<uint32_t>(value);
            value >>= c_digit_size;

            if (value == 0) {
                break;
            }
        }

        return result;
    }

    constexpr void reduce() {
        F_256 s1({m_digits[0],
                  m_digits[1],
                  m_digits[2],
                  m_digits[3],
                  m_digits[4],
                  m_digits[5],
                  m_digits[6],
                  m_digits[7]});
        F_256 s2({0, 0, 0, m_digits[11], m_digits[12], m_digits[13], m_digits[14], m_digits[15]});
        F_256 s3({0, 0, 0, m_digits[12], m_digits[13], m_digits[14], m_digits[15], 0});
        F_256 s4({m_digits[8], m_digits[9], m_digits[10], 0, 0, 0, m_digits[14], m_digits[15]});
        F_256 s5({m_digits[9],
                  m_digits[10],
                  m_digits[11],
                  m_digits[13],
                  m_digits[14],
                  m_digits[15],
                  m_digits[13],
                  m_digits[8]});
        F_256 s6({m_digits[11], m_digits[12], m_digits[13], 0, 0, 0, m_digits[8], m_digits[10]});
        F_256 s7({m_digits[12], m_digits[13], m_digits[14], m_digits[15], 0, 0, m_digits[9], m_digits[11]});
        F_256 s8({m_digits[13],
                  m_digits[14],
                  m_digits[15],
                  m_digits[8],
                  m_digits[9],
                  m_digits[10],
                  0,
                  m_digits[12]});
        F_256 s9({m_digits[14], m_digits[15], 0, m_digits[9], m_digits[10], m_digits[11], 0, m_digits[13]});
        *this = s1 + s2 + s2 + s3 + s3 + s4 + s5 - s6 - s7 - s8 - s9;
        assert(is_valid());
    }

    constexpr void increment() {
        for (size_t i = 0; i < c_digit_number; ++i) {
            m_digits[i] += 1;

            if (m_digits[i] != 0) {
                break;
            }
        }

        if (m_digits == p_values) {
            m_digits = {};
        }

        assert(is_valid());
    }

    static constexpr digits max_mod_p = {4294967294U, 4294967295U, 4294967295U, 0U, 0U, 0U, 1U, 4294967295U};

    constexpr void decrement() {
        if (*this == 0) {
            m_digits = max_mod_p;
            assert(is_valid());
            return;
        }

        for (size_t i = 0; i < c_digit_number; ++i) {
            uint32_t temp = m_digits[i];
            m_digits[i] -= 1;

            if (temp >= m_digits[i]) {
                break;
            }
        }

        assert(is_valid());
    }

    static constexpr digits p_values = {4294967295U, 4294967295U, 4294967295U, 0U, 0U, 0U, 1U, 4294967295U};

    constexpr void add_p_uncheck() {
        uint32_t carry = 0;

        for (size_t i = 0; i < c_digit_number; ++i) {
            uint32_t sum = carry + p_values[i];
            m_digits[i] += sum;
            carry = (m_digits[i] < sum) || (sum < carry);
        }
    }

    constexpr void subtract_p_uncheck() {
        uint32_t remainder = 0;

        for (size_t i = 0; i < c_digit_number; ++i) {
            uint32_t prev = m_digits[i];
            uint32_t sum = p_values[i] + remainder;
            m_digits[i] -= sum;
            remainder = (m_digits[i] > prev) || (sum < remainder);
        }
    }

    constexpr bool is_valid() const {
        constexpr F_256 p(p_values);
        return *this < p;
    }
};

constexpr size_t c_width = 3;

struct Coefficient {
    uint16_t value;
    bool is_negative;
};

static constexpr uint16_t c_mask_modulo_2_pow_w = (1 << c_width) - 1;
static constexpr size_t c_kp_number = static_cast<size_t>(1) << (c_width - 2);

using wnaf_form = std::vector<Coefficient>;

static constexpr wnaf_form get_wnaf(F_256 value) {
    wnaf_form result;

    while (value > 0) {
        if (!value.is_even()) {
            uint16_t coef_value = static_cast<uint16_t>(value.first_digit()) & c_mask_modulo_2_pow_w;

            if (coef_value >= (1 << (c_width - 1))) {
                coef_value = (1 << c_width) - coef_value;
                result.push_back({.value = coef_value, .is_negative = true});
                value += coef_value;
            } else {
                result.push_back({.value = coef_value, .is_negative = false});
                value -= coef_value;
            }
        } else {
            result.push_back({.value = 0, .is_negative = false});
        }

        value >>= 1;
    }

    return result;
}

class EllipticCurvePoint {
public:
    constexpr EllipticCurvePoint(const F_256& x, const F_256& y, bool is_null = false) :
        m_x {x}, m_y {y}, m_is_null(is_null) {}

    friend constexpr EllipticCurvePoint operator+(const EllipticCurvePoint& lhs,
                                                  const EllipticCurvePoint& rhs) {
        EllipticCurvePoint result = lhs;
        result += rhs;
        return result;
    }

    friend constexpr EllipticCurvePoint operator+(EllipticCurvePoint&& lhs, const EllipticCurvePoint& rhs) {
        lhs += rhs;
        return lhs;
    }

    friend constexpr EllipticCurvePoint operator+(const EllipticCurvePoint& lhs, EllipticCurvePoint&& rhs) {
        rhs += lhs;
        return rhs;
    }

    friend constexpr EllipticCurvePoint operator+(EllipticCurvePoint&& lhs, EllipticCurvePoint&& rhs) {
        lhs += rhs;
        return lhs;
    }

    friend constexpr EllipticCurvePoint operator-(const EllipticCurvePoint& lhs,
                                                  const EllipticCurvePoint& rhs) {
        EllipticCurvePoint result = lhs;
        result -= rhs;
        return result;
    }

    friend constexpr EllipticCurvePoint operator-(EllipticCurvePoint&& lhs, const EllipticCurvePoint& rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend constexpr EllipticCurvePoint operator-(const EllipticCurvePoint& lhs, EllipticCurvePoint&& rhs) {
        rhs -= lhs;
        rhs.negative();
        return rhs;
    }

    friend constexpr EllipticCurvePoint operator-(EllipticCurvePoint&& lhs, EllipticCurvePoint&& rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend constexpr EllipticCurvePoint operator*(const EllipticCurvePoint& point, const F_256& value) {
        EllipticCurvePoint result = point;
        result *= value;
        return result;
    }

    friend constexpr EllipticCurvePoint operator*(EllipticCurvePoint&& point, const F_256& value) {
        point *= value;
        return point;
    }

    friend constexpr EllipticCurvePoint operator*(const F_256& value, const EllipticCurvePoint& point) {
        EllipticCurvePoint result = point;
        result *= value;
        return result;
    }

    friend constexpr EllipticCurvePoint operator*(const F_256& value, EllipticCurvePoint&& point) {
        point *= value;
        return point;
    }

    friend constexpr bool operator==(const EllipticCurvePoint& lhs, const EllipticCurvePoint& rhs) {
        return (lhs.m_is_null && rhs.m_is_null) || (lhs.m_x == rhs.m_x && lhs.m_y == rhs.m_y);
    }

    constexpr EllipticCurvePoint operator-() const {
        EllipticCurvePoint result = *this;
        result.negative();
        return result;
    }

    constexpr EllipticCurvePoint& operator+=(const EllipticCurvePoint& other) {
        if (m_is_null) {
            return *this = other;
        } else if (other.m_is_null) {
            return *this;
        }

        if (m_x == other.m_x) {
            if (m_y != other.m_y) {
                m_is_null = true;
            } else {
                twice();
            }

            return *this;
        }

        const F_256 k = (other.m_y - m_y) / (other.m_x - m_x);
        const F_256 x = F_256::pow(k, 2) - m_x - other.m_x;
        m_y = k * (m_x - x) - m_y;
        m_x = x;

        return *this;
    }

    constexpr EllipticCurvePoint& operator-=(const EllipticCurvePoint& other) {
        EllipticCurvePoint temp = other;
        temp.negative();
        return *this += temp;
    }

    constexpr EllipticCurvePoint& operator-=(EllipticCurvePoint&& other) {
        other.negative();
        return *this += other;
    }

    constexpr EllipticCurvePoint& operator*=(const F_256& value) {
        multiply(value);
        return *this;
    }

    constexpr F_256 get_x() const {
        return m_x;
    }

    constexpr F_256 get_y() const {
        return m_y;
    }

    constexpr bool is_zero() const {
        return m_is_null;
    }

    static constexpr EllipticCurvePoint null_point() {
        return EllipticCurvePoint(0, 1, true);
    }

    static std::optional<EllipticCurvePoint> point_with_x_equal_to(const F_256& x) {
        if (!x.is_invertible()) {
            return null_point();
        }

        std::optional<F_256> y = find_y(x);

        if (!y.has_value()) {
            return std::nullopt;
        }

        return EllipticCurvePoint(x, std::move(y.value()));
    }

private:
    static constexpr std::optional<F_256> find_y(const F_256& x) {
        F_256 value = F_256::pow(x, 3) + c_a * x + c_b;

        if (!value.is_invertible()) {
            return std::nullopt;
        }

        if (F_256::pow(value, c_p_1) != F_256(1)) {
            return std::nullopt;
        }

        return F_256::pow(value, c_p_2);
    }

    constexpr void multiply(const F_256& value) {
        auto wnaf_form = get_wnaf(value);
        EllipticCurvePoint two_p = *this + *this;
        std::vector<EllipticCurvePoint> kp = {*this};

        for (size_t i = 1; i < c_kp_number; ++i) {
            kp.emplace_back(kp.back() + two_p);
        }

        m_is_null = true;

        for (size_t i = wnaf_form.size(); i > 0; --i) {
            twice();

            if (wnaf_form[i - 1].value != 0) {
                if (!wnaf_form[i - 1].is_negative) {
                    *this += kp[wnaf_form[i - 1].value >> 1];
                } else {
                    *this -= kp[wnaf_form[i - 1].value >> 1];
                }
            }
        }
    }

    constexpr void negative() {
        m_y = -m_y;
    }

    constexpr void twice() {
        if (m_is_null) {
            return;
        }

        if (!m_y.is_invertible()) {
            m_is_null = true;
            return;
        }

        const F_256 k = (3 * F_256::pow(m_x, 2) + c_a) / (m_y << 1);
        const F_256 x = F_256::pow(k, 2) - (m_x << 1);
        m_y = k * (m_x - x) - m_y;
        m_x = x;
    }

    static constexpr F_256 c_a = "0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc";
    static constexpr F_256 c_b = "0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b";
    static constexpr F_256 c_p_1 = "0x7FFFFFFF800000008000000000000000000000007FFFFFFFFFFFFFFFFFFFFFFF";
    static constexpr F_256 c_p_2 = "0x3FFFFFFFC0000000400000000000000000000000400000000000000000000000";

    F_256 m_x;
    F_256 m_y;
    bool m_is_null;
};

static std::mt19937 gen32(42);

constexpr F_256 generate_random_non_zero_value_modulo(const F_256& modulus) {
    F_256 result;

    for (size_t i = 0; i < 8; ++i) {
        result[i] = gen32();
    }

    if (result >= modulus) {
        result = modulus - 1;
    }

    return result;
}

namespace ElGamal {
    using Point = EllipticCurvePoint;

    namespace {
        static constexpr Point m_generator =
            Point(F_256("0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296"),
                  F_256("0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5"));
        static constexpr F_256 m_generator_order =
            "0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551";

        Point map_to_curve(const F_256& message) {
            for (;;) {
                F_256 x = generate_random_non_zero_value_modulo(m_generator_order);

                for (size_t i = 0; i < 6; ++i) {
                    x[i] = message[i];
                }

                auto opt = Point::point_with_x_equal_to(x);

                if (opt.has_value()) {
                    return opt.value();
                }
            }
        }

        F_256 map_to_uint(const Point& message) {
            F_256 result;
            const F_256& x = message.get_x();

            for (size_t i = 0; i < 6; ++i) {
                result[i] = x[i];
            }

            return result;
        }
    }   // namespace

    struct Keys {
        F_256 private_key;
        Point public_key;
    };

    struct EncryptedMessage {
        Point generator_degree;
        Point message_with_salt;
    };

    constexpr Keys generate_keys() {
        F_256 private_key = generate_random_non_zero_value_modulo(m_generator_order);
        Point public_key = private_key * m_generator;
        return Keys {.private_key = private_key, .public_key = public_key};
    }

    EncryptedMessage encrypt(const F_256& message, const Point& public_key) {
        const F_256 k = generate_random_non_zero_value_modulo(m_generator_order);
        const Point generator_degree = k * m_generator;
        const Point message_with_salt = map_to_curve(message) + k * public_key;
        return {.generator_degree = generator_degree, .message_with_salt = message_with_salt};
    }

    F_256 decrypt(const EncryptedMessage& encrypted_message, const F_256& private_key) {
        Point M = encrypted_message.message_with_salt - private_key * encrypted_message.generator_degree;
        return map_to_uint(M);
    }
};   // namespace ElGamal

int main() {
    std::cout << "Enter hexadecimal message:\n";
    std::string msg;
    std::cin >> msg;
    F_256 message = msg.c_str();
    std::cout << "Generating keys...\n";
    auto keys = ElGamal::generate_keys();
    std::cout << "Encrypting message...\n";
    auto encrypted_message = ElGamal::encrypt(message, keys.public_key);
    std::cout << "Decrypting message...\n";
    auto decrypted_message = ElGamal::decrypt(encrypted_message, keys.private_key);
    std::cout << "Decrypted message is:\n";
    std::string decrypted_msg = decrypted_message.convert_to_hex_string();
    std::cout << decrypted_msg << '\n';
    std::cout << "Checking correctness...\n";

    if (message != decrypted_message) {
        std::cout << "Fail ";

        if (msg.ends_with(decrypted_msg.substr(2))) {
            std::cout << "due to insufficient number of encryption bits\n";
        } else {
            std::cout << "due to implementation problems\n";
        }
    } else {
        std::cout << "Success!\n";
    }
}
