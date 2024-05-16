#include <array>
#include <cassert>
#include <string>

template<typename From, typename To>
concept is_convertible_to = requires(From f) {
    { static_cast<To>(f) } noexcept;
};

template<typename From, typename To>
concept is_upcastable_to = sizeof(From) <= sizeof(To) && is_convertible_to<From, To>;

template<typename From, typename To>
concept is_downcastable_to = sizeof(From) > sizeof(To) && is_convertible_to<From, To>;

template<typename Container, typename T>
concept is_convertible_container = requires(Container t, size_t i) {
    { t[i] } -> is_convertible_to<T>;
    { t.size() } -> std::same_as<size_t>;
};

class F_256 {
    using block_t = uint32_t;
    using double_block_t = uint64_t;

    static constexpr size_t c_bits = 512;
    static constexpr size_t c_bits_in_byte = 8;
    static constexpr size_t c_block_size = sizeof(block_t) * c_bits_in_byte;
    static constexpr size_t c_block_number = c_bits / c_block_size;
    static constexpr size_t c_double_block_size = sizeof(double_block_t) * c_bits_in_byte;
    static constexpr size_t c_double_block_number = c_bits / c_double_block_size;

    using blocks = std::array<block_t, c_block_number>;
    blocks m_blocks = {};

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

    static constexpr F_256 pow(const F_256& element, const uint32_t& power);

    constexpr F_256() = default;

    template<typename T>
    constexpr F_256(const T& value) : m_blocks(split_into_blocks<T>(value)) {}

    constexpr F_256(blocks&& value) : m_blocks(std::move(value)) {}

    constexpr F_256(const blocks& value) : m_blocks(value) {}

    constexpr F_256(const char* str) : m_blocks(parse_into_uint(str).m_blocks) {};

    constexpr F_256& operator=(const F_256& value) = default;

    template<typename T>
    constexpr F_256& operator=(const T& value) {
        m_blocks = split_into_blocks<T>(value);
        return *this;
    }

    constexpr F_256& operator=(const char* str) {
        return *this = parse_into_uint(str);
    }

    friend constexpr std::strong_ordering operator<=>(const F_256& lhs, const F_256& rhs) {
        for (size_t i = c_block_number; i > 0; --i) {
            if (lhs[i - 1] != rhs[i - 1]) {
                return lhs[i - 1] <=> rhs[i - 1];
            }
        }

        return std::strong_ordering::equal;
    }

    friend constexpr bool operator==(const F_256& lhs, const F_256& rhs) {
        return lhs.m_blocks == rhs.m_blocks;
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
    friend constexpr F_256 operator*(const F_256& lhs, const F_256& rhs) {}   // TODO

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
        return F_256(p_values) - *this;
    }

    constexpr F_256& operator+=(const F_256& other) {
        block_t carry = 0;

        for (size_t i = 0; i < c_block_number; ++i) {
            block_t sum = carry + other[i];
            m_blocks[i] += sum;
            carry = (m_blocks[i] < sum) || (sum < carry);
        }

        if (!is_valid()) {
            subtract_p();
        }
        assert(is_valid());
        return *this;
    }

    constexpr F_256& operator-=(const F_256& other) {
        block_t remainder = 0;

        for (size_t i = 0; i < c_block_number; ++i) {
            block_t prev = m_blocks[i];
            block_t sum = other[i] + remainder;
            m_blocks[i] -= sum;
            remainder = (m_blocks[i] > prev) || (sum < remainder);
        }

        if (!is_valid()) {
            add_p();
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
        size_t block_shift = shift_size >> 6;
        auto data = reinterpret_cast<double_block_t*>(m_blocks.data());

        if (block_shift > 0) {
            for (size_t i = 0; i < c_double_block_number; ++i) {
                if (i + block_shift < c_double_block_number) {
                    data[i] = data[i + block_shift];
                } else {
                    data[i] = 0;
                }
            }
        }

        shift_size %= c_double_block_size;

        if (shift_size == 0) {
            return *this;
        }

        for (size_t i = 0; i + block_shift < c_double_block_number; ++i) {
            data[i] >>= shift_size;

            if (i + 1 < c_double_block_number) {
                data[i] |= data[i + 1] << (c_double_block_size - shift_size);
            }
        }

        return *this;
    }

    constexpr F_256& operator<<=(size_t shift_size) {
        size_t block_shift = shift_size >> 6;
        auto data = reinterpret_cast<double_block_t*>(m_blocks.data());

        if (block_shift > 0) {
            for (size_t i = c_double_block_number; i > 0; --i) {
                if (i > block_shift) {
                    data[i - 1] = data[i - block_shift - 1];
                } else {
                    data[i - 1] = 0;
                }
            }
        }

        shift_size %= c_double_block_size;

        if (shift_size == 0) {
            return *this;
        }

        for (size_t i = c_double_block_number; i > block_shift; --i) {
            data[i - 1] <<= shift_size;

            if (i - 1 > 0) {
                data[i - 1] |= data[i - 2] >> (c_double_block_size - shift_size);
            }
        }

        while (!is_valid()) {
            subtract_p();
        }

        return *this;
    }

    [[nodiscard("Optimization")]]
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

    [[nodiscard("Optimization")]]
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
        size_t pos = c_block_number;

        while (pos > 0 && m_blocks[pos - 1] == 0) {
            --pos;
        }

        if (pos == 0) {
            return "0x0";
        }

        block_t up_value = m_blocks[pos - 1] >> 16;
        block_t down_value = m_blocks[pos - 1] & UINT16_MAX;

        auto push = [&](const block_t& value) {
            if (value < 10) {
                result.push_back(value + '0');
            } else {
                result.push_back(value - 10 + 'A');
            }
        };

        if (up_value != 0) {
            push(up_value);
        }

        push(down_value);
        --pos;

        while (pos > 0) {
            up_value = m_blocks[pos - 1] >> 16;
            down_value = m_blocks[pos - 1] & UINT16_MAX;
            push(up_value);
            push(down_value);
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

        constexpr auto is_even = [&](const F_256& value) {
            return (value.m_blocks[0] & 0b1) == 0;
        };

        F_256 u = *this;
        F_256 v(p_values);

        F_256 x_1 = one;
        F_256 x_2;

        while (u != one && v != one) {
            while (is_even(u)) {
                u >>= 1;

                if (is_even(x_1)) {
                    x_1 >>= 1;
                } else {
                    x_1.add_p_uncheck();
                    x_1 >>= 1;
                    assert(x_1.is_valid());
                }
            }

            while (is_even(v)) {
                v >>= 1;

                if (is_even(x_2)) {
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

            block_t symbol_value = 0;

            if (*str >= '0' && *str <= '9') {
                symbol_value = static_cast<block_t>(*str - '0');
            } else if (*str >= 'a' && *str <= 'f') {
                symbol_value = static_cast<block_t>(*str - 'a') + 10;
            } else if (*str >= 'A' && *str <= 'F') {
                symbol_value = static_cast<block_t>(*str - 'A') + 10;
            }

            value += symbol_value;
            ++str;
        }

        return value;
    }

    static constexpr size_t size() {
        return c_block_number;
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

        for (size_t i = 0; i < c_block_number; ++i) {
            result[i] = static_cast<block_t>(value);
            value >>= c_block_size;

            if (value == 0) {
                break;
            }
        }

        return result;
    }

    constexpr void reduce() {
        F_256 s1({m_blocks[0],
                  m_blocks[1],
                  m_blocks[2],
                  m_blocks[3],
                  m_blocks[4],
                  m_blocks[5],
                  m_blocks[6],
                  m_blocks[7]});
        F_256 s2({0, 0, 0, m_blocks[11], m_blocks[12], m_blocks[13], m_blocks[14], m_blocks[15]});
        F_256 s3({0, 0, 0, m_blocks[12], m_blocks[13], m_blocks[14], m_blocks[15], 0});
        F_256 s4({m_blocks[8], m_blocks[9], m_blocks[10], 0, 0, 0, m_blocks[14], m_blocks[15]});
        F_256 s5({m_blocks[9],
                  m_blocks[10],
                  m_blocks[11],
                  m_blocks[13],
                  m_blocks[14],
                  m_blocks[15],
                  m_blocks[13],
                  m_blocks[8]});
        F_256 s6({m_blocks[11], m_blocks[12], m_blocks[13], 0, 0, 0, m_blocks[8], m_blocks[10]});
        F_256 s7({m_blocks[12], m_blocks[13], m_blocks[14], m_blocks[15], 0, 0, m_blocks[9], m_blocks[11]});
        F_256 s8({m_blocks[13],
                  m_blocks[14],
                  m_blocks[15],
                  m_blocks[8],
                  m_blocks[9],
                  m_blocks[10],
                  0,
                  m_blocks[12]});
        F_256 s9({m_blocks[14], m_blocks[15], 0, m_blocks[9], m_blocks[10], m_blocks[11], 0, m_blocks[13]});
        *this = s1 + s2 + s2 + s3 + s3 + s4 + s5 - s6 - s7 - s8 - s9;
        assert(is_valid());
    }

    constexpr void increment() {
        for (size_t i = 0; i < c_block_number; ++i) {
            m_blocks[i] += 1;

            if (m_blocks[i] != 0) {
                break;
            }
        }

        if (m_blocks == p_values) {
            m_blocks = {};
        }

        assert(is_valid());
    }

    static constexpr blocks max_mod_p = {-2, -1, -1, 0, 0, 0, 1, -1};

    constexpr void decrement() {
        if (*this == 0) {
            m_blocks = max_mod_p;
            assert(is_valid());
            return;
        }

        for (size_t i = 0; i < c_block_number; ++i) {
            block_t temp = m_blocks[i];
            m_blocks[i] -= 1;

            if (temp >= m_blocks[i]) {
                break;
            }
        }

        assert(is_valid());
    }

    static constexpr blocks p_values = {-1, -1, -1, 0, 0, 0, 1, -1};

    constexpr void add_p() {
        constexpr F_256 p(p_values);
        *this += p_values;
        assert(is_valid());
    }

    constexpr void add_p_uncheck() {
        block_t carry = 0;

        for (size_t i = 0; i < c_block_number; ++i) {
            block_t sum = carry + p_values[i];
            m_blocks[i] += sum;
            carry = (m_blocks[i] < sum) || (sum < carry);
        }
    }

    constexpr void subtract_p() {
        constexpr F_256 p(p_values);
        *this -= p_values;
        assert(is_valid());
    }

    constexpr bool is_valid() const {
        constexpr F_256 p(p_values);
        return *this < p;
    }

    constexpr const block_t& operator[](size_t pos) const {
        return m_blocks[pos];
    }

    constexpr block_t& operator[](size_t pos) {
        return m_blocks[pos];
    }
};

int main() {}
