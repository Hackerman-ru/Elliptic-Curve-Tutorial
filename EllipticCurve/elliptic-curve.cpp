#include "elliptic-curve.h"

using ECG::ECE;

ECE::ECE(const Field& x) : m_x(x) {
    if (!find_y(x, &m_y)) {
        m_x = Field(1);
        m_y = Field(0);
    }
}

void ECE::set_params(const Field& a, const Field& b) {
    m_a = a;
    m_b = b;
}

ECE ECE::operator+(const ECE& other) const {
    if (is_inf() || other.is_inf()) {
        return (is_inf() ? other : *this);
    }

    if (*this == -other) {
        return ECE(1, 0);
    }

    Field k;

    if (*this == other) {
        k = (Field(3) * m_x * m_x + m_a) / (Field(2) * m_y);   // char(F) > 3 assumption
    } else {
        k = (other.m_y - m_y) / (other.m_x - m_x);
    }

    Field x = k * k - m_x - other.m_x;
    Field y = -k * (x - m_x) - m_y;
    return ECE(x, y);
}

ECE ECE::operator-(const ECE& other) const {
    return *this + (-other);
}

ECE ECE::operator-() const {
    return ECE(m_x, -m_y);
}

bool ECE::is_inf() const {
    return (m_x == Field(1) && m_y == Field(0));
}

ECE ECE::find_point(const std::string& str) {
    while (true) {
        const size_t STR_SIZE = str.size();
        const size_t BYTES = (Field::get_p() >> 4).convert_to<size_t>();
        Field::uint512_t x_ = generate_random_uint();

        for (size_t i = 0; i < STR_SIZE && i < BYTES; ++i) {
            x_ <<= 3;
            x_ |= Field::uint(str[STR_SIZE - i - 1]);
        }

        Field x(x_);
        Field y;

        if (find_y(x, &y)) {
            return ECE(x, y);
        }
    };
}

ECE& ECE::operator+=(const ECE& other) {
    return (*this = *this + other);
}

ECE& ECE::operator-=(const ECE& other) {
    return (*this = *this - other);
}

bool ECE::operator==(const ECE& other) const {
    return (m_x == other.m_x && m_y == other.m_y);
}

bool ECE::find_y(const Field& x, Field* y) {
    Field a = x * x * x + x * m_a + m_b;
    Field a_pow_spec = a.fast_pow((Field::get_p() - Field::uint512_t(1)) >> 1);

    if (a_pow_spec != Field(1)) {
        return false;
    }

    if ((Field::get_p() & Field::uint512_t(0b11)) == Field::uint512_t(0b11)) {
        *y = a.fast_pow((Field::get_p() + Field::uint512_t(1)) >> 2);
        return true;
    }

    // TODO: insert return statement here
}

ECE::Field::uint512_t ECE::generate_random_uint() {
    static constexpr size_t BITS = sizeof(Field::uint512_t) * CHAR_BIT;
    Field::uint512_t result;

    for (size_t i = 0; i < BITS; ++i) {
        result <<= 1;
        result |= Field::uint((rand()) & 0b1);
    }

    return result;
}

ECE::Field ECE::find_n() {   // SEA algorithm
    // TODO: insert return statement here
}
