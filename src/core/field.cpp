#include "field.h"

using ECG::Field;
using ECG::FieldElement;
using ECG::uint;

FieldElement::FieldElement(const uint& value, std::shared_ptr<const uint> modulus) :
    m_value(normalize(value, modulus)), m_modulus(modulus) {
    assert(is_valid() && "FieldElement::FieldElement() : Field element value must be less than p");
}

FieldElement FieldElement::operator-() const {
    return FieldElement(*m_modulus - m_value, m_modulus);
}

// operator+
FieldElement ECG::operator+(const FieldElement& lhs, const FieldElement& rhs) {
    FieldElement result = lhs;
    result += rhs;
    return result;
}

FieldElement ECG::operator+(FieldElement&& lhs, const FieldElement& rhs) {
    lhs += rhs;
    return lhs;
}

FieldElement ECG::operator+(const FieldElement& lhs, FieldElement&& rhs) {
    rhs += lhs;
    return rhs;
}

FieldElement ECG::operator+(FieldElement&& lhs, FieldElement&& rhs) {
    lhs += rhs;
    return lhs;
}

// operator-
FieldElement ECG::operator-(const FieldElement& lhs, const FieldElement& rhs) {
    FieldElement result = lhs;
    result -= rhs;
    return result;
}

FieldElement ECG::operator-(FieldElement&& lhs, const FieldElement& rhs) {
    lhs -= rhs;
    return lhs;
}

FieldElement ECG::operator-(const FieldElement& lhs, FieldElement&& rhs) {
    return -(rhs -= lhs);
}

FieldElement ECG::operator-(FieldElement&& lhs, FieldElement&& rhs) {
    lhs -= rhs;
    return lhs;
}

// operator*
FieldElement ECG::operator*(const FieldElement& lhs, const FieldElement& rhs) {
    FieldElement result = lhs;
    result *= rhs;
    return result;
}

FieldElement ECG::operator*(FieldElement&& lhs, const FieldElement& rhs) {
    lhs *= rhs;
    return lhs;
}

FieldElement ECG::operator*(const FieldElement& lhs, FieldElement&& rhs) {
    rhs *= lhs;
    return rhs;
}

FieldElement ECG::operator*(FieldElement&& lhs, FieldElement&& rhs) {
    lhs *= rhs;
    return lhs;
}

// operator/
FieldElement ECG::operator/(const FieldElement& lhs, const FieldElement& rhs) {
    FieldElement result = lhs;
    result /= rhs;
    return result;
}

FieldElement ECG::operator/(FieldElement&& lhs, const FieldElement& rhs) {
    lhs /= rhs;
    return lhs;
}

FieldElement ECG::operator/(const FieldElement& lhs, FieldElement&& rhs) {
    rhs.inverse();
    return lhs * rhs;
}

FieldElement ECG::operator/(FieldElement&& lhs, FieldElement&& rhs) {
    lhs /= rhs;
    return lhs;
}

FieldElement& FieldElement::operator+=(const FieldElement& other) {
    m_value += other.m_value;

    if (m_value > *m_modulus) {
        m_value -= *m_modulus;
    }

    assert(is_valid() && "FieldElement::operator+= : Field element value must be less than p");
    return *this;
}

FieldElement& FieldElement::operator-=(const FieldElement& other) {
    if (m_value < other.m_value) {
        m_value += *m_modulus;
    }

    m_value -= other.m_value;

    assert(is_valid() && "FieldElement::operator-= : Field element value must be less than p");
    return *this;
}

FieldElement& FieldElement::operator*=(const FieldElement& other) {
    m_value *= other.m_value;

    if (m_value > *m_modulus) {
        m_value %= *m_modulus;
    }

    assert(is_valid() && "FieldElement::operator*= : Field element value must be less than p");
    return *this;
}

FieldElement& FieldElement::operator/=(const FieldElement& other) {
    return (*this *= inverse(other));
}

static uint gcdex(const uint& a, const uint& b, uint& x, uint& y) {
    if (a == 0) {
        x = 0;
        y = 1;
        return b;
    }

    uint x1, y1;
    uint d = gcdex(b % a, a, x1, y1);
    x = y1 - (b / a) * x1;
    y = x1;
    return d;
}

void FieldElement::inverse() {
    uint x, y;
    gcdex(m_value, *m_modulus, x, y);
    m_value = x;

    assert(is_valid() && "FieldElement::inverse : Field element value must be less than p");
}

bool FieldElement::is_invertible() const {
    return m_value != 0;
}

FieldElement FieldElement::pow(const uint& power) const {
    if ((power & 1) != 0) {
        if (power == 1) {
            return *this;
        }

        return *this * pow(power - 1);
    }

    FieldElement temp = pow(power >> 1);
    return temp * temp;
}

FieldElement FieldElement::inverse(const FieldElement& element) {
    FieldElement result = element;
    result.inverse();
    return result;
}

const uint& FieldElement::modulus() const {
    return *m_modulus;
}

const uint& FieldElement::value() const {
    return m_value;
}

uint FieldElement::normalize(const uint& value, std::shared_ptr<const uint> modulus) {
    return value % *modulus;
}

bool FieldElement::is_valid() const {
    return m_value < *m_modulus;
}

Field::Field(const uint& modulus) : m_modulus(std::make_shared<const uint>(modulus)) {};

FieldElement Field::element(const uint& value) {
    return FieldElement(value, m_modulus);
}
