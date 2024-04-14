#include "field.h"

using ECG::Field;
using ECG::FieldElement;

FieldElement::FieldElement(const uint& value, std::shared_ptr<const uint> modulus) :
    m_value(normalize(value, modulus)), m_modulus(modulus) {}

FieldElement FieldElement::operator-() const {
    assert(is_valid() && "Field element value must be less than p");

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
    rhs.m_value -= lhs.m_value;

    if (rhs.m_value < 0) {
        rhs.m_value += *lhs.m_modulus;
    }

    return FieldElement(rhs.m_value, lhs.m_modulus);
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
    FieldElement result = lhs;
    result /= rhs;
    return result;
}

FieldElement ECG::operator/(FieldElement&& lhs, FieldElement&& rhs) {
    lhs /= rhs;
    return lhs;
}

FieldElement& FieldElement::operator+=(const FieldElement& other) {
    assert(is_valid() && "FieldElement::operator+= : Field element value must be less than p");
    assert(other.is_valid() && "FieldElement::operator+= : Field element value must be less than p");

    m_value += other.m_value;

    if (m_value > *m_modulus) {
        m_value -= *m_modulus;
    }

    return *this;
}

FieldElement& FieldElement::operator-=(const FieldElement& other) {
    assert(is_valid() && "FieldElement::operator-= : Field element value must be less than p");
    assert(other.is_valid() && "FieldElement::operator-= : Field element value must be less than p");

    m_value -= other.m_value;

    if (m_value < 0) {
        m_value += *m_modulus;
    }

    return *this;
}

FieldElement& FieldElement::operator*=(const FieldElement& other) {
    assert(is_valid() && "FieldElement::operator*= : Field element value must be less than p");
    assert(other.is_valid() && "FieldElement::operator*= : Field element value must be less than p");

    m_value *= other.m_value;

    if (m_value > *m_modulus) {
        m_value %= *m_modulus;
    }

    return *this;
}

FieldElement& FieldElement::operator/=(const FieldElement& other) {
    assert(is_valid() && "FieldElement::operator/= : Field element value must be less than p");
    assert(other.is_valid() && "FieldElement::operator/= : Field element value must be less than p");

    return (*this *= inverse(other));
}

void FieldElement::inverse() {
    assert(is_valid() && "FieldElement::inverse : Field element value must be less than p");

    return;   // TODO
}

bool ECG::FieldElement::is_invertible() const {
    return m_value != 0;
}

FieldElement FieldElement::pow(const uint& power) const {
    assert(is_valid() && "FieldElement::pow : Field element value must be less than p");

    if (power.convert_to<uint32_t>() & 1) {
        if (power == 1) {
            return *this;
        }

        return *this * pow(power - 1);
    }

    FieldElement temp = pow(power >> 1);
    return temp * temp;
}

FieldElement FieldElement::inverse(const FieldElement& element) {
    assert(element.is_valid() && "static FieldElement::inverse : Field element value must be less than p");

    FieldElement result = element;
    result.inverse();
    return result;
}

const ECG::uint& FieldElement::modulus() const {
    return *m_modulus;
}

const ECG::uint& ECG::FieldElement::value() const {
    return m_value;
}

ECG::uint FieldElement::normalize(const uint& value, std::shared_ptr<const uint> modulus) {
    return value % *modulus;
}

bool ECG::FieldElement::is_valid() const {
    return m_value < *m_modulus;
}

ECG::Field::Field(const uint& modulus) : m_modulus(std::make_shared<const uint>(modulus)) {};

FieldElement Field::operator()(const uint& value) {
    return FieldElement(value, m_modulus);
}
