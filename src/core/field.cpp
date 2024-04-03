#include "Field.h"

using ECG::Field;
using ECG::FieldElement;

FieldElement::FieldElement(uint value, std::shared_ptr<const uint> p) : m_value(value), m_p(p) {
    if (m_value > *m_p) {
        m_value %= *m_p;
    }
}

FieldElement FieldElement::operator-() const {
    assert(m_value < *m_p && "Field element value must be less than p");

    return FieldElement(*m_p - m_value, m_p);
}

// operator+
FieldElement ECG::operator+(const FieldElement& lhs, const FieldElement& rhs) {
    FieldElement result = lhs;
    return result += rhs;
}

FieldElement ECG::operator+(FieldElement&& lhs, const FieldElement& rhs) {
    return lhs += rhs;
}

FieldElement ECG::operator+(const FieldElement& lhs, FieldElement&& rhs) {
    return rhs += lhs;
}

FieldElement ECG::operator+(FieldElement&& lhs, FieldElement&& rhs) {
    return lhs += rhs;
}

// operator-
FieldElement ECG::operator-(const FieldElement& lhs, const FieldElement& rhs) {
    FieldElement result = lhs;
    return result -= rhs;
}

FieldElement ECG::operator-(FieldElement&& lhs, const FieldElement& rhs) {
    return lhs -= rhs;
}

FieldElement ECG::operator-(const FieldElement& lhs, FieldElement&& rhs) {
    return -(rhs -= lhs);
}

FieldElement ECG::operator-(FieldElement&& lhs, FieldElement&& rhs) {
    return lhs -= rhs;
}

// operator*
FieldElement ECG::operator*(const FieldElement& lhs, const FieldElement& rhs) {
    FieldElement result = lhs;
    return result *= rhs;
}

FieldElement ECG::operator*(FieldElement&& lhs, const FieldElement& rhs) {
    return lhs *= rhs;
}

FieldElement ECG::operator*(const FieldElement& lhs, FieldElement&& rhs) {
    return rhs *= lhs;
}

FieldElement ECG::operator*(FieldElement&& lhs, FieldElement&& rhs) {
    return lhs *= rhs;
}

// operator/
FieldElement ECG::operator/(const FieldElement& lhs, const FieldElement& rhs) {
    FieldElement result = lhs;
    return result /= rhs;
}

FieldElement ECG::operator/(FieldElement&& lhs, const FieldElement& rhs) {
    return lhs /= rhs;
}

FieldElement ECG::operator/(const FieldElement& lhs, FieldElement&& rhs) {
    FieldElement result = lhs;
    return result /= rhs;
}

FieldElement ECG::operator/(FieldElement&& lhs, FieldElement&& rhs) {
    return lhs /= rhs;
}

bool ECG::operator==(const FieldElement& lhs, const FieldElement& rhs) {
    return lhs.m_value == rhs.m_value;
}

bool ECG::operator!=(const FieldElement& lhs, const FieldElement& rhs) {
    return lhs.m_value != rhs.m_value;
}

bool ECG::operator>(const FieldElement& lhs, const FieldElement& rhs) {
    return lhs.m_value > rhs.m_value;
}

bool ECG::operator<(const FieldElement& lhs, const FieldElement& rhs) {
    return lhs.m_value < rhs.m_value;
}

bool ECG::operator>=(const FieldElement& lhs, const FieldElement& rhs) {
    return lhs.m_value >= rhs.m_value;
}

bool ECG::operator<=(const FieldElement& lhs, const FieldElement& rhs) {
    return lhs.m_value <= rhs.m_value;
}

FieldElement& FieldElement::operator+=(const FieldElement& other) {
    assert(m_value < *m_p && "Field element value must be less than p");
    assert(other.m_value < *m_p && "Field element value must be less than p");

    m_value += other.m_value;

    if (m_value > *m_p) {
        m_value -= *m_p;
    }

    return *this;
}

FieldElement& FieldElement::operator-=(const FieldElement& other) {
    assert(m_value < *m_p && "Field element value must be less than p");
    assert(other.m_value < *m_p && "Field element value must be less than p");

    if (m_value < other.m_value) {
        m_value += *m_p;
    }

    m_value -= other.m_value;
    return *this;
}

FieldElement& FieldElement::operator*=(const FieldElement& other) {
    assert(m_value < *m_p && "Field element value must be less than p");
    assert(other.m_value < *m_p && "Field element value must be less than p");

    m_value *= other.m_value;

    if (m_value > *m_p) {
        m_value %= *m_p;
    }

    return *this;
}

FieldElement& FieldElement::operator/=(const FieldElement& other) {
    assert(m_value < *m_p && "Field element value must be less than p");
    assert(other.m_value < *m_p && "Field element value must be less than p");

    return (*this *= other.inverse());
}

FieldElement FieldElement::inverse() const {
    assert(m_value < *m_p && "Field element value must be less than p");

    return *this;   // TODO
}

bool ECG::FieldElement::is_inversible() const {
    return m_value != 0;
}

FieldElement FieldElement::fast_pow(const uint& pow) const {
    assert(m_value < *m_p && "Field element value must be less than p");

    if (pow.convert_to<uint32_t>() & 1) {
        if (pow == 1) {
            return *this;
        }

        return *this * fast_pow(pow - 1);
    }

    FieldElement temp = fast_pow(pow >> 1);
    return temp * temp;
}

const ECG::uint& FieldElement::get_p() const {
    return *m_p;
}

ECG::Field::Field(uint p) : m_p(std::make_shared<const uint>(std::move(p))) {};

FieldElement Field::operator()(uint value) {
    return FieldElement(std::move(value), m_p);
}
