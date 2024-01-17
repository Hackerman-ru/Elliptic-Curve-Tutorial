#ifndef ECG_FIELD_HPP
#define ECG_FIELD_HPP
#include "field.h"

using ECG::Field;
using ECG::FieldElement;

template<const Field* field>
bool FieldElement<field>::operator==(const FieldElement& other) const {
    return m_value == other.m_value;
}

template<const Field* field>
FieldElement<field> FieldElement<field>::operator+(const FieldElement& other) const {
    static const uint& p = field->get_p();
    assert(m_value < p && other.m_value < p);

    uint result = m_value + other.m_value;

    if (result > p) {
        result -= p;
    }

    return FieldElement(result);
}

template<const Field* field>
FieldElement<field> FieldElement<field>::operator-() const {
    static const uint& p = field->get_p();
    return FieldElement(p - m_value);
}

template<const Field* field>
FieldElement<field> FieldElement<field>::fast_pow(const ECG::uint& pow) const {
    static const uint& p = field->get_p();
    assert(m_value < p);

    if (pow == uint(1)) {
        return *this;
    }

    if ((pow & uint(0b1)) == 0b1) {
        FieldElement f = fast_pow(pow - uint(1));
        FieldElement fe = *this * f;
        return fe;
    }

    FieldElement temp = fast_pow(pow >> 1);
    FieldElement t = temp * temp;
    return t;
}

template<const Field* field>
FieldElement<field> FieldElement<field>::operator-(const FieldElement& other) const {
    static const uint& p = field->get_p();
    assert(m_value < p && other.m_value < p);

    uint result = m_value;

    if (result < other.m_value) {
        result += p;
    }

    result -= other.m_value;

    if (result > p) {
        result -= p;
    }

    return FieldElement(result);
}

template<const Field* field>
FieldElement<field> FieldElement<field>::operator*(const FieldElement& other) const {
    static const uint& p = field->get_p();
    assert(m_value < p && other.m_value < p);

    uint result = m_value * other.m_value;
    return FieldElement(result);
}

template<const Field* field>
FieldElement<field> FieldElement<field>::operator/(const FieldElement& other) const {
    static const uint& p = field->get_p();
    assert(m_value < p);
    return *this * other.inverse();
}

template<const Field* field>
FieldElement<field> FieldElement<field>::inverse() const {   // Extended Euclidean algorithm
    static const uint& p = field->get_p();
    assert(m_value < p);
    FieldElement t(uint(0));
    uint r(p);
    FieldElement new_t(uint(1));
    uint new_r(m_value);

    while (new_r != 0) {
        uint quotien = r / new_r;
        FieldElement temp(quotien * new_t.m_value, p);

        std::make_pair(t, new_t) = std::make_pair(new_t, t - temp);
        std::make_pair(r, new_r) = std::make_pair(new_r, r % new_r);
    }

    assert(r == 1);
    return t;
}

template<const Field* field>
FieldElement<field> FieldElement<field>::operator+=(const FieldElement& other) {
    static const uint& p = field->get_p();
    assert(m_value < p && other.m_value < p);
    return (*this = *this + other);
}

template<const Field* field>
FieldElement<field> FieldElement<field>::operator-=(const FieldElement& other) {
    static const uint& p = field->get_p();
    assert(m_value < p && other.m_value < p);
    return (*this = *this - other);
}

template<const Field* field>
FieldElement<field> FieldElement<field>::operator*=(const FieldElement& other) {
    static const uint& p = field->get_p();
    assert(m_value < p && other.m_value < p);
    return (*this = *this * other);
}

template<const Field* field>
FieldElement<field> FieldElement<field>::operator/=(const FieldElement& other) {
    static const uint& p = field->get_p();
    assert(m_value < p && other.m_value < p);
    return (*this = *this / other);
}

template<const Field* field>
std::string FieldElement<field>::into_string(StringType str_type) const {
    static const uint& p = field->get_p();
    assert(m_value < p);
    return m_value.into_string(str_type);
}

template<const Field* field>
std::string FieldElement<field>::into_string(std::function<uint32_t(char)> map, size_t shift) const {
    static const uint& p = field->get_p();
    assert(m_value < p);
    return m_value.into_string(map, shift);
}

template<const Field* field>
FieldElement<field>::FieldElement(const ECG::uint& value) : m_value(value) {
    static const uint& p = field->get_p();

    if (m_value >= p) {
        m_value %= p;
    }

    assert(m_value < p);
}

#endif
