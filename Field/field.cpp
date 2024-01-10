#include "Field.h"

using ECG::PFE;

void PFE::set_p(const uint& p) {
    m_p = p;
}

PFE::PFE(const std::string& str, StringType type) : m_value(str, type) {
    m_value %= m_p;
    assert(m_value < m_p);
}

PFE PFE::operator+(const PFE& other) const {
    assert(m_value < m_p && other.m_value < m_p);

    uint result = m_value + other.m_value;

    if (result > m_p) {
        result -= m_p;
    }

    return PFE(result, true);
}

PFE PFE::operator-(const PFE& other) const {
    assert(m_value < m_p && other.m_value < m_p);

    uint result = m_value;

    if (result < other.m_value) {
        result += m_p;
    }

    result -= other.m_value;

    if (result > m_p) {
        result -= m_p;
    }

    return PFE(result, true);
}

PFE PFE::operator*(const PFE& other) const {
    assert(m_value < m_p && other.m_value < m_p);

    uint result = m_value * other.m_value;
    return PFE(result, false);
}

PFE PFE::operator/(const PFE& other) const {
    return *this * other.inverse();
}

PFE PFE::inverse() const {   // Extended Euclidean algorithm
    assert(m_value < m_p);
    PFE t(0);
    uint r(m_p);
    PFE new_t(1);
    uint new_r(m_value);

    while (new_r != uint(0)) {
        uint quotien = r / new_r;
        PFE temp(quotien * new_t.m_value);

        std::make_pair(t, new_t) = std::make_pair(new_t, t - temp);
        std::make_pair(r, new_r) = std::make_pair(new_r, r % new_r);
    }

    assert(r == uint(1));
    return t;
}

PFE PFE::operator+=(const PFE& other) {
    return (*this = *this + other);
}

PFE PFE::operator-=(const PFE& other) {
    return (*this = *this - other);
}

PFE PFE::operator*=(const PFE& other) {
    return (*this = *this * other);
}

PFE PFE::operator/=(const PFE& other) {
    return (*this = *this / other);
}

std::string PFE::into_string(StringType str_type) const {
    return m_value.into_string(str_type);
}

std::string PFE::into_string(auto map, size_t shift) const {
    return m_value.into_string(map, shift);
}
