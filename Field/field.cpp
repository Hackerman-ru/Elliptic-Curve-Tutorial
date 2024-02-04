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

std::string FieldElement::into_string(StringType str_type) const {
    assert(m_value < *m_p && "Field element value must be less than p");

#ifdef ECG_BOOST
    std::string result;
    uint clone = m_value;

    switch (str_type) {
    case ECG::StringType::BINARY :
        while (clone != 0) {
            result += ((clone & 1) != 0) + '0';
            clone >>= 1;
        }

        std::reverse(result.begin(), result.end());
        break;
    case ECG::StringType::DECIMAL :
        result = m_value.convert_to<std::string>();
        break;
    case ECG::StringType::HEXADECIMAL :
        while (clone != 0) {
            auto n = clone.convert_to<uint32_t>() & 0xF;

            if (n < 10) {
                result += n + '0';
            } else {
                n -= 10;
                result += n + 'a';
            }

            clone >>= 4;
        }

        std::reverse(result.begin(), result.end());
        break;
    }

    return result;
#else
    return m_value.into_string(str_type);
#endif
}

std::string FieldElement::into_string(std::function<char(uint32_t)> map, size_t shift) const {
    assert(m_value < *m_p && "Field element value must be less than p");

#ifdef ECG_BOOST
    std::string result;
    uint clone = m_value;

    do {
        result += map(static_cast<uint32_t>(clone));
        clone >>= shift;
    } while (clone > 0);

    std::reverse(result.begin(), result.end());

    return result;
#else
    return m_value.into_string(map, shift);
#endif
}

const ECG::uint& FieldElement::get_p() const {
    return *m_p;
}

ECG::Field::Field(uint p) : m_p(std::make_shared<const uint>(std::move(p))) {};

FieldElement Field::operator()(uint value) {
    return FieldElement(std::move(value), m_p);
}
