#include "Field.h"

using ECG::Field;
using ECG::FieldElement;

FieldElement::FieldElement(std::unique_ptr<ECG::FieldElement::BaseElement>&& ptr) : m_self(std::move(ptr)) {};
FieldElement::FieldElement(uint value, std::shared_ptr<const uint> p) :
    m_self(std::make_unique<PrimeElement>(std::move(value), std::move(p))) {};

FieldElement::FieldElement(std::vector<uint> poly, std::shared_ptr<const std::vector<uint>> reducer) :
    m_self(std::make_unique<PolyElement>(std::move(poly), std::move(reducer))) {};

FieldElement::FieldElement(const FieldElement& other) : m_self(other.m_self->copy_()) {};

FieldElement& FieldElement::operator=(const FieldElement& other) {
    return *this = FieldElement(other);
}

FieldElement FieldElement::operator-() const {
    return {m_self->operator-()};
}

// operator+
FieldElement operator+(const FieldElement& lhs, const FieldElement& rhs) {
    FieldElement result = lhs;
    return result += rhs;
}

FieldElement operator+(FieldElement&& lhs, const FieldElement& rhs) {
    return lhs += rhs;
}

FieldElement operator+(const FieldElement& lhs, FieldElement&& rhs) {
    return rhs += lhs;
}

FieldElement operator+(FieldElement&& lhs, FieldElement&& rhs) {
    return lhs += rhs;
}

// operator-
FieldElement operator-(const FieldElement& lhs, const FieldElement& rhs) {
    FieldElement result = lhs;
    return result -= rhs;
}

FieldElement operator-(FieldElement&& lhs, const FieldElement& rhs) {
    return lhs -= rhs;
}

FieldElement operator-(const FieldElement& lhs, FieldElement&& rhs) {
    return -(rhs -= lhs);
}

FieldElement operator-(FieldElement&& lhs, FieldElement&& rhs) {
    return lhs -= rhs;
}

// operator*
FieldElement operator*(const FieldElement& lhs, const FieldElement& rhs) {
    FieldElement result = lhs;
    return result *= rhs;
}

FieldElement operator*(FieldElement&& lhs, const FieldElement& rhs) {
    return lhs *= rhs;
}

FieldElement operator*(const FieldElement& lhs, FieldElement&& rhs) {
    return rhs *= lhs;
}

FieldElement operator*(FieldElement&& lhs, FieldElement&& rhs) {
    return lhs *= rhs;
}

// operator/
FieldElement operator/(const FieldElement& lhs, const FieldElement& rhs) {
    FieldElement result = lhs;
    return result /= rhs;
}

FieldElement operator/(FieldElement&& lhs, const FieldElement& rhs) {
    return lhs /= rhs;
}

FieldElement operator/(const FieldElement& lhs, FieldElement&& rhs) {
    FieldElement result = lhs;
    return result /= rhs;
}

FieldElement operator/(FieldElement&& lhs, FieldElement&& rhs) {
    return lhs /= rhs;
}

std::unique_ptr<ECG::FieldElement::BaseElement> FieldElement::PrimeElement::operator-() const {
    return std::make_unique<PrimeElement>(*m_p - m_value, m_p);
}

FieldElement& FieldElement::operator+=(const FieldElement& other) {
    m_self->operator+=(*other.m_self);
    return *this;
}

FieldElement& FieldElement::operator-=(const FieldElement& other) {
    m_self->operator-=(*other.m_self);
    return *this;
}

FieldElement& FieldElement::operator*=(const FieldElement& other) {
    m_self->operator*=(*other.m_self);
    return *this;
}

FieldElement& FieldElement::operator/=(const FieldElement& other) {
    m_self->operator/=(*other.m_self);
    return *this;
}

FieldElement FieldElement::inverse() const {
    return {m_self->inverse()};
}

FieldElement FieldElement::fast_pow(const uint& pow) const {
    if ((pow & 1) == 1) {
        if (pow == 1) {
            return *this;
        }

        return *this * fast_pow(pow - 1);
    }

    FieldElement temp = fast_pow(pow >> 1);
    return temp * temp;
}

std::string FieldElement::into_string(StringType str_type) const {
    return m_self->into_string(str_type);
}

std::string FieldElement::into_string(std::function<char(uint32_t)> map, size_t shift) const {
    return m_self->into_string(map, shift);
}

FieldElement::PrimeElement::PrimeElement(const uint& value, const std::shared_ptr<const uint>& p) :
    m_value(value), m_p(p) {
    m_value %= *m_p;
}

bool ECG::FieldElement::PrimeElement::operator==(const BaseElement& other) const {
    return m_value == static_cast<const PrimeElement&>(other).m_value;
}

void FieldElement::PrimeElement::operator+=(const BaseElement& other) {
    m_value += static_cast<const PrimeElement&>(other).m_value;

    if (m_value > *m_p) {
        m_value -= *m_p;
    }
}

void FieldElement::PrimeElement::operator-=(const BaseElement& other) {
    m_value -= static_cast<const PrimeElement&>(other).m_value;

    if (m_value > *m_p) {
        m_value -= *m_p;
    }
}

void FieldElement::PrimeElement::operator*=(const BaseElement& other) {
    m_value *= static_cast<const PrimeElement&>(other).m_value;
    m_value %= *m_p;
}

void FieldElement::PrimeElement::operator/=(const BaseElement& other) {
    *this *= *other.inverse();
}

std::unique_ptr<ECG::FieldElement::BaseElement> FieldElement::PrimeElement::inverse() const {
    return std::make_unique<PrimeElement>(*this);   // TODO
}

std::unique_ptr<ECG::FieldElement::BaseElement> FieldElement::PrimeElement::copy_() const {
    return std::make_unique<PrimeElement>(*this);
}

std::string FieldElement::PrimeElement::into_string(StringType str_type) const {
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

std::string FieldElement::PrimeElement::into_string(std::function<char(uint32_t)> map, size_t shift) const {
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

FieldElement::PolyElement::PolyElement(const std::vector<uint>& poly,
                                       const std::shared_ptr<const std::vector<uint>>& reducer) :
    m_poly(poly), m_reducer(reducer) {}

bool ECG::FieldElement::PolyElement::operator==(const BaseElement& other) const {
    return m_poly == static_cast<const PolyElement&>(other).m_poly;
}

std::unique_ptr<ECG::FieldElement::BaseElement> FieldElement::PolyElement::operator-() const {
    return std::make_unique<PolyElement>(*this);   // TODO
}

void FieldElement::PolyElement::operator+=(const BaseElement& other) {
    // TODO
}

void FieldElement::PolyElement::operator-=(const BaseElement& other) {
    // TODO
}

void FieldElement::PolyElement::operator*=(const BaseElement& other) {
    // TODO
}

void FieldElement::PolyElement::operator/=(const BaseElement& other) {
    *this *= *other.inverse();
}

std::string FieldElement::PolyElement::into_string(StringType str_type) const {
    return "0";   // TODO
}

std::unique_ptr<ECG::FieldElement::BaseElement> FieldElement::PolyElement::inverse() const {
    return std::make_unique<PolyElement>(*this);   // TODO
}

std::unique_ptr<ECG::FieldElement::BaseElement> FieldElement::PolyElement::copy_() const {
    return std::make_unique<PolyElement>(*this);
}

std::string FieldElement::PolyElement::into_string(std::function<char(uint32_t)> map, size_t shift) const {
    return "0";   // TODO
}

Field::Field(uint p, uint n) : m_p(std::make_shared<const uint>(std::move(p))), m_n(n) {
    if (n > 1) {
        // TODO
    }
}

FieldElement Field::create_element(uint value) {
    return FieldElement(std::move(value), m_p);
}

FieldElement Field::create_element(std::vector<uint> value) {
    return FieldElement(std::move(value), m_reducer);
}
