#include "field.h"

#include "utils/fast-pow.h"
#include "utils/modulo_inversion.h"

namespace elliptic_curve_guide::field {
    FieldElement::FieldElement(const uint& value, std::shared_ptr<const uint> modulus) :
        m_value(normalize(value, modulus)), m_modulus(std::move(modulus)) {
        assert(is_valid() && "FieldElement::FieldElement() : Field element value must be less than modulus");
    }

    FieldElement FieldElement::operator-() const {
        if (!is_invertible()) {
            return *this;
        }

        return FieldElement(*m_modulus - m_value, m_modulus);
    }

    // operator+
    FieldElement operator+(const FieldElement& lhs, const FieldElement& rhs) {
        FieldElement result = lhs;
        result += rhs;
        return result;
    }

    FieldElement operator+(FieldElement&& lhs, const FieldElement& rhs) {
        lhs += rhs;
        return lhs;
    }

    FieldElement operator+(const FieldElement& lhs, FieldElement&& rhs) {
        rhs += lhs;
        return rhs;
    }

    FieldElement operator+(FieldElement&& lhs, FieldElement&& rhs) {
        lhs += rhs;
        return lhs;
    }

    // operator-
    FieldElement operator-(const FieldElement& lhs, const FieldElement& rhs) {
        FieldElement result = lhs;
        result -= rhs;
        return result;
    }

    FieldElement operator-(FieldElement&& lhs, const FieldElement& rhs) {
        lhs -= rhs;
        return lhs;
    }

    FieldElement operator-(const FieldElement& lhs, FieldElement&& rhs) {
        return -(rhs -= lhs);
    }

    FieldElement operator-(FieldElement&& lhs, FieldElement&& rhs) {
        lhs -= rhs;
        return lhs;
    }

    // operator*
    FieldElement operator*(const FieldElement& lhs, const FieldElement& rhs) {
        FieldElement result = lhs;
        result *= rhs;
        return result;
    }

    FieldElement operator*(FieldElement&& lhs, const FieldElement& rhs) {
        lhs *= rhs;
        return lhs;
    }

    FieldElement operator*(const FieldElement& lhs, FieldElement&& rhs) {
        rhs *= lhs;
        return rhs;
    }

    FieldElement operator*(FieldElement&& lhs, FieldElement&& rhs) {
        lhs *= rhs;
        return lhs;
    }

    // operator/
    FieldElement operator/(const FieldElement& lhs, const FieldElement& rhs) {
        FieldElement result = lhs;
        result /= rhs;
        return result;
    }

    FieldElement operator/(FieldElement&& lhs, const FieldElement& rhs) {
        lhs /= rhs;
        return lhs;
    }

    FieldElement operator/(const FieldElement& lhs, FieldElement&& rhs) {
        rhs.inverse();
        return lhs * rhs;
    }

    FieldElement operator/(FieldElement&& lhs, FieldElement&& rhs) {
        lhs /= rhs;
        return lhs;
    }

    FieldElement operator<<(const FieldElement& value, const uint& shift) {
        FieldElement result = value;
        result <<= shift;
        return result;
    }

    FieldElement operator<<(FieldElement&& value, const uint& shift) {
        value <<= shift;
        return value;
    }

    FieldElement& FieldElement::operator+=(const FieldElement& other) {
        m_value += other.m_value;

        if (m_value >= *m_modulus) {
            m_value -= *m_modulus;
        }

        assert(is_valid() && "FieldElement::operator+= : Field element value must be less than modulus");
        return *this;
    }

    FieldElement& FieldElement::operator-=(const FieldElement& other) {
        if (m_value < other.m_value) {
            m_value += *m_modulus;
        }

        m_value -= other.m_value;

        assert(is_valid() && "FieldElement::operator-= : Field element value must be less than modulus");
        return *this;
    }

    FieldElement& FieldElement::operator*=(const FieldElement& other) {
        m_value *= other.m_value;

        if (m_value >= *m_modulus) {
            m_value %= *m_modulus;
        }

        assert(is_valid() && "FieldElement::operator*= : Field element value must be less than modulus");
        return *this;
    }

    FieldElement& FieldElement::operator/=(const FieldElement& other) {
        return (*this *= inverse(other));
    }

    FieldElement& FieldElement::operator<<=(const uint& shift) {
        for (size_t i = 0; i < shift; ++i) {
            m_value <<= 1;

            if (m_value >= *m_modulus) {
                m_value %= *m_modulus;
            }
        }

        assert(is_valid() && "FieldElement::operator<<= : Field element value must be less than modulus");
        return *this;
    }

    bool operator==(const FieldElement& lhs, const FieldElement& rhs) {
        return lhs.m_value == rhs.m_value;
    }

    void FieldElement::inverse() {
        m_value = algorithm::inverse_modulo(m_value, *m_modulus);
        assert(is_valid() && "FieldElement::inverse : Field element value must be less than modulus");
    }

    bool FieldElement::is_invertible() const {
        return m_value != 0;
    }

    void FieldElement::pow(const uint& power) {
        *this = pow(*this, power);
    }

    FieldElement FieldElement::inverse(const FieldElement& element) {
        FieldElement result = element;
        result.inverse();
        return result;
    }

    FieldElement FieldElement::inverse(FieldElement&& element) {
        element.inverse();
        return element;
    }

    FieldElement FieldElement::pow(const FieldElement& element, const uint& power) {
        if (power == 0) {
            return FieldElement(1, element.m_modulus);
        }

        return algorithm::fast_pow<FieldElement>(element, power);
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
#ifdef ECG_USE_BOOST
    Field::Field(const char* str) : Field(uint(str)) {}

    FieldElement Field::element(const char* str) const {
        return FieldElement(uint(str), m_modulus);
    }
#endif

    Field::Field(const uint& modulus) : m_modulus(std::make_shared<const uint>(modulus)) {};

    FieldElement Field::element(const uint& value) const {
        return FieldElement(value, m_modulus);
    }

    const uint& Field::modulus() const {
        return *m_modulus;
    }

    bool Field::operator==(const Field& other) const {
        return *m_modulus == *other.m_modulus;
    }
}   // namespace elliptic_curve_guide::field
