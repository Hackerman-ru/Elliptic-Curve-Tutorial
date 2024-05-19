#include "ring.h"

#include "utils/fast-pow.h"

namespace elliptic_curve_guide::ring {
    RingElement::RingElement(const Poly& value, std::shared_ptr<const Poly> modulus) :
        m_value(normalize(value, modulus)), m_modulus(std::move(modulus)) {}

    RingElement RingElement::pow(const RingElement& element, const uint& power) {
        return algorithm::fast_pow<RingElement>(element, power);
    }

    // operator+
    RingElement operator+(const RingElement& lhs, const RingElement& rhs) {
        RingElement result = lhs;
        result += rhs;
        return result;
    }

    RingElement operator+(RingElement&& lhs, const RingElement& rhs) {
        lhs += rhs;
        return lhs;
    }

    RingElement operator+(const RingElement& lhs, RingElement&& rhs) {
        rhs += lhs;
        return rhs;
    }

    RingElement operator+(RingElement&& lhs, RingElement&& rhs) {
        lhs += rhs;
        return lhs;
    }

    // operator-
    RingElement operator-(const RingElement& lhs, const RingElement& rhs) {
        RingElement result = lhs;
        result -= rhs;
        return result;
    }

    RingElement operator-(RingElement&& lhs, const RingElement& rhs) {
        lhs -= rhs;
        return lhs;
    }

    RingElement operator-(const RingElement& lhs, RingElement&& rhs) {
        return -(rhs -= lhs);
    }

    RingElement operator-(RingElement&& lhs, RingElement&& rhs) {
        lhs -= rhs;
        return lhs;
    }

    // operator*
    RingElement operator*(const RingElement& lhs, const RingElement& rhs) {
        RingElement result = lhs;
        result *= rhs;
        return result;
    }

    RingElement operator*(RingElement&& lhs, const RingElement& rhs) {
        lhs *= rhs;
        return lhs;
    }

    RingElement operator*(const RingElement& lhs, RingElement&& rhs) {
        rhs *= lhs;
        return rhs;
    }

    RingElement operator*(RingElement&& lhs, RingElement&& rhs) {
        lhs *= rhs;
        return lhs;
    }

    RingElement operator*(const RingElement& element, const RingElement::Element& value) {
        RingElement result = element;
        result *= value;
        return result;
    }

    RingElement operator*(RingElement&& element, const RingElement::Element& value) {
        element *= value;
        return element;
    }

    RingElement operator*(const RingElement::Element& value, const RingElement& element) {
        RingElement result = element;
        result *= value;
        return result;
    }

    RingElement operator*(const RingElement::Element& value, RingElement&& element) {
        element *= value;
        return element;
    }

    bool operator==(const RingElement& lhs, const RingElement& rhs) {
        return false;
    }

    bool operator!=(const RingElement& lhs, const RingElement& rhs) {
        return false;
    }

    RingElement RingElement::operator-() const {
        return RingElement(-m_value, m_modulus);
    }

    RingElement& RingElement::operator+=(const RingElement& other) {
        m_value += other.m_value;
        assert(is_valid()
               && "RingElement::operator+= : Degree of ring element must be less than that of modulus");
        return *this;
    }

    RingElement& RingElement::operator-=(const RingElement& other) {
        m_value -= other.m_value;
        assert(is_valid()
               && "RingElement::operator-= : Degree of ring element must be less than that of modulus");
        return *this;
    }

    RingElement& RingElement::operator*=(const RingElement& other) {
        m_value *= other.value();
        normalize();
        assert(is_valid()
               && "RingElement::operator*= : Degree of ring element must be less than that of modulus");
        return *this;
    }

    RingElement& RingElement::operator*=(const Element& value) {
        m_value *= value;
        assert(is_valid()
               && "RingElement::operator*= : Degree of ring element must be less than that of modulus");
        return *this;
    }

    void RingElement::pow(const uint& power) {
        *this = pow(*this, power);
    }

    const RingElement::Poly& RingElement::modulus() const {
        return *m_modulus;
    }

    const RingElement::Poly& RingElement::value() const {
        return m_value;
    }

    RingElement::Poly RingElement::normalize(const Poly& value, std::shared_ptr<const Poly> modulus) {
        return value % *modulus;
    }

    void RingElement::normalize() {
        m_value %= *m_modulus;
        assert(is_valid()
               && "RingElement::normalize : Degree of ring element must be less than that of modulus");
    }

    bool RingElement::is_valid() const {
        return m_value.degree() < m_modulus->degree();
    }

    Ring::Ring(const Poly& modulus) : m_modulus(std::make_shared<const Poly>(modulus)) {};

    RingElement Ring::element(const Poly& value) const {
        return RingElement(value, m_modulus);
    }

    const Ring::Poly& Ring::modulus() const {
        return *m_modulus;
    }

    const field::Field& Ring::get_field() const {
        return m_modulus->get_field();
    }

    bool Ring::operator==(const Ring& other) const {
        return *m_modulus == *other.m_modulus;
    }
}   // namespace elliptic_curve_guide::ring
