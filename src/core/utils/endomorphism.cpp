#include "endomorphism.h"

namespace elliptic_curve_guide::endomorphism {
    End::End(const Ring& ring, const Poly& a, const Poly& b) :
        m_ring(ring), m_a(ring.element(a)), m_b(ring.element(b)) {}

    End::End(const Ring& ring, const Element& a, const Element& b) : m_ring(ring), m_a(a), m_b(b) {}

    End::End(const Ring& ring, Element&& a, Element&& b) :
        m_ring(ring), m_a(std::move(a)), m_b(std::move(b)) {}

    End::AdditionResult operator+(const End& lhs, const End& rhs) {}

    End::AdditionResult operator-(const End& lhs, const End& rhs) {
        return lhs + (-rhs);
    }

    End operator*(const End& lhs, const End& rhs) {
        End result = lhs;
        result *= rhs;
        return result;
    }

    End operator*(End&& lhs, const End& rhs) {
        lhs *= rhs;
        return lhs;
    }

    End operator*(const End& lhs, End&& rhs) {
        rhs *= lhs;
        return rhs;
    }

    End operator*(End&& lhs, End&& rhs) {
        lhs *= rhs;
        return lhs;
    }

    End operator*(const End& end, const uint& value) {
        End result = end;
        result *= value;
        return result;
    }

    End operator*(End&& end, const uint& value) {
        end *= value;
        return end;
    }

    End operator*(const uint& value, const End& end) {
        End result = end;
        result *= value;
        return result;
    }

    End operator*(const uint& value, End&& end) {
        end *= value;
        return end;
    }

    End End::operator-() const {
        End result = *this;
        result.m_b = -result.m_b;
        return result;
    }

    End& End::operator*=(const End& other) {
        m_a = m_ring.element(polynomial::Poly::compose(m_a.value(), other.m_a.value()));
        m_b = m_ring.element(polynomial::Poly::compose(m_b.value(), other.m_b.value()));
    }

    End& End::operator*=(const uint& value) {
        return *this = algorithm::wnaf_addition<End>(*this, value);
    }

    void End::twice() {}

    void End::nullify() {
        m_a = m_ring.element(Poly(m_ring.get_field(), {0, 1}));
        m_b = m_ring.element(Poly(m_ring.get_field(), {1}));
    }

    void End::change_modulus(const Poly& modulus) {
        m_ring = Ring(modulus);
        m_a = m_ring.element(m_a.value());
        m_b = m_ring.element(m_b.value());
    }
}   // namespace elliptic_curve_guide::endomorphism
