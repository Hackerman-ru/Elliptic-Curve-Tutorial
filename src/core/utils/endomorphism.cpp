#include "endomorphism.h"

namespace elliptic_curve_guide::endomorphism {
    End::End(const Ring& ring, const Poly& a, const Poly& b,
             const std::shared_ptr<const Element>& curve_function) :
        m_ring(ring), m_a(ring.element(a)), m_b(ring.element(b)), m_curve_function(curve_function) {}

    End::End(const Ring& ring, const Element& a, const Element& b,
             const std::shared_ptr<const Element>& curve_function) :
        m_ring(ring), m_a(a), m_b(b), m_curve_function(curve_function) {}

    End::End(const Ring& ring, Element&& a, Element&& b, std::shared_ptr<const Element>&& curve_function) :
        m_ring(ring), m_a(std::move(a)), m_b(std::move(b)), m_curve_function(std::move(curve_function)) {}

    End::AdditionResult End::twice(const End& end) {
        const ring::Ring& R = end.m_ring;
        const field::Field& F = R.modulus().get_field();
        const Element& f = *end.m_curve_function;

        Element denominator = F.element(2) * end.m_b * f;
        algorithm::ModulusGcdResult gcd_result =
            algorithm::modulus_gcd(denominator.value(), denominator.modulus());

        const polynomial::Poly& g = gcd_result.gcd;

        if (g.degree() > 0) {
            return End::AdditionResult {.end = std::nullopt, .g = g};
        }

        assert(g[0].value() == 1 && "End::operator+ : gcd must be 1");
        const polynomial::Poly& inverse_denominator = gcd_result.value_multiplier;
        const field::FieldElement& A = f.value()[1];
        End::Element r = (F.element(3) * Element::pow(end.m_a, 2) + R.element(Poly(F, {A})) * R.element(inverse_denominator);
        End::Element a = End::Element::pow(r, 2) * f - end.m_a * F.element(2);
        End::Element b = r * (end.m_a - a) - end.m_b;
        End result(R, a, b, end.m_curve_function);
        return End::AdditionResult {.end = result, .g = std::nullopt};
    }

    End::AdditionResult operator+(const End& lhs, const End& rhs) {
        assert(lhs.m_ring == rhs.m_ring && "End::operator+ : endomorphisms from different rings");

        if (lhs.m_a == rhs.m_a) {
            return End::twice(lhs);
        }

        End::Element denominator = lhs.m_a - rhs.m_a;
        algorithm::ModulusGcdResult gcd_result =
            algorithm::modulus_gcd(denominator.value(), denominator.modulus());

        const polynomial::Poly& g = gcd_result.gcd;
        const ring::Ring& R = lhs.m_ring;
        const End::Element& f = *lhs.m_curve_function;

        if (g.degree() > 0) {
            return End::AdditionResult {.end = std::nullopt, .g = g};
        }

        assert(g[0].value() == 1 && "End::operator+ : gcd must be 1");
        const polynomial::Poly& inverse_denominator = gcd_result.value_multiplier;
        End::Element r = (lhs.m_b - rhs.m_b) * R.element(inverse_denominator);
        End::Element a = End::Element::pow(r, 2) * f - lhs.m_a - rhs.m_a;
        End::Element b = r * (lhs.m_a - a) - lhs.m_b;
        End result(R, a, b, lhs.m_curve_function);
        return End::AdditionResult {.end = result, .g = std::nullopt};
    }

    End::AdditionResult operator-(const End& lhs, const End& rhs) {
        assert(lhs.m_ring == rhs.m_ring && "End::operator- : endomorphisms from different rings");
        return lhs + (-rhs);
    }

    End operator*(const End& lhs, const End& rhs) {
        assert(lhs.m_ring == rhs.m_ring && "End::operator* : endomorphisms from different rings");
        End result = lhs;
        result *= rhs;
        return result;
    }

    End operator*(End&& lhs, const End& rhs) {
        assert(lhs.m_ring == rhs.m_ring && "End::operator* : endomorphisms from different rings");
        lhs *= rhs;
        return lhs;
    }

    End operator*(const End& lhs, End&& rhs) {
        assert(lhs.m_ring == rhs.m_ring && "End::operator* : endomorphisms from different rings");
        rhs *= lhs;
        return rhs;
    }

    End operator*(End&& lhs, End&& rhs) {
        assert(lhs.m_ring == rhs.m_ring && "End::operator* : endomorphisms from different rings");
        lhs *= rhs;
        return lhs;
    }

    End::AdditionResult operator*(const End& end, const uint& value) {
        // TODO : double and add
    }

    End::AdditionResult operator*(const uint& value, const End& end) {
        // TODO : double and add
    }

    End End::operator-() const {
        End result = *this;
        result.m_b = -result.m_b;
        return result;
    }

    End& End::operator*=(const End& other) {
        assert(m_ring == other.m_ring && "End::operator*= : endomorphisms from different rings");
        m_a = m_ring.element(polynomial::Poly::compose(m_a.value(), other.m_a.value()));
        m_b = m_ring.element(polynomial::Poly::compose(m_b.value(), other.m_b.value()));
    }

    End& End::operator*=(const uint& value) {
        return *this = algorithm::wnaf_addition<End>(*this, value);
    }

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
