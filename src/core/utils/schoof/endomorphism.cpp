#include "endomorphism.h"

#include "utils/wnaf.h"

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
            return g;
        }

        const polynomial::Poly& inverse_denominator = gcd_result.value_multiplier;
        const field::FieldElement& A = f.value()[1];
        End::Element r = F.element(3) * Element::pow(end.m_a, 2) + R.element(Poly(F, {A}));
        r *= R.element(inverse_denominator);
        End::Element a = End::Element::pow(r, 2) * f - end.m_a * F.element(2);
        End::Element b = r * (end.m_a - a) - end.m_b;
        End result(R, a, b, end.m_curve_function);
        return result;
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
            return g;
        }

        const polynomial::Poly& inverse_denominator =
            gcd_result.value_multiplier * field::FieldElement::inverse(gcd_result.gcd[0]);
        End::Element r = (lhs.m_b - rhs.m_b) * R.element(inverse_denominator);
        End::Element a = End::Element::pow(r, 2) * f - lhs.m_a - rhs.m_a;
        End::Element b = r * (lhs.m_a - a) - lhs.m_b;
        End result(R, a, b, lhs.m_curve_function);
        return result;
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

    static End::AdditionResult multiply(End value, const uint& n) {
        if (n == 0) {
            value.nullify();
            return value;
        }

        if ((n & 0b1) == 1) {
            if (n == 1) {
                return value;
            }

            End::AdditionResult var = multiply(value, n - 1);

            if (std::holds_alternative<polynomial::Poly>(var)) {
                return std::get<polynomial::Poly>(var);
            }

            return value + std::get<End>(var);
        } else {
            End::AdditionResult var = multiply(value, n >> 1);

            if (std::holds_alternative<polynomial::Poly>(var)) {
                return std::get<polynomial::Poly>(var);
            }

            End temp = std::get<End>(var);
            return temp + temp;
        }

        /*algorithm::WnafForm wform = algorithm::get_wnaf(n);
        End::AdditionResult var = value + value;

        if (std::holds_alternative<polynomial::Poly>(var)) {
            return std::get<polynomial::Poly>(var);
        }

        End two_value = std::get<End>(var);
        std::vector<End> k_values = {value};

        for (size_t i = 1; i < algorithm::c_k_number; ++i) {
            var = k_values.back() + two_value;

            if (std::holds_alternative<polynomial::Poly>(var)) {
                return std::get<polynomial::Poly>(var);
            }

            k_values.emplace_back(std::get<End>(var));
        }

        value.nullify();

        for (size_t i = wform.size(); i > 0; --i) {
            var = value + value;

            if (std::holds_alternative<polynomial::Poly>(var)) {
                return std::get<polynomial::Poly>(var);
            }

            value = std::get<End>(var);

            if (wform[i - 1].value != 0) {
                if (!wform[i - 1].is_negative) {
                    var = value + k_values[wform[i - 1].value >> 1];

                    if (std::holds_alternative<polynomial::Poly>(var)) {
                        return std::get<polynomial::Poly>(var);
                    }

                    value = std::get<End>(var);
                } else {
                    var = value - k_values[wform[i - 1].value >> 1];

                    if (std::holds_alternative<polynomial::Poly>(var)) {
                        return std::get<polynomial::Poly>(var);
                    }

                    value = std::get<End>(var);
                }
            }
        }

        return value;*/
    }

    End::AdditionResult operator*(const End& end, const uint& value) {
        return multiply(end, value);
    }

    End::AdditionResult operator*(const uint& value, const End& end) {
        return multiply(end, value);
    }

    bool End::operator==(const End& other) const {
        return m_a == other.m_a && m_b == other.m_b;
    }

    End End::operator-() const {
        End result = *this;
        result.m_b = -result.m_b;
        return result;
    }

    End& End::operator*=(const End& other) {
        assert(m_ring == other.m_ring && "End::operator*= : endomorphisms from different rings");
        m_a.compose(other.m_a);
        m_b.compose(other.m_a);
        m_b *= other.m_b;
        return *this;
    }

    void End::nullify() {
        m_a = m_ring.element(Poly(m_ring.get_field(), {0, 1}));
        m_b = m_ring.element(Poly(m_ring.get_field(), {1}));
    }
}   // namespace elliptic_curve_guide::endomorphism
