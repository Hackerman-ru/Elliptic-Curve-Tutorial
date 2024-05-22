#include "endomorphism.h"

#include "utils/wnaf.h"

namespace elliptic_curve_guide::endomorphism {
    End::End(const Poly& a, const Poly& b, std::shared_ptr<const Info> info) :
        m_a_x((*info).ring.element(a)), m_b_x((*info).ring.element(b)), m_info(std::move(info)) {}

    End::End(const RingElement& a, const RingElement& b, std::shared_ptr<const Info> info) :
        m_a_x(a), m_b_x(b), m_info(std::move(info)) {}

    End::End(RingElement&& a, RingElement&& b, std::shared_ptr<const Info> info) :
        m_a_x(std::move(a)), m_b_x(std::move(b)), m_info(std::move(info)) {}

    End::AdditionResult End::twice(const End& end) {
        using namespace algorithm;

        const field::FieldElement& A = (*end.m_info).a;
        const Ring& R = (*end.m_info).ring;
        const Field& F = R.modulus().get_field();
        const RingElement& curve_function = (*end.m_info).curve_function;

        RingElement denominator = F.element(2) * end.m_b_x * curve_function;
        polynomial::Poly inverse_denominator(F);

        if (denominator.value().degree() == 0) {
            inverse_denominator = denominator.value();
            assert(inverse_denominator[0].value() != 0);
            inverse_denominator[0].inverse();
        } else {
            ModulusGcdResult gcd_result = modulus_gcd(denominator.value(), denominator.modulus());
            const polynomial::Poly& g = gcd_result.gcd;

            if (g.degree() > 0) {
                return g;
            }

            const polynomial::Poly& inverse_denominator = gcd_result.value_multiplier;
        }

        RingElement r = F.element(3) * RingElement::pow(end.m_a_x, 2) + R.element(Poly(F, {A}));
        r *= R.element(inverse_denominator);
        RingElement a = RingElement::pow(r, 2) * curve_function - end.m_a_x * F.element(2);
        RingElement b = r * (end.m_a_x - a) - end.m_b_x;
        End result(a, b, end.m_info);
        return result;
    }

    End::AdditionResult operator+(const End& lhs, const End& rhs) {
        using namespace ring;
        using namespace algorithm;

        if (lhs.m_a_x == rhs.m_a_x) {
            return End::twice(lhs);
        }

        const Ring& R = (*lhs.m_info).ring;
        const field::Field F = R.get_field();
        const RingElement& curve_function = (*lhs.m_info).curve_function;

        RingElement denominator = lhs.m_a_x - rhs.m_a_x;
        polynomial::Poly inverse_denominator(F);

        if (denominator.value().degree() == 0) {
            inverse_denominator = denominator.value();
            assert(inverse_denominator[0].value() != 0);
            inverse_denominator[0].inverse();
        } else {
            ModulusGcdResult gcd_result = modulus_gcd(denominator.value(), denominator.modulus());
            const polynomial::Poly& g = gcd_result.gcd;

            if (g.degree() > 0) {
                return g;
            }

            inverse_denominator = gcd_result.value_multiplier;
        }

        RingElement r = (lhs.m_b_x - rhs.m_b_x) * R.element(inverse_denominator);
        RingElement a = RingElement::pow(r, 2) * curve_function - lhs.m_a_x - rhs.m_a_x;
        RingElement b = r * (lhs.m_a_x - a) - lhs.m_b_x;
        End result(a, b, lhs.m_info);
        return result;
    }

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
    }

    End::AdditionResult operator*(const End& end, const uint& value) {
        return multiply(end, value);
    }

    End::AdditionResult operator*(const uint& value, const End& end) {
        return multiply(end, value);
    }

    bool End::operator==(const End& other) const {
        return m_a_x == other.m_a_x && m_b_x == other.m_b_x;
    }

    End End::operator-() const {
        End result = *this;
        result.m_b_x = -result.m_b_x;
        return result;
    }

    End& End::operator*=(const End& other) {
        m_a_x.compose(other.m_a_x);
        m_b_x.compose(other.m_a_x);
        m_b_x *= other.m_b_x;
        return *this;
    }

    void End::nullify() {
        const Ring& R = (*m_info).ring;
        m_a_x = R.element(Poly(R.get_field(), {0, 1}));
        m_b_x = R.element(Poly(R.get_field(), {1}));
    }
}   // namespace elliptic_curve_guide::endomorphism
