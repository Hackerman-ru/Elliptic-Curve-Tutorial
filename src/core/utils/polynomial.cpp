#include "polynomial.h"

#include "utils/fast-pow.h"
#include "utils/ring.h"

namespace elliptic_curve_guide::polynomial {

    Poly Poly::pow(const Poly& poly, const uint& power) {
        return algorithm::fast_pow<Poly>(poly, power);
    }

    Poly::Poly(const Field& field) : m_field(field), m_coeffs({field.element(0)}) {};
    Poly::Poly(const Field& field, const std::vector<Element>& coeffs) : m_field(field), m_coeffs(coeffs) {};
    Poly::Poly(const Field& field, std::vector<Element>&& coeffs) :
        m_field(field), m_coeffs(std::move(coeffs)) {};

    static std::vector<field::FieldElement> convert_to_field_coeffs(const field::Field& field,
                                                                    const std::vector<uint>& coeffs) {
        std::vector<field::FieldElement> result;
        result.reserve(coeffs.size());

        for (const auto& coef : coeffs) {
            result.emplace_back(field.element(coef));
        }

        return result;
    }

    Poly::Poly(const Field& field, const std::vector<uint>& coeffs) :
        m_field(field), m_coeffs(convert_to_field_coeffs(field, coeffs)) {};

    Poly operator+(const Poly& lhs, const Poly& rhs) {
        Poly result = lhs;
        result += rhs;
        return result;
    }

    Poly operator+(Poly&& lhs, const Poly& rhs) {
        lhs += rhs;
        return lhs;
    }

    Poly operator+(const Poly& lhs, Poly&& rhs) {
        rhs += lhs;
        return rhs;
    }

    Poly operator+(Poly&& lhs, Poly&& rhs) {
        lhs += rhs;
        return lhs;
    }

    Poly operator-(const Poly& lhs, const Poly& rhs) {
        Poly result = lhs;
        result -= rhs;
        return result;
    }

    Poly operator-(Poly&& lhs, const Poly& rhs) {
        lhs -= rhs;
        return lhs;
    }

    Poly operator-(const Poly& lhs, Poly&& rhs) {
        rhs -= lhs;
        rhs.negative();
        return rhs;
    }

    Poly operator-(Poly&& lhs, Poly&& rhs) {
        lhs -= rhs;
        return lhs;
    }

    Poly operator*(const Poly& lhs, const Poly& rhs) {
        const Poly::Field& field = lhs.get_field();

        size_t degree = lhs.degree() + rhs.degree();
        Poly result(field);
        result.m_coeffs.resize(degree + 1, field.element(0));

        for (size_t lhs_pos = 0; lhs_pos < lhs.len(); ++lhs_pos) {
            for (size_t rhs_pos = 0; rhs_pos < rhs.len(); ++rhs_pos) {
                result[lhs_pos + rhs_pos] += lhs[lhs_pos] * rhs[rhs_pos];
            }
        }

        result.clean();
        assert(result.is_valid() && "Poly::operator%= : invalid representation of polynomial");
        return result;
    }

    Poly operator*(const Poly& poly, const Poly::Element& value) {
        Poly result = poly;
        result *= value;
        return result;
    }

    Poly operator*(Poly&& poly, const Poly::Element& value) {
        poly *= value;
        return poly;
    }

    Poly operator*(const Poly::Element& value, const Poly& poly) {
        Poly result = poly;
        result *= value;
        return result;
    }

    Poly operator*(const Poly::Element& value, Poly&& poly) {
        poly *= value;
        return poly;
    }

    Poly operator%(const Poly& lhs, const Poly& rhs) {
        Poly result = lhs;
        result %= rhs;
        return result;
    }

    Poly operator%(Poly&& lhs, const Poly& rhs) {
        lhs %= rhs;
        return lhs;
    }

    Poly Poly::operator-() const {
        Poly result = *this;
        result.negative();
        return result;
    }

    Poly& Poly::operator+=(const Poly& other) {
        size_t other_len = other.len();
        m_coeffs.resize(std::max(len(), other_len), m_field.element(0));

        for (size_t i = 0; i < other_len; ++i) {
            m_coeffs[i] += other[i];
        }

        clean();
        assert(is_valid() && "Poly::operator+= : invalid representation of polynomial");
        return *this;
    }

    Poly& Poly::operator-=(const Poly& other) {
        size_t other_len = other.len();
        m_coeffs.resize(std::max(len(), other_len), m_field.element(0));

        for (size_t i = 0; i < other_len; ++i) {
            m_coeffs[i] -= other.m_coeffs[i];
        }

        clean();
        assert(is_valid() && "Poly::operator-= : invalid representation of polynomial");
        return *this;
    }

    Poly& Poly::operator*=(const Poly& other) {
        return *this = *this * other;
    }

    Poly& Poly::operator*=(const Element& value) {
        for (auto& coeff : m_coeffs) {
            coeff *= value;
        }

        clean();
        assert(is_valid() && "Poly::operator*= : invalid representation of polynomial");
        return *this;
    }

    Poly& Poly::operator%=(const Poly& other) {
        if (other.degree() == 0) {
            return *this;
        }

        while (degree() >= other.degree()) {
            Element factor = top_coef() / other.top_coef();
            *this -= other * factor;
        }

        assert(is_valid() && "Poly::operator%= : invalid representation of polynomial");
        return *this;
    }

    void Poly::pow(const uint& power) {
        *this = pow(*this, power);
        assert(is_valid() && "Poly::pow : invalid representation of polynomial");
    }

    size_t Poly::degree() const {
        return m_coeffs.size() - 1;
    }

    const Poly::Field& Poly::get_field() const {
        return m_field;
    }

    Poly::Element& Poly::operator[](const size_t& pos) {
        return m_coeffs[pos];
    }

    const Poly::Element& Poly::operator[](const size_t& pos) const {
        return m_coeffs[pos];
    }

    size_t Poly::len() const {
        return m_coeffs.size();
    }

    void Poly::clean() {
        size_t size = m_coeffs.size();

        while (size > 0) {
            if (!m_coeffs[size].is_invertible()) {
                m_coeffs.pop_back();
                --size;
            } else {
                break;
            }
        }

        if (size == 0) {
            m_coeffs = {m_field.element(0)};
        }
    }

    void Poly::negative() {
        for (auto& coef : m_coeffs) {
            coef = -coef;
        }

        assert(is_valid() && "Poly::negative : invalid representation of polynomial");
    }

    bool Poly::is_valid() const {
        size_t size = len();

        if (size == 0) {
            return false;
        }

        if (!m_coeffs[size - 1].is_invertible()) {
            return false;
        }

        return true;
    }

    const Poly::Element& Poly::top_coef() const {
        return m_coeffs[degree()];
    }

}   // namespace elliptic_curve_guide::polynomial

namespace elliptic_curve_guide::algorithm {
    polynomial::Poly gcd(const polynomial::Poly& lhs, const polynomial::Poly& rhs) {
        size_t lhs_degree = lhs.degree();

        if (lhs_degree == 0) {
            return lhs;
        }

        size_t rhs_degree = rhs.degree();

        if (rhs_degree == 0) {
            return rhs;
        }

        if (lhs_degree >= rhs_degree) {
            return gcd(lhs % rhs, rhs);
        } else {
            return gcd(lhs, rhs % lhs);
        }
    }

    bool has_root(const polynomial::Poly& poly) {
        const field::Field& F = poly.get_field();
        const ring::Ring R(poly);

        polynomial::Poly variable(F, {0, 1});
        ring::RingElement x = R.element(variable);
        ring::RingElement x_p_minus_x = ring::RingElement::pow(x, F.modulus()) - x;

        polynomial::Poly gcd = algorithm::gcd(poly, x_p_minus_x.value());
        return gcd.degree() > 0;
    }
}   // namespace elliptic_curve_guide::algorithm
