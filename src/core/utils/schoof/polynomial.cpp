#include "polynomial.h"

#include "ring.h"
#include "utils/fast-pow.h"

namespace elliptic_curve_guide::polynomial {
    Poly Poly::pow(const Poly& poly, const uint& power) {
        if (power == 0) {
            return Poly(poly.get_field(), {1});
        }

        Poly result = algorithm::fast_pow<Poly>(poly, power);
        result.clean();
        assert(result.is_valid() && "Poly::pow : invalid representation of polynomial");
        return result;
    }

    Poly Poly::compose(const Poly& outside_poly, const Poly& inside_poly) {
        Poly result(outside_poly.get_field(), {outside_poly[0]});

        for (size_t i = 1; i <= outside_poly.degree(); ++i) {
            result += pow(inside_poly, i) * outside_poly[i];
        }

        return result;
    }

    Poly Poly::get_x_power_n(const Field& field, size_t n) {
        Poly result = Poly(field, {1});
        result.increase_degree_by(n);
        return result;
    }

    Poly Poly::increase_degree_by(const Poly& poly, size_t shift) {
        Poly result = poly;
        result.increase_degree_by(shift);
        return result;
    }

    Poly::Poly(const Field& field) : m_field(field), m_coeffs({field.element(0)}) {};

    Poly::Poly(const Field& field, const std::vector<Element>& coeffs) : m_field(field), m_coeffs(coeffs) {
        clean();
        assert(is_valid() && "Poly::pow : invalid representation of polynomial");
    }

    Poly::Poly(const Field& field, std::vector<Element>&& coeffs) :
        m_field(field), m_coeffs(std::move(coeffs)) {
        clean();
        assert(is_valid() && "Poly::pow : invalid representation of polynomial");
    }

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
        m_field(field), m_coeffs(convert_to_field_coeffs(field, coeffs)) {
        clean();
        assert(is_valid() && "Poly::pow : invalid representation of polynomial");
    }

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
        const Poly::Field& F = lhs.get_field();

        size_t degree = lhs.degree() + rhs.degree();
        Poly result(F);
        result.m_coeffs.resize(degree + 1, F.element(0));

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
        assert(other.degree() > 0 && "Poly::operator%= : invalid modulus");

        while (degree() >= other.degree()) {
            Element factor = top_coef() / other.top_coef();
            *this -= increase_degree_by(other * factor, degree() - other.degree());
        }

        clean();
        assert(is_valid() && "Poly::operator%= : invalid representation of polynomial");
        return *this;
    }

    bool operator==(const Poly& lhs, const Poly& rhs) {
        return lhs.m_coeffs == rhs.m_coeffs;
    }

    void Poly::pow(const uint& power) {
        *this = pow(*this, power);
    }

    void Poly::increase_degree_by(size_t shift) {
        std::vector<Element> result(shift, m_field.element(0));
        result.insert(result.end(), m_coeffs.begin(), m_coeffs.end());
        m_coeffs = result;
        clean();
        assert(is_valid() && "Poly::increase_degree_by : invalid representation of polynomial");
    }

    size_t Poly::degree() const {
        return m_coeffs.size() - 1;
    }

    void Poly::compose(const Poly& inside_poly) {
        *this = compose(*this, inside_poly);
        assert(is_valid() && "Poly::compose : invalid representation of polynomial");
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
            if (!m_coeffs[size - 1].is_invertible()) {
                m_coeffs.pop_back();
                --size;
            } else {
                break;
            }
        }

        if (size == 0) {
            m_coeffs = {m_field.element(0)};
        }

        assert(is_valid() && "Poly::clean : invalid representation of polynomial");
    }

    void Poly::negative() {
        for (auto& coef : m_coeffs) {
            coef = -coef;
        }

        clean();
        assert(is_valid() && "Poly::negative : invalid representation of polynomial");
    }

    bool Poly::is_valid() const {
        size_t size = len();

        if (size == 0) {
            return false;
        }

        if (size == 1) {
            return true;
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
    namespace {
        struct DivisionResult {
            polynomial::Poly quotient;
            polynomial::Poly remainder;
        };
    }   // namespace

    static DivisionResult divide(const polynomial::Poly& lhs, const polynomial::Poly& rhs) {
        using Element = field::FieldElement;
        using polynomial::Poly;
        const field::Field& F = lhs.get_field();

        const size_t lhs_degree = lhs.degree();
        const size_t rhs_degree = rhs.degree();

        if (rhs_degree == 0) {
            return DivisionResult {.quotient = lhs * Element::inverse(rhs[0]), .remainder = Poly(F)};
        }

        if (lhs_degree < rhs_degree) {
            return DivisionResult {.quotient = Poly(F), .remainder = lhs};
        }

        Poly quotient(F, {1});
        Poly remainder = lhs;

        size_t degree_delta = remainder.degree() - rhs_degree;
        quotient.increase_degree_by(degree_delta);
        Element k = remainder.top_coef() / rhs.top_coef();
        remainder -= k * Poly::increase_degree_by(rhs, degree_delta);
        quotient[degree_delta] = k;

        while (remainder.degree() >= rhs_degree) {
            degree_delta = remainder.degree() - rhs_degree;
            k = remainder.top_coef() / rhs.top_coef();
            remainder -= k * Poly::increase_degree_by(rhs, degree_delta);
            quotient[degree_delta] = k;
        }

        return DivisionResult {.quotient = quotient, .remainder = remainder};
    }

    static polynomial::Poly polynomial_extended_modular_gcd(const polynomial::Poly& a,
                                                            const polynomial::Poly& b,
                                                            polynomial::Poly& x,
                                                            polynomial::Poly& y,
                                                            const polynomial::Poly& modulus) {
        using polynomial::Poly;
        const field::Field& F = a.get_field();

        if (b.degree() == 0 && !b[0].is_invertible()) {
            x = Poly(F, {1});
            y = Poly(F);
            return a;
        }

        DivisionResult division_result = divide(a, b);
        Poly x1(F), y1(F);
        Poly d = polynomial_extended_modular_gcd(b, division_result.remainder, x1, y1, modulus);
        x = y1;
        y = x1 - y1 * division_result.quotient;
        return d;
    }

    ModulusGcdResult modulus_gcd(const polynomial::Poly& value, const polynomial::Poly& modulus) {
        const field::Field& F = value.get_field();
        polynomial::Poly x(F), y(F);
        polynomial::Poly d = polynomial_extended_modular_gcd(value, modulus, x, y, modulus);

        if (d.degree() == 0) {
            const field::FieldElement inverse_value = field::FieldElement::inverse(d.top_coef());
            x *= inverse_value;
            y *= inverse_value;
            d = polynomial::Poly(F, {1});
        }

        return ModulusGcdResult {
            .gcd = d,
            .value_multiplier = x,
            .modulus_multiplier = y,
        };
    }

    polynomial::Poly gcd(const polynomial::Poly& lhs, const polynomial::Poly& rhs) {
        if (rhs.degree() == 0 && !rhs[0].is_invertible()) {
            return lhs;
        }

        DivisionResult division_result = divide(lhs, rhs);
        return gcd(rhs, division_result.remainder);
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
