#include "elliptic-curve.h"

#include "map"

namespace elliptic_curve_guide::elliptic_curve {
    EllipticCurve::EllipticCurve(const Element& a, const Element& b, Field F) :
        m_a {std::make_shared<const Element>(a)},
        m_b {std::make_shared<const Element>(b)},
        m_Field {std::make_shared<const Field>(std::move(F))} {}

    EllipticCurve::EllipticCurve(Element&& a, const Element& b, Field F) :
        m_a {std::make_shared<const Element>(std::move(a))},
        m_b {std::make_shared<const Element>(b)},
        m_Field {std::make_shared<const Field>(std::move(F))} {}

    EllipticCurve::EllipticCurve(const Element& a, Element&& b, Field F) :
        m_a {std::make_shared<const Element>(a)},
        m_b {std::make_shared<const Element>(std::move(b))},
        m_Field {std::make_shared<const Field>(std::move(F))} {}

    EllipticCurve::EllipticCurve(Element&& a, Element&& b, Field F) :
        m_a {std::make_shared<const Element>(std::move(a))},
        m_b {std::make_shared<const Element>(std::move(b))},
        m_Field {std::make_shared<const Field>(std::move(F))} {}

    const EllipticCurve::Field& EllipticCurve::get_field() const {
        return *m_Field;
    }

    const EllipticCurve::Element& EllipticCurve::get_a() const {
        return *m_a;
    }

    const EllipticCurve::Element& EllipticCurve::get_b() const {
        return *m_b;
    }

    bool EllipticCurve::is_valid_coordinates(const Element& x, const Element& y) const {
        const Element lhs = Element::pow(y, 2);
        const Element rhs = Element::pow(x, 3) + *m_a * x + *m_b;
        return lhs == rhs;
    }

    bool EllipticCurve::is_null_coordinates(const Element& x, const Element& y) const {
        return x.value() == 0 && y.value() == 1;
    }

    namespace {
        struct Cache {
            size_t power_of_two;
            uint residue;   // p - 1 = 2.pow(power_of_two) * residue
            std::vector<field::FieldElement>
                b_second_powers;   // b.pow(2) ... b.pow(2.pow(power_of_two - 1)), where b is a quadratic nonresidue
            std::vector<field::FieldElement>
                b_second_u_powers;   // b.pow(residue), b.pow(2 * residue), ... , b.pow(2.pow(power_of_two - 1) * residue)
        };

        struct Decomposition {
            size_t power_of_two;
            uint residue;
        };
    }   // namespace

    static std::map<uint, Cache> p_cache;

    static field::FieldElement find_b(const field::FieldElement& one, const uint& degree) {
        field::FieldElement b = one + one;

        while (field::FieldElement::pow(b, degree) == one) {
            b += one;
        }

        return b;
    }

    static Decomposition decompose(const uint& value) {
        size_t power_of_two = 0;
        uint residue = value;

        while (residue != 0 && (residue & 0b1) == 0) {
            ++power_of_two;
            residue >>= 1;
        }

        return {power_of_two, residue};
    }

    std::optional<EllipticCurve::Element> EllipticCurve::find_y(const Element& x) const {
        Element value = Element::pow(x, 3) + *m_a * x + *m_b;

        if (!value.is_invertible()) {
            return std::nullopt;
        }

        const uint& p = m_Field->modulus();
        const Element one = m_Field->element(1);

        if (Element::pow(value, (p - 1) >> 1) != one) {
            return std::nullopt;
        }

        if ((p & 0b11) == 3) {
            return Element::pow(value, (p + 1) >> 2);
        }

        if (!p_cache.contains(p)) {
            Decomposition decomposition = decompose(p - 1);
            Element b = find_b(one, (p - 1) >> 1);
            size_t e = decomposition.power_of_two;

            std::vector<Element> second_powers = {b};
            second_powers.reserve(e - 1);

            for (size_t i = 1; i < e; ++i) {
                second_powers.emplace_back(Element::pow(second_powers[i - 1], 2));
            }

            std::vector<Element> second_u_powers = {Element::pow(b, decomposition.residue)};
            second_u_powers.reserve(e);

            for (size_t i = 1; i < e; ++i) {
                second_u_powers.emplace_back(Element::pow(second_u_powers[i - 1], 2));
            }

            p_cache.insert({
                p,
                {.power_of_two = std::move(decomposition.power_of_two),
                  .residue = std::move(decomposition.residue),
                  .b_second_powers = std::move(second_powers),
                  .b_second_u_powers = std::move(second_u_powers)}
            });
        }

        const Cache& cache = p_cache.at(p);
        std::vector<Element> z_u_powers_of_2 = {Element::pow(value, cache.residue)};
        size_t current_r = 0;

        while (z_u_powers_of_2.back() != one) {
            z_u_powers_of_2.emplace_back(Element::pow(z_u_powers_of_2.back(), 2));
            ++current_r;
        }

        const size_t& e = cache.power_of_two;
        std::vector<size_t> two_orders = {current_r};
        Element current_z = value;

        while (current_r != 0) {
            current_z *= cache.b_second_powers[e - current_r];

            for (size_t next_r = 0; next_r < current_r; ++next_r) {
                z_u_powers_of_2[next_r] *= cache.b_second_u_powers[e - (current_r - next_r)];

                if (z_u_powers_of_2[next_r] == one) {
                    current_r = next_r;
                    break;
                }
            }

            two_orders.emplace_back(current_r);
        }

        Element current_x = Element::pow(current_z, (cache.residue + 1) >> 1);
        const size_t n = two_orders.size();

        for (size_t i = 0; i + 1 < n; ++i) {
            current_x /= cache.b_second_powers[e - two_orders[n - i - 2] - 1];
        }

        return current_x;
    }
}   // namespace elliptic_curve_guide::elliptic_curve
