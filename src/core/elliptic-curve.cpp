#include "elliptic-curve.h"

#include "map"

namespace ECG {
    EllipticCurve::EllipticCurve(const FieldElement& a, const FieldElement& b, Field F) :
        m_a {std::make_shared<const FieldElement>(a)},
        m_b {std::make_shared<const FieldElement>(b)},
        m_Field {std::make_shared<const Field>(std::move(F))} {}

    EllipticCurve::EllipticCurve(FieldElement&& a, const FieldElement& b, Field F) :
        m_a {std::make_shared<const FieldElement>(std::move(a))},
        m_b {std::make_shared<const FieldElement>(b)},
        m_Field {std::make_shared<const Field>(std::move(F))} {}

    EllipticCurve::EllipticCurve(const FieldElement& a, FieldElement&& b, Field F) :
        m_a {std::make_shared<const FieldElement>(a)},
        m_b {std::make_shared<const FieldElement>(std::move(b))},
        m_Field {std::make_shared<const Field>(std::move(F))} {}

    EllipticCurve::EllipticCurve(FieldElement&& a, FieldElement&& b, Field F) :
        m_a {std::make_shared<const FieldElement>(std::move(a))},
        m_b {std::make_shared<const FieldElement>(std::move(b))},
        m_Field {std::make_shared<const Field>(std::move(F))} {}

    uint EllipticCurve::points_number() const {   // TODO
        return uint();
    }

    const Field& EllipticCurve::get_field() const {
        return *m_Field;
    }

    bool EllipticCurve::is_zero(const FieldElement& x, const FieldElement& y) {
        return x.value() == 0 && y.value() == 1;
    }

    namespace {
        struct Cache {
            size_t power_of_two;
            uint residue;   // p - 1 = 2.pow(power_of_two) * residue
            std::vector<FieldElement>
                b_second_powers;   // b.pow(2) ... b.pow(2.pow(power_of_two - 1)), where b is a quadratic nonresidue
            std::vector<FieldElement>
                b_second_u_powers;   // b.pow(residue), b.pow(2 * residue), ... , b.pow(2.pow(power_of_two - 1) * residue)
        };

        struct Decomposition {
            size_t power_of_two;
            uint residue;
        };
    }   // namespace

    static std::map<uint, Cache> p_cache;

    static FieldElement find_b(const FieldElement& one, const uint& degree) {
        FieldElement b = one + one;

        while (FieldElement::pow(b, degree) == one) {
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

    std::optional<FieldElement> EllipticCurve::find_y(const FieldElement& x) const {
        FieldElement value = FieldElement::pow(x, 3) + *m_a * x + *m_b;

        if (!value.is_invertible()) {
            return std::nullopt;
        }

        const uint& p = m_Field->modulus();
        const FieldElement one = m_Field->element(1);

        if (FieldElement::pow(value, (p - 1) >> 1) != one) {
            return std::nullopt;
        }

        if ((p & 0b11) == 3) {
            return FieldElement::pow(value, (p + 1) >> 2);
        }

        if (!p_cache.contains(p)) {
            Decomposition decomposition = decompose(p - 1);
            FieldElement b = find_b(one, (p - 1) >> 1);
            size_t e = decomposition.power_of_two;

            std::vector<FieldElement> second_powers = {b};
            second_powers.reserve(e - 1);

            for (size_t i = 1; i < e; ++i) {
                second_powers.emplace_back(FieldElement::pow(second_powers[i - 1], 2));
            }

            std::vector<FieldElement> second_u_powers = {FieldElement::pow(b, decomposition.residue)};
            second_u_powers.reserve(e);

            for (size_t i = 1; i < e; ++i) {
                second_u_powers.emplace_back(FieldElement::pow(second_u_powers[i - 1], 2));
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
        std::vector<FieldElement> z_u_powers_of_2 = {FieldElement::pow(value, cache.residue)};
        size_t current_r = 0;

        while (z_u_powers_of_2.back() != one) {
            z_u_powers_of_2.emplace_back(FieldElement::pow(z_u_powers_of_2.back(), 2));
            ++current_r;
        }

        const size_t& e = cache.power_of_two;
        std::vector<size_t> two_orders = {current_r};
        FieldElement current_z = value;

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

        FieldElement current_x = FieldElement::pow(current_z, (cache.residue + 1) >> 1);
        const size_t n = two_orders.size();

        for (size_t i = 0; i + 1 < n; ++i) {
            current_x /= cache.b_second_powers[e - two_orders[n - i - 2] - 1];
        }

        return current_x;
    }
}   // namespace ECG
