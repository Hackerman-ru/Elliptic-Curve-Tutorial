#include "elliptic-curve.h"

#include "map"

namespace ECG {
    EllipticCurve::EllipticCurve(FieldElement a, FieldElement b, Field F) :
        m_a {std::make_shared<const FieldElement>(std::move(a))},
        m_b {std::make_shared<const FieldElement>(std::move(b))},
        m_F {std::make_shared<const Field>(std::move(F))} {}

    uint EllipticCurve::points_number() const {   // TODO
        return uint();
    }

    static std::map<uint, FieldElement> b_cache;

    static FieldElement find_b(const FieldElement& one, const uint& degree) {
        if (b_cache.contains(degree)) {
            return b_cache.at(degree);
        }

        FieldElement b = one + one;

        while (b.pow(degree) == one) {
            b += one;
        }

        b_cache.insert({degree, b});
        return b;
    }

    struct Decomposition {
        uint power_of_two;
        uint residue;
    };

    static Decomposition decompose(const uint& value) {
        uint power_of_two = 0;
        uint residue = value;

        while (residue != 0 && (residue & 0b1) == 0) {
            ++power_of_two;
            residue >>= 1;
        }

        return {power_of_two, residue};
    }

    static uint order_of_two(uint value) {
        uint power_of_two = 0;

        while (value != 0 && (value & 0b1) == 0) {
            ++power_of_two;
        }

        return power_of_two;
    }

    std::optional<FieldElement> EllipticCurve::find_y(const FieldElement& x) const {
        FieldElement value = x.pow(3) + *m_a * x + *m_b;

        if (!value.is_invertible()) {
            return std::nullopt;
        }

        const uint& p = m_F->modulus();
        const FieldElement one = m_F->element(1);

        if (!m_cache) {
            uint degree = (p - 1) >> 1;
            Decomposition decomposition = decompose(p - 1);
            FieldElement b = find_b(one, degree);
            size_t e = decomposition.residue.convert_to<size_t>();

            std::vector<FieldElement> second_powers = {b << 1};
            second_powers.reserve(e - 1);

            for (size_t i = 1; i < e - 1; ++i) {
                second_powers.emplace_back(second_powers[i - 1] << 1);
            }

            std::vector<FieldElement> second_u_powers = {b.pow(decomposition.residue)};
            second_u_powers.reserve(e);

            for (size_t i = 1; i < e; ++i) {
                second_u_powers.emplace_back(second_u_powers[i - 1] << 1);
            }

            Cache cache = {.degree = std::move(degree),
                           .power_of_two = std::move(decomposition.power_of_two),
                           .residue = std::move(decomposition.residue),
                           .b = std::move(b),
                           .b_second_powers = std::move(second_powers),
                           .b_second_u_powers = std::move(second_u_powers)};

            m_cache.emplace(std::move(cache));
        }

        const Cache& cache = m_cache.value();

        if (value.pow(cache.degree) != one) {
            return std::nullopt;
        }

        if ((p & 0b11) == 3) {
            return value.pow((cache.degree + 1) >> 1);
        }
    }

}   // namespace ECG
