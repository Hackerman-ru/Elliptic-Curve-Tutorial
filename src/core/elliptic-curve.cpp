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
            FieldElement b = find_b(one, degree);
            Decomposition decomposition = decompose(p - 1);
            Cache cache = {.degree = degree,
                           .power_of_two = decomposition.power_of_two,
                           .residue = decomposition.residue,
                           .b = b};
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
