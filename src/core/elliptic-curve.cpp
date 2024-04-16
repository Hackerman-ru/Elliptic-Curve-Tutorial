#include "elliptic-curve.h"

#include <random>

using ECG::EllipticCurve;
using ECG::uint;

uint ECG::random_uint() {
    static constexpr size_t uint_size = sizeof(uint);
    static constexpr size_t GENERATOR_SIZE = 8;

    std::mt19937_64 generator;
    uint result = 0;

    for (size_t size = 0; size < uint_size; size += GENERATOR_SIZE) {
        result += generator();
        result <<= GENERATOR_SIZE * CHAR_BIT;
    }

    return result;
}

ECG::EllipticCurve::EllipticCurve(FieldElement a, FieldElement b, Field F) :
    m_a {std::make_shared<const FieldElement>(std::move(a))},
    m_b {std::make_shared<const FieldElement>(std::move(b))},
    m_F {std::make_shared<const Field>(std::move(F))} {};

bool ECG::EllipticCurve::find_y(const FieldElement& x, FieldElement& y) const {
    return false;
}
