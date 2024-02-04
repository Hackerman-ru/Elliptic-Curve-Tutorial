#include "elliptic-curve.h"

#include <random>

using ECG::EllipticCurve;
using ECG::uint;

uint ECG::random_uint() {
    static constexpr size_t UINT_SIZE = sizeof(uint);
    static constexpr size_t GENERATOR_SIZE = 8;

    std::mt19937_64 generator;
    uint result = 0;

    for (size_t size = 0; size < UINT_SIZE; size += GENERATOR_SIZE) {
        result += generator();
        result <<= GENERATOR_SIZE * CHAR_BIT;
    }

    return result;
}

ECG::EllipticCurve::EllipticCurve(const FieldElement& a, const FieldElement& b, const Field& F) :
    m_a(a), m_b(b), m_F(F) {};

bool ECG::EllipticCurve::find_y(const FieldElement& x, FieldElement& y) const {
    return false;
}
