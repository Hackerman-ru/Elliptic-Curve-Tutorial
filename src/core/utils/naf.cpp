#include "naf.h"

namespace ECG {
    NAF::NAF(const uint& positive_bits, const uint& negative_bits) :
        m_positive_bits(positive_bits), m_negative_bits(negative_bits) {}

    bool NAF::empty() const {
        return m_positive_bits == 0 && m_negative_bits == 0;
    }

    NAF NAF::operator>>(const size_t& shift) const {
        NAF result = *this;
        result >>= shift;
        return result;
    }

    NAF& NAF::operator>>=(const size_t& shift) {
        m_positive_bits >>= shift;
        m_negative_bits >>= shift;
        return *this;
    }

    bool NAF::negative_bit() const {
        return (m_negative_bits & 0b1) == 1;
    }

    bool NAF::positive_bit() const {
        return (m_positive_bits & 0b1) == 1;
    }

    static uint reverse(uint value) {
        uint result = 0;

        while (value > 0) {
            result <<= 1;
            result += value & 0b1;
            value >>= 1;
        }

        return result;
    }

    NAF get_naf(const uint& value) {
        uint half_value = value >> 1;
        uint one_and_a_half_value = value + half_value;
        uint c = half_value ^ one_and_a_half_value;
        uint positive_bits = one_and_a_half_value & c;
        uint negative_bits = half_value & c;
        return {reverse(positive_bits), reverse(negative_bits)};
    }

}   // namespace ECG
