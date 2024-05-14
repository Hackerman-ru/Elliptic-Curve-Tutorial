#include "bitsize.h"

namespace elliptic_curve_guide::algorithm {
    size_t actual_bit_size(const uint& value) {
        size_t l = 1;
        size_t r = uint_info::uint_bits_number + 1;

        while (r - l > 1) {
            size_t m = (l + r) / 2;
            uint temp = 1;
            temp <<= m;

            if (temp > value) {
                r = m;
            } else {
                l = m;
            }
        }

        return l;
    }
}   // namespace elliptic_curve_guide::algorithm
