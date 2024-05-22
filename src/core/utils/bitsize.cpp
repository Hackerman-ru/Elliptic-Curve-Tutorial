#include "bitsize.h"

namespace elliptic_curve_guide::algorithm {
    size_t actual_bit_size(uint value) {
        size_t result = 0;

        while (value != 0) {
            ++result;
            value >>= 1;
        }

        return result;
    }
}   // namespace elliptic_curve_guide::algorithm
