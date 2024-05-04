#include "bitsize.h"

namespace ECG {
    size_t actual_bit_size(const uint& value) {
        size_t l = 1;
        size_t r = 513;

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
}   // namespace ECG
