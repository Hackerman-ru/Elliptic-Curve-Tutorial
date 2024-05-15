#ifndef ECG_UINT_ALGORITHMS_H
#define ECG_UINT_ALGORITHMS_H

#include "uint.h"

namespace elliptic_curve_guide {
    namespace algorithm {
        namespace non_adjacent_form {
            constexpr size_t c_width = 3;

            struct Coefficient {
                uint16_t value;
                bool is_negative;
            };

            using wnaf_form = std::vector<Coefficient>;

            wnaf_form get_wnaf(uint value);
        }   // namespace non_adjacent_form

        size_t actual_bit_size(uint value);
    }   // namespace algorithm
}   // namespace elliptic_curve_guide
#endif
