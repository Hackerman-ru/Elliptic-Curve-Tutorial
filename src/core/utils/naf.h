#ifndef ECG_NAF_H
#define ECG_NAF_H

#include "uint.h"

namespace ECG {
    namespace NAF {
        constexpr size_t c_width = 3;
        constexpr size_t c_kp_number = 1 << (c_width - 2);
        constexpr uint16_t c_mask_modulo_2_pow_w = (1 << c_width) - 1;

        struct Coefficient {
            uint16_t value;
            bool is_negative;
        };

        using wnaf_form = std::vector<Coefficient>;

        wnaf_form get_wnaf(uint value);
    }   // namespace NAF
}   // namespace ECG
#endif
