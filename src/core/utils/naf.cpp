#include "naf.h"

namespace elliptic_curve_guide::algorithm::non_adjacent_form {
    static constexpr uint16_t c_mask_modulo_2_pow_w = (1 << c_width) - 1;

    wnaf_form get_wnaf(uint value) {
        wnaf_form result;

        while (value > 0) {
            if ((value & 0b1) == 1) {
                uint16_t coef_value = value.convert_to<uint16_t>() & c_mask_modulo_2_pow_w;

                if (coef_value >= (1 << (c_width - 1))) {
                    coef_value = (1 << c_width) - coef_value;
                    result.push_back({.value = coef_value, .is_negative = true});
                    value += coef_value;
                } else {
                    result.push_back({.value = coef_value, .is_negative = false});
                    value -= coef_value;
                }
            } else {
                result.push_back({.value = 0, .is_negative = false});
            }

            value >>= 1;
        }

        return result;
    }
}   // namespace elliptic_curve_guide::algorithm::non_adjacent_form
