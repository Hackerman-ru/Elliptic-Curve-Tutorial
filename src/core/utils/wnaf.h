#ifndef ECG_WNAF_H
#define ECG_WNAF_H

#include "uint.h"

namespace elliptic_curve_guide {
    namespace algorithm {
        namespace {
            constexpr size_t c_width = 3;

            struct Coefficient {
                uint16_t value;
                bool is_negative;
            };

            using WnafForm = std::vector<Coefficient>;
        }   // namespace

        static constexpr uint16_t c_mask_modulo_2_pow_w = (1 << c_width) - 1;

        static WnafForm get_wnaf(uint value) {
            WnafForm result;

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

        static constexpr size_t c_k_number = static_cast<size_t>(1) << (c_width - 2);

        template<typename T>
        T wnaf_addition(T value, const uint& n) {
            WnafForm wnaf_form = get_wnaf(n);
            T two_value = value + value;
            std::vector<T> k_values = {value};

            for (size_t i = 1; i < c_k_number; ++i) {
                k_values.emplace_back(k_values.back() + two_value);
            }

            value.nullify();

            for (size_t i = wnaf_form.size(); i > 0; --i) {
                value.twice();

                if (wnaf_form[i - 1].value != 0) {
                    if (!wnaf_form[i - 1].is_negative) {
                        value += k_values[wnaf_form[i - 1].value >> 1];
                    } else {
                        value -= k_values[wnaf_form[i - 1].value >> 1];
                    }
                }
            }

            return value;
        }
    }   // namespace algorithm
}   // namespace elliptic_curve_guide
#endif
