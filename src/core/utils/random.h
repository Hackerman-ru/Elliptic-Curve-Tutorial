#ifndef ECG_RANDOM_H
#define ECG_RANDOM_H

#include "field.h"
#include "uint.h"

namespace elliptic_curve_guide {
    namespace algorithm {
        namespace random {
            uint generate_random_uint();
            uint generate_random_uint_modulo(const uint& modulus);
            uint generate_random_non_zero_uint_modulo(const uint& modulus);
            field::FieldElement generate_random_field_element(const field::Field& field);
            field::FieldElement generate_random_non_zero_field_element(const field::Field& field);
        }   // namespace random
    }       // namespace algorithm
}   // namespace elliptic_curve_guide
#endif
