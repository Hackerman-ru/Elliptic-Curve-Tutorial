#ifndef ECG_RANDOM_H
#define ECG_RANDOM_H

#include "uint.h"

namespace elliptic_curve_guide {
    namespace algorithm {
        namespace random {
            uint generate_random_uint();
            uint generate_random_uint_modulo(const uint& modulus);
            uint generate_random_non_zero_uint_modulo(const uint& modulus);
        }   // namespace random
    }       // namespace algorithm
}   // namespace elliptic_curve_guide
#endif
