#ifndef ECG_MODULO_INVERSION_H
#define ECG_MODULO_INVERSION_H

#include "uint.h"

namespace elliptic_curve_guide {
    namespace algorithm {
        uint inverse_modulo(const uint& value, const uint& modulus);
        uint64_t inverse_modulo(const uint64_t& value, const uint64_t& modulus);
    }   // namespace algorithm
}   // namespace elliptic_curve_guide
#endif
