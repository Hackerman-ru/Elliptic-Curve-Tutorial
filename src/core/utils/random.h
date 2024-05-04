#ifndef ECG_RANDOM_H
#define ECG_RANDOM_H

#include "uint.h"

namespace ECG {
    uint generate_random_uint();
    uint generate_random_uint_modulo(const uint& modulus);
}   // namespace ECG
#endif
