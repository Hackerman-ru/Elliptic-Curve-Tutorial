#ifndef ECG_FAST_POW_H
#define ECG_FAST_POW_H

#include "uint.h"

namespace elliptic_curve_guide {
    namespace algorithm {
        template<typename T>
        T fast_pow(const T& value, const size_t& power) {
            if ((power & 1) != 0) {
                if (power == 1) {
                    return value;
                }

                return value * fast_pow<T>(value, power - 1);
            }

            T temp = fast_pow<T>(value, power >> 1);
            return temp * temp;
        }

        template<typename T>
        T fast_pow(const T& value, const uint& power) {
            if ((power & 1) != 0) {
                if (power == 1) {
                    return value;
                }

                return value * fast_pow<T>(value, power - 1);
            }

            T temp = fast_pow<T>(value, power >> 1);
            return temp * temp;
        }
    }   // namespace algorithm
}   // namespace elliptic_curve_guide
#endif
