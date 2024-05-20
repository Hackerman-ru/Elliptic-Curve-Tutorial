#ifndef ECG_MODULO_INVERSION_H
#define ECG_MODULO_INVERSION_H

#include "uint.h"

namespace elliptic_curve_guide {
    namespace algorithm {
        template<typename T>
        static T extended_modular_gcd(const T& a, const T& b, T& x, T& y, const T& modulus) {
            if (b == 0) {
                x = 1;
                y = 0;
                return a;
            }

            T x1, y1;
            T d = extended_modular_gcd<T>(b, a % b, x1, y1, modulus);
            x = y1;
            T temp = y1 * (a / b);

            while (x1 < temp) {
                x1 += modulus;
            }

            y = x1 - temp;
            return d;
        }

        template<typename T>
        T inverse_modulo(const T& value, const T& modulus) {
            T result, temp;
            extended_modular_gcd<T>(value, modulus, result, temp, modulus);
            assert(result < modulus && "inverse_modulo : value must be less than modulus");
            assert((result * value) % modulus == 1 && "inverse_modulo : incorrect answer");
            return result;
        }
    }   // namespace algorithm
}   // namespace elliptic_curve_guide
#endif
