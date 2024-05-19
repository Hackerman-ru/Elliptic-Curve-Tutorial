#include "modulo_inversion.h"

namespace elliptic_curve_guide::algorithm {
    template<typename T>
    static T extended_modular_gcd(const T& a, const T& b, T& x, T& y, const T& modulus) {
        if (b == 0) {
            x = 1;
            y = 0;
            return a;
        }

        T x1, y1;
        T d = extended_modular_gcd(b, a % b, x1, y1, modulus);
        x = y1;
        T temp = y1 * (a / b);

        while (x1 < temp) {
            x1 += modulus;
        }

        y = x1 - temp;
        return d;
    }

    uint inverse_modulo(const uint& value, const uint& modulus) {
        uint result, temp;
        extended_modular_gcd<uint>(value, modulus, result, temp, modulus);
        assert(result < modulus && "inverse_modulo : value must be less than modulus");
        assert((result * value) % modulus == 1 && "inverse_modulo : incorrect answer");
        return result;
    }

    uint64_t inverse_modulo(const uint64_t& value, const uint64_t& modulus) {
        uint64_t result, temp;
        extended_modular_gcd<uint64_t>(value, modulus, result, temp, modulus);
        assert(result < modulus && "inverse_modulo : value must be less than modulus");
        assert((result * value) % modulus == 1 && "inverse_modulo : incorrect answer");
        return result;
    }
}   // namespace elliptic_curve_guide::algorithm
