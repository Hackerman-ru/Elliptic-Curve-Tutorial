#ifndef ECG_ECDSA_H
#define ECG_ECDSA_H

#include "../EllipticCurve/elliptic-curve.h"

namespace ECG {
    class ECDSA {
    public:
        ECDSA(const uint& field_order, const uint& security_level) {
            Field F(field_order);
            FieldElement a = F(0);
            FieldElement b = F(0);

            do {
                a = generate_random_element(F);
                b = generate_random_element(F);
            } while (!((F(4) * a.fast_pow(3) + F(27) * b.fast_pow(2)).is_inversible()));

            EllipticCurve E(a, b);
            uint N = E.points_number();
            // TODO
        }

    private:
        FieldElement generate_random_element(const Field& F) const;
    };
}   // namespace ECG
#endif
