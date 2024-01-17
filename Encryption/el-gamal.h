#ifndef ECG_EL_GAMAL_H
#define ECG_EL_GAMAL_H

#include "../EllipticCurve/elliptic-curve.h"

namespace ECG {
    class ElGamal {
    public:
        using EC = EllipticCurvePoint;

        ElGamal();
        ElGamal(const EC& gen);

        std::pair<EC, EC> encrypt(const std::string& message) const;
        std::string decrypt(const EC::Field& key, const std::pair<EC, EC>& value) const;

    private:
        EC m_gen;
    };
}   // namespace ECG

#endif
