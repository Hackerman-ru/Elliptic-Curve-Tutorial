#ifndef ECG_EL_GAMAL_H
#define ECG_EL_GAMAL_H

#include "../EllipticCurve/elliptic-curve.h"

namespace ECG {

    class ElGamal {   // TODO
    public:
        ElGamal() = delete;
        ElGamal(const EllipticCurve& elliptic_curve, const EllipticCurvePoint& public_key,
                const EllipticCurvePoint& generator) :
            m_elliptic_curve(elliptic_curve), m_public_key(public_key), m_generator(generator) {};

        std::pair<std::string, std::string> encrypt(const std::string& message, CoordinatesType type) const;
        std::pair<std::string, std::string> decrypt(const std::string& a, const std::string& b,
                                                    CoordinatesType type) const;

    private:
        const EllipticCurve m_elliptic_curve;
        const EllipticCurvePoint& m_public_key;
        const EllipticCurvePoint& m_generator;
    };
}   // namespace ECG

#endif
