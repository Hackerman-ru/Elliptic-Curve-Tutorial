#ifndef ECG_EL_GAMAL_H
#define ECG_EL_GAMAL_H

#include "../EllipticCurve/elliptic-curve.h"

namespace ECG {
    struct EncryptedInfo {
        const EllipticCurvePoint generator_degree;
        const EllipticCurvePoint message_with_noise;
    };

    class ElGamal {   // TODO
    public:
        using Type = EllipticCurvePoint::CoordinatesType;

        ElGamal() = delete;
        ElGamal(const EllipticCurve& elliptic_curve, const EllipticCurvePoint& generator,
                const uint& generator_order) :
            m_elliptic_curve(elliptic_curve), m_generator(generator), m_generator_order(generator_order) {};

        EncryptedInfo encrypt(const std::string& message, Type type) const;
        std::string decrypt(EncryptedInfo info, Type type) const;

    private:
        const EllipticCurve m_elliptic_curve;
        const EllipticCurvePoint m_generator;
        const uint m_generator_order;
    };
}   // namespace ECG

#endif
