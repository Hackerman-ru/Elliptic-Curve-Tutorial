#ifndef ECG_EL_GAMAL_H
#define ECG_EL_GAMAL_H

#include "../EllipticCurve/elliptic-curve.h"

namespace ECG {

    template<const Field* field>
    class ElGamal {
        using Curve = EllipticCurve<field>;
        static const Curve m_elliptic_curve;

    public:
        ElGamal() = delete;
        ElGamal(const Curve& ec) : m_elliptic_curve(ec) {};

        std::pair<EllipticCurvePointBasic<field, &m_elliptic_curve>,
                  EllipticCurvePointBasic<field, &m_elliptic_curve>>
            encrypt(const std::string& message) const;
        std::string decrypt(const EllipticCurve::Field& key,
                            const std::pair<EllipticCurve, EllipticCurve>& value) const;

    private:
    };
}   // namespace ECG

#endif
