#ifndef ECG_ECDSA_H
#define ECG_ECDSA_H

#include "elliptic-curve.h"
#include "utils/random.h"

namespace ECG {
    class ECDSA {
    public:
        struct Keys {
            EllipticCurvePoint<> public_key;
            uint private_key;
        };

        struct Signature {
            uint r;
            uint s;
        };

        static ECDSA generate(const uint& field_order, const uint& security_level);

        ECDSA(const Field& field, const EllipticCurve& elliptic_curve, const EllipticCurvePoint<>& generator,
              const uint& n, const uint& h) :
            m_field(field), m_elliptic_curve(elliptic_curve), m_generator(generator), m_n(n), m_h(h) {}

        Keys generate_keys() const;
        Signature generate_signature(const uint& private_key, const uint& message) const;
        bool is_correct_signature(const EllipticCurvePoint<>& public_key, const uint& message,
                                  const Signature& signature) const;

    private:
        Field m_field;
        EllipticCurve m_elliptic_curve;
        EllipticCurvePoint<> m_generator;
        uint m_n;
        uint m_h;
    };
}   // namespace ECG
#endif
