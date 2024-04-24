#ifndef ECG_ECDSA_H
#define ECG_ECDSA_H

#include "elliptic-curve.h"
#include "utils/random.h"

namespace ECG {
    namespace {
        struct Parameters {
            Field m_F;
            EllipticCurve m_E;
            EllipticCurvePoint<> m_generator;
            uint m_n;
            uint m_order;
        };
    }   // namespace

    struct Keys {
        EllipticCurvePoint<> public_key;
        FieldElement private_key;
    };

    struct Signature {
        FieldElement r;
        FieldElement s;
    };

    class ECDSA {
    public:
        static std::optional<ECDSA> generate(const uint& field_order, const uint& security_level);

        ECDSA(const Field& F, const EllipticCurve& E, const EllipticCurvePoint<>& generator, uint n, uint h) :
            m_parameters {F, E, generator, n, h} {}

        Keys generate_keys() const;
        Signature generate_signature(const FieldElement& private_key, const uint& message) const;
        bool is_correct_signature(const EllipticCurvePoint<>& public_key, const uint& message,
                                  const Signature& signature) const;

    private:
        Parameters m_parameters;
    };
}   // namespace ECG
#endif
