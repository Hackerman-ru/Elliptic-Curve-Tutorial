#ifndef ECG_ECDSA_H
#define ECG_ECDSA_H

#include "elliptic-curve.h"

namespace elliptic_curve_guide {
    namespace algorithm {
        namespace encryption {
            class ECDSA {
                using Field = field::Field;
                using Element = field::FieldElement;
                using Curve = elliptic_curve::EllipticCurve;

            public:
                static constexpr elliptic_curve::CoordinatesType point_type =
                    elliptic_curve::CoordinatesType::ModifiedJacobi;
                using Point = elliptic_curve::EllipticCurvePoint<point_type>;

                struct Keys {
                    Point public_key;
                    uint private_key;
                };

                struct Signature {
                    uint r;
                    uint s;
                };

                ECDSA(const Field& field, const Curve& elliptic_curve, const Point& generator, const uint& n,
                      const uint& h);

                Keys generate_keys() const;
                Signature generate_signature(const uint& private_key, const uint& message) const;
                bool is_correct_signature(const Point& public_key, const uint& message,
                                          const Signature& signature) const;

            private:
                Field m_field;
                Curve m_elliptic_curve;
                Point m_generator;
                uint m_n;
                uint m_h;
            };
        }   // namespace encryption
    }       // namespace algorithm
}   // namespace elliptic_curve_guide
#endif
