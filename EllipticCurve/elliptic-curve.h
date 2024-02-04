#ifndef ECG_ELLIPTIC_CURVE_H
#define ECG_ELLIPTIC_CURVE_H

#include "../Field/Field.h"

namespace ECG {
    template<CoordinatesType type = CoordinatesType::Normal>
    class EllipticCurvePoint;

    class EllipticCurve {
    public:
        EllipticCurve(const FieldElement& a, const FieldElement& b, const Field& F);
        uint points_number() const;   // SEA algorithm

        template<CoordinatesType type = CoordinatesType::Normal>
        EllipticCurvePoint<type> operator()(FieldElement x) const {
            if (!find_y(x, &y)) {
            }
        }

        template<CoordinatesType type = CoordinatesType::Normal>
        EllipticCurvePoint<type> operator()(FieldElement x, FieldElement y) const {
            return EllipticCurvePoint<type>(x, y);
        }

        template<CoordinatesType type = CoordinatesType::Normal>
        EllipticCurvePoint<type> create_random_point() const;

    private:
        const FieldElement m_a;
        const FieldElement m_b;
        const Field m_F;

        bool find_y(const FieldElement& x, FieldElement& y) const;
    };

    template<>
    class EllipticCurvePoint<CoordinatesType::Normal> {
    public:

    private:
        FieldElement x;
        FieldElement y;
    };

    uint random_uint();
}   // namespace ECG

#endif
