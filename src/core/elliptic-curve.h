#ifndef ECG_ELLIPTIC_CURVE_H
#define ECG_ELLIPTIC_CURVE_H

#include "Field.h"

namespace ECG {
    enum class CoordinatesType {
        Normal,
        Projective,
        Jacobi,
        ModifiedJacobi,
        JacobiChudnovski,
        SimplifiedJacobiChudnovski,
    };

    namespace {
        class EllipticCurvePointConcept {
        public:
            EllipticCurvePointConcept(std::shared_ptr<const FieldElement> a,
                                      std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F) :
                m_a {std::move(a)}, m_b {std::move(b)}, m_F {std::move(F)} {};

        protected:
            const std::shared_ptr<const FieldElement> m_a;
            const std::shared_ptr<const FieldElement> m_b;
            const std::shared_ptr<const Field> m_F;
        };
    }   // namespace

    template<CoordinatesType type = CoordinatesType::Normal>
    class EllipticCurvePoint;

    class EllipticCurve {
    public:
        EllipticCurve(FieldElement a, FieldElement b, Field F);
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
        std::shared_ptr<const FieldElement> m_a;
        std::shared_ptr<const FieldElement> m_b;
        std::shared_ptr<const Field> m_F;

        bool find_y(const FieldElement& x, FieldElement& y) const;
    };

    template<>
    class EllipticCurvePoint<CoordinatesType::Normal> : private EllipticCurvePointConcept {
    public:
        EllipticCurvePoint(FieldElement x, FieldElement y, std::shared_ptr<const FieldElement> a,
                           std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F) :
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F)),
            m_x {std::move(x)},
            m_y {std::move(y)} {};

    private:
        FieldElement m_x;
        FieldElement m_y;
    };

    uint random_uint();
}   // namespace ECG

#endif
