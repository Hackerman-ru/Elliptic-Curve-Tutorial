#ifndef ECG_ELLIPTIC_CURVE_H
#define ECG_ELLIPTIC_CURVE_H

#include "field.h"

#include <optional>

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
            std::shared_ptr<const FieldElement> m_a;
            std::shared_ptr<const FieldElement> m_b;
            std::shared_ptr<const Field> m_F;
        };
    }   // namespace

    template<CoordinatesType type = CoordinatesType::Normal>
    class EllipticCurvePoint;

    class EllipticCurve {
    public:
        EllipticCurve(FieldElement a, FieldElement b, Field F);
        uint points_number() const;   // SEA algorithm

        template<CoordinatesType type = CoordinatesType::Normal>
        EllipticCurvePoint<type> operator()(FieldElement x) const {}

        template<CoordinatesType type = CoordinatesType::Normal>
        EllipticCurvePoint<type> operator()(FieldElement x, FieldElement y) const {
            return EllipticCurvePoint<type>(x, y);
        }

        template<CoordinatesType type = CoordinatesType::Normal>
        EllipticCurvePoint<type> random_point() const;

    private:
        std::optional<FieldElement> find_y(const FieldElement& x) const;

        struct Cache {
            uint degree;   // (p - 1) / 2
            size_t power_of_two;
            uint residue;                                // p - 1 = 2.pow(power_of_two) * residue
            FieldElement b;                              // quadratic nonresidue
            FieldElement inverse_two;                    // 1/2 mod p
            std::vector<FieldElement> b_second_powers;   // b.pow(2) ... b.pow(2.pow(power_of_two - 1))
            std::vector<FieldElement>
                b_second_u_powers;   // b.pow(residue), b.pow(2 * residue), ... , b.pow(2.pow(power_of_two - 1) * residue)
        };

        mutable std::optional<Cache> m_cache;
        std::shared_ptr<const FieldElement> m_a;
        std::shared_ptr<const FieldElement> m_b;
        std::shared_ptr<const Field> m_F;
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
}   // namespace ECG

#endif
