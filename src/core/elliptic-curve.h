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
        protected:
            EllipticCurvePointConcept(std::shared_ptr<const FieldElement> a,
                                      std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F,
                                      bool is_null = false) :
                m_a {std::move(a)}, m_b {std::move(b)}, m_F {std::move(F)}, m_is_null(is_null) {};

            std::shared_ptr<const FieldElement> m_a;
            std::shared_ptr<const FieldElement> m_b;
            std::shared_ptr<const Field> m_F;
            bool m_is_null;
        };
    }   // namespace

    template<CoordinatesType type = CoordinatesType::Normal>
    class EllipticCurvePoint;

    class EllipticCurve {
    public:
        EllipticCurve(const FieldElement& a, const FieldElement& b, Field F);
        EllipticCurve(FieldElement&& a, const FieldElement& b, Field F);
        EllipticCurve(const FieldElement& a, FieldElement&& b, Field F);
        EllipticCurve(FieldElement&& a, FieldElement&& b, Field F);

        uint points_number() const;   // SEA algorithm

        template<CoordinatesType type = CoordinatesType::Normal>
        std::optional<EllipticCurvePoint<type>> point_with_x_equal_to(FieldElement x) const {
            std::optional<FieldElement> y = find_y(x);

            if (!y.has_value()) {
                return std::nullopt;
            }

            return EllipticCurvePoint<type>(std::move(x), std::move(y), m_a, m_b, m_F);
        }

        template<CoordinatesType type = CoordinatesType::Normal>
        EllipticCurvePoint<type> null_point() const {
            return EllipticCurvePoint<type>::null_point(m_a, m_b, m_F);
        }

        template<CoordinatesType type = CoordinatesType::Normal>
        EllipticCurvePoint<type> random_point() const;   // TODO

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
    private:
        friend class EllipticCurve;

        static EllipticCurvePoint null_point(std::shared_ptr<const FieldElement> a,
                                             std::shared_ptr<const FieldElement>
                                                 b,
                                             std::shared_ptr<const Field>
                                                 F) {
            return EllipticCurvePoint(F->element(0), F->element(1), a, b, F);
        }

        EllipticCurvePoint(FieldElement x, FieldElement y, std::shared_ptr<const FieldElement> a,
                           std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F) :
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F)),
            m_x {std::move(x)},
            m_y {std::move(y)} {};

    public:
        friend EllipticCurvePoint operator+(const EllipticCurvePoint& lhs, const EllipticCurvePoint& rhs);
        friend EllipticCurvePoint operator+(EllipticCurvePoint&& lhs, const EllipticCurvePoint& rhs);
        friend EllipticCurvePoint operator+(const EllipticCurvePoint& lhs, EllipticCurvePoint&& rhs);
        friend EllipticCurvePoint operator+(EllipticCurvePoint&& lhs, EllipticCurvePoint&& rhs);

        friend EllipticCurvePoint operator+(const EllipticCurvePoint& lhs, const EllipticCurvePoint& rhs);
        friend EllipticCurvePoint operator+(EllipticCurvePoint&& lhs, const EllipticCurvePoint& rhs);
        friend EllipticCurvePoint operator+(const EllipticCurvePoint& lhs, EllipticCurvePoint&& rhs);
        friend EllipticCurvePoint operator+(EllipticCurvePoint&& lhs, EllipticCurvePoint&& rhs);

        friend EllipticCurvePoint operator*(const EllipticCurvePoint& point, const uint& value);
        friend EllipticCurvePoint operator*(EllipticCurvePoint&& point, const uint& value);
        friend EllipticCurvePoint operator*(const uint& value, const EllipticCurvePoint& point);
        friend EllipticCurvePoint operator*(const uint& value, EllipticCurvePoint&& point);

        EllipticCurvePoint& operator+=(const EllipticCurvePoint& other);
        EllipticCurvePoint& operator-=(const EllipticCurvePoint& other);
        EllipticCurvePoint& operator*=(const EllipticCurvePoint& other);

    private:
        FieldElement m_x;
        FieldElement m_y;
    };
}   // namespace ECG

#endif
