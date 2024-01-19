#ifndef ECG_ELLIPTIC_CURVE_H
#define ECG_ELLIPTIC_CURVE_H

#include "../Field/Field.h"

namespace ECG {

    class EllipticCurve {
    public:
        constexpr EllipticCurve(const FieldElement& a, const FieldElement& b);
        uint find_points_number() const;   // SEA algorithm

        const FieldElement& get_a() const;
        const FieldElement& get_b() const;

    private:
        const FieldElement m_a;
        const FieldElement m_b;

        FieldElement find_y(const FieldElement& x) const;
        FieldElement find_n(const FieldElement& x) const;
        static uint generate_random_uint();
    };

    class EllipticCurvePoint {
        enum class CoordinatesType {
            Normal,
            Projective,
            Jacobi,
            ModifiedJacobi,
            JacobiChudnovski,
            SimplifiedJacobiChudnovski,
        };
        friend class EllipticCurve;
        class BasePoint;

        EllipticCurvePoint(std::unique_ptr<BasePoint>&& ptr);
        EllipticCurvePoint(FieldElement x, FieldElement y, CoordinatesType type);

    public:
        EllipticCurvePoint(const EllipticCurvePoint& other);
        EllipticCurvePoint(EllipticCurvePoint&& other) noexcept = default;

        EllipticCurvePoint& operator+=(const EllipticCurvePoint& other) {
            m_self->operator+=(*other.m_self);
            return *this;
        }

        EllipticCurvePoint& operator-=(const EllipticCurvePoint& other) {
            m_self->operator-=(*other.m_self);
            return *this;
        }

        EllipticCurvePoint& operator*=(const uint& other) {
            // TODO
            return *this;
        }

        EllipticCurvePoint operator-() const {
            return {m_self->operator-()};
        }

    private:
        class BasePoint {
        public:
            virtual ~BasePoint() = default;

            virtual void operator+=(const BasePoint& other) const = 0;
            virtual void operator-=(const BasePoint& other) const = 0;
            virtual void operator*=(const uint& other) const = 0;

            virtual std::unique_ptr<BasePoint> operator-() const = 0;

            virtual bool operator==(const BasePoint& other) const = 0;
        };

        class NormalPoint final : public BasePoint {};

        class ProjectivePoint final : public BasePoint {};

        class JacobiPoint final : public BasePoint {};

        class ModifiedJacobiPoint final : public BasePoint {};

        class JacobiChudnovskyPoint final : public BasePoint {};

        class SimplifiedJacobiChudnovskyPoint final : public BasePoint {};

        std::unique_ptr<BasePoint> m_self;
    };
}   // namespace ECG

#endif
