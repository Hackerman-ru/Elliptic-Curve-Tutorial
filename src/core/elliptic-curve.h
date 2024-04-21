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
        std::optional<EllipticCurvePoint<type>> point_with_x_equal_to(const FieldElement& x) const {
            std::optional<FieldElement> y = find_y(x);

            if (!y.has_value()) {
                return std::nullopt;
            }

            return EllipticCurvePoint<type>(x, std::move(y), m_a, m_b, m_F);
        }

        template<CoordinatesType type = CoordinatesType::Normal>
        std::optional<EllipticCurvePoint<type>> point_with_x_equal_to(FieldElement&& x) const {
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

        EllipticCurvePoint null_point() {
            return EllipticCurvePoint(m_F->element(0), m_F->element(1), m_a, m_b, m_F);
        }

        EllipticCurvePoint(FieldElement x, FieldElement y, std::shared_ptr<const FieldElement> a,
                           std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F) :
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F)),
            m_x {std::move(x)},
            m_y {std::move(y)} {};

    public:
        friend EllipticCurvePoint operator+(const EllipticCurvePoint& lhs, const EllipticCurvePoint& rhs) {
            EllipticCurvePoint result = lhs;
            result += rhs;
            return result;
        }

        friend EllipticCurvePoint operator+(EllipticCurvePoint&& lhs, const EllipticCurvePoint& rhs) {
            lhs += rhs;
            return lhs;
        }

        friend EllipticCurvePoint operator+(const EllipticCurvePoint& lhs, EllipticCurvePoint&& rhs) {
            rhs += lhs;
            return rhs;
        }

        friend EllipticCurvePoint operator+(EllipticCurvePoint&& lhs, EllipticCurvePoint&& rhs) {
            lhs += rhs;
            return lhs;
        }

        friend EllipticCurvePoint operator-(const EllipticCurvePoint& lhs, const EllipticCurvePoint& rhs) {
            EllipticCurvePoint result = lhs;
            result -= rhs;
            return result;
        }

        friend EllipticCurvePoint operator-(EllipticCurvePoint&& lhs, const EllipticCurvePoint& rhs) {
            lhs -= rhs;
            return lhs;
        }

        friend EllipticCurvePoint operator-(const EllipticCurvePoint& lhs, EllipticCurvePoint&& rhs) {
            rhs -= lhs;
            rhs.negative();
            return rhs;
        }

        friend EllipticCurvePoint operator-(EllipticCurvePoint&& lhs, EllipticCurvePoint&& rhs) {
            lhs -= rhs;
            return lhs;
        }

        friend EllipticCurvePoint operator*(const EllipticCurvePoint& point, const uint& value) {
            EllipticCurvePoint result = point;
            result *= value;
            return result;
        }

        friend EllipticCurvePoint operator*(EllipticCurvePoint&& point, const uint& value) {
            point *= value;
            return point;
        }

        friend EllipticCurvePoint operator*(const uint& value, const EllipticCurvePoint& point) {
            EllipticCurvePoint result = point;
            result *= value;
            return result;
        }

        friend EllipticCurvePoint operator*(const uint& value, EllipticCurvePoint&& point) {
            point *= value;
            return point;
        }

        EllipticCurvePoint operator-() const {
            EllipticCurvePoint result = *this;
            result.negative();
            return result;
        }

        EllipticCurvePoint& operator+=(const EllipticCurvePoint& other) {
            if (m_is_null) {
                return *this = other;
            } else if (other.m_is_null) {
                return *this;
            }

            FieldElement k = m_F->element(0);

            if (m_x == other.m_x) {
                if (m_y == -other.m_y) {
                    m_x = m_F->element(0);
                    m_y = m_F->element(1);
                    m_is_null = true;
                    return *this;
                }

                k = (m_F->element(3) * m_x.pow(2) + *m_a) / (m_y << 1);
            } else {
                k = (other.m_y - m_y) / (other.m_x - m_x);
            }

            FieldElement x = k.pow(2) - m_x - other.m_x;
            m_y = k * (m_x - x) - m_y;
            m_x = x;

            return *this;
        }

        EllipticCurvePoint& operator-=(const EllipticCurvePoint& other) {
            EllipticCurvePoint temp = other;
            temp.negative();
            return *this += temp;
        }

        EllipticCurvePoint& operator-=(EllipticCurvePoint&& other) {
            other.negative();
            return *this += other;
        }

        EllipticCurvePoint& operator*=(const uint& value);

        friend bool operator==(const EllipticCurvePoint& lhs, const EllipticCurvePoint& rhs) {
            return (lhs.m_is_null && rhs.m_is_null) || (lhs.m_x == rhs.m_x && lhs.m_y == rhs.m_y);
        }

    private:
        void negative() {
            m_y = -m_y;
        }

        FieldElement m_x;
        FieldElement m_y;

        static EllipticCurvePoint two_p_plus_q(const EllipticCurvePoint& P, const EllipticCurvePoint& Q) {
            assert(P != Q && "EllipticCurvePoint<CoordinatesType::Normal>::two_p_plus_q : wrong arguments");
        }
    };
}   // namespace ECG

#endif
