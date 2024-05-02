#ifndef ECG_ELLIPTIC_CURVE_H
#define ECG_ELLIPTIC_CURVE_H

#include "field.h"
#include "utils/naf.h"

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

        uint points_number() const;   // Scoof's algorithm

        template<CoordinatesType type = CoordinatesType::Normal>
        std::optional<EllipticCurvePoint<type>> point_with_x_equal_to(const FieldElement& x) const {
            std::optional<FieldElement> y = find_y(x);

            if (!y.has_value()) {
                return std::nullopt;
            }

            return EllipticCurvePoint<type>(x, std::move(y.value()), m_a, m_b, m_F);
        }

        template<CoordinatesType type = CoordinatesType::Normal>
        std::optional<EllipticCurvePoint<type>> point_with_x_equal_to(FieldElement&& x) const {
            std::optional<FieldElement> y = find_y(x);

            if (!y.has_value()) {
                return std::nullopt;
            }

            return EllipticCurvePoint<type>(std::move(x), std::move(y.value()), m_a, m_b, m_F);
        }

        template<CoordinatesType type = CoordinatesType::Normal>
        std::optional<EllipticCurvePoint<type>> point(const FieldElement& x, const FieldElement& y) const {
            EllipticCurvePoint<type> point(x, y, m_a, m_b, m_F, true);
            return point.is_valid() ? point : std::nullopt;
        }

        template<CoordinatesType type = CoordinatesType::Normal>
        std::optional<EllipticCurvePoint<type>> point(FieldElement&& x, const FieldElement& y) const {
            EllipticCurvePoint<type> point(std::move(x), y, m_a, m_b, m_F, true);
            return point.is_valid() ? point : std::nullopt;
        }

        template<CoordinatesType type = CoordinatesType::Normal>
        std::optional<EllipticCurvePoint<type>> point(const FieldElement& x, FieldElement&& y) const {
            EllipticCurvePoint<type> point(x, std::move(y), m_a, m_b, m_F, true);
            return point.is_valid() ? point : std::nullopt;
        }

        template<CoordinatesType type = CoordinatesType::Normal>
        std::optional<EllipticCurvePoint<type>> point(FieldElement&& x, FieldElement&& y) const {
            EllipticCurvePoint<type> point(std::move(x), std::move(y), m_a, m_b, m_F, true);
            return point.is_valid() ? point : std::nullopt;
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

        friend bool operator==(const EllipticCurvePoint& lhs, const EllipticCurvePoint& rhs) {
            return (lhs.m_is_null && rhs.m_is_null) || (lhs.m_x == rhs.m_x && lhs.m_y == rhs.m_y);
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

                k = (m_F->element(3) * FieldElement::pow(m_x, 2) + *m_a) / (m_y << 1);
            } else {
                k = (other.m_y - m_y) / (other.m_x - m_x);
            }

            FieldElement x = FieldElement::pow(k, 2) - m_x - other.m_x;
            m_y = k * (m_x - x) - m_y;
            m_x = x;

            assert(is_valid()
                   && "EllipticCurvePoint<CoordinatesType::Normal>::operator+= : invalid coordinates");
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

        EllipticCurvePoint& operator*=(const uint& value) {
            NAF::wnaf_form wnaf_form = NAF::get_wnaf(value);
            EllipticCurvePoint two_p = *this + *this;
            std::vector<EllipticCurvePoint> kp = {*this};

            for (size_t i = 1; i < NAF::c_kp_number; ++i) {
                kp.emplace_back(kp.back() + two_p);
            }

            m_is_null = true;

            for (size_t i = wnaf_form.size(); i > 0; --i) {
                *this += *this;

                if (wnaf_form[i - 1].value != 0) {
                    if (!wnaf_form[i - 1].is_negative) {
                        *this += kp[wnaf_form[i - 1].value >> 1];
                    } else {
                        *this -= kp[wnaf_form[i - 1].value >> 1];
                    }
                }
            }

            return *this;
        }

        const FieldElement& get_x() const {
            return m_x;
        }

        const FieldElement& get_y() const {
            return m_y;
        }

        bool is_zero() const {
            return m_is_null;
        }

    private:
        static EllipticCurvePoint null_point(std::shared_ptr<const FieldElement> a,
                                             std::shared_ptr<const FieldElement>
                                                 b,
                                             std::shared_ptr<const Field>
                                                 F) {
            return EllipticCurvePoint(F->element(0), F->element(1), a, b, F, true);
        }

        EllipticCurvePoint(const FieldElement& x, const FieldElement& y,
                           std::shared_ptr<const FieldElement> a, std::shared_ptr<const FieldElement> b,
                           std::shared_ptr<const Field> F, bool is_null = false) :
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null), m_x {x}, m_y {y} {
            assert(
                is_valid()
                && "EllipticCurvePoint<CoordinatesType::Normal>::EllipticCurvePoint : invalid coordinates");
        }

        EllipticCurvePoint(FieldElement&& x, FieldElement&& y, std::shared_ptr<const FieldElement> a,
                           std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F,
                           bool is_null = false) :
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
            m_x {std::move(x)},
            m_y {std::move(y)} {
            assert(
                is_valid()
                && "EllipticCurvePoint<CoordinatesType::Normal>::EllipticCurvePoint : invalid coordinates");
        }

        EllipticCurvePoint null_point() const {
            return EllipticCurvePoint(m_F->element(0), m_F->element(1), m_a, m_b, m_F, true);
        }

        void negative() {
            m_y = -m_y;
        }

        bool is_valid() const {
            if (m_is_null) {
                return true;
            }

            FieldElement value = FieldElement::pow(m_x, 3) + *m_a * m_x + *m_b;
            return FieldElement::pow(m_y, 2) == value;
        }

        FieldElement m_x;
        FieldElement m_y;
    };
}   // namespace ECG

#endif
