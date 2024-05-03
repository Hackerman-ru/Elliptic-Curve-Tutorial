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
        public:
            virtual FieldElement get_x() const = 0;
            virtual FieldElement get_y() const = 0;

            bool is_zero() const {
                return m_is_null;
            }

        protected:
            EllipticCurvePointConcept(std::shared_ptr<const FieldElement> a,
                                      std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F,
                                      bool is_null = false) :
                m_a {std::move(a)}, m_b {std::move(b)}, m_F {std::move(F)}, m_is_null(is_null) {};

            virtual void negative() = 0;
            virtual void twice() = 0;
            virtual bool is_valid() const = 0;

            std::shared_ptr<const FieldElement> m_a;
            std::shared_ptr<const FieldElement> m_b;
            std::shared_ptr<const Field> m_F;
            bool m_is_null;
        };
    }   // namespace

    template<CoordinatesType type = CoordinatesType::Normal>
    class EllipticCurvePoint;

    template<CoordinatesType type>
    EllipticCurvePoint<type> operator+(const EllipticCurvePoint<type>& lhs,
                                       const EllipticCurvePoint<type>& rhs) {
        EllipticCurvePoint result = lhs;
        result += rhs;
        return result;
    }

    template<CoordinatesType type>
    EllipticCurvePoint<type> operator+(EllipticCurvePoint<type>&& lhs, const EllipticCurvePoint<type>& rhs) {
        lhs += rhs;
        return lhs;
    }

    template<CoordinatesType type>
    EllipticCurvePoint<type> operator+(const EllipticCurvePoint<type>& lhs, EllipticCurvePoint<type>&& rhs) {
        rhs += lhs;
        return rhs;
    }

    template<CoordinatesType type>
    EllipticCurvePoint<type> operator+(EllipticCurvePoint<type>&& lhs, EllipticCurvePoint<type>&& rhs) {
        lhs += rhs;
        return lhs;
    }

    template<CoordinatesType type>
    EllipticCurvePoint<type> operator-(const EllipticCurvePoint<type>& lhs,
                                       const EllipticCurvePoint<type>& rhs) {
        EllipticCurvePoint result = lhs;
        result -= rhs;
        return result;
    }

    template<CoordinatesType type>
    EllipticCurvePoint<type> operator-(EllipticCurvePoint<type>&& lhs, const EllipticCurvePoint<type>& rhs) {
        lhs -= rhs;
        return lhs;
    }

    template<CoordinatesType type>
    EllipticCurvePoint<type> operator-(const EllipticCurvePoint<type>& lhs, EllipticCurvePoint<type>&& rhs) {
        rhs -= lhs;
        rhs.negative();
        return rhs;
    }

    template<CoordinatesType type>
    EllipticCurvePoint<type> operator-(EllipticCurvePoint<type>&& lhs, EllipticCurvePoint<type>&& rhs) {
        lhs -= rhs;
        return lhs;
    }

    template<CoordinatesType type>
    static void multiply(EllipticCurvePoint<type>& point, const uint& value) {
        NAF::wnaf_form wnaf_form = NAF::get_wnaf(value);
        EllipticCurvePoint<type> two_p = point + point;
        std::vector<EllipticCurvePoint<type>> kp = {point};

        for (size_t i = 1; i < NAF::c_kp_number; ++i) {
            kp.emplace_back(kp.back() + two_p);
        }

        point.m_is_null = true;

        for (size_t i = wnaf_form.size(); i > 0; --i) {
            point.twice();

            if (wnaf_form[i - 1].value != 0) {
                if (!wnaf_form[i - 1].is_negative) {
                    point += kp[wnaf_form[i - 1].value >> 1];
                } else {
                    point -= kp[wnaf_form[i - 1].value >> 1];
                }
            }
        }
    }

    template<CoordinatesType type>
    EllipticCurvePoint<type> operator*(const EllipticCurvePoint<type>& point, const uint& value) {
        EllipticCurvePoint<type> result = point;
        result *= value;
        return result;
    }

    template<CoordinatesType type>
    EllipticCurvePoint<type> operator*(EllipticCurvePoint<type>&& point, const uint& value) {
        point *= value;
        return point;
    }

    template<CoordinatesType type>
    EllipticCurvePoint<type> operator*(const uint& value, const EllipticCurvePoint<type>& point) {
        EllipticCurvePoint<type> result = point;
        result *= value;
        return result;
    }

    template<CoordinatesType type>
    EllipticCurvePoint<type> operator*(const uint& value, EllipticCurvePoint<type>&& point) {
        point *= value;
        return point;
    }

    template<>
    class EllipticCurvePoint<CoordinatesType::Normal> : public EllipticCurvePointConcept {
    private:
        friend class EllipticCurve;
        friend EllipticCurvePoint operator-(const EllipticCurvePoint& lhs, EllipticCurvePoint&& rhs);
        friend void multiply<CoordinatesType::Normal>(EllipticCurvePoint& point, const uint& value);

    public:
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

            if (m_x == other.m_x) {
                if (m_y != other.m_y) {
                    m_is_null = true;
                } else {
                    twice();
                }

                return *this;
            }

            FieldElement k = (other.m_y - m_y) / (other.m_x - m_x);
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
            multiply<CoordinatesType::Normal>(*this, value);
            return *this;
        }

        FieldElement get_x() const final {
            return m_x;
        }

        FieldElement get_y() const final {
            return m_y;
        }

    private:
        static EllipticCurvePoint null_point(std::shared_ptr<const FieldElement> a,
                                             std::shared_ptr<const FieldElement>
                                                 b,
                                             std::shared_ptr<const Field>
                                                 F) {
            return EllipticCurvePoint(F->element(0), F->element(1), a, b, F, true);
        }

        static EllipticCurvePoint null_point_from(const EllipticCurvePoint& point) {
            return point.null_point();
        }

        EllipticCurvePoint(const FieldElement& x, const FieldElement& y,
                           std::shared_ptr<const FieldElement> a, std::shared_ptr<const FieldElement> b,
                           std::shared_ptr<const Field> F, bool is_null = false) :
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null), m_x {x}, m_y {y} {
            assert(
                is_valid()
                && "EllipticCurvePoint<CoordinatesType::Normal>::EllipticCurvePoint : invalid coordinates");
        }

        EllipticCurvePoint(FieldElement&& x, const FieldElement& y, std::shared_ptr<const FieldElement> a,
                           std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F,
                           bool is_null = false) :
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
            m_x {std::move(x)},
            m_y {y} {
            assert(
                is_valid()
                && "EllipticCurvePoint<CoordinatesType::Normal>::EllipticCurvePoint : invalid coordinates");
        }

        EllipticCurvePoint(const FieldElement& x, FieldElement&& y, std::shared_ptr<const FieldElement> a,
                           std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F,
                           bool is_null = false) :
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
            m_x {x},
            m_y {std::move(y)} {
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

        void negative() final {
            m_y = -m_y;
        }

        void twice() final {
            if (m_is_null) {
                return;
            }

            if (!m_y.is_invertible()) {
                m_is_null = true;
                return;
            }

            FieldElement k = (m_F->element(3) * FieldElement::pow(m_x, 2) + *m_a) / (m_y << 1);
            FieldElement x = FieldElement::pow(k, 2) - (m_x << 1);
            m_y = k * (m_x - x) - m_y;
            m_x = x;
            assert(is_valid() && "EllipticCurvePoint<CoordinatesType::Normal>::twice : invalid coordinates");
        }

        bool is_valid() const final {
            if (m_is_null) {
                return true;
            }

            FieldElement lhs = FieldElement::pow(m_y, 2);
            FieldElement rhs = FieldElement::pow(m_x, 3) + *m_a * m_x + *m_b;
            return lhs == rhs;
        }

        FieldElement m_x;
        FieldElement m_y;
    };

    template<>
    class EllipticCurvePoint<CoordinatesType::Projective> : public EllipticCurvePointConcept {
    private:
        friend class EllipticCurve;
        friend EllipticCurvePoint operator-(const EllipticCurvePoint& lhs, EllipticCurvePoint&& rhs);
        friend void multiply<CoordinatesType::Projective>(EllipticCurvePoint& point, const uint& value);

    public:
        friend bool operator==(const EllipticCurvePoint& lhs, const EllipticCurvePoint& rhs) {
            FieldElement X1Z2 = lhs.m_X * rhs.m_Z;
            FieldElement X2Z1 = rhs.m_X * lhs.m_Z;
            FieldElement Y1Z2 = lhs.m_Y * rhs.m_Z;
            FieldElement Y2Z1 = rhs.m_Y * lhs.m_Z;
            return (lhs.m_is_null && rhs.m_is_null) || (X1Z2 == X2Z1 && Y1Z2 == Y2Z1);
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

            const FieldElement X1Z2 = m_X * other.m_Z;
            const FieldElement X2Z1 = other.m_X * m_Z;
            const FieldElement Y1Z2 = m_Y * other.m_Z;
            const FieldElement Y2Z1 = other.m_Y * m_Z;

            if (X1Z2 == X2Z1) {
                if (Y1Z2 != Y2Z1) {
                    m_is_null = true;
                } else {
                    twice();
                }

                return *this;
            }

            FieldElement u = Y2Z1 - Y1Z2;
            FieldElement v = X2Z1 - X1Z2;
            FieldElement v2 = FieldElement::pow(v, 2);
            FieldElement v3 = v2 * v;
            const FieldElement Z1Z2 = m_Z * other.m_Z;
            FieldElement A = FieldElement::pow(u, 2) * Z1Z2 - v3 - ((v2 * X1Z2) << 1);

            m_X = v * A;
            m_Y = u * (v2 * X1Z2 - A) - v3 * Y1Z2;
            m_Z = v3 * Z1Z2;

            assert(is_valid()
                   && "EllipticCurvePoint<CoordinatesType::Projective>::operator+= : invalid coordinates");
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
            multiply<CoordinatesType::Projective>(*this, value);
            return *this;
        }

        FieldElement get_x() const final {
            return m_X / m_Z;
        }

        FieldElement get_y() const final {
            return m_Y / m_Z;
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
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
            m_X {x},
            m_Y {y},
            m_Z {m_F->element(1)} {
            assert(
                is_valid()
                && "EllipticCurvePoint<CoordinatesType::Projective>::EllipticCurvePoint : invalid coordinates");
        }

        EllipticCurvePoint(FieldElement&& x, const FieldElement& y, std::shared_ptr<const FieldElement> a,
                           std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F,
                           bool is_null = false) :
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
            m_X {std::move(x)},
            m_Y {y},
            m_Z {m_F->element(1)} {
            assert(
                is_valid()
                && "EllipticCurvePoint<CoordinatesType::Projective>::EllipticCurvePoint : invalid coordinates");
        }

        EllipticCurvePoint(const FieldElement& x, FieldElement&& y, std::shared_ptr<const FieldElement> a,
                           std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F,
                           bool is_null = false) :
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
            m_X {x},
            m_Y {std::move(y)},
            m_Z {m_F->element(1)} {
            assert(
                is_valid()
                && "EllipticCurvePoint<CoordinatesType::Projective>::EllipticCurvePoint : invalid coordinates");
        }

        EllipticCurvePoint(FieldElement&& x, FieldElement&& y, std::shared_ptr<const FieldElement> a,
                           std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F,
                           bool is_null = false) :
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
            m_X {std::move(x)},
            m_Y {std::move(y)},
            m_Z {m_F->element(1)} {
            assert(
                is_valid()
                && "EllipticCurvePoint<CoordinatesType::Projective>::EllipticCurvePoint : invalid coordinates");
        }

        EllipticCurvePoint null_point() const {
            return EllipticCurvePoint(m_F->element(0), m_F->element(1), m_a, m_b, m_F, true);
        }

        void negative() final {
            m_Y = -m_Y;
        }

        void twice() final {
            if (m_is_null) {
                return;
            }

            if (!m_Y.is_invertible()) {
                m_is_null = true;
                return;
            }

            FieldElement w = *m_a * FieldElement::pow(m_Z, 2) + m_F->element(3) * FieldElement::pow(m_X, 2);
            FieldElement s = m_Y * m_Z;
            FieldElement s2 = FieldElement::pow(s, 2);
            FieldElement s3 = s2 * s;
            FieldElement B = m_X * m_Y * s;
            FieldElement h = FieldElement::pow(w, 2) - (B << 3);
            m_X = (h * s) << 1;
            m_Y = w * ((B << 2) - h) - ((FieldElement::pow(m_Y, 2) * s2) << 3);
            m_Z = s3 << 3;
            assert(is_valid()
                   && "EllipticCurvePoint<CoordinatesType::Projective>::twice : invalid coordinates");
        }

        bool is_valid() const final {
            if (m_is_null) {
                return true;
            }

            FieldElement Z2 = FieldElement::pow(m_Z, 2);
            FieldElement Z3 = m_Z * Z2;
            FieldElement lhs = FieldElement::pow(m_Y, 2) * m_Z;
            FieldElement rhs = FieldElement::pow(m_X, 3) + *m_a * m_X * Z2 + *m_b * Z3;
            return lhs == rhs;
        }

        FieldElement m_X;
        FieldElement m_Y;
        FieldElement m_Z;
    };

    template<>
    class EllipticCurvePoint<CoordinatesType::Jacobi> : public EllipticCurvePointConcept {
    private:
        friend class EllipticCurve;
        friend EllipticCurvePoint operator-(const EllipticCurvePoint& lhs, EllipticCurvePoint&& rhs);
        friend void multiply<CoordinatesType::Jacobi>(EllipticCurvePoint& point, const uint& value);

    public:
        friend bool operator==(const EllipticCurvePoint& lhs, const EllipticCurvePoint& rhs) {
            FieldElement X1Z2 = lhs.m_X * FieldElement::pow(rhs.m_Z, 2);
            FieldElement X2Z1 = rhs.m_X * FieldElement::pow(lhs.m_Z, 2);
            FieldElement Y1Z2 = lhs.m_Y * FieldElement::pow(rhs.m_Z, 3);
            FieldElement Y2Z1 = rhs.m_Y * FieldElement::pow(lhs.m_Z, 3);
            return (lhs.m_is_null && rhs.m_is_null) || (X1Z2 == X2Z1 && Y1Z2 == Y2Z1);
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

            const FieldElement X1Z2 = m_X * FieldElement::pow(other.m_Z, 2);
            const FieldElement X2Z1 = other.m_X * FieldElement::pow(m_Z, 2);
            const FieldElement Y1Z2 = m_Y * FieldElement::pow(other.m_Z, 3);
            const FieldElement Y2Z1 = other.m_Y * FieldElement::pow(m_Z, 3);

            if (X1Z2 == X2Z1) {
                if (Y1Z2 != Y2Z1) {
                    m_is_null = true;
                } else {
                    twice();
                }

                return *this;
            }

            FieldElement H = X2Z1 - X1Z2;
            FieldElement H2 = FieldElement::pow(H, 2);
            FieldElement H3 = H2 * H;
            FieldElement r = Y2Z1 - Y1Z2;

            m_X = -H3 - ((X1Z2 * H2) << 1) + FieldElement::pow(r, 2);
            m_Y = -Y1Z2 * H3 + r * (X1Z2 * H2 - m_X);
            m_Z = m_Z * other.m_Z * H;

            assert(is_valid()
                   && "EllipticCurvePoint<CoordinatesType::Jacobi>::operator+= : invalid coordinates");
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
            multiply<CoordinatesType::Jacobi>(*this, value);
            return *this;
        }

        FieldElement get_x() const final {
            return m_X / FieldElement::pow(m_Z, 2);
        }

        FieldElement get_y() const final {
            return m_Y / FieldElement::pow(m_Z, 3);
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
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
            m_X {x},
            m_Y {y},
            m_Z {m_F->element(1)} {
            assert(
                is_valid()
                && "EllipticCurvePoint<CoordinatesType::Jacobi>::EllipticCurvePoint : invalid coordinates");
        }

        EllipticCurvePoint(FieldElement&& x, const FieldElement& y, std::shared_ptr<const FieldElement> a,
                           std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F,
                           bool is_null = false) :
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
            m_X {std::move(x)},
            m_Y {y},
            m_Z {m_F->element(1)} {
            assert(
                is_valid()
                && "EllipticCurvePoint<CoordinatesType::Jacobi>::EllipticCurvePoint : invalid coordinates");
        }

        EllipticCurvePoint(const FieldElement& x, FieldElement&& y, std::shared_ptr<const FieldElement> a,
                           std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F,
                           bool is_null = false) :
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
            m_X {x},
            m_Y {std::move(y)},
            m_Z {m_F->element(1)} {
            assert(
                is_valid()
                && "EllipticCurvePoint<CoordinatesType::Jacobi>::EllipticCurvePoint : invalid coordinates");
        }

        EllipticCurvePoint(FieldElement&& x, FieldElement&& y, std::shared_ptr<const FieldElement> a,
                           std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F,
                           bool is_null = false) :
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
            m_X {std::move(x)},
            m_Y {std::move(y)},
            m_Z {m_F->element(1)} {
            assert(
                is_valid()
                && "EllipticCurvePoint<CoordinatesType::Jacobi>::EllipticCurvePoint : invalid coordinates");
        }

        EllipticCurvePoint null_point() const {
            return EllipticCurvePoint(m_F->element(0), m_F->element(1), m_a, m_b, m_F, true);
        }

        void negative() final {
            m_Y = -m_Y;
        }

        void twice() final {
            if (m_is_null) {
                return;
            }

            if (!m_Y.is_invertible()) {
                m_is_null = true;
                return;
            }

            FieldElement Y2 = FieldElement::pow(m_Y, 2);
            FieldElement Y4 = FieldElement::pow(Y2, 2);
            FieldElement V = (m_X * Y2) << 2;
            FieldElement W = m_F->element(3) * FieldElement::pow(m_X, 2) + *m_a * FieldElement::pow(m_Z, 4);
            m_X = -(V << 1) + FieldElement::pow(W, 2);
            m_Z = (m_Y * m_Z) << 1;
            m_Y = -(Y4 << 3) + W * (V - m_X);
            assert(is_valid() && "EllipticCurvePoint<CoordinatesType::Jacobi>::twice : invalid coordinates");
        }

        bool is_valid() const final {
            if (m_is_null) {
                return true;
            }

            FieldElement Z2 = FieldElement::pow(m_Z, 2);
            FieldElement Z4 = FieldElement::pow(Z2, 2);
            FieldElement Z6 = Z4 * Z2;
            FieldElement value = FieldElement::pow(m_X, 3) + *m_a * m_X * Z4 + *m_b * Z6;
            return FieldElement::pow(m_Y, 2) == value;
        }

        FieldElement m_X;
        FieldElement m_Y;
        FieldElement m_Z;
    };

    template<>
    class EllipticCurvePoint<CoordinatesType::ModifiedJacobi> : public EllipticCurvePointConcept {
    private:
        friend class EllipticCurve;
        friend EllipticCurvePoint operator-(const EllipticCurvePoint& lhs, EllipticCurvePoint&& rhs);
        friend void multiply<CoordinatesType::ModifiedJacobi>(EllipticCurvePoint& point, const uint& value);

    public:
        friend bool operator==(const EllipticCurvePoint& lhs, const EllipticCurvePoint& rhs) {
            FieldElement X1Z2 = lhs.m_X * FieldElement::pow(rhs.m_Z, 2);
            FieldElement X2Z1 = rhs.m_X * FieldElement::pow(lhs.m_Z, 2);
            FieldElement Y1Z2 = lhs.m_Y * FieldElement::pow(rhs.m_Z, 3);
            FieldElement Y2Z1 = rhs.m_Y * FieldElement::pow(lhs.m_Z, 3);
            return (lhs.m_is_null && rhs.m_is_null) || (X1Z2 == X2Z1 && Y1Z2 == Y2Z1);
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

            const FieldElement X1Z2 = m_X * FieldElement::pow(other.m_Z, 2);
            const FieldElement X2Z1 = other.m_X * FieldElement::pow(m_Z, 2);
            const FieldElement Y1Z2 = m_Y * FieldElement::pow(other.m_Z, 3);
            const FieldElement Y2Z1 = other.m_Y * FieldElement::pow(m_Z, 3);

            if (X1Z2 == X2Z1) {
                if (Y1Z2 != Y2Z1) {
                    m_is_null = true;
                } else {
                    twice();
                }

                return *this;
            }

            FieldElement H = X2Z1 - X1Z2;
            FieldElement H2 = FieldElement::pow(H, 2);
            FieldElement H3 = H2 * H;
            FieldElement r = Y2Z1 - Y1Z2;

            m_X = -H3 - ((X1Z2 * H2) << 1) + FieldElement::pow(r, 2);
            m_Y = -Y1Z2 * H3 + r * (X1Z2 * H2 - m_X);
            m_Z = m_Z * other.m_Z * H;
            m_aZ4 = *m_a * FieldElement::pow(m_Z, 4);

            assert(
                is_valid()
                && "EllipticCurvePoint<CoordinatesType::ModifiedJacobi>::operator+= : invalid coordinates");
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
            multiply<CoordinatesType::ModifiedJacobi>(*this, value);
            return *this;
        }

        FieldElement get_x() const final {
            return m_X / FieldElement::pow(m_Z, 2);
        }

        FieldElement get_y() const final {
            return m_Y / FieldElement::pow(m_Z, 3);
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
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
            m_X {x},
            m_Y {y},
            m_Z {m_F->element(1)},
            m_aZ4 {*m_a} {
            assert(
                is_valid()
                && "EllipticCurvePoint<CoordinatesType::ModifiedJacobi>::EllipticCurvePoint : invalid coordinates");
        }

        EllipticCurvePoint(FieldElement&& x, const FieldElement& y, std::shared_ptr<const FieldElement> a,
                           std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F,
                           bool is_null = false) :
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
            m_X {std::move(x)},
            m_Y {y},
            m_Z {m_F->element(1)},
            m_aZ4 {*m_a} {
            assert(
                is_valid()
                && "EllipticCurvePoint<CoordinatesType::ModifiedJacobi>::EllipticCurvePoint : invalid coordinates");
        }

        EllipticCurvePoint(const FieldElement& x, FieldElement&& y, std::shared_ptr<const FieldElement> a,
                           std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F,
                           bool is_null = false) :
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
            m_X {x},
            m_Y {std::move(y)},
            m_Z {m_F->element(1)},
            m_aZ4 {*m_a} {
            assert(
                is_valid()
                && "EllipticCurvePoint<CoordinatesType::ModifiedJacobi>::EllipticCurvePoint : invalid coordinates");
        }

        EllipticCurvePoint(FieldElement&& x, FieldElement&& y, std::shared_ptr<const FieldElement> a,
                           std::shared_ptr<const FieldElement> b, std::shared_ptr<const Field> F,
                           bool is_null = false) :
            EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
            m_X {std::move(x)},
            m_Y {std::move(y)},
            m_Z {m_F->element(1)},
            m_aZ4 {*m_a} {
            assert(
                is_valid()
                && "EllipticCurvePoint<CoordinatesType::ModifiedJacobi>::EllipticCurvePoint : invalid coordinates");
        }

        EllipticCurvePoint null_point() const {
            return EllipticCurvePoint(m_F->element(0), m_F->element(1), m_a, m_b, m_F, true);
        }

        void negative() final {
            m_Y = -m_Y;
        }

        void twice() final {
            if (m_is_null) {
                return;
            }

            if (!m_Y.is_invertible()) {
                m_is_null = true;
                return;
            }

            FieldElement Y2 = FieldElement::pow(m_Y, 2);
            FieldElement V = (m_X * Y2) << 2;
            FieldElement U = FieldElement::pow(Y2, 2) << 3;
            FieldElement W = m_F->element(3) * FieldElement::pow(m_X, 2) + m_aZ4;
            m_X = -(V << 1) + FieldElement::pow(W, 2);
            m_Z = (m_Y * m_Z) << 1;
            m_Y = W * (V - m_X) - U;
            m_aZ4 = (U * m_aZ4) << 1;
            assert(is_valid()
                   && "EllipticCurvePoint<CoordinatesType::ModifiedJacobi>::twice : invalid coordinates");
        }

        bool is_valid() const final {
            if (m_is_null) {
                return true;
            }

            FieldElement Z2 = FieldElement::pow(m_Z, 2);
            FieldElement Z4 = FieldElement::pow(Z2, 2);
            FieldElement Z6 = Z4 * Z2;
            FieldElement value = FieldElement::pow(m_X, 3) + m_X * m_aZ4 + *m_b * Z6;
            return m_aZ4 == (*m_a * Z4) && FieldElement::pow(m_Y, 2) == value;
        }

        FieldElement m_X;
        FieldElement m_Y;
        FieldElement m_Z;
        FieldElement m_aZ4;
    };

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
        EllipticCurvePoint<type> random_point() const {}

    private:
        std::optional<FieldElement> find_y(const FieldElement& x) const;

        std::shared_ptr<const FieldElement> m_a;
        std::shared_ptr<const FieldElement> m_b;
        std::shared_ptr<const Field> m_F;
    };
}   // namespace ECG

#endif
