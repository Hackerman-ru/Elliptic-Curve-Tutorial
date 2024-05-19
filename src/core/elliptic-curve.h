#ifndef ECG_ELLIPTIC_CURVE_H
#define ECG_ELLIPTIC_CURVE_H

#include "field.h"
#include "utils/random.h"
#include "utils/wnaf.h"
//#include "utils/schoof.h"

#include <optional>

namespace elliptic_curve_guide {
    namespace elliptic_curve {
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
                using Field = field::Field;
                using Element = field::FieldElement;

            public:
                virtual Element get_x() const = 0;
                virtual Element get_y() const = 0;

                bool is_zero() const {
                    return m_is_null;
                }

            protected:
                EllipticCurvePointConcept(std::shared_ptr<const Element> a, std::shared_ptr<const Element> b,
                                          std::shared_ptr<const Field> F, bool is_null = false) :
                    m_a {std::move(a)}, m_b {std::move(b)}, m_field {std::move(F)}, m_is_null(is_null) {};

                virtual void negative() = 0;
                virtual void twice() = 0;
                virtual bool is_valid() const = 0;

                void nullify() {
                    m_is_null = true;
                }

                std::shared_ptr<const Element> m_a;
                std::shared_ptr<const Element> m_b;
                std::shared_ptr<const Field> m_field;
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
        EllipticCurvePoint<type> operator+(EllipticCurvePoint<type>&& lhs,
                                           const EllipticCurvePoint<type>& rhs) {
            lhs += rhs;
            return lhs;
        }

        template<CoordinatesType type>
        EllipticCurvePoint<type> operator+(const EllipticCurvePoint<type>& lhs,
                                           EllipticCurvePoint<type>&& rhs) {
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
        EllipticCurvePoint<type> operator-(EllipticCurvePoint<type>&& lhs,
                                           const EllipticCurvePoint<type>& rhs) {
            lhs -= rhs;
            return lhs;
        }

        template<CoordinatesType type>
        EllipticCurvePoint<type> operator-(const EllipticCurvePoint<type>& lhs,
                                           EllipticCurvePoint<type>&& rhs) {
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

        template<CoordinatesType type>
        EllipticCurvePoint<type> operator*(const EllipticCurvePoint<type>& point,
                                           const field::FieldElement& value) {
            EllipticCurvePoint<type> result = point;
            result *= value;
            return result;
        }

        template<CoordinatesType type>
        EllipticCurvePoint<type> operator*(EllipticCurvePoint<type>&& point,
                                           const field::FieldElement& value) {
            point *= value;
            return point;
        }

        template<CoordinatesType type>
        EllipticCurvePoint<type> operator*(const field::FieldElement& value,
                                           const EllipticCurvePoint<type>& point) {
            EllipticCurvePoint<type> result = point;
            result *= value;
            return result;
        }

        template<CoordinatesType type>
        EllipticCurvePoint<type> operator*(const field::FieldElement& value,
                                           EllipticCurvePoint<type>&& point) {
            point *= value;
            return point;
        }

        template<>
        class EllipticCurvePoint<CoordinatesType::Normal> : public EllipticCurvePointConcept {
        private:
            friend class EllipticCurve;
            friend EllipticCurvePoint operator-(const EllipticCurvePoint& lhs, EllipticCurvePoint&& rhs);
            friend EllipticCurvePoint algorithm::wnaf_addition<EllipticCurvePoint>(EllipticCurvePoint value,
                                                                                   const uint& n);

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

                const Element k = (other.m_y - m_y) / (other.m_x - m_x);
                const Element x = Element::pow(k, 2) - m_x - other.m_x;
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
                *this = algorithm::wnaf_addition<EllipticCurvePoint>(*this, value);
                return *this;
            }

            EllipticCurvePoint& operator*=(const field::FieldElement& element) {
                *this = algorithm::wnaf_addition<EllipticCurvePoint>(*this, element.value());
                return *this;
            }

            Element get_x() const final {
                return m_x;
            }

            Element get_y() const final {
                return m_y;
            }

        private:
            static EllipticCurvePoint null_point(std::shared_ptr<const Element> a,
                                                 std::shared_ptr<const Element>
                                                     b,
                                                 std::shared_ptr<const Field>
                                                     F) {
                return EllipticCurvePoint(F->element(0), F->element(1), a, b, F, true);
            }

            static EllipticCurvePoint null_point_from(const EllipticCurvePoint& point) {
                return point.null_point();
            }

            EllipticCurvePoint(const Element& x, const Element& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_x {x},
                m_y {y} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::Normal>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(Element&& x, const Element& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_x {std::move(x)},
                m_y {y} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::Normal>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(const Element& x, Element&& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_x {x},
                m_y {std::move(y)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::Normal>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(Element&& x, Element&& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_x {std::move(x)},
                m_y {std::move(y)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::Normal>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint null_point() const {
                return EllipticCurvePoint(m_field->element(0), m_field->element(1), m_a, m_b, m_field, true);
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

                const Element k = (m_field->element(3) * Element::pow(m_x, 2) + *m_a) / (m_y << 1);
                const Element x = Element::pow(k, 2) - (m_x << 1);
                m_y = k * (m_x - x) - m_y;
                m_x = x;
                assert(is_valid()
                       && "EllipticCurvePoint<CoordinatesType::Normal>::twice : invalid coordinates");
            }

            bool is_valid() const final {
                if (m_is_null) {
                    return true;
                }

                const Element lhs = Element::pow(m_y, 2);
                const Element rhs = Element::pow(m_x, 3) + *m_a * m_x + *m_b;
                return lhs == rhs;
            }

            Element m_x;
            Element m_y;
        };

        template<>
        class EllipticCurvePoint<CoordinatesType::Projective> : public EllipticCurvePointConcept {
        private:
            friend class EllipticCurve;
            friend EllipticCurvePoint operator-(const EllipticCurvePoint& lhs, EllipticCurvePoint&& rhs);
            friend EllipticCurvePoint algorithm::wnaf_addition<EllipticCurvePoint>(EllipticCurvePoint value,
                                                                                   const uint& n);

        public:
            friend bool operator==(const EllipticCurvePoint& lhs, const EllipticCurvePoint& rhs) {
                const Element X1Z2 = lhs.m_X * rhs.m_Z;
                const Element X2Z1 = rhs.m_X * lhs.m_Z;
                const Element Y1Z2 = lhs.m_Y * rhs.m_Z;
                const Element Y2Z1 = rhs.m_Y * lhs.m_Z;
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

                const Element X1Z2 = m_X * other.m_Z;
                const Element X2Z1 = other.m_X * m_Z;
                const Element Y1Z2 = m_Y * other.m_Z;
                const Element Y2Z1 = other.m_Y * m_Z;

                if (X1Z2 == X2Z1) {
                    if (Y1Z2 != Y2Z1) {
                        m_is_null = true;
                    } else {
                        twice();
                    }

                    return *this;
                }

                const Element u = Y2Z1 - Y1Z2;
                const Element v = X2Z1 - X1Z2;
                const Element v2 = Element::pow(v, 2);
                const Element v3 = v2 * v;
                const Element Z1Z2 = m_Z * other.m_Z;
                const Element A = Element::pow(u, 2) * Z1Z2 - v3 - ((v2 * X1Z2) << 1);

                m_X = v * A;
                m_Y = u * (v2 * X1Z2 - A) - v3 * Y1Z2;
                m_Z = v3 * Z1Z2;

                assert(
                    is_valid()
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
                *this = algorithm::wnaf_addition<EllipticCurvePoint>(*this, value);
                return *this;
            }

            EllipticCurvePoint& operator*=(const field::FieldElement& element) {
                *this = algorithm::wnaf_addition<EllipticCurvePoint>(*this, element.value());
                return *this;
            }

            Element get_x() const final {
                return m_X / m_Z;
            }

            Element get_y() const final {
                return m_Y / m_Z;
            }

        private:
            static EllipticCurvePoint null_point(std::shared_ptr<const Element> a,
                                                 std::shared_ptr<const Element>
                                                     b,
                                                 std::shared_ptr<const Field>
                                                     F) {
                return EllipticCurvePoint(F->element(0), F->element(1), a, b, F, true);
            }

            EllipticCurvePoint(const Element& x, const Element& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {x},
                m_Y {y},
                m_Z {m_field->element(1)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::Projective>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(Element&& x, const Element& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {std::move(x)},
                m_Y {y},
                m_Z {m_field->element(1)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::Projective>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(const Element& x, Element&& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {x},
                m_Y {std::move(y)},
                m_Z {m_field->element(1)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::Projective>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(Element&& x, Element&& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {std::move(x)},
                m_Y {std::move(y)},
                m_Z {m_field->element(1)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::Projective>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint null_point() const {
                return EllipticCurvePoint(m_field->element(0), m_field->element(1), m_a, m_b, m_field, true);
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

                const Element w = *m_a * Element::pow(m_Z, 2) + m_field->element(3) * Element::pow(m_X, 2);
                const Element s = m_Y * m_Z;
                const Element s2 = Element::pow(s, 2);
                const Element s3 = s2 * s;
                const Element B = m_X * m_Y * s;
                const Element h = Element::pow(w, 2) - (B << 3);
                m_X = (h * s) << 1;
                m_Y = w * ((B << 2) - h) - ((Element::pow(m_Y, 2) * s2) << 3);
                m_Z = s3 << 3;
                assert(is_valid()
                       && "EllipticCurvePoint<CoordinatesType::Projective>::twice : invalid coordinates");
            }

            bool is_valid() const final {
                if (m_is_null) {
                    return true;
                }

                const Element Z2 = Element::pow(m_Z, 2);
                const Element Z3 = m_Z * Z2;
                const Element lhs = Element::pow(m_Y, 2) * m_Z;
                const Element rhs = Element::pow(m_X, 3) + *m_a * m_X * Z2 + *m_b * Z3;
                return lhs == rhs;
            }

            Element m_X;
            Element m_Y;
            Element m_Z;
        };

        template<>
        class EllipticCurvePoint<CoordinatesType::Jacobi> : public EllipticCurvePointConcept {
        private:
            friend class EllipticCurve;
            friend EllipticCurvePoint operator-(const EllipticCurvePoint& lhs, EllipticCurvePoint&& rhs);
            friend EllipticCurvePoint algorithm::wnaf_addition<EllipticCurvePoint>(EllipticCurvePoint value,
                                                                                   const uint& n);

        public:
            friend bool operator==(const EllipticCurvePoint& lhs, const EllipticCurvePoint& rhs) {
                const Element X1Z2 = lhs.m_X * Element::pow(rhs.m_Z, 2);
                const Element X2Z1 = rhs.m_X * Element::pow(lhs.m_Z, 2);
                const Element Y1Z2 = lhs.m_Y * Element::pow(rhs.m_Z, 3);
                const Element Y2Z1 = rhs.m_Y * Element::pow(lhs.m_Z, 3);
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

                const Element X1Z2 = m_X * Element::pow(other.m_Z, 2);
                const Element X2Z1 = other.m_X * Element::pow(m_Z, 2);
                const Element Y1Z2 = m_Y * Element::pow(other.m_Z, 3);
                const Element Y2Z1 = other.m_Y * Element::pow(m_Z, 3);

                if (X1Z2 == X2Z1) {
                    if (Y1Z2 != Y2Z1) {
                        m_is_null = true;
                    } else {
                        twice();
                    }

                    return *this;
                }

                const Element H = X2Z1 - X1Z2;
                const Element H2 = Element::pow(H, 2);
                const Element H3 = H2 * H;
                const Element r = Y2Z1 - Y1Z2;

                m_X = -H3 - ((X1Z2 * H2) << 1) + Element::pow(r, 2);
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
                *this = algorithm::wnaf_addition<EllipticCurvePoint>(*this, value);
                return *this;
            }

            EllipticCurvePoint& operator*=(const field::FieldElement& element) {
                *this = algorithm::wnaf_addition<EllipticCurvePoint>(*this, element.value());
                return *this;
            }

            Element get_x() const final {
                return m_X / Element::pow(m_Z, 2);
            }

            Element get_y() const final {
                return m_Y / Element::pow(m_Z, 3);
            }

        private:
            static EllipticCurvePoint null_point(std::shared_ptr<const Element> a,
                                                 std::shared_ptr<const Element>
                                                     b,
                                                 std::shared_ptr<const Field>
                                                     F) {
                return EllipticCurvePoint(F->element(0), F->element(1), a, b, F, true);
            }

            EllipticCurvePoint(const Element& x, const Element& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {x},
                m_Y {y},
                m_Z {m_field->element(1)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::Jacobi>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(Element&& x, const Element& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {std::move(x)},
                m_Y {y},
                m_Z {m_field->element(1)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::Jacobi>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(const Element& x, Element&& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {x},
                m_Y {std::move(y)},
                m_Z {m_field->element(1)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::Jacobi>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(Element&& x, Element&& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {std::move(x)},
                m_Y {std::move(y)},
                m_Z {m_field->element(1)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::Jacobi>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint null_point() const {
                return EllipticCurvePoint(m_field->element(0), m_field->element(1), m_a, m_b, m_field, true);
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

                const Element Y2 = Element::pow(m_Y, 2);
                const Element Y4 = Element::pow(Y2, 2);
                const Element V = (m_X * Y2) << 2;
                const Element W = m_field->element(3) * Element::pow(m_X, 2) + *m_a * Element::pow(m_Z, 4);
                m_X = -(V << 1) + Element::pow(W, 2);
                m_Z = (m_Y * m_Z) << 1;
                m_Y = -(Y4 << 3) + W * (V - m_X);
                assert(is_valid()
                       && "EllipticCurvePoint<CoordinatesType::Jacobi>::twice : invalid coordinates");
            }

            bool is_valid() const final {
                if (m_is_null) {
                    return true;
                }

                const Element Z2 = Element::pow(m_Z, 2);
                const Element Z4 = Element::pow(Z2, 2);
                const Element Z6 = Z4 * Z2;
                const Element value = Element::pow(m_X, 3) + *m_a * m_X * Z4 + *m_b * Z6;
                return Element::pow(m_Y, 2) == value;
            }

            Element m_X;
            Element m_Y;
            Element m_Z;
        };

        template<>
        class EllipticCurvePoint<CoordinatesType::JacobiChudnovski> : public EllipticCurvePointConcept {
        private:
            friend class EllipticCurve;
            friend EllipticCurvePoint operator-(const EllipticCurvePoint& lhs, EllipticCurvePoint&& rhs);
            friend EllipticCurvePoint algorithm::wnaf_addition<EllipticCurvePoint>(EllipticCurvePoint value,
                                                                                   const uint& n);

        public:
            friend bool operator==(const EllipticCurvePoint& lhs, const EllipticCurvePoint& rhs) {
                const Element X1Z2 = lhs.m_X * rhs.m_Z2;
                const Element X2Z1 = rhs.m_X * lhs.m_Z2;
                const Element Y1Z2 = lhs.m_Y * rhs.m_Z3;
                const Element Y2Z1 = rhs.m_Y * lhs.m_Z3;
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

                const Element X1Z2 = m_X * other.m_Z2;
                const Element X2Z1 = other.m_X * m_Z2;
                const Element Y1Z2 = m_Y * other.m_Z3;
                const Element Y2Z1 = other.m_Y * m_Z3;

                if (X1Z2 == X2Z1) {
                    if (Y1Z2 != Y2Z1) {
                        m_is_null = true;
                    } else {
                        twice();
                    }

                    return *this;
                }

                const Element H = X2Z1 - X1Z2;
                const Element H2 = Element::pow(H, 2);
                const Element H3 = H2 * H;
                const Element r = Y2Z1 - Y1Z2;

                m_X = -H3 - ((X1Z2 * H2) << 1) + Element::pow(r, 2);
                m_Y = -Y1Z2 * H3 + r * (X1Z2 * H2 - m_X);
                m_Z = m_Z * other.m_Z * H;
                m_Z2 = Element::pow(m_Z, 2);
                m_Z3 = m_Z * m_Z2;

                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::JacobiChudnovski>::operator+= : invalid coordinates");
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
                *this = algorithm::wnaf_addition<EllipticCurvePoint>(*this, value);
                return *this;
            }

            EllipticCurvePoint& operator*=(const field::FieldElement& element) {
                *this = algorithm::wnaf_addition<EllipticCurvePoint>(*this, element.value());
                return *this;
            }

            Element get_x() const final {
                return m_X / m_Z2;
            }

            Element get_y() const final {
                return m_Y / m_Z3;
            }

        private:
            static EllipticCurvePoint null_point(std::shared_ptr<const Element> a,
                                                 std::shared_ptr<const Element>
                                                     b,
                                                 std::shared_ptr<const Field>
                                                     F) {
                return EllipticCurvePoint(F->element(0), F->element(1), a, b, F, true);
            }

            EllipticCurvePoint(const Element& x, const Element& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {x},
                m_Y {y},
                m_Z {m_field->element(1)},
                m_Z2 {m_field->element(1)},
                m_Z3 {m_field->element(1)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::JacobiChudnovski>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(Element&& x, const Element& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {std::move(x)},
                m_Y {y},
                m_Z {m_field->element(1)},
                m_Z2 {m_field->element(1)},
                m_Z3 {m_field->element(1)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::JacobiChudnovski>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(const Element& x, Element&& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {x},
                m_Y {std::move(y)},
                m_Z {m_field->element(1)},
                m_Z2 {m_field->element(1)},
                m_Z3 {m_field->element(1)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::JacobiChudnovski>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(Element&& x, Element&& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {std::move(x)},
                m_Y {std::move(y)},
                m_Z {m_field->element(1)},
                m_Z2 {m_field->element(1)},
                m_Z3 {m_field->element(1)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::JacobiChudnovski>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint null_point() const {
                return EllipticCurvePoint(m_field->element(0), m_field->element(1), m_a, m_b, m_field, true);
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

                const Element Y2 = Element::pow(m_Y, 2);
                const Element Y4 = Element::pow(Y2, 2);
                const Element V = (m_X * Y2) << 2;
                const Element W = m_field->element(3) * Element::pow(m_X, 2) + *m_a * Element::pow(m_Z2, 2);
                m_X = -(V << 1) + Element::pow(W, 2);
                m_Z = (m_Y * m_Z) << 1;
                m_Y = -(Y4 << 3) + W * (V - m_X);
                m_Z2 = Element::pow(m_Z, 2);
                m_Z3 = m_Z * m_Z2;
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::JacobiChudnovski>::twice : invalid coordinates");
            }

            bool is_valid() const final {
                if (m_is_null) {
                    return true;
                }

                const Element Z4 = Element::pow(m_Z2, 2);
                const Element Z6 = Element::pow(m_Z3, 2);
                const Element value = Element::pow(m_X, 3) + *m_a * m_X * Z4 + *m_b * Z6;
                return Element::pow(m_Y, 2) == value;
            }

            Element m_X;
            Element m_Y;
            Element m_Z;
            Element m_Z2;
            Element m_Z3;
        };

        template<>
        class EllipticCurvePoint<CoordinatesType::ModifiedJacobi> : public EllipticCurvePointConcept {
        private:
            friend class EllipticCurve;
            friend EllipticCurvePoint operator-(const EllipticCurvePoint& lhs, EllipticCurvePoint&& rhs);
            friend EllipticCurvePoint algorithm::wnaf_addition<EllipticCurvePoint>(EllipticCurvePoint value,
                                                                                   const uint& n);

        public:
            friend bool operator==(const EllipticCurvePoint& lhs, const EllipticCurvePoint& rhs) {
                const Element X1Z2 = lhs.m_X * Element::pow(rhs.m_Z, 2);
                const Element X2Z1 = rhs.m_X * Element::pow(lhs.m_Z, 2);
                const Element Y1Z2 = lhs.m_Y * Element::pow(rhs.m_Z, 3);
                const Element Y2Z1 = rhs.m_Y * Element::pow(lhs.m_Z, 3);
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

                const Element X1Z2 = m_X * Element::pow(other.m_Z, 2);
                const Element X2Z1 = other.m_X * Element::pow(m_Z, 2);
                const Element Y1Z2 = m_Y * Element::pow(other.m_Z, 3);
                const Element Y2Z1 = other.m_Y * Element::pow(m_Z, 3);

                if (X1Z2 == X2Z1) {
                    if (Y1Z2 != Y2Z1) {
                        m_is_null = true;
                    } else {
                        twice();
                    }

                    return *this;
                }

                const Element H = X2Z1 - X1Z2;
                const Element H2 = Element::pow(H, 2);
                const Element H3 = H2 * H;
                const Element r = Y2Z1 - Y1Z2;

                m_X = -H3 - ((X1Z2 * H2) << 1) + Element::pow(r, 2);
                m_Y = -Y1Z2 * H3 + r * (X1Z2 * H2 - m_X);
                m_Z = m_Z * other.m_Z * H;
                m_aZ4 = *m_a * Element::pow(m_Z, 4);

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
                *this = algorithm::wnaf_addition<EllipticCurvePoint>(*this, value);
                return *this;
            }

            EllipticCurvePoint& operator*=(const field::FieldElement& element) {
                *this = algorithm::wnaf_addition<EllipticCurvePoint>(*this, element.value());
                return *this;
            }

            Element get_x() const final {
                return m_X / Element::pow(m_Z, 2);
            }

            Element get_y() const final {
                return m_Y / Element::pow(m_Z, 3);
            }

        private:
            static EllipticCurvePoint null_point(std::shared_ptr<const Element> a,
                                                 std::shared_ptr<const Element>
                                                     b,
                                                 std::shared_ptr<const Field>
                                                     F) {
                return EllipticCurvePoint(F->element(0), F->element(1), a, b, F, true);
            }

            EllipticCurvePoint(const Element& x, const Element& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {x},
                m_Y {y},
                m_Z {m_field->element(1)},
                m_aZ4 {*m_a} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::ModifiedJacobi>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(Element&& x, const Element& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {std::move(x)},
                m_Y {y},
                m_Z {m_field->element(1)},
                m_aZ4 {*m_a} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::ModifiedJacobi>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(const Element& x, Element&& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {x},
                m_Y {std::move(y)},
                m_Z {m_field->element(1)},
                m_aZ4 {*m_a} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::ModifiedJacobi>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(Element&& x, Element&& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {std::move(x)},
                m_Y {std::move(y)},
                m_Z {m_field->element(1)},
                m_aZ4 {*m_a} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::ModifiedJacobi>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint null_point() const {
                return EllipticCurvePoint(m_field->element(0), m_field->element(1), m_a, m_b, m_field, true);
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

                const Element Y2 = Element::pow(m_Y, 2);
                const Element V = (m_X * Y2) << 2;
                const Element U = Element::pow(Y2, 2) << 3;
                const Element W = m_field->element(3) * Element::pow(m_X, 2) + m_aZ4;
                m_X = -(V << 1) + Element::pow(W, 2);
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

                const Element Z2 = Element::pow(m_Z, 2);
                const Element Z4 = Element::pow(Z2, 2);
                const Element Z6 = Z4 * Z2;
                const Element value = Element::pow(m_X, 3) + m_X * m_aZ4 + *m_b * Z6;
                return m_aZ4 == (*m_a * Z4) && Element::pow(m_Y, 2) == value;
            }

            Element m_X;
            Element m_Y;
            Element m_Z;
            Element m_aZ4;
        };

        template<>
        class EllipticCurvePoint<CoordinatesType::SimplifiedJacobiChudnovski>
            : public EllipticCurvePointConcept {
        private:
            friend class EllipticCurve;
            friend EllipticCurvePoint operator-(const EllipticCurvePoint& lhs, EllipticCurvePoint&& rhs);
            friend EllipticCurvePoint algorithm::wnaf_addition<EllipticCurvePoint>(EllipticCurvePoint value,
                                                                                   const uint& n);

        public:
            friend bool operator==(const EllipticCurvePoint& lhs, const EllipticCurvePoint& rhs) {
                const Element X1Z2 = lhs.m_X * rhs.m_Z2;
                const Element X2Z1 = rhs.m_X * lhs.m_Z2;
                const Element Y1Z2 = lhs.m_Y * rhs.m_Z2 * rhs.m_Z;
                const Element Y2Z1 = rhs.m_Y * lhs.m_Z2 * lhs.m_Z;
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

                const Element X1Z2 = m_X * other.m_Z2;
                const Element X2Z1 = other.m_X * m_Z2;
                const Element Y1Z2 = m_Y * other.m_Z2 * other.m_Z;
                const Element Y2Z1 = other.m_Y * m_Z2 * m_Z;

                if (X1Z2 == X2Z1) {
                    if (Y1Z2 != Y2Z1) {
                        m_is_null = true;
                    } else {
                        twice();
                    }

                    return *this;
                }

                const Element H = X2Z1 - X1Z2;
                const Element H2 = Element::pow(H, 2);
                const Element H3 = H2 * H;
                const Element r = Y2Z1 - Y1Z2;

                m_X = -H3 - ((X1Z2 * H2) << 1) + Element::pow(r, 2);
                m_Y = -Y1Z2 * H3 + r * (X1Z2 * H2 - m_X);
                m_Z = m_Z * other.m_Z * H;
                m_Z2 = Element::pow(m_Z, 2);

                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::SimplifiedJacobiChudnovski>::operator+= : invalid coordinates");
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
                *this = algorithm::wnaf_addition<EllipticCurvePoint>(*this, value);
                return *this;
            }

            EllipticCurvePoint& operator*=(const field::FieldElement& element) {
                *this = algorithm::wnaf_addition<EllipticCurvePoint>(*this, element.value());
                return *this;
            }

            Element get_x() const final {
                return m_X / m_Z2;
            }

            Element get_y() const final {
                return m_Y / (m_Z * m_Z2);
            }

        private:
            static EllipticCurvePoint null_point(std::shared_ptr<const Element> a,
                                                 std::shared_ptr<const Element>
                                                     b,
                                                 std::shared_ptr<const Field>
                                                     F) {
                return EllipticCurvePoint(F->element(0), F->element(1), a, b, F, true);
            }

            EllipticCurvePoint(const Element& x, const Element& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {x},
                m_Y {y},
                m_Z {m_field->element(1)},
                m_Z2 {m_field->element(1)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::SimplifiedJacobiChudnovski>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(Element&& x, const Element& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {std::move(x)},
                m_Y {y},
                m_Z {m_field->element(1)},
                m_Z2 {m_field->element(1)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::SimplifiedJacobiChudnovski>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(const Element& x, Element&& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {x},
                m_Y {std::move(y)},
                m_Z {m_field->element(1)},
                m_Z2 {m_field->element(1)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::SimplifiedJacobiChudnovski>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint(Element&& x, Element&& y, std::shared_ptr<const Element> a,
                               std::shared_ptr<const Element> b, std::shared_ptr<const Field> F,
                               bool is_null = false) :
                EllipticCurvePointConcept(std::move(a), std::move(b), std::move(F), is_null),
                m_X {std::move(x)},
                m_Y {std::move(y)},
                m_Z {m_field->element(1)},
                m_Z2 {m_field->element(1)} {
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::SimplifiedJacobiChudnovski>::EllipticCurvePoint : invalid coordinates");
            }

            EllipticCurvePoint null_point() const {
                return EllipticCurvePoint(m_field->element(0), m_field->element(1), m_a, m_b, m_field, true);
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

                const Element Y2 = Element::pow(m_Y, 2);
                const Element Y4 = Element::pow(Y2, 2);
                const Element V = (m_X * Y2) << 2;
                const Element W = m_field->element(3) * Element::pow(m_X, 2) + *m_a * Element::pow(m_Z2, 2);
                m_X = -(V << 1) + Element::pow(W, 2);
                m_Z = (m_Y * m_Z) << 1;
                m_Y = -(Y4 << 3) + W * (V - m_X);
                m_Z2 = Element::pow(m_Z, 2);
                assert(
                    is_valid()
                    && "EllipticCurvePoint<CoordinatesType::SimplifiedJacobiChudnovski>::twice : invalid coordinates");
            }

            bool is_valid() const final {
                if (m_is_null) {
                    return true;
                }

                const Element Z4 = Element::pow(m_Z2, 2);
                const Element Z6 = Z4 * m_Z2;
                const Element value = Element::pow(m_X, 3) + *m_a * m_X * Z4 + *m_b * Z6;
                return Element::pow(m_Y, 2) == value;
            }

            Element m_X;
            Element m_Y;
            Element m_Z;
            Element m_Z2;
        };

        class EllipticCurve {
            using Field = field::Field;
            using Element = field::FieldElement;

        public:
            EllipticCurve(const Element& a, const Element& b, Field F);
            EllipticCurve(Element&& a, const Element& b, Field F);
            EllipticCurve(const Element& a, Element&& b, Field F);
            EllipticCurve(Element&& a, Element&& b, Field F);

            const Field& get_field() const;
            const Element& get_a() const;
            const Element& get_b() const;

            template<CoordinatesType type = CoordinatesType::Normal>
            std::optional<EllipticCurvePoint<type>> point_with_x_equal_to(const Element& x) const {
                if (!x.is_invertible()) {
                    return null_point<type>();
                }

                std::optional<Element> y = find_y(x);

                if (!y.has_value()) {
                    return std::nullopt;
                }

                return EllipticCurvePoint<type>(x, std::move(y.value()), m_a, m_b, m_field);
            }

            template<CoordinatesType type = CoordinatesType::Normal>
            std::optional<EllipticCurvePoint<type>> point_with_x_equal_to(Element&& x) const {
                if (!x.is_invertible()) {
                    return null_point<type>();
                }

                std::optional<Element> y = find_y(x);

                if (!y.has_value()) {
                    return std::nullopt;
                }

                return EllipticCurvePoint<type>(std::move(x), std::move(y.value()), m_a, m_b, m_field);
            }

            template<CoordinatesType type = CoordinatesType::Normal>
            std::optional<EllipticCurvePoint<type>> point(const Element& x, const Element& y) const {
                if (is_null_coordinates(x, y)) {
                    return null_point();
                }

                if (!is_valid_coordinates(x, y)) {
                    return std::nullopt;
                }

                return EllipticCurvePoint<type>(x, y, m_a, m_b, m_field);
            }

            template<CoordinatesType type = CoordinatesType::Normal>
            std::optional<EllipticCurvePoint<type>> point(Element&& x, const Element& y) const {
                if (is_null_coordinates(x, y)) {
                    return null_point();
                }

                if (!is_valid_coordinates(x, y)) {
                    return std::nullopt;
                }

                return EllipticCurvePoint<type>(std::move(x), y, m_a, m_b, m_field);
            }

            template<CoordinatesType type = CoordinatesType::Normal>
            std::optional<EllipticCurvePoint<type>> point(const Element& x, Element&& y) const {
                if (is_null_coordinates(x, y)) {
                    return null_point();
                }

                if (!is_valid_coordinates(x, y)) {
                    return std::nullopt;
                }

                return EllipticCurvePoint<type>(x, std::move(y), m_a, m_b, m_field);
            }

            template<CoordinatesType type = CoordinatesType::Normal>
            std::optional<EllipticCurvePoint<type>> point(Element&& x, Element&& y) const {
                if (is_null_coordinates(x, y)) {
                    return null_point();
                }

                if (!is_valid_coordinates(x, y)) {
                    return std::nullopt;
                }

                return EllipticCurvePoint<type>(std::move(x), std::move(y), m_a, m_b, m_field);
            }

            template<CoordinatesType type = CoordinatesType::Normal>
            EllipticCurvePoint<type> null_point() const {
                return EllipticCurvePoint<type>::null_point(m_a, m_b, m_field);
            }

            template<CoordinatesType type = CoordinatesType::Normal>
            EllipticCurvePoint<type> random_point() const {
                static constexpr size_t c_repeat_number = 1000;
                const uint& p = m_field->modulus();

                for (size_t i = 0; i < c_repeat_number; ++i) {
                    uint x = algorithm::random::generate_random_uint_modulo(p);
                    auto opt = point_with_x_equal_to<type>(m_field->element(x));

                    if (opt.has_value()) {
                        return opt.value();
                    }
                }

                return null_point<type>();
            }

        private:
            bool is_null_coordinates(const Element& x, const Element& y) const;
            bool is_valid_coordinates(const Element& x, const Element& y) const;
            std::optional<Element> find_y(const Element& x) const;

            std::shared_ptr<const Element> m_a;
            std::shared_ptr<const Element> m_b;
            std::shared_ptr<const Field> m_field;
        };
    }   // namespace elliptic_curve
}   // namespace elliptic_curve_guide

#endif
