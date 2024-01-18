#include "elliptic-curve.h"

using ECG::EllipticCurve;
using ECG::EllipticCurvePointBasic;

template<const Field* field>
constexpr EllipticCurve<field>::EllipticCurve(const EllipticCurve<field>::Element& a,
                                              const EllipticCurve<field>::Element& b) :
    m_a(a), m_b(b) {};

template<const Field* field>
ECG::uint EllipticCurve<field>::find_points_number() const {
    // SEA algorithm
}

template<const Field* field>
const typename EllipticCurve<field>::Element& ECG::EllipticCurve<field>::get_a() const {
    return m_a;
}

template<const Field* field>
const typename EllipticCurve<field>::Element& ECG::EllipticCurve<field>::get_b() const {
    return m_b;
}

template<const Field* field>
typename EllipticCurve<field>::Element EllipticCurve<field>::find_y(const Element& x) const {
    return Element();
}

template<const Field* field>
ECG::uint EllipticCurve<field>::generate_random_uint() {
    return uint();
}

template<const Field* field, EllipticCurve<field>* elliptic_curve>
ECG::EllipticCurvePointBasic<field, elliptic_curve>::EllipticCurvePointBasic(const Element& x,
                                                                             const Element& y) :
    m_x(x), m_y(y) {}

template<const Field* field, EllipticCurve<field>* elliptic_curve>
EllipticCurvePointBasic<field, elliptic_curve> ECG::EllipticCurvePointBasic<field, elliptic_curve>::operator+(
    const EllipticCurvePointBasic& other) const {
    if (is_inf() || other.is_inf()) {
        return (is_inf() ? other : *this);
    }

    if (*this == -other) {
        return EllipticCurvePointBasic(1, 0);
    }

    Element k;

    if (*this == other) {
        k = (Element(3) * m_x * m_x + elliptic_curve->get_a())
          / (Element(2) * m_y);   // char(F) > 3 assumption
    } else {
        k = (other.m_y - m_y) / (other.m_x - m_x);
    }

    Element x = k * k - m_x - other.m_x;
    Element y = -k * (x - m_x) - m_y;
    return EllipticCurvePointBasic(x, y);
}

template<const Field* field, EllipticCurve<field>* elliptic_curve>
EllipticCurvePointBasic<field, elliptic_curve> ECG::EllipticCurvePointBasic<field, elliptic_curve>::operator-(
    const EllipticCurvePointBasic& other) const {
    return *this + (-other);
}

template<const Field* field, EllipticCurve<field>* elliptic_curve>
EllipticCurvePointBasic<field, elliptic_curve>
    ECG::EllipticCurvePointBasic<field, elliptic_curve>::operator*(const uint& other) const {
    return EllipticCurvePointBasic();
}

template<const Field* field, EllipticCurve<field>* elliptic_curve>
EllipticCurvePointBasic<field, elliptic_curve>
    ECG::EllipticCurvePointBasic<field, elliptic_curve>::operator-() const {
    return EllipticCurvePointBasic(m_x, -m_y);
}

template<const Field* field, EllipticCurve<field>* elliptic_curve>
EllipticCurvePointBasic<field, elliptic_curve>&
    ECG::EllipticCurvePointBasic<field, elliptic_curve>::operator+=(const EllipticCurvePointBasic& other) {
    return (*this = *this + other);
}

template<const Field* field, EllipticCurve<field>* elliptic_curve>
EllipticCurvePointBasic<field, elliptic_curve>&
    ECG::EllipticCurvePointBasic<field, elliptic_curve>::operator-=(const EllipticCurvePointBasic& other) {
    return (*this = *this - other);
}

template<const Field* field, EllipticCurve<field>* elliptic_curve>
bool ECG::EllipticCurvePointBasic<field, elliptic_curve>::operator==(
    const EllipticCurvePointBasic& other) const {
    return false;
}

template<const Field* field, EllipticCurve<field>* elliptic_curve>
bool ECG::EllipticCurvePointBasic<field, elliptic_curve>::is_inf() const {
    return (m_x == Element(1) && m_y == Element(0));
}

//EllipticCurvePoint EllipticCurvePoint::operator+(const EllipticCurvePoint& other) const {
//    if (is_inf() || other.is_inf()) {
//        return (is_inf() ? other : *this);
//    }
//
//    if (*this == -other) {
//        return EllipticCurvePoint(1, 0);
//    }
//
//    Field k;
//
//    if (*this == other) {
//        k = (Field(3) * m_x * m_x + m_a) / (Field(2) * m_y);   // char(F) > 3 assumption
//    } else {
//        k = (other.m_y - m_y) / (other.m_x - m_x);
//    }
//
//    Field x = k * k - m_x - other.m_x;
//    Field y = -k * (x - m_x) - m_y;
//    return EllipticCurvePoint(x, y);
//}
//
//EllipticCurvePoint EllipticCurvePoint::operator-(const EllipticCurvePoint& other) const {
//    return *this + (-other);
//}
//
//EllipticCurvePoint ECG::EllipticCurvePoint::operator*(const uint& other) const {
//    return EllipticCurvePoint();
//}
//
//EllipticCurvePoint EllipticCurvePoint::operator-() const {
//    return EllipticCurvePoint(m_x, -m_y);
//}
//
//bool EllipticCurvePoint::is_inf() const {
//    return (m_x == Field(1) && m_y == Field(0));
//}
//
//EllipticCurvePoint EllipticCurvePoint::find_point(const std::string& str) {
//    while (true) {
//        const size_t STR_SIZE = str.size();
//        const size_t BYTES = (Field::get_p() >> 4).convert_to<size_t>();
//        uint512_t x_ = generate_random_uint();
//
//        for (size_t i = 0; i < STR_SIZE && i < BYTES; ++i) {
//            x_ <<= 3;
//            x_ |= Field::uint(str[STR_SIZE - i - 1]);
//        }
//
//        Field x(x_);
//        Field y;
//
//        if (find_y(x, &y)) {
//            return EllipticCurvePoint(x, y);
//        }
//    };
//}
//
//EllipticCurvePoint& EllipticCurvePoint::operator+=(const EllipticCurvePoint& other) {
//    return (*this = *this + other);
//}
//
//EllipticCurvePoint& EllipticCurvePoint::operator-=(const EllipticCurvePoint& other) {
//    return (*this = *this - other);
//}
//
//bool EllipticCurvePoint::operator==(const EllipticCurvePoint& other) const {
//    return (m_x == other.m_x && m_y == other.m_y);
//}
