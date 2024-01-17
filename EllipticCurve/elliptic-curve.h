#ifndef ECG_ELLIPTIC_CURVE_H
#define ECG_ELLIPTIC_CURVE_H

#include "../Field/field.h"

namespace ECG {

    template<const Field* field>
    class EllipticCurve {
        using Element = FieldElement<field>;

    public:
        constexpr EllipticCurve(const Element& a, const Element& b);
        uint find_points_number() const;   // SEA algorithm

        const Element& get_a() const;
        const Element& get_b() const;

    private:
        const Element m_a;
        const Element m_b;

        Element find_y(const Element& x) const;
        static uint generate_random_uint();
    };

    template<const Field* field, EllipticCurve<field>* elliptic_curve>
    class EllipticCurvePoint {
        using Element = FieldElement<field>;

    public:
        EllipticCurvePoint(const Element& x, const Element& y);

        EllipticCurvePoint operator+(const EllipticCurvePoint& other) const;
        EllipticCurvePoint operator-(const EllipticCurvePoint& other) const;
        EllipticCurvePoint operator*(const uint& other) const;
        EllipticCurvePoint operator-() const;

        EllipticCurvePoint& operator+=(const EllipticCurvePoint& other);
        EllipticCurvePoint& operator-=(const EllipticCurvePoint& other);

        bool operator==(const EllipticCurvePoint& other) const;
        bool is_inf() const;

    private:
        Element m_x;
        Element m_y;
    };
}   // namespace ECG

#endif
