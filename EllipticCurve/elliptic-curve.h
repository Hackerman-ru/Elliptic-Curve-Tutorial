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
    class EllipticCurvePointBasic {
        using Element = FieldElement<field>;

    public:
        EllipticCurvePointBasic(const Element& x, const Element& y);

        EllipticCurvePointBasic operator+(const EllipticCurvePointBasic& other) const;
        EllipticCurvePointBasic operator-(const EllipticCurvePointBasic& other) const;
        EllipticCurvePointBasic operator*(const uint& other) const;
        EllipticCurvePointBasic operator-() const;

        EllipticCurvePointBasic& operator+=(const EllipticCurvePointBasic& other);
        EllipticCurvePointBasic& operator-=(const EllipticCurvePointBasic& other);

        bool operator==(const EllipticCurvePointBasic& other) const;
        bool is_inf() const;

    private:
        Element m_x;
        Element m_y;
    };

    template<const Field* field, EllipticCurve<field>* elliptic_curve>
    class EllipticCurvePointProjective {
        using Element = FieldElement<field>;

    public:
        EllipticCurvePointProjective(const Element& x, const Element& y);

    private:
        Element m_X;
        Element m_Y;
        Element m_Z;
    };

    template<const Field* field, EllipticCurve<field>* elliptic_curve>
    class EllipticCurvePointJacobi {
        using Element = FieldElement<field>;

    public:
        EllipticCurvePointJacobi(const Element& x, const Element& y);

    private:
        Element m_X;
        Element m_Y;
        Element m_Z;
    };

    template<const Field* field, EllipticCurve<field>* elliptic_curve>
    class EllipticCurvePointJacobiChudnovski {
        using Element = FieldElement<field>;

    public:
        EllipticCurvePointJacobiChudnovski(const Element& x, const Element& y);

    private:
        Element m_X;
        Element m_Y;
        Element m_Z;
        Element m_Z2;
        Element m_Z3;
    };

    template<const Field* field, EllipticCurve<field>* elliptic_curve>
    class EllipticCurvePointModifiedJacobi {
        using Element = FieldElement<field>;

    public:
        EllipticCurvePointModifiedJacobi(const Element& x, const Element& y);

    private:
        Element m_X;
        Element m_Y;
        Element m_Z;
        Element m_aZ4;
    };

    template<const Field* field, EllipticCurve<field>* elliptic_curve>
    class EllipticCurvePointSimplifiedJacobiChudnovski {
        using Element = FieldElement<field>;

    public:
        EllipticCurvePointSimplifiedJacobiChudnovski(const Element& x, const Element& y);

    private:
        Element m_X;
        Element m_Y;
        Element m_Z;
        Element m_Z2;
    };
}   // namespace ECG

#endif
