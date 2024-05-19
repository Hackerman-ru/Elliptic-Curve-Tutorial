#ifndef ECG_DIVISION_POLYNOMIAL_H
#define ECG_DIVISION_POLYNOMIAL_H

#include "field.h"
#include "polynomial.h"

namespace elliptic_curve_guide {
    namespace polynomial {
        class DivisionPoly {
            using Field = field::Field;
            using Element = field::FieldElement;

        public:
            static DivisionPoly pow(const DivisionPoly& division_poly, const uint& power);

            DivisionPoly(const Poly& x_poly, const Poly& curve_poly, const Element& y_coef,
                         const size_t& y_power = 1);
            DivisionPoly(Poly&& x_poly, Poly&& curve_poly, Element&& y_coef, const size_t& y_power = 1);

            friend DivisionPoly operator+(const DivisionPoly& lhs, const DivisionPoly& rhs);
            friend DivisionPoly operator+(DivisionPoly&& lhs, const DivisionPoly& rhs);
            friend DivisionPoly operator+(const DivisionPoly& lhs, DivisionPoly&& rhs);
            friend DivisionPoly operator+(DivisionPoly&& lhs, DivisionPoly&& rhs);

            friend DivisionPoly operator-(const DivisionPoly& lhs, const DivisionPoly& rhs);
            friend DivisionPoly operator-(DivisionPoly&& lhs, const DivisionPoly& rhs);
            friend DivisionPoly operator-(const DivisionPoly& lhs, DivisionPoly&& rhs);
            friend DivisionPoly operator-(DivisionPoly&& lhs, DivisionPoly&& rhs);

            friend DivisionPoly operator*(const DivisionPoly& division_poly, const Element& value);
            friend DivisionPoly operator*(DivisionPoly&& division_poly, const Element& value);
            friend DivisionPoly operator*(const Element& value, const DivisionPoly& division_poly);
            friend DivisionPoly operator*(const Element& value, DivisionPoly&& division_poly);

            friend DivisionPoly operator*(const DivisionPoly& lhs, const DivisionPoly& rhs);
            friend DivisionPoly operator*(DivisionPoly&& lhs, const DivisionPoly& rhs);
            friend DivisionPoly operator*(const DivisionPoly& lhs, DivisionPoly&& rhs);
            friend DivisionPoly operator*(DivisionPoly&& lhs, DivisionPoly&& rhs);

            DivisionPoly operator-() const;

            DivisionPoly& operator+=(const DivisionPoly& other);
            DivisionPoly& operator-=(const DivisionPoly& other);
            DivisionPoly& operator*=(const DivisionPoly& other);
            DivisionPoly& operator*=(const Element& value);

            void pow(const uint& power);
            void divide_by_2_y();
            void reduce_y();
            const Poly& get_x_poly() const;
            const Element& get_y_coef() const;

        private:
            Poly m_x_poly;
            Poly m_curve_poly;
            Element m_y_coef;
            size_t m_y_power;
        };
    }   // namespace polynomial
}   // namespace elliptic_curve_guide
#endif
