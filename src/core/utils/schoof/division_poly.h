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

            DivisionPoly(const Poly& x_poly, const std::shared_ptr<const Poly>& curve_poly,
                         const uint& y_power);
            DivisionPoly(Poly&& x_poly, const std::shared_ptr<const Poly>& curve_poly, const uint& y_power);

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
            void divide_by_y();
            void reduce_y();
            const Poly& get_x_poly() const;

        private:
            bool is_x_poly_null() const;

            Poly m_x_poly;
            std::shared_ptr<const Poly> m_curve_poly;
            uint m_y_power;
        };
    }   // namespace polynomial
}   // namespace elliptic_curve_guide
#endif
