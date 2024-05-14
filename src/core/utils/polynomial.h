#ifndef ECG_POLYNOMIAL_H
#define ECG_POLYNOMIAL_H

#include "field.h"

namespace elliptic_curve_guide {
    namespace polynomial {
        class Poly {
            using Field = field::Field;
            using Element = field::FieldElement;

        public:
            static Poly pow(const Poly& poly, const uint& power);

            Poly(const Field& field);
            Poly(const Field& field, const std::vector<Element>& coeffs);
            Poly(const Field& field, std::vector<Element>&& coeffs);

            friend Poly operator+(const Poly& lhs, const Poly& rhs);
            friend Poly operator+(Poly&& lhs, const Poly& rhs);
            friend Poly operator+(const Poly& lhs, Poly&& rhs);
            friend Poly operator+(Poly&& lhs, Poly&& rhs);

            friend Poly operator-(const Poly& lhs, const Poly& rhs);
            friend Poly operator-(Poly&& lhs, const Poly& rhs);
            friend Poly operator-(const Poly& lhs, Poly&& rhs);
            friend Poly operator-(Poly&& lhs, Poly&& rhs);

            friend Poly operator*(const Poly& poly, const Element& value);
            friend Poly operator*(Poly&& poly, const Element& value);
            friend Poly operator*(const Element& value, const Poly& poly);
            friend Poly operator*(const Element& value, Poly&& poly);

            friend Poly operator*(const Poly& lhs, const Poly& rhs);
            friend Poly operator%(const Poly& lhs, const Poly& rhs);
            friend Poly operator%(Poly&& lhs, const Poly& rhs);

            Poly operator-() const;

            Poly& operator+=(const Poly& other);
            Poly& operator-=(const Poly& other);
            Poly& operator*=(const Poly& other);
            Poly& operator*=(const Element& value);
            Poly& operator%=(const Poly& other);

            void pow(const uint& power);
            size_t degree() const;
            const Field& get_field() const;

            Element& operator[](const size_t& pos);
            const Element& operator[](const size_t& pos) const;

        private:
            size_t len() const;
            void clean();
            void negative();
            bool is_valid() const;
            const Element& top_coef() const;

            Field m_field;
            std::vector<Element> m_coeffs;
        };
    }   // namespace polynomial
}   // namespace elliptic_curve_guide

#endif
