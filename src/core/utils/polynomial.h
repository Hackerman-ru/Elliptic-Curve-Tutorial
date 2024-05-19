#ifndef ECG_POLYNOMIAL_H
#define ECG_POLYNOMIAL_H

#include "field.h"
#include "optional"

namespace elliptic_curve_guide {
    namespace polynomial {
        class Poly {
            using Field = field::Field;
            using Element = field::FieldElement;

        public:
            static Poly pow(const Poly& poly, const uint& power);
            static Poly compose(const Poly& outside_poly, const Poly& inside_poly);

            Poly(const Field& field);
            Poly(const Field& field, const std::vector<Element>& coeffs);
            Poly(const Field& field, std::vector<Element>&& coeffs);
            Poly(const Field& field, const std::vector<uint>& coeffs);

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
            void reduce_degree();
            size_t degree() const;
            void compose(const Poly& inside_poly);
            const Field& get_field() const;

            Element& operator[](const size_t& pos);
            const Element& operator[](const size_t& pos) const;

        private:
            size_t len() const;
            void multiply_by_x();
            void clean();
            void negative();
            bool is_valid() const;
            const Element& top_coef() const;

            Field m_field;
            std::vector<Element> m_coeffs;
        };
    }   // namespace polynomial

    namespace algorithm {
        static polynomial::Poly gcd(const polynomial::Poly& lhs, const polynomial::Poly& rhs);
        bool has_root(const polynomial::Poly& poly);
    }   // namespace algorithm
}   // namespace elliptic_curve_guide

#endif
