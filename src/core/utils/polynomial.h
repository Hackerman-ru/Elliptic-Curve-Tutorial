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
            static Poly decrease_degree_by(const Poly& poly, size_t shift);
            static Poly increase_degree_by(const Poly& poly, size_t shift);

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

            bool operator==(const Poly& other) const = default;

            void pow(const uint& power);
            void decrease_degree();
            void increase_degree();
            void decrease_degree_by(size_t shift);
            void increase_degree_by(size_t shift);
            size_t degree() const;
            void compose(const Poly& inside_poly);
            const Field& get_field() const;
            const Element& top_coef() const;

            Element& operator[](const size_t& pos);
            const Element& operator[](const size_t& pos) const;

        private:
            size_t len() const;
            void clean();
            void negative();
            bool is_valid() const;

            Field m_field;
            std::vector<Element> m_coeffs;
        };
    }   // namespace polynomial

    namespace algorithm {
        struct ModulusGcdResult {
            polynomial::Poly gcd;
            polynomial::Poly value_multiplier;
            polynomial::Poly modulus_multiplier;
        };

        ModulusGcdResult modulus_gcd(const polynomial::Poly& value, const polynomial::Poly& modulus);
        polynomial::Poly gcd(const polynomial::Poly& lhs, const polynomial::Poly& rhs);
        bool has_root(const polynomial::Poly& poly);
    }   // namespace algorithm
}   // namespace elliptic_curve_guide

#endif
