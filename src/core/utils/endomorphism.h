#ifndef ECG_ENDOMORPHISM_H
#define ECG_ENDOMORPHISM_H

#include "polynomial.h"
#include "ring.h"
#include "uint.h"
#include "wnaf.h"

namespace elliptic_curve_guide {
    namespace endomorphism {
        class End {
            using Poly = polynomial::Poly;
            using Ring = ring::Ring;
            using Element = ring::RingElement;

            friend End algorithm::wnaf_addition<End>(End end, const uint& value);

        public:
            End(const Ring& ring, const Poly& a, const Poly& b,
                const std::shared_ptr<const Element>& curve_function);
            End(const Ring& ring, const Element& a, const Element& b,
                const std::shared_ptr<const Element>& curve_function);
            End(const Ring& ring, Element&& a, Element&& b, std::shared_ptr<const Element>&& curve_function);

            struct AdditionResult {
                std::optional<End> end;
                std::optional<Poly> g;
            };

            static AdditionResult twice(const End& end);

            friend AdditionResult operator+(const End& lhs, const End& rhs);

            friend AdditionResult operator-(const End& lhs, const End& rhs);

            friend End operator*(const End& lhs, const End& rhs);
            friend End operator*(End&& lhs, const End& rhs);
            friend End operator*(const End& lhs, End&& rhs);
            friend End operator*(End&& lhs, End&& rhs);

            friend AdditionResult operator*(const End& end, const uint& value);
            friend AdditionResult operator*(const uint& value, const End& end);

            bool operator==(const End& other) const = default;

            End operator-() const;

            End& operator*=(const End& other);
            End& operator*=(const uint& value);

        private:
            void nullify();
            void change_modulus(const Poly& modulus);

            Ring m_ring;
            Element m_a;
            Element m_b;
            std::shared_ptr<const Element> m_curve_function;
        };
    }   // namespace endomorphism
}   // namespace elliptic_curve_guide
#endif
