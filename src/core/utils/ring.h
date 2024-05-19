#ifndef ECG_RING_H
#define ECG_RING_H

#include "polynomial.h"

namespace elliptic_curve_guide {
    namespace ring {
        class RingElement {
            using Element = field::FieldElement;
            using Poly = polynomial::Poly;
            friend class Ring;

            RingElement(const Poly& value, std::shared_ptr<const Poly> modulus);

        public:
            static RingElement pow(const RingElement& element, const uint& power);

            friend RingElement operator+(const RingElement& lhs, const RingElement& rhs);
            friend RingElement operator+(RingElement&& lhs, const RingElement& rhs);
            friend RingElement operator+(const RingElement& lhs, RingElement&& rhs);
            friend RingElement operator+(RingElement&& lhs, RingElement&& rhs);

            friend RingElement operator-(const RingElement& lhs, const RingElement& rhs);
            friend RingElement operator-(RingElement&& lhs, const RingElement& rhs);
            friend RingElement operator-(const RingElement& lhs, RingElement&& rhs);
            friend RingElement operator-(RingElement&& lhs, RingElement&& rhs);

            friend RingElement operator*(const RingElement& lhs, const RingElement& rhs);
            friend RingElement operator*(RingElement&& lhs, const RingElement& rhs);
            friend RingElement operator*(const RingElement& lhs, RingElement&& rhs);
            friend RingElement operator*(RingElement&& lhs, RingElement&& rhs);

            friend RingElement operator*(const RingElement& element, const Element& value);
            friend RingElement operator*(RingElement&& element, const Element& value);
            friend RingElement operator*(const Element& value, const RingElement& element);
            friend RingElement operator*(const Element& value, RingElement&& element);

            RingElement operator-() const;

            RingElement& operator+=(const RingElement& other);
            RingElement& operator-=(const RingElement& other);
            RingElement& operator*=(const RingElement& other);
            RingElement& operator*=(const Element& value);

            friend bool operator==(const RingElement& lhs, const RingElement& rhs);
            friend bool operator!=(const RingElement& lhs, const RingElement& rhs);

            void pow(const uint& power);
            const Poly& modulus() const;
            const Poly& value() const;

        private:
            static Poly normalize(const Poly& value, std::shared_ptr<const Poly> modulus);
            void normalize();
            bool is_valid() const;

            Poly m_value;
            std::shared_ptr<const Poly> m_modulus;
        };

        class Ring {
            using Poly = polynomial::Poly;

        public:
            Ring(const Poly& modulus);
            RingElement element(const Poly& value) const;
            const Poly& modulus() const;
            const field::Field& get_field() const;

        private:
            std::shared_ptr<const Poly> m_modulus;
        };
    }   // namespace ring
}   // namespace elliptic_curve_guide

#endif
