#ifndef ECG_FIELD_H
#define ECG_FIELD_H

#include "uint.h"

namespace elliptic_curve_guide {
    namespace field {
        class FieldElement {
            friend class Field;

            FieldElement(const uint& value, std::shared_ptr<const uint> modulus);

        public:
            static FieldElement inverse(const FieldElement& element);
            static FieldElement inverse(FieldElement&& element);
            static FieldElement pow(const FieldElement& element, const uint& power);

            friend FieldElement operator+(const FieldElement& lhs, const FieldElement& rhs);
            friend FieldElement operator+(FieldElement&& lhs, const FieldElement& rhs);
            friend FieldElement operator+(const FieldElement& lhs, FieldElement&& rhs);
            friend FieldElement operator+(FieldElement&& lhs, FieldElement&& rhs);

            friend FieldElement operator-(const FieldElement& lhs, const FieldElement& rhs);
            friend FieldElement operator-(FieldElement&& lhs, const FieldElement& rhs);
            friend FieldElement operator-(const FieldElement& lhs, FieldElement&& rhs);
            friend FieldElement operator-(FieldElement&& lhs, FieldElement&& rhs);

            friend FieldElement operator*(const FieldElement& lhs, const FieldElement& rhs);
            friend FieldElement operator*(FieldElement&& lhs, const FieldElement& rhs);
            friend FieldElement operator*(const FieldElement& lhs, FieldElement&& rhs);
            friend FieldElement operator*(FieldElement&& lhs, FieldElement&& rhs);

            friend FieldElement operator/(const FieldElement& lhs, const FieldElement& rhs);
            friend FieldElement operator/(FieldElement&& lhs, const FieldElement& rhs);
            friend FieldElement operator/(const FieldElement& lhs, FieldElement&& rhs);
            friend FieldElement operator/(FieldElement&& lhs, FieldElement&& rhs);

            friend FieldElement operator<<(const FieldElement& value, const uint& shift);
            friend FieldElement operator<<(FieldElement&& value, const uint& shift);

            FieldElement operator-() const;

            FieldElement& operator+=(const FieldElement& other);
            FieldElement& operator-=(const FieldElement& other);
            FieldElement& operator*=(const FieldElement& other);
            FieldElement& operator/=(const FieldElement& other);
            FieldElement& operator<<=(const uint& shift);

            friend bool operator==(const FieldElement& lhs, const FieldElement& rhs);

#ifdef ECG_USE_BOOST
            friend bool operator<(const FieldElement& lhs, const FieldElement& rhs) {
                return lhs.m_value < rhs.m_value;
            }

            friend bool operator>(const FieldElement& lhs, const FieldElement& rhs) {
                return lhs.m_value > rhs.m_value;
            }

            friend bool operator<=(const FieldElement& lhs, const FieldElement& rhs) {
                return lhs.m_value <= rhs.m_value;
            }

            friend bool operator>=(const FieldElement& lhs, const FieldElement& rhs) {
                return lhs.m_value >= rhs.m_value;
            }

            friend bool operator!=(const FieldElement& lhs, const FieldElement& rhs) {
                return lhs.m_value != rhs.m_value;
            }
#else
            friend std::strong_ordering operator<=>(const FieldElement& lhs, const FieldElement& rhs) {
                return lhs.m_value <=> rhs.m_value;
            }
#endif

            bool is_invertible() const;
            void pow(const uint& power);
            void inverse();
            const uint& modulus() const;
            const uint& value() const;

        private:
            static uint normalize(const uint& value, std::shared_ptr<const uint> modulus);
            bool is_valid() const;

            uint m_value;
            std::shared_ptr<const uint> m_modulus;
        };

        class Field {
        public:
#ifdef ECG_USE_BOOST
            Field(const char* str);
            FieldElement element(const char* str) const;
#endif
            Field(const uint& modulus);
            FieldElement element(const uint& value) const;
            const uint& modulus() const;
            bool operator==(const Field& other) const;

        private:
            std::shared_ptr<const uint> m_modulus;
        };
    }   // namespace field
}   // namespace elliptic_curve_guide

#endif
