#ifndef ECG_FIELD_H
#define ECG_FIELD_H

#include "uint.h"

namespace ECG {
    class FieldElement {
        friend class Field;

        FieldElement(const uint& value, std::shared_ptr<const uint> modulus);

    public:
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

        friend auto operator<=>(const FieldElement& lhs, const FieldElement& rhs);
        friend bool operator==(const FieldElement& lhs, const FieldElement& rhs);

        bool is_invertible() const;
        FieldElement pow(const uint& power) const;
        static FieldElement inverse(const FieldElement& element);
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
        Field(const uint& modulus);
        FieldElement element(const uint& value) const;
        const uint& modulus() const;

    private:
        const std::shared_ptr<const uint> m_modulus;
    };
}   // namespace ECG

#endif
