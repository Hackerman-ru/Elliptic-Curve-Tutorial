#ifndef ECG_FIELD_H
#define ECG_FIELD_H

#ifdef ECG_USE_BOOST
    #include "boost/multiprecision/cpp_int.hpp"
    #include "boost/multiprecision/fwd.hpp"

namespace ECG {
    using uint = boost::multiprecision::uint512_t;
}   // namespace ECG
#else
    #include "uint.h"

namespace ECG {
    using uint = uint_t<512>;
}   // namespace ECG
#endif

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

        FieldElement operator-() const;

        FieldElement& operator+=(const FieldElement& other);
        FieldElement& operator-=(const FieldElement& other);
        FieldElement& operator*=(const FieldElement& other);
        FieldElement& operator/=(const FieldElement& other);

        friend auto operator<=>(const FieldElement& lhs, const FieldElement& rhs) = default;

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
        FieldElement operator()(const uint& value);

    private:
        const std::shared_ptr<const uint> m_modulus;
    };
}   // namespace ECG

#endif
