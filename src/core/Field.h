#ifndef ECG_FIELD_H
#define ECG_FIELD_H

#include "util.h"

#ifdef ECG_BOOST
    #include "boost/multiprecision/cpp_int.hpp"
#else
    #include "uint.h"
#endif

#include <string>
#include <vector>

namespace ECG {
#ifdef ECG_BOOST
    using uint = boost::multiprecision::uint512_t;   // should be uint512_t or bigger
#else
    using uint = uint_t<512>;   // should be uint512_t or bigger
#endif
    class FieldElement {
        friend class Field;

        FieldElement(uint value, std::shared_ptr<const uint> p);

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

        friend bool operator==(const FieldElement& lhs, const FieldElement& rhs);
        friend bool operator!=(const FieldElement& lhs, const FieldElement& rhs);
        friend bool operator>(const FieldElement& lhs, const FieldElement& rhs);
        friend bool operator<(const FieldElement& lhs, const FieldElement& rhs);
        friend bool operator>=(const FieldElement& lhs, const FieldElement& rhs);
        friend bool operator<=(const FieldElement& lhs, const FieldElement& rhs);

        bool is_inversible() const;
        FieldElement fast_pow(const uint& pow) const;
        FieldElement inverse() const;

        const uint& get_p() const;

    private:
        uint m_value;
        std::shared_ptr<const uint> m_p;
    };

    class Field {
    public:
        Field(uint p);
        FieldElement operator()(uint value);

    private:
        const std::shared_ptr<const uint> m_p;
    };
}   // namespace ECG

#endif
