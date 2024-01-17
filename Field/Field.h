#ifndef ECG_FIELD_H
#define ECG_FIELD_H

#include "../Uint/uint.h"

#include <string>
#include <vector>

namespace ECG {
    using uint = uint_t<512>;   // should be uint512_t or bigger

    class Field {
    public:
        template<is_convertible<uint> T>
        constexpr explicit Field(const T& p) : p(p) {};

        const uint& get_p() const {
            return p;
        }

    private:
        const uint p;
    };

    template<const Field* field>
    class FieldElement {
        friend class Field;

    public:
        FieldElement(const uint& value);

        auto operator<=>(const FieldElement& other) const = default;
        bool operator==(const FieldElement& other) const;

        FieldElement operator+(const FieldElement& other) const;
        FieldElement operator-(const FieldElement& other) const;
        FieldElement operator*(const FieldElement& other) const;
        FieldElement operator/(const FieldElement& other) const;
        FieldElement operator-() const;
        FieldElement inverse() const;
        FieldElement fast_pow(const uint& pow) const;

        FieldElement operator+=(const FieldElement& other);
        FieldElement operator-=(const FieldElement& other);
        FieldElement operator*=(const FieldElement& other);
        FieldElement operator/=(const FieldElement& other);

        std::string into_string(StringType str_type) const;
        std::string into_string(std::function<uint32_t(char)> map, size_t shift) const;

    private:
        uint m_value;
    };
}   // namespace ECG

#include "field.tpp"

#endif
