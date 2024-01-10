#ifndef ECG_FIELD_H
#define ECG_FIELD_H

#include "../LongInt/long-int.h"

#include <string>
#include <vector>

namespace ECG {
    using uint = uint_t<512>;   // should be uint512_t or bigger

    class PFE {
    public:
        static void set_p(const uint& p);

        template<is_convertible<uint> T>
        explicit PFE(const T& value, bool is_normalized = false);
        PFE(const std::string& str, StringType type);

        PFE operator+(const PFE& other) const;
        PFE operator-(const PFE& other) const;
        PFE operator*(const PFE& other) const;
        PFE operator/(const PFE& other) const;
        PFE inverse() const;

        PFE operator+=(const PFE& other);
        PFE operator-=(const PFE& other);
        PFE operator*=(const PFE& other);
        PFE operator/=(const PFE& other);

        std::string into_string(StringType str_type) const;
        std::string into_string(auto map, size_t shift) const;

    private:
        uint m_value;
        static uint m_p;
    };

    uint PFE::m_p(0);

    template<is_convertible<uint> T>
    inline PFE::PFE(const T& value, bool is_normalized) : m_value(value) {
        if (!is_normalized) {
            m_value %= m_p;
        }

        assert(m_value < m_p);
    }
}   // namespace ECG

#endif
