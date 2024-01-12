#ifndef ECG_FIELD_H
#define ECG_FIELD_H

#include "../LongInt/long-int.h"

#include <string>
#include <vector>

namespace ECG {
    class PFE {
    public:
        using uint = uint_t<512>;   // should be uint512_t or bigger

        static void set_p(const uint& p);
        static uint get_p();

        PFE() = default;

        template<is_convertible<uint> T>
        explicit PFE(const T& value, bool is_normalized = false) : m_value(value) {
            if (!is_normalized && m_value >= m_p) {
                m_value %= m_p;
            }

            assert(m_value < m_p);
        }

        PFE(const std::string& str, StringType type);

        auto operator<=>(const PFE& other) const = default;
        bool operator==(const PFE& other) const;

        PFE operator+(const PFE& other) const;
        PFE operator-(const PFE& other) const;
        PFE operator*(const PFE& other) const;
        PFE operator/(const PFE& other) const;
        PFE operator-() const;
        PFE inverse() const;
        PFE fast_pow(const uint& pow) const;

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

    PFE::uint PFE::m_p(0);
}   // namespace ECG

#endif
