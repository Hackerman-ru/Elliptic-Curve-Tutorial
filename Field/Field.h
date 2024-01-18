#ifndef ECG_FIELD_H
#define ECG_FIELD_H

#define ECG_BOOST

#include "../util.h"

#ifdef ECG_BOOST
    #include "boost/multiprecision/cpp_int.hpp"
#else
    #include "../Uint/uint.h"
#endif

#include <string>
#include <vector>

namespace ECG {
#ifdef ECG_BOOST
    using uint = boost::multiprecision::uint512_t;   // should be uint512_t or bigger
#else
    using uint = uint_t<512>;   // should be uint512_t or bigger
#endif

    class Field {
    public:
        template<is_convertible<uint> T>
        constexpr explicit Field(const T& m_p) : m_p(m_p) {};

        const uint& get_p() const {
            return m_p;
        }

    private:
        const uint m_p;
    };

    template<const Field* field>
    class FieldElement {
    public:
        FieldElement(const uint& value) : m_value(value) {
            static const uint& m_p = field->get_p();

            if (m_value >= m_p) {
                m_value %= m_p;
            }

            assert(m_value < m_p);
        }

        auto operator<=>(const FieldElement& other) const = default;

        bool operator==(const FieldElement& other) const {
            return m_value == other.m_value;
        }

        // operator+
        friend FieldElement operator+(const FieldElement& lhs, const FieldElement& rhs) {
            Field result = lhs;
            return result += rhs;
        }

        friend FieldElement operator+(FieldElement&& lhs, const FieldElement& rhs) {
            return lhs += rhs;
        }

        friend FieldElement operator+(const FieldElement& lhs, FieldElement&& rhs) {
            return rhs += lhs;
        }

        friend FieldElement operator+(FieldElement&& lhs, FieldElement&& rhs) {
            return lhs += rhs;
        }

        // operator-
        friend FieldElement operator-(const FieldElement& lhs, const FieldElement& rhs) {
            Field result = lhs;
            return result -= rhs;
        }

        friend FieldElement operator-(FieldElement&& lhs, const FieldElement& rhs) {
            return lhs -= rhs;
        }

        friend FieldElement operator-(const FieldElement& lhs, FieldElement&& rhs) {
            return -(rhs -= lhs);
        }

        friend FieldElement operator-(FieldElement&& lhs, FieldElement&& rhs) {
            return lhs -= rhs;
        }

        // operator*
        friend FieldElement operator*(const FieldElement& lhs, const FieldElement& rhs) {
            Field result = lhs;
            return result *= rhs;
        }

        friend FieldElement operator*(FieldElement&& lhs, const FieldElement& rhs) {
            return lhs *= rhs;
        }

        friend FieldElement operator*(const FieldElement& lhs, FieldElement&& rhs) {
            return rhs *= lhs;
        }

        friend FieldElement operator*(FieldElement&& lhs, FieldElement&& rhs) {
            return lhs *= rhs;
        }

        // operator*
        friend FieldElement operator/(const FieldElement& lhs, const FieldElement& rhs) {
            Field result = lhs;
            return result /= rhs;
        }

        friend FieldElement operator/(FieldElement&& lhs, const FieldElement& rhs) {
            return lhs /= rhs;
        }

        friend FieldElement operator/(const FieldElement& lhs, FieldElement&& rhs) {
            return (rhs.inverse()) * lhs;
        }

        friend FieldElement operator/(FieldElement&& lhs, FieldElement&& rhs) {
            return lhs /= rhs;
        }

        FieldElement operator-() const {
            static const uint& m_p = field->get_p();
            assert(m_value < m_p);
            return FieldElement(m_p - m_value);
        }

        FieldElement inverse() const {
            static const uint& m_p = field->get_p();
            assert(m_value < m_p);
            FieldElement t(0);
            uint r(m_p);
            FieldElement new_t(1);
            uint new_r(m_value);

            while (new_r != 0) {
                uint quotien = r / new_r;
                FieldElement temp(quotien * new_t.m_value, m_p);

                std::make_pair(t, new_t) = std::make_pair(new_t, t - temp);
                std::make_pair(r, new_r) = std::make_pair(new_r, r % new_r);
            }

            assert(r == 1);
            return t;
        }

        FieldElement fast_pow(const uint& pow) const {
            static const uint& m_p = field->get_p();
            assert(m_value < m_p);

            if (pow == 1) {
                return *this;
            }

            if ((pow & 0b1) == 0b1) {
                FieldElement f = fast_pow(pow - 1);
                FieldElement fe = *this * f;
                return fe;
            }

            FieldElement temp = fast_pow(pow >> 1);
            FieldElement t = temp * temp;
            return t;
        }

        FieldElement& operator+=(const FieldElement& other) {
            static const uint& m_p = field->get_p();
            assert(m_value < m_p && other.m_value < m_p);
            m_value += other.m_value;

            if (m_value > m_p) {
                m_value -= m_p;
            }

            return *this;
        }

        FieldElement& operator-=(const FieldElement& other) {
            static const uint& m_p = field->get_p();
            assert(m_value < m_p && other.m_value < m_p);
            uint result;

            if (m_value < other.m_value) {
                m_value += m_p;
            }

            m_value -= other.m_value;

            return *this;
        }

        FieldElement& operator*=(const FieldElement& other) {
            static const uint& m_p = field->get_p();
            assert(m_value < m_p && other.m_value < m_p);
            m_value *= other.m_value;
            m_value %= m_p;
            return *this;
        }

        FieldElement& operator/=(const FieldElement& other) {
            return *this *= other.inverse();
        }

        std::string into_string(StringType str_type = StringType::DECIMAL) const {
#ifdef ECG_BOOST
            std::string str;
            uint clone = m_value;

            switch (str_type) {
            case ECG::StringType::BINARY :
                while (clone != 0) {
                    str += ((clone & 1) != 0) + '0';
                    clone >>= 1;
                }

                std::reverse(str.begin(), str.end());
                break;
            case ECG::StringType::DECIMAL :
                str = m_value.convert_to<std::string>();
                break;
            case ECG::StringType::HEXADECIMAL :
                while (clone != 0) {
                    auto n = clone.convert_to<uint32_t>() & 0xF;

                    if (n < 10) {
                        str += n + '0';
                    } else {
                        n -= 10;
                        str += n + 'a';
                    }

                    clone >>= 4;
                }

                std::reverse(str.begin(), str.end());
                break;
            }

            return str;
#else
            return m_value.into_string(str_type);
#endif
        }

#ifndef ECG_BOOST
        std::string into_string(std::function<uint32_t(char)> map, size_t shift) const {
            return m_value.into_string(map, shift);
        }
#endif

    private:
        uint m_value;
    };
}   // namespace ECG

#endif
