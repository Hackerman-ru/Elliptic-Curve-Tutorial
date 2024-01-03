#ifndef ECG_LONGINT_H
#define ECG_LONGINT_H

#include <array>
#include <functional>
#include <string>

namespace ECG {
    class uint512_t {
        using T = uint64_t;

    public:
        enum class StringType {
            BINARY,
            DECIMAL,
            HEXADECIMAL,
        };

        constexpr explicit uint512_t(T value);
        explicit uint512_t(const std::string& str, StringType str_type);
        explicit uint512_t(const std::string& str, std::function<T(char)> map);

        uint512_t operator+(const uint512_t& other) const;
        uint512_t operator-(const uint512_t& other) const;
        uint512_t operator*(const uint512_t& other) const;
        uint512_t operator/(const uint512_t& other) const;
        uint512_t operator%(const uint512_t& other) const;
        uint512_t operator>>(size_t shift) const;
        uint512_t operator<<(size_t shift) const;
        uint512_t operator^(const uint512_t& other) const;
        uint512_t operator|(const uint512_t& other) const;
        uint512_t operator&(const uint512_t& other) const;

        uint512_t& operator+=(const uint512_t& other);
        uint512_t& operator-=(const uint512_t& other);
        uint512_t& operator*=(const uint512_t& other);
        uint512_t& operator/=(const uint512_t& other);
        uint512_t& operator%=(const uint512_t& other);
        uint512_t& operator>>=(size_t shift);
        uint512_t& operator<<=(size_t shift);
        uint512_t& operator^=(const uint512_t& other);
        uint512_t& operator|=(const uint512_t& other);
        uint512_t& operator&=(const uint512_t& other);

        auto operator<=>(const uint512_t& other) const;
        bool operator==(const uint512_t& other) const;

        explicit operator T() const;

        std::string into_string(StringType str_type) const;
        std::string into_string(std::function<char(T)> map) const;

    private:
        static constexpr size_t Bytes = 64;
        std::array<T, Bytes / sizeof(T)> m_value;
    };
}   // namespace ECG

#endif
