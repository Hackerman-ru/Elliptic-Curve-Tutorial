#ifndef ECG_LONGINT_H
#define ECG_LONGINT_H

#ifndef CHAR_BIT
    #define CHAR_BIT 8
#endif

#include <array>
#include <functional>
#include <string>

template<typename T>
concept UnsignedIntegral = std::is_integral_v<T> && !std::is_signed_v<T>;

template<typename T, typename W>
concept Convertable = std::is_convertible_v<T, W>;

namespace ECG {
    class uint512_t {
        using T = uint32_t;   // must be unsigned integer

    public:
        enum class StringType {
            BINARY = 0b1,
            DECIMAL = 10,
            HEXADECIMAL = 0xF,
        };

        constexpr uint512_t() = default;
        constexpr explicit uint512_t(T value);
        uint512_t(const std::string& str, StringType str_type);
        uint512_t(const std::string& str, std::function<T(char)> map, size_t shift);

        auto operator<=>(const uint512_t& other) const;
        bool operator==(const uint512_t& other) const = default;

        uint512_t operator+(const uint512_t& other) const;
        uint512_t operator-(const uint512_t& other) const;
        uint512_t operator*(const uint512_t& other) const;
        uint512_t operator/(T value) const;
        uint512_t operator%(T value) const;
        uint512_t operator>>(size_t shift) const;
        uint512_t operator<<(size_t shift) const;
        uint512_t operator^(const uint512_t& other) const;
        uint512_t operator|(const uint512_t& other) const;
        uint512_t operator&(const uint512_t& other) const;

        uint512_t operator-() const;

        uint512_t& operator+=(const uint512_t& other);
        uint512_t& operator-=(const uint512_t& other);
        uint512_t& operator*=(const uint512_t& other);
        uint512_t& operator/=(T value);
        uint512_t& operator%=(T value);
        uint512_t& operator>>=(size_t shift);
        uint512_t& operator<<=(size_t shift);
        uint512_t& operator^=(const uint512_t& other);
        uint512_t& operator|=(const uint512_t& other);
        uint512_t& operator&=(const uint512_t& other);

        uint512_t& operator++(int);
        uint512_t& operator++();
        uint512_t& operator--();

        explicit operator T() const;

        // Experimental feature //
        template<typename W>
        requires std::is_nothrow_convertible_v<T, W>
        W convert_to() const {
            size_t shift = sizeof(W) * CHAR_BIT;
            size_t bucket_number = shift / c_BUCKET_SIZE;
            W result = {};

            for (size_t i = 0; i < m_value.size() && i < bucket_number; ++i) {
                result |= W(m_value[i]) << (i * c_BUCKET_SIZE);
            }

            shift %= c_BUCKET_SIZE;

            result |= ((W(m_value[bucket_number]) << (c_BUCKET_SIZE - shift)) >> (c_BUCKET_SIZE - shift))
                   << (bucket_number * c_BUCKET_SIZE);

            return result;
        }

        //

        std::string into_string(StringType str_type) const;
        std::string into_string(std::function<char(T)> map, size_t shift) const;

    private:
        static constexpr size_t c_BITS = 512;
        static constexpr size_t c_BUCKET_SIZE = sizeof(T) * CHAR_BIT;
        std::array<T, c_BITS / c_BUCKET_SIZE> m_value = {};
    };
}   // namespace ECG

#endif