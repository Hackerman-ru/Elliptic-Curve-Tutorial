#pragma once

#include <array>

class uint512_t {
public:
    uint512_t(uint64_t value) : m_value({value}) {};

    uint512_t operator+(const uint512_t& other) const;
    uint512_t operator-(const uint512_t& other) const;
    uint512_t operator*(const uint512_t& other) const;
    uint512_t operator/(const uint512_t& other) const;

    uint512_t& operator+=(const uint512_t& other);
    uint512_t& operator-=(const uint512_t& other);
    uint512_t& operator*=(const uint512_t& other);
    uint512_t& operator/=(const uint512_t& other);

private:
    std::array<uint64_t, 8> m_value;
};
