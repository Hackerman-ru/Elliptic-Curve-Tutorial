#pragma once

#include "../LongInt/LongInt.h"

#include <string>
#include <vector>

namespace ECG {
    using _Ty = uint512_t;   // should be uint256_t or uint512_t

    template<_Ty p>
    class PrimeField {
    public:
        PrimeField() : m_value({0}) {};
        PrimeField(const _Ty& value) : m_value(value) {};

        PrimeField operator+(const PrimeField& other) const {

        };
        PrimeField operator-(const PrimeField& other) const;
        PrimeField operator-() const;
        PrimeField operator*(const PrimeField& other) const;
        PrimeField operator/(const PrimeField& other) const;

        PrimeField& operator+=(const PrimeField& other);
        PrimeField& operator-=(const PrimeField& other);
        PrimeField& operator*=(const PrimeField& other);
        PrimeField& operator/=(const PrimeField& other);

        PrimeField inverse() const;

    private:
        std::vector<uint64_t> m_value;
    };

    template<_Ty p, _Ty n>
    class Field {};
}   // namespace ECG

class PNumber {
public:
    PNumber() {};
    PNumber(int64_t base) : BASE(base) {};

    PNumber(int64_t number, int64_t base) : BASE(base) {
        int64_t remainder;
        remainder = number % base;

        do {
            emplace_back(remainder);
            number /= base;
            remainder = number % base;
        } while (number != 0);
    }

    PNumber(const std::string& str) : BASE(STRING_BASE) {
        num.reserve(str.size());

        for (const auto& symbol : str) {
            emplace_back(char_to_int(symbol));
        }
    }

    bool operator==(const PNumber& other) const {
        size_t n = num.size(), m = other.num.size();

        return !(n != m || num != other.num);
    }

    bool operator<(const PNumber& other) const {
        size_t n = num.size(), m = other.num.size();

        if (n != m) {
            return n < m;
        }

        for (size_t i = n; i > 0; --i) {
            if (num[i - 1] != other.num[i - 1]) {
                return num[i - 1] < other.num[i - 1];
            }
        }

        return false;
    }

    bool operator>(const PNumber& other) const {
        return other < *this;
    }

    PNumber& operator+=(const PNumber& other) {
        size_t n = num.size(), m = other.num.size();

        if (n < m) {
            num.resize(m);
            n = m;
        }

        int64_t remainder = 0, sum;

        for (size_t i = 0; i < n; ++i) {
            if (i < m) {
                sum = num[i] + other.num[i] + remainder;
            } else {
                sum = num[i] + remainder;
            }

            remainder = sum / BASE;
            num[i] = sum - (remainder * BASE);
        }

        if (remainder) {
            emplace_back(remainder);
        }

        return *this;
    }

    PNumber& operator-=(const PNumber& other) {
        if (*this < other) {
            throw("Wrong subtraction");
        }

        size_t t = std::min(num.size(), other.num.size());
        int64_t diff;

        for (size_t i = 0; i < t; ++i) {
            diff = num[i] - other.num[i];

            if (diff < 0) {
                num[i] = diff + BASE;
                num[i + 1] -= 1;
            }
        }

        return *this;
    }

    const PNumber& ConvertToOtherBase(const int64_t& new_base) {
        if (num.empty()) {
            BASE = new_base;

            return *this;
        }

        std::vector<int64_t> res;
        std::vector<int64_t> a = num;
        std::vector<int64_t> b;

        size_t pos = a.size() - 1;
        int64_t curr = a[pos];
        bool shift = false;
        bool first_del = false;

        while (true) {
            if (curr >= new_base) {
                b.emplace_back(curr / new_base);
                curr %= new_base;
                shift = false;
                first_del = true;
            }

            if (pos == 0) {
                res.emplace_back(curr);

                if (b.empty()) {
                    break;
                }

                if (shift) {
                    b.emplace_back(0);
                }

                std::reverse(b.begin(), b.end());
                a = b;
                b.clear();

                pos = a.size() - 1;
                curr = a[pos];
                shift = false;
                first_del = false;

                continue;
            }

            if (shift && first_del) {
                b.emplace_back(0);
            }

            shift = true;
            --pos;
            curr *= BASE;
            curr += a[pos];
        }

        num = res;
        BASE = new_base;

        return *this;
    }

    const std::vector<int64_t>& GetPresentation() const {
        return num;
    }

    void emplace_back(int64_t p_symbol) {
        num.emplace_back(p_symbol);
    }

    std::string str() const {
        PNumber temp = *this;

        if (BASE != STRING_BASE) {
            temp.ConvertToOtherBase(STRING_BASE);
        }

        std::string res;

        for (const auto& symbol : temp.num) {
            res.push_back(int_to_char(symbol));
        }

        return res;
    }

private:
    std::vector<int64_t> num;
    int64_t BASE;
    static constexpr int64_t STRING_BASE = 64;

    static int64_t char_to_int(char symbol) {
        if (symbol >= 48 && symbol <= 57) {
            return symbol - 48;
        }
        if (symbol >= 65 && symbol <= 90) {
            return symbol - 55;
        }
        if (symbol >= 97 && symbol <= 122) {
            return symbol - 61;
        }
        if (symbol == 32) {
            return 62;
        }
        if (symbol == 46) {
            return 63;
        }
        return 64;
    }

    static char int_to_char(int64_t num) {
        if (num >= 0 && num <= 9) {
            return num + 48;
        }
        if (num >= 10 && num <= 35) {
            return num + 55;
        }
        if (num >= 36 && num <= 61) {
            return num + 61;
        }
        if (num == 62) {
            return 32;
        }
        if (num == 63) {
            return 46;
        }
        return 64;
    }
};
