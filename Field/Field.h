#ifndef ECG_FIELD_H
#define ECG_FIELD_H

#include "../LongInt/long-int.h"

#include <string>
#include <vector>

namespace ECG {
    using FieldElement = uint512_t;   // should be uint512_t or bigger

    class PrimeField {
    public:
        explicit PrimeField(const FieldElement& p);

        template<typename T>
        requires std::is_nothrow_convertible_v<T, FieldElement>
        FieldElement normalize(const T& value) const;

        FieldElement add(const FieldElement& lhs, const FieldElement& rhs) const;
        FieldElement sub(const FieldElement& lhs, const FieldElement& rhs) const;
        FieldElement neg(const FieldElement& element) const;
        FieldElement mul(const FieldElement& lhs, const FieldElement& rhs) const;
        FieldElement div(const FieldElement& lhs, const FieldElement& rhs) const;

    private:
        FieldElement m_p;
    };

    class GaloisField {};
}   // namespace ECG

class PNumber {
public:
    PNumber(int64_t number, int64_t base) : BASE(base) {
        int64_t remainder;
        remainder = number % base;

        do {
            emplace_back(remainder);
            number /= base;
            remainder = number % base;
        } while (number != 0);
    }

    bool operator==(const PNumber& other) const {
        size_t n = num.size(), m = other.num.size();

        return !(n != m || num != other.num);
    }

    auto operator<(const PNumber& other) const -> bool {
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

private:
    std::vector<int64_t> num;
    int64_t BASE;
    static constexpr int64_t STRING_BASE = 64;
};

#endif
