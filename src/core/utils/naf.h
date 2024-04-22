#ifndef ECG_NAF_H
#define ECG_NAF_H

#include "uint.h"

namespace ECG {
    enum class Bit {
        Positive,
        Negative,
        Zero,
    };

    class NAF {
    public:
        NAF(const uint& positive_bits, const uint& negative_bits);

        NAF operator>>(const size_t& shift) const;
        NAF& operator>>=(const size_t& shift);

        bool empty() const;
        Bit bit() const;

    private:
        uint m_positive_bits;
        uint m_negative_bits;
    };

    NAF get_naf(const uint& value);
}   // namespace ECG
#endif
