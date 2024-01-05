#ifndef ECG_ELLIPTIC_CURVE_H
#define ECG_ELLIPTIC_CURVE_H

#include "../Field/field.h"

namespace ECG {
    template<typename Field>
    class EllipticCurve {
    public:
        Field m_field;
        FieldElement m_a;
        FieldElement m_b;
        FieldElement m_p;
    };
}   // namespace ECG

#endif
