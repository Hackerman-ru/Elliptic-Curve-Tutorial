#ifndef ECG_SCHOOF_H
#define ECG_SCHOOF_H

#include "elliptic-curve.h"
#include "uint.h"

namespace elliptic_curve_guide {
    namespace algorithm {
        namespace schoof {
            uint points_number(const elliptic_curve::EllipticCurve& curve);
        }   // namespace schoof
    }       // namespace algorithm
}   // namespace elliptic_curve_guide

#endif
