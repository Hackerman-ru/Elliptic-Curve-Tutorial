#include "elliptic-curve.h"

#include "utils/field_root.h"

namespace elliptic_curve_guide::elliptic_curve {
    EllipticCurve::EllipticCurve(const Element& a, const Element& b, Field F) :
        m_a {std::make_shared<const Element>(a)},
        m_b {std::make_shared<const Element>(b)},
        m_field {std::make_shared<const Field>(std::move(F))} {}

    EllipticCurve::EllipticCurve(Element&& a, const Element& b, Field F) :
        m_a {std::make_shared<const Element>(std::move(a))},
        m_b {std::make_shared<const Element>(b)},
        m_field {std::make_shared<const Field>(std::move(F))} {}

    EllipticCurve::EllipticCurve(const Element& a, Element&& b, Field F) :
        m_a {std::make_shared<const Element>(a)},
        m_b {std::make_shared<const Element>(std::move(b))},
        m_field {std::make_shared<const Field>(std::move(F))} {}

    EllipticCurve::EllipticCurve(Element&& a, Element&& b, Field F) :
        m_a {std::make_shared<const Element>(std::move(a))},
        m_b {std::make_shared<const Element>(std::move(b))},
        m_field {std::make_shared<const Field>(std::move(F))} {}

    const EllipticCurve::Field& EllipticCurve::get_field() const {
        return *m_field;
    }

    const EllipticCurve::Element& EllipticCurve::get_a() const {
        return *m_a;
    }

    const EllipticCurve::Element& EllipticCurve::get_b() const {
        return *m_b;
    }

    bool EllipticCurve::is_valid_coordinates(const Element& x, const Element& y) const {
        const Element lhs = Element::pow(y, 2);
        const Element rhs = Element::pow(x, 3) + *m_a * x + *m_b;
        return lhs == rhs;
    }

    bool EllipticCurve::is_null_coordinates(const Element& x, const Element& y) const {
        return x.value() == 0 && y.value() == 1;
    }

    std::optional<EllipticCurve::Element> EllipticCurve::find_y(const Element& x) const {
        Element value = Element::pow(x, 3) + *m_a * x + *m_b;
        return algorithm::find_root(value, *m_field);
    }
}   // namespace elliptic_curve_guide::elliptic_curve
