#include "division_poly.h"

namespace elliptic_curve_guide::polynomial {
    DivisionPoly DivisionPoly::pow(const DivisionPoly& division_poly, const uint& power) {
        DivisionPoly result = division_poly;
        result.pow(power);
        return result;
    }

    DivisionPoly::DivisionPoly(const Poly& x_poly, const std::shared_ptr<const Poly>& curve_poly,
                               const uint& y_power) :
        m_x_poly(x_poly), m_curve_poly(curve_poly), m_y_power(y_power) {}

    DivisionPoly::DivisionPoly(Poly&& x_poly, const std::shared_ptr<const Poly>& curve_poly,
                               const uint& y_power) :
        m_x_poly(std::move(x_poly)), m_curve_poly(curve_poly), m_y_power(y_power) {}

    DivisionPoly operator+(const DivisionPoly& lhs, const DivisionPoly& rhs) {
        DivisionPoly result = lhs;
        result += rhs;
        return result;
    }

    DivisionPoly operator+(DivisionPoly&& lhs, const DivisionPoly& rhs) {
        lhs += rhs;
        return lhs;
    }

    DivisionPoly operator+(const DivisionPoly& lhs, DivisionPoly&& rhs) {
        rhs += lhs;
        return rhs;
    }

    DivisionPoly operator+(DivisionPoly&& lhs, DivisionPoly&& rhs) {
        lhs += rhs;
        return lhs;
    }

    DivisionPoly operator-(const DivisionPoly& lhs, const DivisionPoly& rhs) {
        DivisionPoly result = lhs;
        result -= rhs;
        return result;
    }

    DivisionPoly operator-(DivisionPoly&& lhs, const DivisionPoly& rhs) {
        lhs -= rhs;
        return lhs;
    }

    DivisionPoly operator-(const DivisionPoly& lhs, DivisionPoly&& rhs) {
        rhs -= lhs;
        rhs = -rhs;
        return rhs;
    }

    DivisionPoly operator-(DivisionPoly&& lhs, DivisionPoly&& rhs) {
        lhs -= rhs;
        return lhs;
    }

    DivisionPoly operator*(const DivisionPoly& division_poly, const DivisionPoly::Element& value) {
        DivisionPoly result = division_poly;
        result *= value;
        return result;
    }

    DivisionPoly operator*(DivisionPoly&& division_poly, const DivisionPoly::Element& value) {
        division_poly *= value;
        return division_poly;
    }

    DivisionPoly operator*(const DivisionPoly::Element& value, const DivisionPoly& division_poly) {
        DivisionPoly result = division_poly;
        result *= value;
        return result;
    }

    DivisionPoly operator*(const DivisionPoly::Element& value, DivisionPoly&& division_poly) {
        division_poly *= value;
        return division_poly;
    }

    DivisionPoly operator*(const DivisionPoly& lhs, const DivisionPoly& rhs) {
        DivisionPoly result = lhs;
        result *= rhs;
        return result;
    }

    DivisionPoly operator*(DivisionPoly&& lhs, const DivisionPoly& rhs) {
        lhs *= rhs;
        return lhs;
    }

    DivisionPoly operator*(const DivisionPoly& lhs, DivisionPoly&& rhs) {
        rhs *= lhs;
        return rhs;
    }

    DivisionPoly operator*(DivisionPoly&& lhs, DivisionPoly&& rhs) {
        lhs *= rhs;
        return lhs;
    }

    DivisionPoly DivisionPoly::operator-() const {
        DivisionPoly result = *this;
        result.m_x_poly = -m_x_poly;
        return result;
    }

    DivisionPoly& DivisionPoly::operator+=(const DivisionPoly& other) {
        assert(m_y_power == other.m_y_power && "DivisionPoly::operator+= : incorrect operation");

        if (is_x_poly_null()) {
            m_y_power = other.m_y_power;
        }

        m_x_poly += other.m_x_poly;

        if (is_x_poly_null()) {
            m_y_power = 0;
        }

        return *this;
    }

    DivisionPoly& DivisionPoly::operator-=(const DivisionPoly& other) {
        assert(m_y_power == other.m_y_power && "DivisionPoly::operator-= : incorrect operation");

        if (is_x_poly_null()) {
            m_y_power = other.m_y_power;
        }

        m_x_poly -= other.m_x_poly;

        if (is_x_poly_null()) {
            m_y_power = 0;
        }

        return *this;
    }

    DivisionPoly& DivisionPoly::operator*=(const DivisionPoly& other) {
        m_x_poly *= other.m_x_poly;

        if (is_x_poly_null()) {
            m_y_power = 0;
        } else {
            m_y_power += other.m_y_power;
        }

        return *this;
    }

    DivisionPoly& DivisionPoly::operator*=(const Element& value) {
        m_x_poly *= value;

        if (is_x_poly_null()) {
            m_y_power = 0;
        }

        return *this;
    }

    void DivisionPoly::pow(const uint& power) {
        m_x_poly.pow(power);
        m_y_power *= power;
    }

    void DivisionPoly::divide_by_y() {
        assert(m_y_power > 0 && "DivisionPoly::divide_by_2_y : no y to divide");
        m_y_power -= 1;
    }

    const Poly& DivisionPoly::get_x_poly() const {
        return m_x_poly;
    }

    bool DivisionPoly::is_x_poly_null() const {
        return m_x_poly.degree() == 0 && !m_x_poly[0].is_invertible();
    }

    void DivisionPoly::reduce_y() {
        while (m_y_power > 1) {
            m_x_poly *= *m_curve_poly;
            m_y_power -= 2;
        }
    }
}   // namespace elliptic_curve_guide::polynomial
