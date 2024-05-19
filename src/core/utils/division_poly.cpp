#include "division_poly.h"

namespace elliptic_curve_guide::polynomial {
    DivisionPoly DivisionPoly::pow(const DivisionPoly& division_poly, const uint& power) {
        DivisionPoly result = division_poly;
        result.pow(power);
        return result;
    }

    DivisionPoly::DivisionPoly(const Poly& x_poly, const std::shared_ptr<const Poly>& curve_poly,
                               const Element& y_coef, const size_t& y_power) :
        m_x_poly(x_poly), m_curve_poly(curve_poly), m_y_coef(y_coef), m_y_power(y_power) {}

    DivisionPoly::DivisionPoly(Poly&& x_poly, const std::shared_ptr<const Poly>& curve_poly, Element&& y_coef,
                               const size_t& y_power) :
        m_x_poly(std::move(x_poly)),
        m_curve_poly(curve_poly),
        m_y_coef(std::move(y_coef)),
        m_y_power(y_power) {}

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

        if (m_y_coef != other.m_y_coef) {
            m_x_poly *= m_y_coef / other.m_y_coef;
            m_y_coef = other.m_y_coef;
        }

        m_x_poly += other.m_x_poly;
    }

    DivisionPoly& DivisionPoly::operator-=(const DivisionPoly& other) {
        assert(m_y_power == other.m_y_power && "DivisionPoly::operator-= : incorrect operation");

        if (m_y_coef != other.m_y_coef) {
            m_x_poly *= m_y_coef / other.m_y_coef;
            m_y_coef = other.m_y_coef;
        }

        m_x_poly -= other.m_x_poly;
    }

    DivisionPoly& DivisionPoly::operator*=(const DivisionPoly& other) {
        m_x_poly *= other.m_x_poly;
        m_y_power += other.m_y_power;
        m_y_coef *= other.m_y_coef;
    }

    DivisionPoly& DivisionPoly::operator*=(const Element& value) {
        m_x_poly *= value;
    }

    void DivisionPoly::pow(const uint& power) {
        m_x_poly.pow(power);
        m_y_coef.pow(2);
        m_y_power <<= 1;
    }

    void DivisionPoly::divide_by_2_y() {
        assert(m_y_power > 0 && "DivisionPoly::divide_by_2_y : no y to divide");
        m_y_power -= 1;
        m_y_coef /= m_x_poly.get_field().element(2);
    }

    const Poly& DivisionPoly::get_x_poly() const {
        return m_x_poly;
    }

    const DivisionPoly::Element& DivisionPoly::get_y_coef() const {
        return m_y_coef;
    }

    void DivisionPoly::reduce_y() {
        while (m_y_power > 1) {
            m_x_poly *= m_curve_poly;
            m_y_power -= 2;
        }

        if (m_y_power == 0) {
            m_x_poly *= m_y_coef;
            m_y_coef = m_x_poly.get_field().element(1);
        } else {
            Element two = m_x_poly.get_field().element(2);
            m_x_poly *= m_y_coef / two;
            m_y_coef = two;
        }
    }
}   // namespace elliptic_curve_guide::polynomial
