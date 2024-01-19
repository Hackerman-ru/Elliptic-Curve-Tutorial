#ifndef ECG_FIELD_H
#define ECG_FIELD_H

#include "../util.h"

#define ECG_BOOST

#ifdef ECG_BOOST
    #include "boost/multiprecision/cpp_int.hpp"
#else
    #include "../Uint/uint.h"
#endif

#include <string>
#include <vector>

namespace ECG {
#ifdef ECG_BOOST
    using uint = boost::multiprecision::uint512_t;   // should be uint512_t or bigger
#else
    using uint = uint_t<512>;   // should be uint512_t or bigger
#endif
    class FieldElement {
        friend class Field;
        class BaseElement;

        FieldElement(std::unique_ptr<BaseElement>&& ptr);
        FieldElement(uint value, std::shared_ptr<const uint> p);
        FieldElement(std::vector<uint> poly, std::shared_ptr<const std::vector<uint>> reducer);

    public:
        FieldElement(const FieldElement& other);
        FieldElement(FieldElement&& other) noexcept = default;

        FieldElement& operator=(const FieldElement& other);
        FieldElement& operator=(FieldElement&& other) noexcept = default;

        FieldElement operator-() const;

        FieldElement& operator+=(const FieldElement& other);
        FieldElement& operator-=(const FieldElement& other);
        FieldElement& operator*=(const FieldElement& other);
        FieldElement& operator/=(const FieldElement& other);

        FieldElement fast_pow(const uint& pow) const;
        FieldElement inverse() const;

        std::string into_string(StringType str_type = StringType::DECIMAL) const;
        std::string into_string(std::function<char(uint32_t)> map, size_t shift) const;

    private:
        class BaseElement {
        public:
            virtual ~BaseElement() = default;

            virtual bool operator==(const BaseElement& other) const = 0;

            virtual std::unique_ptr<BaseElement> operator-() const = 0;

            virtual void operator+=(const BaseElement& other) = 0;
            virtual void operator-=(const BaseElement& other) = 0;
            virtual void operator*=(const BaseElement& other) = 0;
            virtual void operator/=(const BaseElement& other) = 0;

            virtual std::unique_ptr<BaseElement> inverse() const = 0;
            virtual std::unique_ptr<BaseElement> copy_() const = 0;

            virtual std::string into_string(StringType str_type = StringType::DECIMAL) const = 0;
            virtual std::string into_string(std::function<char(uint32_t)> map, size_t shift) const = 0;
        };

        class PrimeElement final : public BaseElement {
        public:
            PrimeElement(const uint& value, const std::shared_ptr<const uint>& p);

            bool operator==(const BaseElement& other) const override;

            std::unique_ptr<BaseElement> operator-() const override;

            void operator+=(const BaseElement& other) override;
            void operator-=(const BaseElement& other) override;
            void operator*=(const BaseElement& other) override;
            void operator/=(const BaseElement& other) override;

            std::unique_ptr<BaseElement> inverse() const override;
            std::unique_ptr<BaseElement> copy_() const override;

            std::string into_string(StringType str_type = StringType::DECIMAL) const override;
            std::string into_string(std::function<char(uint32_t)> map, size_t shift) const override;

        private:
            uint m_value;
            const std::shared_ptr<const uint> m_p;
        };

        class PolyElement final : public BaseElement {
        public:
            PolyElement(const std::vector<uint>& poly,
                        const std::shared_ptr<const std::vector<uint>>& reducer);

            bool operator==(const BaseElement& other) const override;

            std::unique_ptr<BaseElement> operator-() const override;

            void operator+=(const BaseElement& other) override;
            void operator-=(const BaseElement& other) override;
            void operator*=(const BaseElement& other) override;
            void operator/=(const BaseElement& other) override;

            std::unique_ptr<BaseElement> inverse() const override;
            std::unique_ptr<BaseElement> copy_() const override;

            std::string into_string(StringType str_type = StringType::DECIMAL) const override;
            std::string into_string(std::function<char(uint32_t)> map, size_t shift) const override;

        private:
            std::vector<uint> m_poly;
            const std::shared_ptr<const std::vector<uint>> m_reducer;
        };

        std::unique_ptr<BaseElement> m_self;
    };

    class Field {
    public:
        Field(uint p, uint n = 1);

        FieldElement create_element(uint value);
        FieldElement create_element(std::vector<uint> value);

    private:
        const std::shared_ptr<const uint> m_p;
        const uint m_n;
        const std::shared_ptr<const std::vector<uint>> m_reducer = nullptr;
    };
}   // namespace ECG

#endif
