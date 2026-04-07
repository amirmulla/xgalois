#ifndef XGALOIS_FIELD_GF_BASE_HPP
#define XGALOIS_FIELD_GF_BASE_HPP

#include <cstdint>
#include <iostream>
#include <vector>

namespace xg {

enum class FieldRepresentation : std::uint8_t {
  INT,
  HEX,
  POW,
  LOG,
  POLY
};

template <typename ElementType>
class GaloisFieldBase {
 public:
  using element_type = ElementType;

  virtual ~GaloisFieldBase() = default;

  virtual uint32_t Characteristic() const = 0;

  virtual uint32_t Order() const = 0;

  inline virtual ElementType Add(const ElementType &a,
                                 const ElementType &b) const = 0;
  inline virtual ElementType Sub(const ElementType &a,
                                 const ElementType &b) const = 0;
  inline virtual ElementType Mul(const ElementType &a,
                                 const ElementType &b) const = 0;
  inline virtual ElementType Div(const ElementType &a,
                                 const ElementType &b) const = 0;
  inline virtual ElementType Neg(const ElementType &a) const = 0;
  inline virtual ElementType Inv(const ElementType &a) const = 0;
  inline virtual ElementType Pow(const ElementType &a, uint32_t exp) const = 0;
  inline virtual ElementType Sqrt(const ElementType &a) const = 0;
  inline virtual uint32_t Log(const ElementType &a) const = 0;
  inline virtual uint32_t Log(const ElementType &a,
                              const ElementType &generator) const = 0;

  virtual ElementType Random() const = 0;

  inline virtual ElementType MultiplicativeIdentity() const = 0;
  inline virtual ElementType AdditiveIdentity() const = 0;

  virtual ElementType MultiplicativeGenerator() const = 0;
  virtual std::vector<ElementType> MultiplicativeGenerators() const = 0;

  virtual ElementType GetElementValue(const ElementType &value) const = 0;
  virtual ElementType SetElementValue(const ElementType &value) const = 0;
  virtual ElementType SetElementValue(const std::string &value_str) const = 0;

  virtual void Print(std::ostream &os) const = 0;
  virtual void Print(const ElementType &a, std::ostream &os) const = 0;
  virtual std::string ToString(const ElementType &a) const = 0;

  virtual FieldRepresentation GetRepresentation() const = 0;
  virtual void SetRepresentation(FieldRepresentation rep) = 0;
};

template <typename ElementType>
inline std::ostream &operator<<(std::ostream &os,
                                const GaloisFieldBase<ElementType> &field) {
  field.Print(os);
  return os;
};

}

#endif
