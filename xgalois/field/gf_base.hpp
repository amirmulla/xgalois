#ifndef XGALOIS_FIELD_GF_BASE_HPP
#define XGALOIS_FIELD_GF_BASE_HPP

#include <cstdint>
#include <iostream>  // Added for std::ostream
#include <memory>
#include <type_traits>
#include <vector>

namespace xg {

enum class FieldRepresentation {
  INT,  // E.g., 5 (default for prime fields)
  HEX,  // E.g., 0x5
  POW,  // E.g., g^k (useful for fields with known generator)
  LOG,  // E.g., log_g(x) (useful for fields with known generator)
  POLY  // E.g., x+1 (default for extension fields)
};

template <typename ElementType>
class GaloisFieldBase {
 public:
  using element_type = ElementType;

  virtual ~GaloisFieldBase() = default;

  // Returns the characteristic of the field (e.g., p for GF(p^n))
  virtual uint32_t Characteristic() const = 0;

  // Returns the order (number of elements) of the field
  virtual uint32_t Order() const = 0;

  // Field operations
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

}  // namespace xg

#endif  // XGALOIS_FIELD_GF_BASE_HPP