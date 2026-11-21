#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <utility>

#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_extension.hpp"

using namespace xg;

namespace {

template <typename Field>
typename Field::ElementType EvaluateOverBase(
    const Field& field, const typename Field::ElementType& poly_over_fp,
    const typename Field::ElementType& point) {
  auto acc = field.AdditiveIdentity();
  for (int d = poly_over_fp.Degree(); d >= 0; --d) {
    acc = field.Mul(acc, point);
    const auto c = poly_over_fp[static_cast<size_t>(d)];
    if (c != 0) {
      acc = field.Add(acc, field.FromCoeffs({c}));
    }
  }
  return acc;
}

}  // namespace

int main() {
  std::cout << "XGalois Extension Field Examples\n";
  std::cout << "================================\n";

  using ExtField = GFPX<uint8_t>;
  using ExtElem = GaloisFieldElement<ExtField>;

  std::cout << "\n=== Creating Extension Fields GF(p^m) ===\n";
  auto gf9 = std::make_shared<ExtField>(std::pair<uint64_t, uint64_t>{3, 2});
  auto gf25 = std::make_shared<ExtField>(std::pair<uint64_t, uint64_t>{5, 2});
  auto gf49 = std::make_shared<ExtField>(std::pair<uint64_t, uint64_t>{7, 2});

  std::cout << "  GF(3^2): order=" << gf9->Order()
            << ", modulus=" << gf9->ToString(gf9->Modulus()) << "\n";
  std::cout << "  GF(5^2): order=" << gf25->Order()
            << ", modulus=" << gf25->ToString(gf25->Modulus()) << "\n";
  std::cout << "  GF(7^2): order=" << gf49->Order()
            << ", modulus=" << gf49->ToString(gf49->Modulus()) << "\n";

  std::cout << "\n=== Arithmetic in GF(3^2) ===\n";
  ExtElem a("α + 1", gf9);
  ExtElem b("2*α + 1", gf9);
  std::cout << "  a = " << a << ", b = " << b << "\n";
  std::cout << "  a + b  = " << (a + b) << "\n";
  std::cout << "  a - b  = " << (a - b) << "\n";
  std::cout << "  a * b  = " << (a * b) << "\n";
  std::cout << "  a / b  = " << (a / b) << "\n";
  std::cout << "  a^(-1) = " << a.Inv()
            << " (check a * a^(-1) = " << (a * a.Inv()) << ")\n";

  std::cout << "\n=== MultiplicativeGenerator & Discrete Log in GF(3^2) ===\n";
  ExtElem gen9(gf9->MultiplicativeGenerator(), gf9);
  std::cout << "  Multiplicative generator g = " << gen9
            << " (order = " << gen9.MultiplicativeOrder() << ")\n";
  ExtElem power(gf9->MultiplicativeIdentity(), gf9);
  for (uint64_t i = 0; i < gf9->Order() - 1; ++i) {
    std::cout << "    g^" << i << " = " << power << "  (Log = " << power.Log()
              << ")\n";
    power *= gen9;
  }

  const uint64_t log_a = a.Log();
  std::cout << "  Log_g(a) = " << log_a
            << ", check g^Log(a) == a: " << std::boolalpha
            << (gen9.Pow(log_a) == a) << "\n";

  std::cout << "\n=== Sqrt, Trace, Norm & MinimalPolynomial in GF(5^2) ===\n";
  ExtElem x("2*α + 3", gf25);
  ExtElem x_sq = x * x;
  ExtElem root = x_sq.Sqrt();
  std::cout << "  x        = " << x << "\n";
  std::cout << "  x^2      = " << x_sq << "\n";
  std::cout << "  Sqrt(x^2)= " << root
            << " (check root^2 == x^2: " << (root * root == x_sq) << ")\n";
  std::cout << "  Frobenius(x) = " << x.Frobenius() << "\n";
  std::cout << "  Trace(x)     = " << x.Trace() << "\n";
  std::cout << "  Norm(x)      = " << x.Norm() << "\n";

  const auto min_poly = gf25->MinimalPolynomial(x.RawValue());
  const auto eval_min = EvaluateOverBase(*gf25, min_poly, x.RawValue());
  std::cout << "  MinimalPolynomial(x) = " << gf25->ToString(min_poly)
            << " (degree " << min_poly.Degree()
            << ", m_x(x) == 0: " << gf25->IsZero(eval_min) << ")\n";

  std::cout << "\n=== Seeded Random Elements in GF(7^2) ===\n";
  std::mt19937_64 rng(42);
  ExtElem rand_u(gf49->Random(rng), gf49);
  ExtElem rand_v(gf49->Random(rng), gf49);
  while (rand_v.IsZero()) {
    rand_v = ExtElem(gf49->Random(rng), gf49);
  }
  std::cout << "  u = " << rand_u << ", v = " << rand_v << "\n";
  std::cout << "  u + v = " << (rand_u + rand_v) << "\n";
  std::cout << "  u * v = " << (rand_u * rand_v) << "\n";
  std::cout << "  u / v = " << (rand_u / rand_v) << "\n";
  std::cout << "  Log_g(v) = " << rand_v.Log() << "\n";

  std::cout << "\nExtension field examples completed successfully!\n";
  return 0;
}
