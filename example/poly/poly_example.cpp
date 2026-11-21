#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <vector>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/poly.hpp"

using namespace xg;

int main() {
  std::cout << "XGalois Polynomial Examples\n";
  std::cout << "===========================\n";

  using PrimeField = GFP<uint8_t>;
  using Elem7 = GaloisFieldElement<PrimeField>;
  using Poly7 = PolynomialDense<PrimeField>;

  auto gf7 = std::make_shared<PrimeField>(7);
  auto gf2 = std::make_shared<GF2>();

  std::cout << "\n=== Polynomial Arithmetic & Division over GF(7) ===\n";
  // p1(x) = 3x^2 + 2x + 1, p2(x) = 5x + 4
  Poly7 p1({Elem7(1, gf7), Elem7(2, gf7), Elem7(3, gf7)});
  Poly7 p2({Elem7(4, gf7), Elem7(5, gf7)});

  const auto sum = p1 + p2;
  const auto diff = p1 - p2;
  const auto prod = p1 * p2;
  const auto [quot, rem] = p1.DivRem(p2);
  const auto dp1 = p1.Derivative();

  std::cout << "  p1(x)      = " << p1.ToString() << " (deg " << p1.Degree()
            << ", monic: " << std::boolalpha << p1.IsMonic() << ")\n";
  std::cout << "  p2(x)      = " << p2.ToString() << "\n";
  std::cout << "  p1 + p2    = " << sum.ToString() << "\n";
  std::cout << "  p1 - p2    = " << diff.ToString() << "\n";
  std::cout << "  p1 * p2    = " << prod.ToString() << "\n";
  std::cout << "  p1 / p2    = " << quot.ToString()
            << ", rem = " << rem.ToString()
            << " (check q*p2 + r == p1: " << (quot * p2 + rem == p1) << ")\n";
  std::cout << "  p1'(x)     = " << dp1.ToString() << "\n";
  std::cout << "  p1(2)      = " << p1(Elem7(2, gf7)) << "\n";

  std::cout << "\n=== GCD & Extended GCD (Bézout Identity) over GF(7) ===\n";
  // a(x) = (x + 1)(x + 2), b(x) = (x + 1)(x + 4)
  Poly7 lin1({Elem7(1, gf7), Elem7(1, gf7)});
  Poly7 lin2({Elem7(2, gf7), Elem7(1, gf7)});
  Poly7 lin4({Elem7(4, gf7), Elem7(1, gf7)});
  Poly7 poly_a = lin1 * lin2;
  Poly7 poly_b = lin1 * lin4;

  const auto gcd_ab = utils::PolynomialGcd(poly_a, poly_b);
  const auto lcm_ab = utils::PolynomialLcm(poly_a, poly_b);
  const auto [ext_g, bezout] =
      utils::PolynomialDenseExtendedGcd(poly_a, poly_b);
  const auto& [s_poly, t_poly] = bezout;

  std::cout << "  a(x)       = " << poly_a.ToString() << "\n";
  std::cout << "  b(x)       = " << poly_b.ToString() << "\n";
  std::cout << "  gcd(a, b)  = " << gcd_ab.ToString() << "\n";
  std::cout << "  lcm(a, b)  = " << lcm_ab.ToString() << "\n";
  std::cout << "  ExtendedGcd: s(x) = " << s_poly.ToString()
            << ", t(x) = " << t_poly.ToString() << "\n";
  std::cout << "  Verify s*a + t*b == g: "
            << (s_poly * poly_a + t_poly * poly_b == ext_g) << "\n";

  std::cout << "\n=== Irreducibility & Primitivity Tests over GF(2) ===\n";
  using Elem2 = GaloisFieldElement<GF2>;
  using Poly2 = PolynomialDense<GF2>;

  // x^4 + x + 1 (primitive & irreducible)
  Poly2 prim_poly({Elem2(1, gf2), Elem2(1, gf2), Elem2(0, gf2), Elem2(0, gf2),
                   Elem2(1, gf2)});
  // x^4 + x^3 + x^2 + x + 1 (irreducible, NOT primitive: order 5 != 15)
  Poly2 irred_non_prim({Elem2(1, gf2), Elem2(1, gf2), Elem2(1, gf2),
                        Elem2(1, gf2), Elem2(1, gf2)});
  // x^4 + 1 = (x + 1)^4 (reducible)
  Poly2 red_poly({Elem2(1, gf2), Elem2(0, gf2), Elem2(0, gf2), Elem2(0, gf2),
                  Elem2(1, gf2)});

  std::cout << "  " << prim_poly.ToString()
            << " -> IsIrreducible=" << utils::IsIrreducible(prim_poly)
            << ", IsPrimitive=" << utils::IsPrimitive(prim_poly) << "\n";
  std::cout << "  " << irred_non_prim.ToString()
            << " -> IsIrreducible=" << utils::IsIrreducible(irred_non_prim)
            << ", IsPrimitive=" << utils::IsPrimitive(irred_non_prim) << "\n";
  std::cout << "  " << red_poly.ToString()
            << " -> IsIrreducible=" << utils::IsIrreducible(red_poly)
            << ", IsPrimitive=" << utils::IsPrimitive(red_poly) << "\n";

  std::cout << "\n=== Complete Factorization (Cantor-Zassenhaus) over "
               "GF(7) ===\n";
  // f(x) = 3 * (x + 1)^2 * (x^2 + 1)  (x^2 + 1 is irreducible over GF(7))
  Poly7 quad_irred({Elem7(1, gf7), Elem7(0, gf7), Elem7(1, gf7)});
  Poly7 composite = lin1 * lin1 * quad_irred * Elem7(3, gf7);

  const auto [lead_coeff, factors] = utils::Factor(composite);
  std::cout << "  f(x) = " << composite.ToString() << "\n";
  std::cout << "  Leading coefficient: " << lead_coeff << "\n";
  Poly7 reconstructed = Poly7::One(gf7) * lead_coeff;
  for (const auto& [fac, mult] : factors) {
    std::cout << "    (" << fac.ToString() << ")^" << mult
              << "  [irreducible: " << utils::IsIrreducible(fac) << "]\n";
    reconstructed *= fac.Pow(mult);
  }
  std::cout << "  Product of factors == f(x): " << (reconstructed == composite)
            << "\n";

  std::cout << "\n=== Root Finding over GF(7) ===\n";
  // h(x) = (x - 1)(x - 3)(x - 5)(x^2 + 1)
  Poly7 f_r1({Elem7(6, gf7), Elem7(1, gf7)});  // x - 1
  Poly7 f_r3({Elem7(4, gf7), Elem7(1, gf7)});  // x - 3
  Poly7 f_r5({Elem7(2, gf7), Elem7(1, gf7)});  // x - 5
  Poly7 h = f_r1 * f_r3 * f_r5 * quad_irred;

  const auto roots = utils::FindRoots(h);
  std::cout << "  h(x) = " << h.ToString() << "\n";
  std::cout << "  Distinct roots in GF(7): [";
  bool all_roots_zero = true;
  for (size_t i = 0; i < roots.size(); ++i) {
    std::cout << roots[i];
    if (i + 1 < roots.size()) std::cout << ", ";
    if (!h(roots[i]).IsZero()) all_roots_zero = false;
  }
  std::cout << "]  (all satisfy h(r) == 0: " << all_roots_zero << ")\n";

  std::cout << "\n=== Newton Polynomial Interpolation over GF(7) ===\n";
  std::vector<Elem7> xs = {Elem7(0, gf7), Elem7(1, gf7), Elem7(2, gf7)};
  std::vector<Elem7> ys = {p1(xs[0]), p1(xs[1]), p1(xs[2])};
  const auto interpolated = utils::Interpolate(xs, ys);
  std::cout << "  Original p1(x)     = " << p1.ToString() << "\n";
  std::cout << "  Interpolated q(x)  = " << interpolated.ToString() << "\n";
  std::cout << "  Recovered exact p1: " << (interpolated == p1) << "\n";

  std::cout << "\nPolynomial examples completed successfully!\n";
  return 0;
}
