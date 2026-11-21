#include <cstdint>
#include <iostream>
#include <memory>

#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"

using namespace xg;

int main() {
  std::cout << "XGalois Prime Field Examples\n";
  std::cout << "============================\n";

  std::cout << "\n=== Creating Prime Fields: Barrett GFP & Table GFPLOG ===\n";
  auto gf11 = std::make_shared<GFP<uint8_t>>(11);
  auto gf11_log = std::make_shared<GFPLOG<uint8_t>>(11);
  auto gf17 = std::make_shared<GFP<uint8_t>>(17);
  auto gf_mersenne = std::make_shared<GFP<uint32_t>>(2147483647U);  // 2^31 - 1

  std::cout << "  GFP<uint8_t>(11)        -> Order: " << gf11->Order()
            << ", Generator: "
            << static_cast<int>(gf11->MultiplicativeGenerator()) << "\n";
  std::cout << "  GFPLOG<uint8_t>(11)     -> Order: " << gf11_log->Order()
            << ", Generator: "
            << static_cast<int>(gf11_log->MultiplicativeGenerator()) << "\n";
  std::cout << "  GFP<uint8_t>(17)        -> Order: " << gf17->Order()
            << ", Generator: "
            << static_cast<int>(gf17->MultiplicativeGenerator()) << "\n";
  std::cout << "  GFP<uint32_t>(2^31 - 1) -> Order: " << gf_mersenne->Order()
            << ", Generator: " << gf_mersenne->MultiplicativeGenerator()
            << "\n";

  std::cout << "\n=== Arithmetic & PowSigned in GF(11) ===\n";
  using Elem11 = GaloisFieldElement<GFP<uint8_t>>;
  Elem11 a(3, gf11);
  Elem11 b(5, gf11);
  Elem11 c(7, gf11);

  std::cout << "  a = " << a << ", b = " << b << ", c = " << c << "\n";
  std::cout << "  a + b = " << (a + b) << ", a - b = " << (a - b) << "\n";
  std::cout << "  a * b = " << (a * b) << ", a / b = " << (a / b) << "\n";
  std::cout << "  (a + b) * c - a * b = " << ((a + b) * c - a * b) << "\n";

  const Elem11 a_pow3 = a.Pow(3);
  const Elem11 a_pow_neg3 = a.PowSigned(-3);
  std::cout << "  a^3 = " << a_pow3 << ", a^(-3) = " << a_pow_neg3
            << " (check a^3 * a^(-3) = " << (a_pow3 * a_pow_neg3) << ")\n";

  std::cout << "\n=== Comparing Barrett GFP vs Table GFPLOG in GF(11) ===\n";
  using Elem11Log = GaloisFieldElement<GFPLOG<uint8_t>>;
  Elem11Log a_tbl(3, gf11_log);
  Elem11Log b_tbl(5, gf11_log);
  const bool backend_match =
      ((a + b).Value() == (a_tbl + b_tbl).Value()) &&
      ((a * b).Value() == (a_tbl * b_tbl).Value()) &&
      ((a / b).Value() == (a_tbl / b_tbl).Value()) &&
      (a.PowSigned(-4).Value() == a_tbl.PowSigned(-4).Value());
  std::cout << "  GFP and GFPLOG agree on +, *, /, PowSigned: "
            << std::boolalpha << backend_match << "\n";

  std::cout << "\n=== Square Roots (p == 3 mod 4 and Tonelli-Shanks p == 1 mod "
               "4) ===\n";
  // GF(11): p = 11 == 3 (mod 4)
  Elem11 q11(4, gf11);
  Elem11 r11 = q11.Sqrt();
  std::cout << "  GF(11): Sqrt(4) = " << r11
            << " (check r^2 == 4: " << (r11 * r11 == q11) << ")\n";

  // GF(17): p = 17 == 1 (mod 4) -> Tonelli-Shanks
  using Elem17 = GaloisFieldElement<GFP<uint8_t>>;
  Elem17 q17(2, gf17);
  Elem17 r17 = q17.Sqrt();
  std::cout << "  GF(17): Sqrt(2) = " << r17
            << " (check r^2 == 2: " << (r17 * r17 == q17) << ")\n";

  std::cout << "\n=== Discrete Logarithms (BSGS / Table) ===\n";
  Elem17 g17(gf17->MultiplicativeGenerator(), gf17);
  Elem17 target17(13, gf17);
  const uint64_t log17 = target17.Log();
  std::cout << "  GF(17): Log_" << g17 << "(" << target17 << ") = " << log17
            << " (check g^log == target: " << (g17.Pow(log17) == target17)
            << ")\n";

  Elem11Log g11_tbl(gf11_log->MultiplicativeGenerator(), gf11_log);
  const uint64_t log11 = a_tbl.Log();
  std::cout << "  GF(11) [GFPLOG]: Log_" << g11_tbl << "(" << a_tbl
            << ") = " << log11
            << " (check g^log == a: " << (g11_tbl.Pow(log11) == a_tbl) << ")\n";

  std::cout << "\n=== Large Prime Barrett Arithmetic in GF(2^31 - 1) ===\n";
  using ElemM = GaloisFieldElement<GFP<uint32_t>>;
  ElemM u(123456789U, gf_mersenne);
  ElemM v(987654321U, gf_mersenne);
  const ElemM u_inv = u.Inv();
  std::cout << "  u * v     = " << (u * v) << "\n";
  std::cout << "  u * u^(-1)= " << (u * u_inv) << "\n";
  std::cout << "  u^(-5) * u^5 = " << (u.PowSigned(-5) * u.Pow(5)) << "\n";

  std::cout << "\n=== All Non-Zero Inverses in GF(11) ===\n";
  for (uint8_t i = 1; i < 11; ++i) {
    Elem11 elem(i, gf11);
    Elem11 inv = elem.Inv();
    std::cout << "  " << static_cast<int>(i) << "^(-1) = " << inv
              << " (check: " << (elem * inv) << ")\n";
  }

  std::cout << "\nPrime field examples completed successfully!\n";
  return 0;
}
