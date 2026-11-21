#include <bitset>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>

#include "xgalois/field/gf128.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"

using namespace xg;

int main() {
  std::cout << "XGalois Binary Field Examples\n";
  std::cout << "=============================\n";

  std::cout << "\n=== Base Binary Field GF(2) ===\n";
  auto gf2 = std::make_shared<GF2>();
  const uint8_t a = 1;
  const uint8_t b = 1;
  const uint8_t c = 0;
  std::cout << "Order: " << gf2->Order()
            << ", Characteristic: " << gf2->Characteristic() << "\n";
  std::cout << "  1 + 1 = " << static_cast<int>(gf2->Add(a, b)) << "\n";
  std::cout << "  1 * 1 = " << static_cast<int>(gf2->Mul(a, b)) << "\n";
  std::cout << "  1 + 0 = " << static_cast<int>(gf2->Add(a, c)) << "\n";

  std::cout << "\n=== Binary Extension GF(2^4) with GF2X<uint8_t> ===\n";
  auto gf16 = std::make_shared<GF2X<uint8_t>>(4);
  const uint8_t x = 3;  // α + 1
  const uint8_t y = 5;  // α^2 + 1

  const uint8_t sum = gf16->Add(x, y);
  const uint8_t prod = gf16->Mul(x, y);
  const uint8_t quot = gf16->Div(x, y);
  const uint8_t x_inv = gf16->Inv(x);

  std::cout << "  x = " << static_cast<int>(x) << " (" << std::bitset<4>(x)
            << "), y = " << static_cast<int>(y) << " (" << std::bitset<4>(y)
            << ")\n";
  std::cout << "  x + y  = " << static_cast<int>(sum) << " ("
            << std::bitset<4>(sum) << ")\n";
  std::cout << "  x * y  = " << static_cast<int>(prod) << " ("
            << std::bitset<4>(prod) << ")\n";
  std::cout << "  x / y  = " << static_cast<int>(quot) << " ("
            << std::bitset<4>(quot) << ")\n";
  std::cout << "  x^(-1) = " << static_cast<int>(x_inv)
            << " (check x * x^(-1) = " << static_cast<int>(gf16->Mul(x, x_inv))
            << ")\n";

  std::cout << "\n=== Sqrt, Frobenius, Trace, and Norm in GF(2^4) ===\n";
  const uint8_t sx = gf16->Sqrt(x);
  const uint8_t tr_x = gf16->Trace(x);
  const uint8_t nm_x = gf16->Norm(x);
  std::cout << "  Sqrt(x)   = " << static_cast<int>(sx)
            << " (check Sqrt(x)^2 = " << static_cast<int>(gf16->Mul(sx, sx))
            << ")\n";
  std::cout << "  Trace(x)  = " << static_cast<int>(tr_x) << " (in GF(2))\n";
  std::cout << "  Norm(x)   = " << static_cast<int>(nm_x)
            << " (1 for all nonzero elements in GF(2^m))\n";

  std::cout << "\n=== GF(2^8) in uint8_t: GF2X vs GF2XLOG vs GF2XZECH ===\n";
  auto gf256_std = std::make_shared<GF2X<uint8_t>>(8);
  auto gf256_log = std::make_shared<GF2XLOG<uint8_t>>(8);
  auto gf256_zech = std::make_shared<GF2XZECH<uint8_t>>(8);

  std::cout << "  GF(2^8) Order = " << gf256_std->Order()
            << ", Degree = " << gf256_std->Degree()
            << ", sizeof(element_type) = "
            << sizeof(GF2X<uint8_t>::element_type) << " byte\n";

  const uint8_t u_val = 0x57;
  const uint8_t v_val = 0x83;

  GaloisFieldElement<GF2X<uint8_t>> u_std(u_val, gf256_std);
  GaloisFieldElement<GF2X<uint8_t>> v_std(v_val, gf256_std);
  GaloisFieldElement<GF2XLOG<uint8_t>> u_log(u_val, gf256_log);
  GaloisFieldElement<GF2XLOG<uint8_t>> v_log(v_val, gf256_log);
  GaloisFieldElement<GF2XZECH<uint8_t>> u_zech(u_val, gf256_zech);
  GaloisFieldElement<GF2XZECH<uint8_t>> v_zech(v_val, gf256_zech);

  auto print_hex = [](uint8_t val) {
    std::cout << "0x" << std::hex << std::setw(2) << std::setfill('0')
              << static_cast<int>(val) << std::dec;
  };

  std::cout << "  u = ";
  print_hex(u_val);
  std::cout << ", v = ";
  print_hex(v_val);
  std::cout << "\n";

  std::cout << "  GF2X     : u+v=";
  print_hex((u_std + v_std).Value());
  std::cout << ", u*v=";
  print_hex((u_std * v_std).Value());
  std::cout << ", u/v=";
  print_hex((u_std / v_std).Value());
  std::cout << ", sqrt(u)=";
  print_hex(u_std.Sqrt().Value());
  std::cout << ", Tr(u)=" << static_cast<int>(u_std.Trace().Value())
            << ", N(u)=" << static_cast<int>(u_std.Norm().Value()) << "\n";

  std::cout << "  GF2XLOG  : u+v=";
  print_hex((u_log + v_log).Value());
  std::cout << ", u*v=";
  print_hex((u_log * v_log).Value());
  std::cout << ", u/v=";
  print_hex((u_log / v_log).Value());
  std::cout << ", sqrt(u)=";
  print_hex(u_log.Sqrt().Value());
  std::cout << ", Tr(u)=" << static_cast<int>(u_log.Trace().Value())
            << ", N(u)=" << static_cast<int>(u_log.Norm().Value()) << "\n";

  std::cout << "  GF2XZECH : u+v=";
  print_hex((u_zech + v_zech).Value());
  std::cout << ", u*v=";
  print_hex((u_zech * v_zech).Value());
  std::cout << ", u/v=";
  print_hex((u_zech / v_zech).Value());
  std::cout << ", sqrt(u)=";
  print_hex(u_zech.Sqrt().Value());
  std::cout << ", Tr(u)=" << static_cast<int>(u_zech.Trace().Value())
            << ", N(u)=" << static_cast<int>(u_zech.Norm().Value()) << "\n";

  const bool all_match =
      ((u_std + v_std).Value() == (u_log + v_log).Value()) &&
      ((u_std + v_std).Value() == (u_zech + v_zech).Value()) &&
      ((u_std * v_std).Value() == (u_log * v_log).Value()) &&
      ((u_std * v_std).Value() == (u_zech * v_zech).Value()) &&
      ((u_std / v_std).Value() == (u_log / v_log).Value()) &&
      ((u_std / v_std).Value() == (u_zech / v_zech).Value()) &&
      (u_std.Sqrt().Value() == u_log.Sqrt().Value()) &&
      (u_std.Sqrt().Value() == u_zech.Sqrt().Value());
  std::cout << "  All three GF(2^8) backends agree: " << std::boolalpha
            << all_match << "\n";

  std::cout << "\n=== Cryptographic GF(2^128) & AES-GCM GHASH (GF128) ===\n";
  auto gf128 = std::make_shared<GF128>();  // x^128 + x^7 + x^2 + x + 1
  const UInt128 h = UInt128::FromHex("0x66e94bd4ef8a2c3b884cfa59ca342b2e");
  const UInt128 c1 = UInt128::FromHex("0x0388dace60b6a392f328c2b971b2fe78");
  const UInt128 expected_x1 =
      UInt128::FromHex("0x5e2ec746917062882c85b0685353deb7");

  const UInt128 ghash_x1 = gf128->GhashMul(c1, h);
  GaloisFieldElement<GF128> h_elem(h, gf128);
  const auto h_inv = h_elem.Inv();
  const auto h_sqrt = h_elem.Sqrt();

  std::cout << "  GF(2^128) Degree = " << gf128->Degree()
            << ", modulus_low = " << gf128->ModulusLow().ToHex() << "\n";
  std::cout << "  NIST SP 800-38D GHASH(C1, H) = "
            << ghash_x1.ToHex(/*zero_pad=*/true)
            << " (matches test vector: " << (ghash_x1 == expected_x1) << ")\n";
  std::cout << "  H * H^(-1) == 1 : " << (h_elem * h_inv).RawValue().IsOne()
            << "\n";
  std::cout << "  Sqrt(H)^2 == H  : " << (h_sqrt * h_sqrt == h_elem) << "\n";
  std::cout << "  Trace(H) = " << h_elem.Trace().RawValue().ToHex()
            << ", Norm(H) = " << h_elem.Norm().RawValue().ToHex() << "\n";

  std::cout << "\n=== Multiplicative Group of GF(2^4) ===\n";
  const uint8_t gen = gf16->MultiplicativeGenerator();
  uint8_t cur = 1;
  std::cout << "Powers of generator " << static_cast<int>(gen) << ":\n";
  for (int i = 0; i < 15; ++i) {
    std::cout << "  g^" << std::setw(2) << std::setfill(' ') << i << " = "
              << std::setw(2) << static_cast<int>(cur) << " ("
              << std::bitset<4>(cur) << ")\n";
    cur = gf16->Mul(cur, gen);
  }

  std::cout << "\nBinary field examples completed successfully!\n";
  return 0;
}
