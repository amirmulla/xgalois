#ifndef XGALOIS_SIMD_DISPATCH_HPP_
#define XGALOIS_SIMD_DISPATCH_HPP_

#include <atomic>
#include <cstdint>
#include <optional>
#include <string_view>

#if defined(__x86_64__) || defined(_M_X64)
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) || defined(__clang__)
#include <cpuid.h>
#endif
#elif defined(__aarch64__) || defined(_M_ARM64)
#if defined(__linux__)
#include <sys/auxv.h>
#if __has_include(<asm/hwcap.h>)
#include <asm/hwcap.h>
#endif
#endif
#endif

namespace xg {
namespace simd {

// Supported instruction-set architecture tiers for xgalois SIMD kernels.
enum class IsaLevel : uint8_t {
  kScalar = 0,
  kSse42Pclmul = 1,
  kAvx2 = 2,
  kAvx512Gfni = 3,
  kNeon = 4,
  kNeonPmull = 5,
};

// Returns a stable human-readable name for the given ISA level.
inline std::string_view IsaLevelName(IsaLevel level) noexcept {
  switch (level) {
    case IsaLevel::kScalar:
      return "Scalar";
    case IsaLevel::kSse42Pclmul:
      return "SSE4.2+PCLMULQDQ";
    case IsaLevel::kAvx2:
      return "AVX2";
    case IsaLevel::kAvx512Gfni:
      return "AVX-512+GFNI+VPCLMULQDQ";
    case IsaLevel::kNeon:
      return "NEON";
    case IsaLevel::kNeonPmull:
      return "NEON+PMULL";
  }
  return "Unknown";
}

namespace detail {

#if defined(__x86_64__) || defined(_M_X64)

inline void CpuidCount(uint32_t leaf, uint32_t subleaf, uint32_t& eax,
                       uint32_t& ebx, uint32_t& ecx, uint32_t& edx) noexcept {
#if defined(_MSC_VER)
  int regs[4] = {0, 0, 0, 0};
  __cpuidex(regs, static_cast<int>(leaf), static_cast<int>(subleaf));
  eax = static_cast<uint32_t>(regs[0]);
  ebx = static_cast<uint32_t>(regs[1]);
  ecx = static_cast<uint32_t>(regs[2]);
  edx = static_cast<uint32_t>(regs[3]);
#elif defined(__GNUC__) || defined(__clang__)
  __cpuid_count(leaf, subleaf, eax, ebx, ecx, edx);
#else
  (void)leaf;
  (void)subleaf;
  eax = ebx = ecx = edx = 0;
#endif
}

inline uint64_t ReadXcr0() noexcept {
#if defined(_MSC_VER)
  return _xgetbv(0);
#elif defined(__GNUC__) || defined(__clang__)
  uint32_t eax = 0;
  uint32_t edx = 0;
  __asm__ __volatile__("xgetbv" : "=a"(eax), "=d"(edx) : "c"(0u));
  return (static_cast<uint64_t>(edx) << 32) | static_cast<uint64_t>(eax);
#else
  return 0;
#endif
}

#endif  // defined(__x86_64__) || defined(_M_X64)

inline IsaLevel DetectIsaLevelUncached() noexcept {
#if defined(__x86_64__) || defined(_M_X64)
  uint32_t eax = 0;
  uint32_t ebx = 0;
  uint32_t ecx = 0;
  uint32_t edx = 0;

  CpuidCount(0u, 0u, eax, ebx, ecx, edx);
  const uint32_t max_leaf = eax;
  if (max_leaf < 1u) {
    return IsaLevel::kScalar;
  }

  CpuidCount(1u, 0u, eax, ebx, ecx, edx);
  const bool has_pclmulqdq = (ecx & (1u << 1)) != 0u;
  const bool has_sse42 = (ecx & (1u << 20)) != 0u;
  const bool has_osxsave = (ecx & (1u << 27)) != 0u;
  const bool has_avx = (ecx & (1u << 28)) != 0u;

  uint64_t xcr0 = 0;
  if (has_osxsave) {
    xcr0 = ReadXcr0();
  }
  // XCR0 bits: XMM (bit 1) + YMM (bit 2) = 0x06.
  const bool os_avx_enabled = has_osxsave && ((xcr0 & 0x6ULL) == 0x6ULL);
  // XCR0 bits: XMM (1) + YMM (2) + OPMASK (5) + ZMM_Hi256 (6) + Hi16_ZMM (7) =
  // 0xE6.
  const bool os_avx512_enabled = has_osxsave && ((xcr0 & 0xE6ULL) == 0xE6ULL);

  if (max_leaf >= 7u) {
    CpuidCount(7u, 0u, eax, ebx, ecx, edx);
    const bool has_avx2 = (ebx & (1u << 5)) != 0u;
    const bool has_avx512f = (ebx & (1u << 16)) != 0u;
    const bool has_avx512bw = (ebx & (1u << 30)) != 0u;
    const bool has_avx512vl = (ebx & (1u << 31)) != 0u;
    const bool has_gfni = (ecx & (1u << 8)) != 0u;
    const bool has_vpclmulqdq = (ecx & (1u << 10)) != 0u;

    if (os_avx512_enabled && has_avx512f && has_avx512bw && has_avx512vl &&
        has_gfni && has_vpclmulqdq) {
      return IsaLevel::kAvx512Gfni;
    }
    if (os_avx_enabled && has_avx && has_avx2) {
      return IsaLevel::kAvx2;
    }
  }

  if (has_sse42 && has_pclmulqdq) {
    return IsaLevel::kSse42Pclmul;
  }
  return IsaLevel::kScalar;

#elif defined(__aarch64__) || defined(_M_ARM64)
#if defined(__linux__)
#if defined(HWCAP_PMULL)
  constexpr unsigned long kHwcapPmull = HWCAP_PMULL;
#else
  constexpr unsigned long kHwcapPmull = (1UL << 4);
#endif
  const unsigned long hwcap = ::getauxval(AT_HWCAP);
  if ((hwcap & kHwcapPmull) != 0UL) {
    return IsaLevel::kNeonPmull;
  }
#elif defined(__ARM_FEATURE_CRYPTO)
  return IsaLevel::kNeonPmull;
#endif
#if defined(__ARM_NEON) || defined(__aarch64__) || defined(_M_ARM64)
  return IsaLevel::kNeon;
#else
  return IsaLevel::kScalar;
#endif

#else
  return IsaLevel::kScalar;
#endif
}

inline std::atomic<int>& ForcedIsaStorage() noexcept {
  static std::atomic<int> forced{-1};
  return forced;
}

}  // namespace detail

// Detects the highest supported ISA level on the current CPU and OS.
// Result is computed once and cached thread-safely.
inline IsaLevel DetectIsaLevel() noexcept {
  static const IsaLevel kDetected = detail::DetectIsaLevelUncached();
  return kDetected;
}

// Forces a specific ISA level for runtime dispatch, or clears any override when
// passed std::nullopt.
inline void SetForcedIsaLevel(std::optional<IsaLevel> forced) noexcept {
  const int value = forced.has_value() ? static_cast<int>(*forced) : -1;
  detail::ForcedIsaStorage().store(value, std::memory_order_release);
}

// Returns the currently forced ISA level, or std::nullopt if auto-detect is
// active.
inline std::optional<IsaLevel> GetForcedIsaLevel() noexcept {
  const int value = detail::ForcedIsaStorage().load(std::memory_order_acquire);
  if (value < 0) {
    return std::nullopt;
  }
  return static_cast<IsaLevel>(value);
}

// Returns the active ISA level used for kernel dispatch (forced override if
// set, otherwise hardware auto-detected level).
inline IsaLevel GetActiveIsaLevel() noexcept {
  const int value = detail::ForcedIsaStorage().load(std::memory_order_acquire);
  if (value >= 0) {
    return static_cast<IsaLevel>(value);
  }
  return DetectIsaLevel();
}

// RAII guard that temporarily sets the forced ISA level and restores the
// previous override state upon destruction.
class ScopedIsaOverride {
 public:
  explicit ScopedIsaOverride(std::optional<IsaLevel> forced) noexcept
      : previous_(GetForcedIsaLevel()) {
    SetForcedIsaLevel(forced);
  }

  ~ScopedIsaOverride() noexcept { SetForcedIsaLevel(previous_); }

  ScopedIsaOverride(const ScopedIsaOverride&) = delete;
  ScopedIsaOverride& operator=(const ScopedIsaOverride&) = delete;
  ScopedIsaOverride(ScopedIsaOverride&&) = delete;
  ScopedIsaOverride& operator=(ScopedIsaOverride&&) = delete;

 private:
  std::optional<IsaLevel> previous_;
};

}  // namespace simd
}  // namespace xg

#endif  // XGALOIS_SIMD_DISPATCH_HPP_
