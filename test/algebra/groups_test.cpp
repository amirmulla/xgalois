#include "xgalois/algebra/groups.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "xgalois/algebra/concepts.hpp"
#include "xgalois/algebra/group_algorithms.hpp"
#include "xgalois/algebra/group_base.hpp"
#include "xgalois/algebra/permutation.hpp"
#include "xgalois/algebra/schreier_sims.hpp"
#include "xgalois/utils/bigint.hpp"

namespace {

using xg::AbelianGroup;
using xg::BigInt;
using xg::Group;
using xg::Magma;
using xg::Monoid;
using xg::Semigroup;
using xg::algebra::AlternatingGroup;
using xg::algebra::BsgsChain;
using xg::algebra::BuildBsgsDeterministic;
using xg::algebra::BuildBsgsRandomized;
using xg::algebra::CayleyTable;
using xg::algebra::Center;
using xg::algebra::CommutatorSubgroup;
using xg::algebra::CyclicGroup;
using xg::algebra::DecomposeAbelianGroupBigInt;
using xg::algebra::DerivedSeries;
using xg::algebra::DihedralElement;
using xg::algebra::DihedralGroup;
using xg::algebra::DirectProductGroup;
using xg::algebra::ElementOrder;
using xg::algebra::FiniteGroup;
using xg::algebra::GroupElement;
using xg::algebra::GroupHomomorphism;
using xg::algebra::IsNormalSubgroup;
using xg::algebra::IsSolvable;
using xg::algebra::LeftCosets;
using xg::algebra::Permutation;
using xg::algebra::PermutationGroup;
using xg::algebra::QuotientGroup;
using xg::algebra::RightCosets;
using xg::algebra::SemidirectProductGroup;
using xg::algebra::Subgroup;
using xg::algebra::SylowSubgroup;
using xg::algebra::SymmetricGroup;
using xg::algebra::VerifyFirstIsomorphismTheorem;

// ============================================================================
// Compile-time Concept Checks
// ============================================================================

struct InvertActionZ4 {
  uint64_t operator()(uint64_t h, uint64_t n) const noexcept {
    return (h % 2 == 0) ? (n % 4) : ((4 - (n % 4)) % 4);
  }
};

static_assert(Magma<CyclicGroup>);
static_assert(Semigroup<CyclicGroup>);
static_assert(Monoid<CyclicGroup>);
static_assert(Group<CyclicGroup>);
static_assert(AbelianGroup<CyclicGroup>);
static_assert(FiniteGroup<CyclicGroup>);

static_assert(Group<SymmetricGroup>);
static_assert(!AbelianGroup<SymmetricGroup>);
static_assert(FiniteGroup<SymmetricGroup>);

static_assert(Group<AlternatingGroup>);
static_assert(!AbelianGroup<AlternatingGroup>);
static_assert(FiniteGroup<AlternatingGroup>);

static_assert(Group<DihedralGroup>);
static_assert(!AbelianGroup<DihedralGroup>);
static_assert(FiniteGroup<DihedralGroup>);

static_assert(Group<PermutationGroup>);
static_assert(!AbelianGroup<PermutationGroup>);
static_assert(FiniteGroup<PermutationGroup>);

using AbelianDirect = DirectProductGroup<CyclicGroup, CyclicGroup>;
using NonAbelianDirect = DirectProductGroup<CyclicGroup, SymmetricGroup>;
static_assert(Group<AbelianDirect>);
static_assert(AbelianGroup<AbelianDirect>);
static_assert(FiniteGroup<AbelianDirect>);
static_assert(Group<NonAbelianDirect>);
static_assert(!AbelianGroup<NonAbelianDirect>);
static_assert(FiniteGroup<NonAbelianDirect>);

using SemidirectD4 =
    SemidirectProductGroup<CyclicGroup, CyclicGroup, InvertActionZ4>;
static_assert(Group<SemidirectD4>);
static_assert(!AbelianGroup<SemidirectD4>);
static_assert(FiniteGroup<SemidirectD4>);

static_assert(Group<QuotientGroup<SymmetricGroup>>);
static_assert(!AbelianGroup<QuotientGroup<SymmetricGroup>>);
static_assert(FiniteGroup<QuotientGroup<SymmetricGroup>>);
static_assert(AbelianGroup<QuotientGroup<CyclicGroup>>);

// ============================================================================
// GroupElement Wrapper Tests
// ============================================================================

TEST(GroupElementTest, ArithmeticNegativePowAndOrder) {
  auto z12 = std::make_shared<CyclicGroup>(12);
  GroupElement<CyclicGroup> a(5, z12);
  GroupElement<CyclicGroup> b(9, z12);

  EXPECT_EQ((a * b).RawValue(), 2u);
  EXPECT_EQ((a / b).RawValue(), 8u);
  EXPECT_EQ(a.Inv().RawValue(), 7u);
  EXPECT_EQ((a * a.Inv()), GroupElement<CyclicGroup>::Identity(z12));
  EXPECT_TRUE((a * a.Inv()).IsIdentity());

  EXPECT_EQ(a.Pow(0).RawValue(), 0u);
  EXPECT_EQ(a.Pow(3).RawValue(), 3u);
  EXPECT_EQ(a.Pow(-1).RawValue(), 7u);
  EXPECT_EQ(a.Pow(-3).RawValue(), 9u);
  EXPECT_EQ((a.Pow(5) * a.Pow(-5)), GroupElement<CyclicGroup>::Identity(z12));

  EXPECT_EQ(a.Order(), 12u);
  EXPECT_EQ(GroupElement<CyclicGroup>(4, z12).Order(), 3u);
  EXPECT_EQ(GroupElement<CyclicGroup>::Identity(z12).Order(), 1u);

  std::ostringstream oss;
  oss << a;
  EXPECT_EQ(oss.str(), "5");
  EXPECT_EQ(a.ToString(), "5");

  auto other_z12 = std::make_shared<CyclicGroup>(12);
  GroupElement<CyclicGroup> c(5, other_z12);
  EXPECT_NE(a, c);
  EXPECT_THROW((void)(a * c), std::invalid_argument);
  EXPECT_THROW((GroupElement<CyclicGroup>(1u, nullptr)), std::invalid_argument);

  // Non-abelian GroupElement on S_4
  auto s4 = std::make_shared<SymmetricGroup>(4);
  GroupElement<SymmetricGroup> p1(s4, Permutation::FromCycles(4, {{0, 1, 2}}));
  GroupElement<SymmetricGroup> p2(s4, Permutation::FromCycles(4, {{0, 1}}));
  EXPECT_NE(p1 * p2, p2 * p1);
  EXPECT_EQ(p1 / p2, p1 * p2.Inv());
  EXPECT_EQ(p1.Order(), 3u);
  EXPECT_EQ(p2.Order(), 2u);
  EXPECT_EQ(p1.Pow(-1), p1.Inv());
  EXPECT_EQ(p1.Pow(-2), p1);
}

// ============================================================================
// Permutation Tests (Left-Action Composition)
// ============================================================================

TEST(PermutationTest, LeftActionCompositionCyclesInverseOrderSign) {
  // a = (0 1 2), b = (0 1) in S_4
  const Permutation a = Permutation::FromCycles(4, {{0, 1, 2}});
  const Permutation b = Permutation::FromCycles(4, {{0, 1}});

  // Left-action: (a * b)(i) = a(b(i))
  // 0 -> b(0)=1 -> a(1)=2
  // 1 -> b(1)=0 -> a(0)=1
  // 2 -> b(2)=2 -> a(2)=0
  // 3 -> b(3)=3 -> a(3)=3
  // So a * b = (0 2)
  const Permutation ab = a * b;
  EXPECT_EQ(ab(0), 2u);
  EXPECT_EQ(ab(1), 1u);
  EXPECT_EQ(ab(2), 0u);
  EXPECT_EQ(ab(3), 3u);
  EXPECT_EQ(ab, Permutation::FromCycles(4, {{0, 2}}));

  // Verify non-commutativity: b * a = (1 2) != (0 2)
  const Permutation ba = b * a;
  EXPECT_EQ(ba, Permutation::FromCycles(4, {{1, 2}}));
  EXPECT_NE(ab, ba);

  // Inverse and Pow (including negative exponents)
  EXPECT_EQ(a * a.Inverse(), Permutation::Identity(4));
  EXPECT_EQ(a.Inverse() * a, Permutation::Identity(4));
  EXPECT_EQ(a.Pow(0), Permutation::Identity(4));
  EXPECT_EQ(a.Pow(3), Permutation::Identity(4));
  EXPECT_EQ(a.Pow(-1), a.Inverse());
  EXPECT_EQ(a.Pow(-2), a);
  EXPECT_EQ(a.Pow(100), a);
  EXPECT_EQ(a.Pow(-100), a.Inverse());

  // Cycles, CycleType, Order, Sign, ToString
  const Permutation p = Permutation::FromCycles(7, {{0, 2, 4}, {1, 5}});
  EXPECT_EQ(p.Cycles(),
            (std::vector<std::vector<uint32_t>>{{0, 2, 4}, {1, 5}}));
  EXPECT_EQ(p.CycleType(), (std::vector<uint32_t>{3, 2, 1, 1}));
  EXPECT_EQ(p.Order(), BigInt(6));
  EXPECT_EQ(p.Sign(), -1);
  EXPECT_FALSE(p.IsEven());
  EXPECT_EQ(p.ToString(), "(0 2 4)(1 5)");
  EXPECT_EQ(Permutation::Identity(5).ToString(), "()");
  EXPECT_EQ(Permutation::Identity(5).Sign(), 1);
  EXPECT_TRUE(Permutation::Identity(5).IsEven());

  // BigInt Order via LCM of disjoint prime cycles in S_50:
  // 2 + 3 + 5 + 7 + 11 + 13 = 41 <= 50; LCM = 30030
  std::vector<std::vector<uint32_t>> prime_cycles;
  uint32_t cursor = 0;
  for (uint32_t len : {2u, 3u, 5u, 7u, 11u, 13u}) {
    std::vector<uint32_t> c(len);
    std::iota(c.begin(), c.end(), cursor);
    cursor += len;
    prime_cycles.push_back(std::move(c));
  }
  const Permutation big_ord_perm = Permutation::FromCycles(50, prime_cycles);
  EXPECT_EQ(big_ord_perm.Order(), BigInt(30030));

  // Validation errors
  EXPECT_THROW(Permutation::FromMapping({0, 1, 1}), std::invalid_argument);
  EXPECT_THROW(Permutation::FromMapping({0, 3, 1}), std::invalid_argument);
  EXPECT_THROW(Permutation::FromCycles(3, {{0, 3}}), std::invalid_argument);
  EXPECT_THROW(Permutation::FromCycles(3, {{0, 1, 0}}), std::invalid_argument);

  // Hash & ordering
  std::unordered_set<Permutation> pset{a, b, ab, ba, a};
  EXPECT_EQ(pset.size(), 4u);
  EXPECT_LT(Permutation::Identity(4), b);
}

// ============================================================================
// BSGS (Schreier-Sims) Tests: S_4, S_5, A_5, S_50, A_50
// ============================================================================

TEST(BsgsTest, SmallGroupsExactOrderAndMembership) {
  // S_4
  const BsgsChain s4 =
      BuildBsgsDeterministic(4, {Permutation::FromCycles(4, {{0, 1}}),
                                 Permutation::FromCycles(4, {{0, 1, 2, 3}})});
  EXPECT_EQ(s4.Order(), BigInt(24));
  EXPECT_TRUE(s4.Contains(Permutation::FromCycles(4, {{1, 3}})));
  const auto s4_elems = s4.EnumerateElements();
  EXPECT_EQ(s4_elems.size(), 24u);
  std::unordered_set<Permutation> s4_unique(s4_elems.begin(), s4_elems.end());
  EXPECT_EQ(s4_unique.size(), 24u);
  const PermutationGroup s4_from_bsgs = PermutationGroup::FromBsgs(s4);
  EXPECT_EQ(s4_from_bsgs.Order(), BigInt(24));

  // S_5
  const BsgsChain s5 = BuildBsgsDeterministic(
      5, {Permutation::FromCycles(5, {{0, 1}}),
          Permutation::FromCycles(5, {{0, 1, 2, 3, 4}})});
  EXPECT_EQ(s5.Order(), BigInt(120));

  // A_5 (3-cycles)
  const BsgsChain a5 =
      BuildBsgsDeterministic(5, {Permutation::FromCycles(5, {{0, 1, 2}}),
                                 Permutation::FromCycles(5, {{0, 1, 3}}),
                                 Permutation::FromCycles(5, {{0, 1, 4}})});
  EXPECT_EQ(a5.Order(), BigInt(60));
  EXPECT_TRUE(a5.Contains(Permutation::FromCycles(5, {{0, 1, 2, 3, 4}})));
  EXPECT_TRUE(a5.Contains(Permutation::FromCycles(5, {{0, 1}, {2, 3}})));
  EXPECT_FALSE(a5.Contains(Permutation::FromCycles(5, {{0, 1}})));
  EXPECT_FALSE(a5.Sift(Permutation::FromCycles(5, {{0, 1}})).IsIdentity());

  // Trivial group
  const BsgsChain triv = BuildBsgsDeterministic(5, {});
  EXPECT_EQ(triv.Order(), BigInt(1));
  EXPECT_TRUE(triv.Contains(Permutation::Identity(5)));
  EXPECT_FALSE(triv.Contains(Permutation::FromCycles(5, {{0, 1}})));
}

TEST(BsgsTest, SymmetricGroup50DeterministicAndPerformance) {
  constexpr uint32_t kN = 50;
  std::vector<uint32_t> n_cycle(kN);
  std::iota(n_cycle.begin(), n_cycle.end(), uint32_t{0});
  const std::vector<Permutation> gens = {
      Permutation::FromCycles(kN, {{0, 1}}),
      Permutation::FromCycles(kN, {n_cycle}),
  };

  const BsgsChain bsgs = BuildBsgsDeterministic(kN, gens);
  EXPECT_EQ(bsgs.Order(), BigInt::Factorial(kN));

  // Fast membership and random element sifting (< 100 ms)
  std::mt19937_64 rng(12345);
  const auto t0 = std::chrono::steady_clock::now();
  for (int trial = 0; trial < 200; ++trial) {
    const Permutation r = bsgs.RandomElement(rng);
    EXPECT_TRUE(bsgs.Sift(r).IsIdentity());
    EXPECT_TRUE(bsgs.Contains(r));
  }
  const auto t1 = std::chrono::steady_clock::now();
  const auto elapsed_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
#if defined(__OPTIMIZE__) && !defined(__SANITIZE_ADDRESS__) && \
    !defined(__SANITIZE_THREAD__)
  EXPECT_LT(elapsed_ms, 500);
#else
  EXPECT_LT(elapsed_ms, 5000);
#endif
}

TEST(BsgsTest, AlternatingGroup50AndRandomizedS50) {
  constexpr uint32_t kN = 50;
  std::vector<Permutation> a50_gens;
  for (uint32_t k = 2; k < kN; ++k) {
    a50_gens.push_back(Permutation::FromCycles(kN, {{0, 1, k}}));
  }
  const BsgsChain a50 = BuildBsgsDeterministic(kN, a50_gens);
  EXPECT_EQ(a50.Order(), BigInt::Factorial(kN) / BigInt(2));
  EXPECT_FALSE(a50.Contains(Permutation::FromCycles(kN, {{0, 1}})));
  EXPECT_TRUE(a50.Contains(Permutation::FromCycles(kN, {{10, 20}, {30, 40}})));

  std::vector<uint32_t> n_cycle(kN);
  std::iota(n_cycle.begin(), n_cycle.end(), uint32_t{0});
  const BsgsChain s50_rand =
      BuildBsgsRandomized(kN,
                          {Permutation::FromCycles(kN, {{0, 1}}),
                           Permutation::FromCycles(kN, {n_cycle})},
                          /*seed=*/2026);
  EXPECT_EQ(s50_rand.Order(), BigInt::Factorial(kN));
  std::mt19937_64 rng(777);
  for (int i = 0; i < 50; ++i) {
    EXPECT_TRUE(s50_rand.Contains(s50_rand.RandomElement(rng)));
  }
}

// ============================================================================
// Concrete Group Tests
// ============================================================================

TEST(ConcreteGroupsTest, CyclicDihedralSymmetricAlternatingPermutation) {
  // CyclicGroup
  CyclicGroup z7(7);
  EXPECT_EQ(z7.Order(), BigInt(7));
  EXPECT_EQ(z7.Op(5, 4), 2u);
  EXPECT_EQ(z7.Div(2, 5), 4u);
  EXPECT_EQ(z7.Pow(3, 5), 1u);

  // DihedralGroup D_5: verify s * r^k = r^{-k} * s
  DihedralGroup d5(5);
  EXPECT_EQ(d5.Order(), BigInt(10));
  const DihedralElement r = d5.Rotation(1);
  const DihedralElement s = d5.Reflection(0);
  for (uint64_t k = 0; k < 5; ++k) {
    const DihedralElement rk = d5.Pow(r, k);
    const DihedralElement r_neg_k = d5.Inv(rk);
    EXPECT_TRUE(d5.Equals(d5.Op(s, rk), d5.Op(r_neg_k, s)));
    EXPECT_TRUE(d5.IsIdentity(d5.Pow(d5.Reflection(k), 2)));
  }

  // SymmetricGroup & AlternatingGroup
  SymmetricGroup s4(4);
  EXPECT_EQ(s4.Order(), BigInt(24));
  EXPECT_EQ(s4.Bsgs().Order(), BigInt(24));

  AlternatingGroup a4(4);
  EXPECT_EQ(a4.Order(), BigInt(12));
  EXPECT_EQ(a4.Bsgs().Order(), BigInt(12));
  EXPECT_TRUE(a4.Contains(Permutation::FromCycles(4, {{0, 1, 2}})));
  EXPECT_FALSE(a4.Contains(Permutation::FromCycles(4, {{0, 1}})));
  EXPECT_THROW(a4.Inv(Permutation::FromCycles(4, {{0, 1}})), std::domain_error);

  // PermutationGroup: Klein four-group V_4 <= S_4
  PermutationGroup v4(4, {Permutation::FromCycles(4, {{0, 1}, {2, 3}}),
                          Permutation::FromCycles(4, {{0, 2}, {1, 3}})});
  EXPECT_EQ(v4.Order(), BigInt(4));
  EXPECT_TRUE(v4.Contains(Permutation::FromCycles(4, {{0, 3}, {1, 2}})));
  EXPECT_FALSE(v4.Contains(Permutation::FromCycles(4, {{0, 1}})));
}

TEST(ConcreteGroupsTest, DirectAndSemidirectProducts) {
  // Z_2 x Z_3 ~= Z_6
  DirectProductGroup<CyclicGroup, CyclicGroup> z2xz3(CyclicGroup(2),
                                                     CyclicGroup(3));
  EXPECT_EQ(z2xz3.Order(), BigInt(6));
  EXPECT_EQ(ElementOrder(z2xz3, {1u, 1u}), 6u);

  // D_4 ~= Z_4 |x| Z_2 where 1 in Z_2 acts by inversion on Z_4
  auto semidirect = std::make_shared<SemidirectD4>(
      CyclicGroup(4), CyclicGroup(2), InvertActionZ4{});
  auto d4 = std::make_shared<DihedralGroup>(4);
  EXPECT_EQ(semidirect->Order(), BigInt(8));

  // Verify group axioms (especially inverse) on all elements of semidirect
  const auto tbl = CayleyTable(*semidirect);
  ASSERT_EQ(tbl.elements.size(), 8u);
  for (const auto& x : tbl.elements) {
    EXPECT_TRUE(semidirect->IsIdentity(semidirect->Op(x, semidirect->Inv(x))));
    EXPECT_TRUE(semidirect->IsIdentity(semidirect->Op(semidirect->Inv(x), x)));
  }

  // Construct explicit isomorphism D_4 -> Z_4 |x| Z_2: r^k s^b |-> (k, b)
  GroupHomomorphism<DihedralGroup, SemidirectD4> iso(
      d4, semidirect, [](const DihedralElement& e) {
        return std::pair<uint64_t, uint64_t>{e.rot % 4, e.flip ? 1u : 0u};
      });
  EXPECT_TRUE(iso.VerifyHomomorphism());
  EXPECT_TRUE(iso.IsIsomorphism());
  EXPECT_TRUE(VerifyFirstIsomorphismTheorem(iso));
}

// ============================================================================
// Subgroup, Cosets, Normality, QuotientGroup Tests
// ============================================================================

TEST(GroupAlgorithmsTest, SubgroupCosetsNormalityAndQuotientS4ModV4) {
  auto s4 = std::make_shared<SymmetricGroup>(4);
  Subgroup<SymmetricGroup> v4(s4,
                              {Permutation::FromCycles(4, {{0, 1}, {2, 3}}),
                               Permutation::FromCycles(4, {{0, 2}, {1, 3}})});

  EXPECT_EQ(v4.Order(), 4u);
  EXPECT_EQ(v4.Index(), BigInt(6));
  EXPECT_TRUE(v4.Contains(Permutation::FromCycles(4, {{0, 3}, {1, 2}})));
  EXPECT_TRUE(IsNormalSubgroup(*s4, v4));

  // Left and right cosets of V_4 in S_4 have size 4 and count 6
  const auto lcosets = LeftCosets(*s4, v4);
  const auto rcosets = RightCosets(*s4, v4);
  EXPECT_EQ(lcosets.size(), 6u);
  EXPECT_EQ(rcosets.size(), 6u);
  for (const auto& c : lcosets) {
    EXPECT_EQ(c.size(), 4u);
  }

  // QuotientGroup S_4 / V_4 ~= S_3 (order 6, non-abelian)
  QuotientGroup<SymmetricGroup> q(s4, v4);
  EXPECT_EQ(q.Order(), BigInt(6));
  const auto q_tbl = CayleyTable(q);
  ASSERT_EQ(q_tbl.elements.size(), 6u);

  // Order profile of S_3: one order-1, three order-2, two order-3
  std::vector<uint64_t> orders;
  for (const auto& rep : q_tbl.elements) {
    orders.push_back(ElementOrder(q, rep));
  }
  std::sort(orders.begin(), orders.end());
  EXPECT_EQ(orders, (std::vector<uint64_t>{1, 2, 2, 2, 3, 3}));

  // Non-normal subgroup in S_3: <(0 1)>
  auto s3 = std::make_shared<SymmetricGroup>(3);
  Subgroup<SymmetricGroup> h2(s3, {Permutation::FromCycles(3, {{0, 1}})});
  EXPECT_EQ(h2.Order(), 2u);
  EXPECT_EQ(h2.Index(), BigInt(3));
  EXPECT_FALSE(IsNormalSubgroup(*s3, h2));
  EXPECT_THROW((QuotientGroup<SymmetricGroup>(s3, h2)), std::invalid_argument);

  // Verify left cosets != right cosets as partitions for non-normal subgroup
  auto to_sorted_partition = [](std::vector<std::vector<Permutation>> cosets) {
    for (auto& c : cosets) {
      std::sort(c.begin(), c.end());
    }
    std::sort(cosets.begin(), cosets.end());
    return cosets;
  };
  EXPECT_NE(to_sorted_partition(LeftCosets(*s3, h2)),
            to_sorted_partition(RightCosets(*s3, h2)));
  EXPECT_EQ(to_sorted_partition(lcosets), to_sorted_partition(rcosets));
}

// ============================================================================
// Center, CommutatorSubgroup, DerivedSeries, IsSolvable Tests
// ============================================================================

TEST(GroupAlgorithmsTest, CenterCommutatorAndSolvability) {
  // Center
  EXPECT_EQ(Center(DihedralGroup(4)).Order(), 2u);
  EXPECT_EQ(Center(DihedralGroup(5)).Order(), 1u);
  EXPECT_EQ(Center(SymmetricGroup(4)).Order(), 1u);
  EXPECT_EQ(Center(CyclicGroup(6)).Order(), 6u);

  // CommutatorSubgroup: [S_4, S_4] = A_4 (order 12)
  auto s4 = std::make_shared<SymmetricGroup>(4);
  const Subgroup<SymmetricGroup> comm_s4 = CommutatorSubgroup(s4);
  EXPECT_EQ(comm_s4.Order(), 12u);
  for (const auto& p : comm_s4.Elements()) {
    EXPECT_TRUE(p.IsEven());
  }

  // DerivedSeries & IsSolvable
  const auto s4_series = DerivedSeries(*s4);
  ASSERT_EQ(s4_series.size(), 4u);
  EXPECT_EQ(s4_series[0].Order(), 24u);
  EXPECT_EQ(s4_series[1].Order(), 12u);
  EXPECT_EQ(s4_series[2].Order(), 4u);
  EXPECT_EQ(s4_series[3].Order(), 1u);
  EXPECT_TRUE(IsSolvable(*s4));

  // A_5 is simple non-abelian -> [A_5, A_5] = A_5, not solvable!
  AlternatingGroup a5(5);
  EXPECT_EQ(CommutatorSubgroup(a5).Order(), 60u);
  const auto a5_series = DerivedSeries(a5);
  ASSERT_EQ(a5_series.size(), 1u);
  EXPECT_EQ(a5_series[0].Order(), 60u);
  EXPECT_FALSE(IsSolvable(a5));

  // S_5 is not solvable; D_6 is solvable
  EXPECT_FALSE(IsSolvable(SymmetricGroup(5)));
  EXPECT_TRUE(IsSolvable(DihedralGroup(6)));
}

// ============================================================================
// Exact SylowSubgroup Tests
// ============================================================================

TEST(GroupAlgorithmsTest, SylowSubgroupsExact) {
  // S_4 (order 24 = 2^3 * 3)
  SymmetricGroup s4(4);
  const auto s4_syl2 = SylowSubgroup(s4, 2);
  const auto s4_syl3 = SylowSubgroup(s4, 3);
  const auto s4_syl5 = SylowSubgroup(s4, 5);
  EXPECT_EQ(s4_syl2.Order(), 8u);
  EXPECT_EQ(s4_syl3.Order(), 3u);
  EXPECT_EQ(s4_syl5.Order(), 1u);

  // Verify closure of Sylow subgroups
  for (const auto& a : s4_syl2.Elements()) {
    for (const auto& b : s4_syl2.Elements()) {
      EXPECT_TRUE(s4_syl2.Contains(s4.Op(a, b)));
    }
  }

  // A_5 (order 60 = 2^2 * 3 * 5)
  AlternatingGroup a5(5);
  EXPECT_EQ(SylowSubgroup(a5, 2).Order(), 4u);
  EXPECT_EQ(SylowSubgroup(a5, 3).Order(), 3u);
  EXPECT_EQ(SylowSubgroup(a5, 5).Order(), 5u);

  // Abelian group Z_12 x Z_6 (order 72 = 2^3 * 3^2) exercises center quotient!
  DirectProductGroup<CyclicGroup, CyclicGroup> z12xz6(CyclicGroup(12),
                                                      CyclicGroup(6));
  EXPECT_EQ(SylowSubgroup(z12xz6, 2).Order(), 8u);
  EXPECT_EQ(SylowSubgroup(z12xz6, 3).Order(), 9u);

  // D_6 (order 12 = 2^2 * 3)
  DihedralGroup d6(6);
  EXPECT_EQ(SylowSubgroup(d6, 2).Order(), 4u);
  EXPECT_EQ(SylowSubgroup(d6, 3).Order(), 3u);

  // S_5 (order 120 = 2^3 * 3 * 5)
  SymmetricGroup s5(5);
  const auto s5_syl2 = SylowSubgroup(s5, 2);
  EXPECT_EQ(s5_syl2.Order(), 8u);
  for (const auto& elem : s5_syl2.Elements()) {
    const uint64_t ord = ElementOrder(s5, elem);
    EXPECT_EQ(8u % ord, 0u);
  }
  EXPECT_EQ(SylowSubgroup(s5, 3).Order(), 3u);
  EXPECT_EQ(SylowSubgroup(s5, 5).Order(), 5u);

  // Non-prime p throws std::invalid_argument
  EXPECT_THROW(SylowSubgroup(s4, 1), std::invalid_argument);
  EXPECT_THROW(SylowSubgroup(s4, 4), std::invalid_argument);
  EXPECT_THROW(SylowSubgroup(s4, 6), std::invalid_argument);
}

// ============================================================================
// Homomorphism & First Isomorphism Theorem Tests
// ============================================================================

TEST(GroupAlgorithmsTest, HomomorphismSignMapAndFirstIsomorphismTheorem) {
  auto s4 = std::make_shared<SymmetricGroup>(4);
  auto z2 = std::make_shared<CyclicGroup>(2);

  // Sign map: sgn: S_4 -> Z_2 (even -> 0, odd -> 1)
  GroupHomomorphism<SymmetricGroup, CyclicGroup> sgn(
      s4, z2,
      [](const Permutation& p) -> uint64_t { return p.IsEven() ? 0u : 1u; });
  EXPECT_TRUE(sgn.VerifyHomomorphism());
  EXPECT_EQ(sgn.Kernel().Order(), 12u);
  EXPECT_EQ(sgn.Image().Order(), 2u);
  EXPECT_FALSE(sgn.IsInjective());
  EXPECT_TRUE(sgn.IsSurjective());
  EXPECT_FALSE(sgn.IsIsomorphism());
  EXPECT_TRUE(VerifyFirstIsomorphismTheorem(sgn));

  // Same sign map constructed via FromGeneratorImages:
  // (0 1) -> 1 (odd), (0 1 2 3) -> 1 (odd, 4-cycle is odd!)
  auto sgn_from_gens =
      GroupHomomorphism<SymmetricGroup, CyclicGroup>::FromGeneratorImages(
          s4, z2, std::vector<uint64_t>{1u, 1u});
  EXPECT_TRUE(sgn_from_gens.VerifyHomomorphism());
  EXPECT_EQ(sgn_from_gens.Kernel(), sgn.Kernel());

  // Invalid generator mapping: (0 1) -> 0 and (0 1 2 3) -> 1 fails because
  // (0 1 2 3) is conjugate to products of transpositions.
  EXPECT_THROW(
      (GroupHomomorphism<SymmetricGroup, CyclicGroup>::FromGeneratorImages(
          s4, z2, std::vector<uint64_t>{0u, 1u})),
      std::invalid_argument);

  // Invalid callable (not a homomorphism)
  GroupHomomorphism<CyclicGroup, CyclicGroup> bad_map(
      z2, z2, [](uint64_t x) -> uint64_t { return (x + 1) % 2; });
  EXPECT_FALSE(bad_map.VerifyHomomorphism());
  EXPECT_FALSE(VerifyFirstIsomorphismTheorem(bad_map));
}

// ============================================================================
// CayleyTable & BigInt SNF Abelian Group Decomposition Tests
// ============================================================================

TEST(GroupAlgorithmsTest, CayleyTableAndBigIntSnfDecomposition) {
  DihedralGroup d3(3);
  const auto res = CayleyTable(d3);
  ASSERT_EQ(res.elements.size(), 6u);
  ASSERT_EQ(res.table.size(), 6u);
  EXPECT_TRUE(d3.IsIdentity(res.elements[0]));

  // Verify Latin square property and identity row/column
  for (std::size_t i = 0; i < 6; ++i) {
    EXPECT_EQ(res.table[0][i], i);
    EXPECT_EQ(res.table[i][0], i);
    std::set<std::size_t> row_vals(res.table[i].begin(), res.table[i].end());
    EXPECT_EQ(row_vals.size(), 6u);
    std::set<std::size_t> col_vals;
    for (std::size_t r = 0; r < 6; ++r) {
      col_vals.insert(res.table[r][i]);
    }
    EXPECT_EQ(col_vals.size(), 6u);
  }

  // BigInt SNF decomposition: Z^3 / rowspace([[2, 4, 4], [-6, 6, 12], [10, -4,
  // -16]]) Known SNF diagonal: (2, 6, 12) -> Z/2Z (+) Z/6Z (+) Z/12Z
  const std::vector<std::vector<BigInt>> rel = {
      {BigInt(2), BigInt(4), BigInt(4)},
      {BigInt(-6), BigInt(6), BigInt(12)},
      {BigInt(10), BigInt(-4), BigInt(-16)},
  };
  const auto decomp = DecomposeAbelianGroupBigInt(rel);
  EXPECT_EQ(decomp.free_rank, 0u);
  EXPECT_EQ(decomp.invariant_factors,
            (std::vector<BigInt>{BigInt(2), BigInt(6), BigInt(12)}));
  EXPECT_EQ(decomp.elementary_divisors,
            (std::vector<std::pair<BigInt, uint64_t>>{
                {BigInt(2), 1},
                {BigInt(2), 1},
                {BigInt(2), 2},
                {BigInt(3), 1},
                {BigInt(3), 1},
            }));

  // BigInt SNF with free rank and > 64-bit entries:
  // diag(2^80, 0) on Z^2 -> Z (+) Z/(2^80)Z
  const BigInt big_pow2 = BigInt(1) << 80;
  const std::vector<std::vector<BigInt>> big_rel = {
      {big_pow2, BigInt(0)},
  };
  const auto big_decomp = DecomposeAbelianGroupBigInt(big_rel);
  EXPECT_EQ(big_decomp.free_rank, 1u);
  ASSERT_EQ(big_decomp.invariant_factors.size(), 1u);
  EXPECT_EQ(big_decomp.invariant_factors[0], big_pow2);
  EXPECT_EQ(big_decomp.elementary_divisors,
            (std::vector<std::pair<BigInt, uint64_t>>{{BigInt(2), 80}}));
}

}  // namespace
