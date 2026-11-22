#ifndef XGALOIS_CHANNEL_SIMULATION_HPP_
#define XGALOIS_CHANNEL_SIMULATION_HPP_

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <iomanip>
#include <locale>
#include <numbers>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/channel/rng.hpp"
#include "xgalois/simd/thread_pool.hpp"

namespace xg {
namespace channel {

// Acklam's rational approximation to the standard normal quantile function
// Phi^{-1}(p) for p in (0, 1), achieving relative error < 1.15e-9 everywhere
// (and refined to near machine precision via one Halley iteration in the
// non-extreme domain).
inline double InverseNormalCdf(double p) {
  if (std::isnan(p) || p <= 0.0 || p >= 1.0) {
    throw std::domain_error("InverseNormalCdf: p must lie strictly in (0, 1)");
  }

  constexpr double kA1 = -3.969683028665376e+01;
  constexpr double kA2 = 2.209460984245205e+02;
  constexpr double kA3 = -2.759285104469687e+02;
  constexpr double kA4 = 1.383577518672690e+02;
  constexpr double kA5 = -3.066479806614716e+01;
  constexpr double kA6 = 2.506628277459239e+00;

  constexpr double kB1 = -5.447609879822406e+01;
  constexpr double kB2 = 1.615858368580409e+02;
  constexpr double kB3 = -1.556989798598866e+02;
  constexpr double kB4 = 6.680131188771972e+01;
  constexpr double kB5 = -1.328068155288572e+01;

  constexpr double kC1 = -7.784894002430293e-03;
  constexpr double kC2 = -3.223964580411365e-01;
  constexpr double kC3 = -2.400758277161838e+00;
  constexpr double kC4 = -2.549732539343734e+00;
  constexpr double kC5 = 4.374664141464968e+00;
  constexpr double kC6 = 2.938163982698783e+00;

  constexpr double kD1 = 7.784695709041462e-03;
  constexpr double kD2 = 3.224671290700398e-01;
  constexpr double kD3 = 2.445134137142996e+00;
  constexpr double kD4 = 3.754408661907416e+00;

  constexpr double kPLow = 0.02425;
  constexpr double kPHigh = 1.0 - kPLow;

  double x = 0.0;
  if (p < kPLow) {
    const double q = std::sqrt(-2.0 * std::log(p));
    x = (((((kC1 * q + kC2) * q + kC3) * q + kC4) * q + kC5) * q + kC6) /
        ((((kD1 * q + kD2) * q + kD3) * q + kD4) * q + 1.0);
  } else if (p <= kPHigh) {
    const double q = p - 0.5;
    const double r = q * q;
    x = (((((kA1 * r + kA2) * r + kA3) * r + kA4) * r + kA5) * r + kA6) * q /
        (((((kB1 * r + kB2) * r + kB3) * r + kB4) * r + kB5) * r + 1.0);
  } else {
    const double q = std::sqrt(-2.0 * std::log(1.0 - p));
    x = -(((((kC1 * q + kC2) * q + kC3) * q + kC4) * q + kC5) * q + kC6) /
        ((((kD1 * q + kD2) * q + kD3) * q + kD4) * q + 1.0);
  }

  // One Halley refinement step where exp(x^2 / 2) cannot overflow.
  if (std::abs(x) < 26.0) {
    const double err = 0.5 * std::erfc(-x / std::numbers::sqrt2) - p;
    const double u =
        err * std::sqrt(2.0 * std::numbers::pi) * std::exp(0.5 * x * x);
    x -= u / (1.0 + 0.5 * x * u);
  }

  return x;
}

// Wilson score confidence interval for a binomial proportion successes/trials
// at two-sided confidence level `confidence` (default 95%).
inline std::pair<double, double> WilsonScoreInterval(uint64_t successes,
                                                     uint64_t trials,
                                                     double confidence = 0.95) {
  if (std::isnan(confidence) || confidence <= 0.0 || confidence >= 1.0) {
    throw std::invalid_argument(
        "WilsonScoreInterval: confidence must be strictly in (0, 1)");
  }
  if (successes > trials) {
    throw std::invalid_argument(
        "WilsonScoreInterval: successes cannot exceed trials");
  }
  if (trials == 0) {
    return {0.0, 1.0};
  }

  const double z = InverseNormalCdf(0.5 + 0.5 * confidence);
  const double z2 = z * z;
  const double n = static_cast<double>(trials);
  const double p_hat = static_cast<double>(successes) / n;

  const double denom = 1.0 + z2 / n;
  const double center = (p_hat + 0.5 * z2 / n) / denom;
  const double half_width =
      (z / denom) *
      std::sqrt((p_hat * (1.0 - p_hat)) / n + (0.25 * z2) / (n * n));

  return {std::clamp(center - half_width, 0.0, 1.0),
          std::clamp(center + half_width, 0.0, 1.0)};
}

// Configuration for deterministic parallel Monte Carlo simulations.
// Note: Named ParallelSimulationConfig so this header can be freely
// co-included alongside xgalois/channel/channel.hpp (which defines the legacy
// serial SimulationConfig).
struct ParallelSimulationConfig {
  uint64_t max_frames = 10000;
  uint64_t target_frame_errors = 100;       // 0 disables error-count stop.
  double min_ci_relative_half_width = 0.0;  // 0 disables Wilson-CI stop.
  uint64_t chunk_size = 64;
  uint64_t base_seed = 42;
  std::size_t num_threads = 0;  // 0 = GlobalThreadPool().NumWorkers()
};

struct FrameTrialOutcome {
  bool frame_error = false;
  uint64_t symbols = 0;
  uint64_t symbol_errors = 0;
  double iterations = 0.0;
  double repair_symbols_read = 0.0;
  bool used_local_repair = false;
};

struct SimulationPointResult {
  double param_value = 0.0;
  uint64_t total_frames = 0;
  uint64_t frame_errors = 0;
  uint64_t total_symbols = 0;
  uint64_t symbol_errors = 0;
  double fer = 0.0;
  double ser = 0.0;
  std::pair<double, double> fer_ci95{0.0, 1.0};
  double avg_decoder_iterations = 0.0;
  double avg_repair_symbols_read = 0.0;
  double local_repair_fraction = 0.0;
  double throughput_mbps = 0.0;
  double elapsed_seconds = 0.0;
};

struct SweepReport {
  std::string schema_version = "1.0";
  std::string param_name;
  std::vector<SimulationPointResult> points;
};

namespace detail {

struct ChunkAggregate {
  uint64_t frames = 0;
  uint64_t frame_errors = 0;
  uint64_t total_symbols = 0;
  uint64_t symbol_errors = 0;
  uint64_t local_repairs = 0;
  double sum_iterations = 0.0;
  double sum_repair_symbols_read = 0.0;
};

inline std::string CsvEscape(const std::string& s) {
  bool needs_quotes = false;
  for (char c : s) {
    if (c == ',' || c == '"' || c == '\n' || c == '\r') {
      needs_quotes = true;
      break;
    }
  }
  if (!needs_quotes) {
    return s;
  }
  std::string out;
  out.reserve(s.size() + 2);
  out.push_back('"');
  for (char c : s) {
    if (c == '"') {
      out.push_back('"');
    }
    out.push_back(c);
  }
  out.push_back('"');
  return out;
}

inline std::string JsonEscape(const std::string& s) {
  std::string out;
  out.reserve(s.size() + 4);
  for (unsigned char c : s) {
    switch (c) {
      case '"':
        out += "\\\"";
        break;
      case '\\':
        out += "\\\\";
        break;
      case '\b':
        out += "\\b";
        break;
      case '\f':
        out += "\\f";
        break;
      case '\n':
        out += "\\n";
        break;
      case '\r':
        out += "\\r";
        break;
      case '\t':
        out += "\\t";
        break;
      default:
        if (c < 0x20U) {
          static constexpr char kHex[] = "0123456789abcdef";
          out += "\\u00";
          out.push_back(kHex[(c >> 4) & 0x0fU]);
          out.push_back(kHex[c & 0x0fU]);
        } else {
          out.push_back(static_cast<char>(c));
        }
        break;
    }
  }
  return out;
}

inline void WriteJsonDouble(std::ostringstream& oss, double v) {
  if (std::isfinite(v)) {
    oss << std::setprecision(12) << v;
  } else {
    oss << "null";
  }
}

}  // namespace detail

// Executes deterministic parallel Monte Carlo frame trials using
// xg::simd::GlobalThreadPool(). Each chunk `chunk_idx` receives an independent
// Pcg64 stream seeded with DeriveStreamSeed(cfg.base_seed, sweep_idx,
// chunk_idx), and chunks are reduced in strict ascending chunk index order with
// per-chunk adaptive stopping checks. Consequently, all statistical fields
// (total_frames, frame_errors, total_symbols, symbol_errors, fer, ser,
// fer_ci95, avg_decoder_iterations, avg_repair_symbols_read,
// local_repair_fraction) are 100% bit-identical regardless of
// `cfg.num_threads`.
template <typename WorkerFn>
SimulationPointResult SimulateCodePerformanceParallel(
    const ParallelSimulationConfig& cfg, WorkerFn&& run_frame,
    double param_value = 0.0, uint64_t sweep_idx = 0,
    double bits_per_symbol = 1.0) {
  if (cfg.chunk_size == 0) {
    throw std::invalid_argument(
        "SimulateCodePerformanceParallel: chunk_size must be >= 1");
  }
  if (std::isnan(cfg.min_ci_relative_half_width) ||
      cfg.min_ci_relative_half_width < 0.0) {
    throw std::invalid_argument(
        "SimulateCodePerformanceParallel: min_ci_relative_half_width must be "
        ">= 0");
  }
  if (std::isnan(bits_per_symbol) || bits_per_symbol <= 0.0) {
    throw std::invalid_argument(
        "SimulateCodePerformanceParallel: bits_per_symbol must be > 0");
  }

  const auto start_time = std::chrono::steady_clock::now();

  SimulationPointResult result;
  result.param_value = param_value;
  if (cfg.max_frames == 0) {
    result.fer_ci95 = {0.0, 1.0};
    return result;
  }

  auto& pool = simd::GlobalThreadPool();
  const std::size_t effective_threads =
      (cfg.num_threads == 0) ? pool.NumWorkers() : cfg.num_threads;

  const uint64_t total_possible_chunks =
      (cfg.max_frames + cfg.chunk_size - 1) / cfg.chunk_size;

  uint64_t next_chunk_idx = 0;
  uint64_t local_repairs = 0;
  double sum_iterations = 0.0;
  double sum_repair_symbols_read = 0.0;
  bool stopped = false;

  auto* run_frame_ptr = std::addressof(run_frame);

  auto execute_chunk = [run_frame_ptr, &cfg, sweep_idx](
                           uint64_t chunk_idx) -> detail::ChunkAggregate {
    const uint64_t frame_start = chunk_idx * cfg.chunk_size;
    const uint64_t frames_in_chunk =
        std::min<uint64_t>(cfg.chunk_size, cfg.max_frames - frame_start);

    Pcg64 rng(DeriveStreamSeed(cfg.base_seed, sweep_idx, chunk_idx));
    detail::ChunkAggregate agg;
    agg.frames = frames_in_chunk;
    for (uint64_t f = 0; f < frames_in_chunk; ++f) {
      const FrameTrialOutcome outcome = (*run_frame_ptr)(rng);
      if (outcome.frame_error) {
        ++agg.frame_errors;
      }
      agg.total_symbols += outcome.symbols;
      agg.symbol_errors += outcome.symbol_errors;
      agg.sum_iterations += outcome.iterations;
      agg.sum_repair_symbols_read += outcome.repair_symbols_read;
      if (outcome.used_local_repair) {
        ++agg.local_repairs;
      }
    }
    return agg;
  };

  while (!stopped && next_chunk_idx < total_possible_chunks) {
    const uint64_t remaining_chunks = total_possible_chunks - next_chunk_idx;
    const std::size_t wave_chunks = static_cast<std::size_t>(
        std::min<uint64_t>(remaining_chunks, effective_threads));

    std::vector<detail::ChunkAggregate> wave_results(wave_chunks);

    if (effective_threads <= 1 || wave_chunks <= 1) {
      for (std::size_t k = 0; k < wave_chunks; ++k) {
        wave_results[k] = execute_chunk(next_chunk_idx + k);
      }
    } else {
      pool.ParallelFor(
          0, wave_chunks,
          [&wave_results, &execute_chunk, next_chunk_idx](std::size_t k) {
            wave_results[k] = execute_chunk(next_chunk_idx + k);
          },
          /*grain_size=*/1);
    }

    // Strict chunk-index-ordered reduction with per-chunk stopping checks.
    // Any speculative chunks beyond the exact stopping chunk are discarded so
    // results are invariant to `effective_threads`.
    for (std::size_t k = 0; k < wave_chunks; ++k) {
      const auto& agg = wave_results[k];
      result.total_frames += agg.frames;
      result.frame_errors += agg.frame_errors;
      result.total_symbols += agg.total_symbols;
      result.symbol_errors += agg.symbol_errors;
      local_repairs += agg.local_repairs;
      sum_iterations += agg.sum_iterations;
      sum_repair_symbols_read += agg.sum_repair_symbols_read;

      if (result.total_frames >= cfg.max_frames) {
        stopped = true;
        break;
      }
      if (cfg.target_frame_errors > 0 &&
          result.frame_errors >= cfg.target_frame_errors) {
        stopped = true;
        break;
      }
      if (cfg.min_ci_relative_half_width > 0.0 && result.frame_errors > 0) {
        const double fer_now = static_cast<double>(result.frame_errors) /
                               static_cast<double>(result.total_frames);
        const auto [ci_lo, ci_hi] =
            WilsonScoreInterval(result.frame_errors, result.total_frames, 0.95);
        const double rel_half_width = 0.5 * (ci_hi - ci_lo) / fer_now;
        if (rel_half_width <= cfg.min_ci_relative_half_width) {
          stopped = true;
          break;
        }
      }
    }

    next_chunk_idx += wave_chunks;
  }

  const auto end_time = std::chrono::steady_clock::now();
  result.elapsed_seconds =
      std::chrono::duration<double>(end_time - start_time).count();

  if (result.total_frames > 0) {
    const double frames_d = static_cast<double>(result.total_frames);
    result.fer = static_cast<double>(result.frame_errors) / frames_d;
    result.fer_ci95 =
        WilsonScoreInterval(result.frame_errors, result.total_frames, 0.95);
    result.avg_decoder_iterations = sum_iterations / frames_d;
    result.avg_repair_symbols_read = sum_repair_symbols_read / frames_d;
    result.local_repair_fraction =
        static_cast<double>(local_repairs) / frames_d;
  }
  if (result.total_symbols > 0) {
    result.ser = static_cast<double>(result.symbol_errors) /
                 static_cast<double>(result.total_symbols);
  }
  if (result.elapsed_seconds > 0.0) {
    result.throughput_mbps =
        (static_cast<double>(result.total_symbols) * bits_per_symbol) /
        (result.elapsed_seconds * 1e6);
  }

  return result;
}

template <typename MakeWorkerFn>
SweepReport RunParameterSweep(const std::string& param_name,
                              std::span<const double> param_values,
                              const ParallelSimulationConfig& cfg,
                              MakeWorkerFn&& make_worker,
                              double bits_per_symbol = 1.0) {
  SweepReport report;
  report.schema_version = "1.0";
  report.param_name = param_name;
  report.points.reserve(param_values.size());

  for (std::size_t idx = 0; idx < param_values.size(); ++idx) {
    const double param_val = param_values[idx];
    auto worker = make_worker(param_val);
    report.points.push_back(SimulateCodePerformanceParallel(
        cfg, worker, param_val, static_cast<uint64_t>(idx), bits_per_symbol));
  }
  return report;
}

template <typename MakeWorkerFn>
SweepReport RunParameterSweep(const std::string& param_name,
                              std::initializer_list<double> param_values,
                              const ParallelSimulationConfig& cfg,
                              MakeWorkerFn&& make_worker,
                              double bits_per_symbol = 1.0) {
  return RunParameterSweep(
      param_name,
      std::span<const double>(param_values.begin(), param_values.size()), cfg,
      std::forward<MakeWorkerFn>(make_worker), bits_per_symbol);
}

inline std::string ExportSweepCsv(const SweepReport& report) {
  std::ostringstream oss;
  oss.imbue(std::locale::classic());

  const std::string param_col = report.param_name.empty()
                                    ? "param"
                                    : detail::CsvEscape(report.param_name);
  oss << param_col
      << ",total_frames,frame_errors,total_symbols,symbol_errors,fer,ser,"
         "fer_ci95_low,fer_ci95_high,avg_decoder_iterations,"
         "avg_repair_symbols_read,local_repair_fraction,throughput_mbps,"
         "elapsed_seconds\n";

  oss << std::setprecision(12);
  for (const auto& pt : report.points) {
    oss << pt.param_value << ',' << pt.total_frames << ',' << pt.frame_errors
        << ',' << pt.total_symbols << ',' << pt.symbol_errors << ',' << pt.fer
        << ',' << pt.ser << ',' << pt.fer_ci95.first << ','
        << pt.fer_ci95.second << ',' << pt.avg_decoder_iterations << ','
        << pt.avg_repair_symbols_read << ',' << pt.local_repair_fraction << ','
        << pt.throughput_mbps << ',' << pt.elapsed_seconds << '\n';
  }
  return oss.str();
}

inline std::string ExportSweepJson(const SweepReport& report) {
  std::ostringstream oss;
  oss.imbue(std::locale::classic());

  oss << "{\n";
  oss << "  \"schema_version\": \"" << detail::JsonEscape(report.schema_version)
      << "\",\n";
  oss << "  \"param_name\": \"" << detail::JsonEscape(report.param_name)
      << "\",\n";
  oss << "  \"points\": [\n";
  for (std::size_t i = 0; i < report.points.size(); ++i) {
    const auto& pt = report.points[i];
    oss << "    {\n";
    oss << "      \"param_value\": ";
    detail::WriteJsonDouble(oss, pt.param_value);
    oss << ",\n";
    oss << "      \"total_frames\": " << pt.total_frames << ",\n";
    oss << "      \"frame_errors\": " << pt.frame_errors << ",\n";
    oss << "      \"total_symbols\": " << pt.total_symbols << ",\n";
    oss << "      \"symbol_errors\": " << pt.symbol_errors << ",\n";
    oss << "      \"fer\": ";
    detail::WriteJsonDouble(oss, pt.fer);
    oss << ",\n";
    oss << "      \"ser\": ";
    detail::WriteJsonDouble(oss, pt.ser);
    oss << ",\n";
    oss << "      \"fer_ci95\": [";
    detail::WriteJsonDouble(oss, pt.fer_ci95.first);
    oss << ", ";
    detail::WriteJsonDouble(oss, pt.fer_ci95.second);
    oss << "],\n";
    oss << "      \"avg_decoder_iterations\": ";
    detail::WriteJsonDouble(oss, pt.avg_decoder_iterations);
    oss << ",\n";
    oss << "      \"avg_repair_symbols_read\": ";
    detail::WriteJsonDouble(oss, pt.avg_repair_symbols_read);
    oss << ",\n";
    oss << "      \"local_repair_fraction\": ";
    detail::WriteJsonDouble(oss, pt.local_repair_fraction);
    oss << ",\n";
    oss << "      \"throughput_mbps\": ";
    detail::WriteJsonDouble(oss, pt.throughput_mbps);
    oss << ",\n";
    oss << "      \"elapsed_seconds\": ";
    detail::WriteJsonDouble(oss, pt.elapsed_seconds);
    oss << "\n    }";
    if (i + 1 < report.points.size()) {
      oss << ',';
    }
    oss << '\n';
  }
  oss << "  ]\n";
  oss << "}\n";
  return oss.str();
}

}  // namespace channel
}  // namespace xg

#endif  // XGALOIS_CHANNEL_SIMULATION_HPP_
