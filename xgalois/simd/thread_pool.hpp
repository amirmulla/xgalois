#ifndef XGALOIS_SIMD_THREAD_POOL_HPP_
#define XGALOIS_SIMD_THREAD_POOL_HPP_

#include <algorithm>
#include <atomic>
#include <concepts>
#include <condition_variable>
#include <cstddef>
#include <deque>
#include <exception>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <span>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace xg {
namespace simd {

// Thread pool backed by std::jthread workers with deadlock-free nested
// ParallelFor execution and clean exception propagation.
class ThreadPool {
 public:
  static std::size_t DefaultNumWorkers() noexcept {
    return std::max<std::size_t>(
        1, static_cast<std::size_t>(std::thread::hardware_concurrency()));
  }

  explicit ThreadPool(std::size_t num_threads = DefaultNumWorkers()) {
    const std::size_t count = std::max<std::size_t>(1, num_threads);
    workers_.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
      workers_.emplace_back([this](std::stop_token) { WorkerLoop(); });
    }
  }

  ~ThreadPool() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stopping_ = true;
    }
    cv_.notify_all();
    workers_.clear();
  }

  ThreadPool(const ThreadPool&) = delete;
  ThreadPool& operator=(const ThreadPool&) = delete;
  ThreadPool(ThreadPool&&) = delete;
  ThreadPool& operator=(ThreadPool&&) = delete;

  // Returns the number of worker threads owned by this pool.
  std::size_t NumWorkers() const noexcept { return workers_.size(); }

  // Submits a callable and arguments for asynchronous execution and returns a
  // std::future carrying the result or any thrown exception.
  template <typename Fn, typename... Args>
  auto Submit(Fn&& fn, Args&&... args)
      -> std::future<
          std::invoke_result_t<std::decay_t<Fn>, std::decay_t<Args>...>> {
    using Result =
        std::invoke_result_t<std::decay_t<Fn>, std::decay_t<Args>...>;

    auto task = std::make_shared<std::packaged_task<Result()>>(
        [callable = std::forward<Fn>(fn),
         bound_args =
             std::make_tuple(std::forward<Args>(args)...)]() mutable -> Result {
          return std::apply(std::move(callable), std::move(bound_args));
        });

    std::future<Result> future = task->get_future();
    {
      std::lock_guard<std::mutex> lock(mutex_);
      if (stopping_) {
        throw std::runtime_error("Cannot submit task to a stopped ThreadPool");
      }
      tasks_.emplace_back([task = std::move(task)]() { (*task)(); });
    }
    cv_.notify_one();
    return future;
  }

  // Splits [begin, end) into balanced chunks and executes them across worker
  // threads plus the calling thread. Supports both range callbacks
  // fn(chunk_begin, chunk_end) and per-index callbacks fn(index).
  // Deadlock-free even when called recursively from inside worker threads, and
  // propagates the first exception thrown by any chunk after all chunks join.
  template <typename Fn>
  void ParallelFor(std::size_t begin, std::size_t end, Fn&& fn,
                   std::size_t grain_size = 0) {
    static_assert(std::invocable<Fn&, std::size_t, std::size_t> ||
                      std::invocable<Fn&, std::size_t>,
                  "ParallelFor callback must be invocable with either "
                  "(size_t, size_t) or (size_t)");

    if (begin >= end) {
      return;
    }

    const std::size_t total = end - begin;
    const std::size_t num_chunks = (grain_size > 0)
                                       ? (1 + (total - 1) / grain_size)
                                       : std::min(total, workers_.size() + 1);

    if (num_chunks <= 1) {
      InvokeChunk(fn, begin, end);
      return;
    }

    struct ParallelState {
      std::atomic<std::size_t> next_chunk{0};
      std::atomic<std::size_t> completed_chunks{0};
      std::mutex exc_mutex;
      std::exception_ptr first_exception;
      std::mutex done_mutex;
      std::condition_variable done_cv;
    };

    auto state = std::make_shared<ParallelState>();
    const std::size_t base = total / num_chunks;
    const std::size_t rem = total % num_chunks;
    auto* fn_ptr = std::addressof(fn);

    auto run_chunks = [state, fn_ptr, begin, num_chunks, base, rem]() {
      while (true) {
        if (state->next_chunk.load(std::memory_order_relaxed) >= num_chunks) {
          return;
        }
        const std::size_t chunk_idx =
            state->next_chunk.fetch_add(1, std::memory_order_relaxed);
        if (chunk_idx >= num_chunks) {
          return;
        }

        const std::size_t c_begin =
            begin + chunk_idx * base + std::min(chunk_idx, rem);
        const std::size_t c_end =
            c_begin + base +
            (chunk_idx < rem ? std::size_t{1} : std::size_t{0});

        try {
          InvokeChunk(*fn_ptr, c_begin, c_end);
        } catch (...) {
          std::lock_guard<std::mutex> exc_lock(state->exc_mutex);
          if (!state->first_exception) {
            state->first_exception = std::current_exception();
          }
        }

        if (state->completed_chunks.fetch_add(1, std::memory_order_acq_rel) +
                1 ==
            num_chunks) {
          std::lock_guard<std::mutex> done_lock(state->done_mutex);
          state->done_cv.notify_one();
        }
      }
    };

    const std::size_t helper_count = std::min(num_chunks - 1, workers_.size());
    std::size_t queued_helpers = 0;
    {
      std::lock_guard<std::mutex> lock(mutex_);
      if (!stopping_) {
        for (std::size_t i = 0; i < helper_count; ++i) {
          tasks_.emplace_back(run_chunks);
          ++queued_helpers;
        }
      }
    }
    if (queued_helpers == 1) {
      cv_.notify_one();
    } else if (queued_helpers > 1) {
      cv_.notify_all();
    }

    // Calling thread participates directly and drains any unclaimed chunks,
    // guaranteeing forward progress even under nested ParallelFor calls.
    run_chunks();

    if (state->completed_chunks.load(std::memory_order_acquire) < num_chunks) {
      std::unique_lock<std::mutex> done_lock(state->done_mutex);
      state->done_cv.wait(done_lock, [&state, num_chunks] {
        return state->completed_chunks.load(std::memory_order_acquire) ==
               num_chunks;
      });
    }

    if (state->first_exception) {
      std::rethrow_exception(state->first_exception);
    }
  }

  // Applies fn element-wise from input into output in parallel.
  template <typename In, typename Out, typename Fn>
  void TransformBatch(std::span<const In> input, std::span<Out> output, Fn&& fn,
                      std::size_t grain_size = 0) {
    if (input.size() != output.size()) {
      throw std::invalid_argument(
          "TransformBatch: input and output spans must have equal size");
    }
    if (input.empty()) {
      return;
    }
    auto* fn_ptr = std::addressof(fn);
    ParallelFor(
        0, input.size(),
        [input, output, fn_ptr](std::size_t lo, std::size_t hi) {
          for (std::size_t i = lo; i < hi; ++i) {
            output[i] = (*fn_ptr)(input[i]);
          }
        },
        grain_size);
  }

 private:
  template <typename Fn>
  static void InvokeChunk(Fn& fn, std::size_t c_begin, std::size_t c_end) {
    if constexpr (std::invocable<Fn&, std::size_t, std::size_t>) {
      fn(c_begin, c_end);
    } else {
      for (std::size_t i = c_begin; i < c_end; ++i) {
        fn(i);
      }
    }
  }

  void WorkerLoop() {
    while (true) {
      std::function<void()> task;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_.wait(lock, [this] { return stopping_ || !tasks_.empty(); });
        if (stopping_ && tasks_.empty()) {
          return;
        }
        task = std::move(tasks_.front());
        tasks_.pop_front();
      }
      task();
    }
  }

  std::mutex mutex_;
  std::condition_variable cv_;
  std::deque<std::function<void()>> tasks_;
  bool stopping_ = false;
  std::vector<std::jthread> workers_;
};

// Meyer's singleton returning the process-wide default ThreadPool.
inline ThreadPool& GlobalThreadPool() {
  static ThreadPool pool;
  return pool;
}

// Convenience wrapper forwarding ParallelFor to GlobalThreadPool().
template <typename Fn>
inline void ParallelFor(std::size_t begin, std::size_t end, Fn&& fn,
                        std::size_t grain_size = 0) {
  GlobalThreadPool().ParallelFor(begin, end, std::forward<Fn>(fn), grain_size);
}

// Convenience wrapper forwarding TransformBatch to GlobalThreadPool().
template <typename In, typename Out, typename Fn>
inline void TransformBatch(std::span<const In> input, std::span<Out> output,
                           Fn&& fn, std::size_t grain_size = 0) {
  GlobalThreadPool().TransformBatch(input, output, std::forward<Fn>(fn),
                                    grain_size);
}

}  // namespace simd
}  // namespace xg

#endif  // XGALOIS_SIMD_THREAD_POOL_HPP_
