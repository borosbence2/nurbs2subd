#include "n2s/metrics/statistics.hpp"

#include <algorithm>
#include <cmath>

namespace n2s::metrics {

void Accumulator::add(double value) {
    if (!std::isfinite(value)) {
        // A non-finite sample is a failed measurement wearing a number. Letting
        // one into the sums turns every statistic into NaN and hides which
        // sample was at fault.
        ++unmeasured_;
        return;
    }

    max_ = std::max(max_, value);
    sum_ += value;
    sum_of_squares_ += value * value;
    ++count_;
}

Stats Accumulator::result() const {
    Stats stats;
    stats.count = count_;
    stats.unmeasured = unmeasured_;
    if (count_ == 0) {
        return stats;
    }

    const auto n = static_cast<double>(count_);
    stats.max = max_;
    stats.mean = sum_ / n;
    stats.rms = std::sqrt(sum_of_squares_ / n);
    return stats;
}

Stats summarise(const std::vector<double>& values) {
    Accumulator accumulator;
    for (const double value : values) {
        accumulator.add(value);
    }
    return accumulator.result();
}

} // namespace n2s::metrics
