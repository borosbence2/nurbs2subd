#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace n2s::metrics {

/// Summary of a set of non-negative samples: a distance field, an angle field,
/// a curvature difference.
///
/// `max` is reported alongside `rms` and `mean` deliberately. A fitting method
/// is easy to make look good on RMS while leaving one region badly wrong, and
/// the maximum is the number a watertightness or tolerance claim actually rests
/// on. Reporting only the average is how a local defect disappears into a
/// good-looking table.
struct Stats {
    double max = 0.0;
    double mean = 0.0;
    double rms = 0.0;
    std::size_t count = 0;

    /// Samples that could not be measured -- a projection that did not
    /// converge, a point where the normal is undefined. Counted separately and
    /// *excluded* from the statistics above, because folding a failed
    /// measurement in as zero flatters the result and folding it in as
    /// something large is invented data. A non-zero value here means the
    /// numbers describe fewer samples than were asked for, and any report has
    /// to say so.
    std::size_t unmeasured = 0;

    bool complete() const { return unmeasured == 0; }
};

/// Accumulates samples one at a time without storing them.
class Accumulator {
public:
    void add(double value);

    /// Records that a sample could not be measured at all.
    void add_unmeasured() { ++unmeasured_; }

    Stats result() const;

private:
    double max_ = 0.0;
    double sum_ = 0.0;
    double sum_of_squares_ = 0.0;
    std::size_t count_ = 0;
    std::size_t unmeasured_ = 0;
};

/// Summarises a vector of samples in one call.
Stats summarise(const std::vector<double>& values);

} // namespace n2s::metrics
