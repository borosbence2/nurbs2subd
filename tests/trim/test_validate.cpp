#include "n2s/trim/cases.hpp"
#include "n2s/trim/validate.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <string>
#include <vector>

using Catch::Approx;
using n2s::NurbsCurve2;
using n2s::Orientation;
using n2s::TrimLoop;
using n2s::TrimRegion;

namespace {

NurbsCurve2 segment(const Eigen::Vector2d& from, const Eigen::Vector2d& to) {
    return NurbsCurve2::bspline(n2s::KnotVector{1, {0, 0, 1, 1}},
                                std::vector<Eigen::Vector2d>{from, to});
}

/// A unit square whose fourth corner misses the first by `gap`.
TrimLoop square_with_gap(double gap) {
    return TrimLoop{{segment({0, 0}, {1, 0}),
                     segment({1, 0}, {1, 1}),
                     segment({1, 1}, {0, 1}),
                     segment({0, 1}, {gap, 0})}};
}

bool any_contains(const std::vector<std::string>& lines, const std::string& needle) {
    for (const std::string& line : lines) {
        if (line.find(needle) != std::string::npos) {
            return true;
        }
    }
    return false;
}

} // namespace

TEST_CASE("a clean region needs no repair", "[trim][validate]") {
    n2s::cases::TrimmedCase testcase =
        n2s::cases::square_with_circular_hole(n2s::cases::plane(), 0.25);

    const n2s::TrimReport report = n2s::validate_and_repair(testcase.region);

    INFO(report.to_string());
    CHECK(report.ok());
    CHECK(report.repairs.empty());
    CHECK(report.errors.empty());
}

TEST_CASE("a hole wound the wrong way is reversed, not rejected", "[trim][validate]") {
    // This is the defect the pass exists for. A counter-clockwise hole
    // triangulates into a region with the hole filled and the surround hollow,
    // which looks perfectly plausible until an error metric is computed on it.
    TrimLoop hole = n2s::cases::circle_loop({0.5, 0.5}, 0.25); // counter-clockwise
    REQUIRE(hole.orientation() == Orientation::CounterClockwise);

    TrimRegion region{n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0}), {std::move(hole)}};
    const n2s::TrimReport report = n2s::validate_and_repair(region);

    INFO(report.to_string());
    CHECK(report.ok());
    REQUIRE(report.repairs.size() == 1);
    CHECK(any_contains(report.repairs, "reversed"));
    CHECK(region.holes().front().orientation() == Orientation::Clockwise);
    CHECK(region.outer().orientation() == Orientation::CounterClockwise);
}

TEST_CASE("an outer loop wound the wrong way is reversed", "[trim][validate]") {
    TrimLoop outer = n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0});
    outer.reverse();
    REQUIRE(outer.orientation() == Orientation::Clockwise);

    TrimRegion region{std::move(outer)};
    const n2s::TrimReport report = n2s::validate_and_repair(region);

    CHECK(report.ok());
    CHECK(any_contains(report.repairs, "reversed"));
    CHECK(region.outer().orientation() == Orientation::CounterClockwise);
}

TEST_CASE("closure gaps are snapped below tolerance and reported above it", "[trim][validate]") {
    SECTION("within tolerance: snapped shut") {
        TrimRegion region{square_with_gap(1e-12)};
        REQUIRE(region.outer().worst_closure_gap() > 0.0);

        const n2s::TrimReport report = n2s::validate_and_repair(region);

        INFO(report.to_string());
        CHECK(report.ok());
        CHECK(any_contains(report.repairs, "snapped"));
        CHECK(region.outer().worst_closure_gap() == Approx(0.0).margin(1e-15));
    }

    SECTION("above tolerance: reported as an error, not silently closed") {
        TrimRegion region{square_with_gap(1e-3)};
        const n2s::TrimReport report = n2s::validate_and_repair(region);

        INFO(report.to_string());
        CHECK_FALSE(report.ok());
        CHECK(any_contains(report.errors, "gap"));
    }

    SECTION("report-only mode leaves the geometry alone") {
        TrimRegion region{square_with_gap(1e-12)};
        n2s::TrimValidationOptions options;
        options.snap_closure_gaps = false;

        const n2s::TrimReport report = n2s::validate_and_repair(region, options);

        CHECK(report.repairs.empty());
        CHECK(region.outer().worst_closure_gap() > 0.0);
    }
}

TEST_CASE("closure tolerance scales with the size of the loop", "[trim][validate]") {
    // The defect real thesis data exposed. Every synthetic case in this project
    // is unit-sized, so an absolute 1e-9 tolerance looked fine. On a loop
    // spanning 450 units the joins of genuinely clean CAD output agree to about
    // a part in 10^12, which is 5e-10 absolute -- close enough to the old fixed
    // floor that whether a loop validated depended on where it sat in space.
    const double span = 450.0;
    const double rounding_scale_gap = 5e-10 * span; // ~2e-7, far above 1e-9

    const auto big_square_with_gap = [&](double gap) {
        return TrimLoop{{segment({0, 0}, {span, 0}),
                         segment({span, 0}, {span, span}),
                         segment({span, span}, {0, span}),
                         segment({0, span}, {gap, 0})}};
    };

    SECTION("a gap at the rounding scale of a large model is snapped") {
        TrimRegion region{big_square_with_gap(rounding_scale_gap)};
        const n2s::TrimReport report = n2s::validate_and_repair(region);

        INFO(report.to_string());
        CHECK(report.ok());
        CHECK(region.outer().worst_closure_gap() < 1e-9);
    }

    SECTION("a gap that is a real fraction of the model is still an error") {
        // 0.2% of the span, which is the size of the four bad joins in the
        // thesis trim loop. Scaling the tolerance must not turn a genuine
        // defect into a silent repair.
        TrimRegion region{big_square_with_gap(0.002 * span)};
        const n2s::TrimReport report = n2s::validate_and_repair(region);

        INFO(report.to_string());
        CHECK_FALSE(report.ok());
        // The message quotes the gap as a fraction of the loop, because an
        // absolute figure means nothing without the scale beside it.
        CHECK(any_contains(report.errors, "% of it"));
    }

    SECTION("the same relative gap behaves identically at unit scale") {
        TrimRegion small{TrimLoop{{segment({0, 0}, {1, 0}),
                                   segment({1, 0}, {1, 1}),
                                   segment({1, 1}, {0, 1}),
                                   segment({0, 1}, {0.002, 0})}}};
        CHECK_FALSE(n2s::validate_and_repair(small).ok());
    }
}

TEST_CASE("degenerate-area detection scales too", "[trim][validate]") {
    // Area grows as the square of the loop, so a fixed threshold calls a small
    // loop degenerate and lets a large sliver through.
    const auto tiny_square = [](double side) {
        return TrimLoop{{segment({0, 0}, {side, 0}),
                         segment({side, 0}, {side, side}),
                         segment({side, side}, {0, side}),
                         segment({0, side}, {0, 0})}};
    };

    // A 1e-4 square has an area of 1e-8, which a fixed 1e-14 threshold would
    // accept and a scaled one also accepts -- it is small, not degenerate.
    TrimRegion small{tiny_square(1e-4)};
    CHECK(n2s::validate_and_repair(small).ok());
}

TEST_CASE("a self-intersecting loop is rejected", "[trim][validate]") {
    // A symmetric bowtie: closed, but crossing itself. Its signed area is
    // exactly zero, which is why the self-intersection test has to run before
    // the orientation test -- otherwise this is reported as "degenerate area"
    // and the real defect is buried.
    TrimLoop bowtie{{segment({0, 0}, {1, 1}),
                     segment({1, 1}, {1, 0}),
                     segment({1, 0}, {0, 1}),
                     segment({0, 1}, {0, 0})}};
    REQUIRE(bowtie.signed_area() == Approx(0.0).margin(1e-12));

    TrimRegion region{std::move(bowtie)};
    const n2s::TrimReport report = n2s::validate_and_repair(region);

    INFO(report.to_string());
    CHECK_FALSE(report.ok());
    CHECK(any_contains(report.errors, "self-intersects"));
}

TEST_CASE("a hole outside the outer loop is rejected", "[trim][validate]") {
    TrimLoop hole = n2s::cases::circle_loop({5.0, 5.0}, 0.25);
    hole.reverse();

    TrimRegion region{n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0}), {std::move(hole)}};
    const n2s::TrimReport report = n2s::validate_and_repair(region);

    INFO(report.to_string());
    CHECK_FALSE(report.ok());
    CHECK(any_contains(report.errors, "not contained"));
}

TEST_CASE("overlapping holes are rejected", "[trim][validate]") {
    TrimLoop first = n2s::cases::circle_loop({0.45, 0.5}, 0.2);
    TrimLoop second = n2s::cases::circle_loop({0.55, 0.5}, 0.2);
    first.reverse();
    second.reverse();

    TrimRegion region{n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0}),
                      {std::move(first), std::move(second)}};
    const n2s::TrimReport report = n2s::validate_and_repair(region);

    INFO(report.to_string());
    CHECK_FALSE(report.ok());
    CHECK(any_contains(report.errors, "overlap"));
}

TEST_CASE("every synthetic case validates cleanly", "[trim][validate][cases]") {
    std::vector<n2s::cases::TrimmedCase> all = n2s::cases::all_synthetic_cases();
    REQUIRE(all.size() == 12); // 4 surfaces x 3 trims

    for (n2s::cases::TrimmedCase& one : all) {
        const n2s::TrimReport report = n2s::validate_and_repair(one.region);
        INFO(one.name << "\n" << report.to_string());
        CHECK(report.ok());
        CHECK(report.repairs.empty());
    }
}

TEST_CASE("every synthetic case reports its analytic domain area", "[trim][cases]") {
    for (const n2s::cases::TrimmedCase& one : n2s::cases::all_synthetic_cases()) {
        INFO(one.name);
        REQUIRE(one.analytic_domain_area > 0.0);
        CHECK(one.region.area(2048) == Approx(one.analytic_domain_area).epsilon(1e-5));
    }
}
