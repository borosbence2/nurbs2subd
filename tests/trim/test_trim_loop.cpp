#include "n2s/trim/cases.hpp"
#include "n2s/trim/trim_loop.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

using Catch::Approx;
using n2s::NurbsCurve2;
using n2s::Orientation;
using n2s::TrimLoop;
using n2s::TrimRegion;

TEST_CASE("reversing a curve keeps its image and flips its parameterisation", "[trim][curve]") {
    const TrimLoop circle = n2s::cases::circle_loop({0.0, 0.0}, 1.0);
    const NurbsCurve2& arc = circle.curves().front();
    const NurbsCurve2 back = n2s::reversed(arc);

    const double a = arc.knots().domain_start();
    const double b = arc.knots().domain_end();

    REQUIRE(back.knots().domain_start() == Approx(a));
    REQUIRE(back.knots().domain_end() == Approx(b));

    for (int i = 0; i <= 50; ++i) {
        const double t = a + (b - a) * static_cast<double>(i) / 50.0;
        INFO("t = " << t);
        // The mirrored parameter of t is a + b - t.
        CHECK((arc.evaluate(t) - back.evaluate(a + b - t)).norm() < 1e-13);
    }
}

TEST_CASE("a rectangle loop closes, encloses the right area and runs counter-clockwise",
          "[trim][loop]") {
    const TrimLoop rect = n2s::cases::rectangle_loop({0.0, 0.0}, {2.0, 3.0});

    CHECK(rect.size() == 4);
    CHECK(rect.worst_closure_gap() == Approx(0.0).margin(1e-15));
    CHECK(rect.orientation() == Orientation::CounterClockwise);
    // A polygon through a straight-sided loop is exact, whatever the sampling.
    CHECK(rect.signed_area() == Approx(6.0).margin(1e-12));
}

TEST_CASE("a circle loop converges to pi r squared", "[trim][loop]") {
    const double radius = 0.3;
    const TrimLoop circle = n2s::cases::circle_loop({0.5, 0.5}, radius);
    const double exact = std::numbers::pi * radius * radius;

    CHECK(circle.worst_closure_gap() == Approx(0.0).margin(1e-14));
    CHECK(circle.orientation() == Orientation::CounterClockwise);

    // An inscribed polygon always understates the area of a convex region, and
    // the deficit falls as the square of the sampling step.
    double previous_error = 1.0;
    for (const int samples : {8, 16, 32, 64}) {
        const double area = circle.signed_area(samples);
        const double error = std::abs(area - exact);
        INFO("samples per curve " << samples << ", area " << area << ", error " << error);
        CHECK(area < exact);
        CHECK(error < previous_error / 3.0);
        previous_error = error;
    }
    CHECK(circle.signed_area(4096) == Approx(exact).epsilon(1e-6));
}

TEST_CASE("reversing a loop flips its orientation and nothing else", "[trim][loop]") {
    TrimLoop loop = n2s::cases::circle_loop({0.0, 0.0}, 1.0);
    const double area = loop.signed_area();

    loop.reverse();
    CHECK(loop.orientation() == Orientation::Clockwise);
    CHECK(loop.signed_area() == Approx(-area));
    CHECK(loop.worst_closure_gap() == Approx(0.0).margin(1e-14));

    loop.reverse();
    CHECK(loop.orientation() == Orientation::CounterClockwise);
    CHECK(loop.signed_area() == Approx(area));
}

TEST_CASE("an empty loop is rejected", "[trim][loop]") {
    CHECK_THROWS_AS(TrimLoop{std::vector<NurbsCurve2>{}}, std::invalid_argument);
}

TEST_CASE("point containment respects holes", "[trim][region]") {
    TrimLoop hole = n2s::cases::circle_loop({0.5, 0.5}, 0.25);
    hole.reverse();
    const TrimRegion region{n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0}), {std::move(hole)}};

    CHECK(region.contains({0.1, 0.1}));
    CHECK(region.contains({0.9, 0.5}));

    CHECK_FALSE(region.contains({0.5, 0.5}));  // dead centre of the hole
    CHECK_FALSE(region.contains({0.5, 0.65})); // still inside the hole
    CHECK_FALSE(region.contains({-0.1, 0.5})); // outside the outer loop
    CHECK_FALSE(region.contains({1.5, 1.5}));
}

TEST_CASE("region area subtracts the holes", "[trim][region]") {
    const double radius = 0.25;
    TrimLoop hole = n2s::cases::circle_loop({0.5, 0.5}, radius);
    hole.reverse();
    const TrimRegion region{n2s::cases::rectangle_loop({0.0, 0.0}, {1.0, 1.0}), {std::move(hole)}};

    const double exact = 1.0 - std::numbers::pi * radius * radius;
    CHECK(region.area(2048) == Approx(exact).epsilon(1e-6));
}

TEST_CASE("the L-shape and arc-corner loops close and enclose their analytic area",
          "[trim][loop]") {
    SECTION("L-shape") {
        const TrimLoop loop = n2s::cases::l_shape_loop({0.0, 0.0}, {1.0, 1.0}, {0.5, 0.5});
        CHECK(loop.size() == 6);
        CHECK(loop.worst_closure_gap() == Approx(0.0).margin(1e-15));
        CHECK(loop.orientation() == Orientation::CounterClockwise);
        CHECK(loop.signed_area() == Approx(0.75).margin(1e-12));
    }

    SECTION("quarter-arc corner, the Shen et al. figure 1 domain") {
        const double radius = 0.5;
        const TrimLoop loop = n2s::cases::quarter_arc_corner_loop({0.0, 0.0}, {1.0, 1.0}, radius);
        CHECK(loop.size() == 5);
        CHECK(loop.worst_closure_gap() == Approx(0.0).margin(1e-14));
        CHECK(loop.orientation() == Orientation::CounterClockwise);

        // Unit square less a quarter disc. The polygonised arc cuts the corner,
        // so the sampled area slightly overstates the true one.
        const double exact = 1.0 - 0.25 * std::numbers::pi * radius * radius;
        CHECK(loop.signed_area(2048) == Approx(exact).epsilon(1e-6));
    }
}
