#include "n2s/nurbs/knot_vector.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <stdexcept>
#include <vector>

using n2s::KnotVector;

TEST_CASE("a clamped cubic knot vector reports its basic quantities", "[nurbs][knots]") {
    // p = 3, m = 9, so n + 1 = m - p = 6 control points.
    const KnotVector kv{3, {0.0, 0.0, 0.0, 0.0, 0.4, 0.7, 1.0, 1.0, 1.0, 1.0}};

    CHECK(kv.degree() == 3);
    CHECK(kv.size() == 10);
    CHECK(kv.num_control_points() == 6);
    CHECK(kv.domain_start() == 0.0);
    CHECK(kv.domain_end() == 1.0);
    CHECK(kv.is_clamped());
    CHECK(kv.num_spans() == 3);
}

TEST_CASE("find_span matches the Piegl & Tiller definition", "[nurbs][knots]") {
    // The worked example from The NURBS Book, section 2.5: p = 2,
    // U = {0,0,0,1,2,3,4,4,5,5,5}.
    const KnotVector kv{2, {0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}};

    SECTION("interior parameters land in the span containing them") {
        CHECK(kv.find_span(0.0) == 2);
        CHECK(kv.find_span(0.5) == 2);
        CHECK(kv.find_span(1.0) == 3);
        CHECK(kv.find_span(2.5) == 4);
        CHECK(kv.find_span(3.0) == 5);
        CHECK(kv.find_span(4.0) == 7); // U[6] == U[7] == 4, so the span is [4,5)
    }

    SECTION("the end of the domain resolves to the last non-empty span") {
        // The textbook special case: u == U[m-p] would otherwise index past the
        // last basis function.
        CHECK(kv.find_span(5.0) == 7);
    }

    SECTION("parameters outside the domain clamp instead of reading out of bounds") {
        CHECK(kv.find_span(-1.0) == kv.find_span(kv.domain_start()));
        CHECK(kv.find_span(6.0) == kv.find_span(kv.domain_end()));
    }
}

TEST_CASE("multiplicity counts repeated knots", "[nurbs][knots]") {
    const KnotVector kv{2, {0, 0, 0, 1, 2, 2, 3, 3, 3}};

    CHECK(kv.multiplicity(0.0) == 3);
    CHECK(kv.multiplicity(1.0) == 1);
    CHECK(kv.multiplicity(2.0) == 2);
    CHECK(kv.multiplicity(3.0) == 3);
    CHECK(kv.multiplicity(1.5) == 0);
}

TEST_CASE("uniform_clamped builds a valid vector", "[nurbs][knots]") {
    const KnotVector kv = KnotVector::uniform_clamped(3, 7);

    REQUIRE(kv.num_control_points() == 7);
    CHECK(kv.size() == 11);
    CHECK(kv.is_clamped());
    CHECK(kv.domain_start() == 0.0);
    CHECK(kv.domain_end() == 1.0);

    // Interior knots are evenly spaced: 3 interior knots at 1/4, 2/4, 3/4.
    CHECK(kv[4] == 0.25);
    CHECK(kv[5] == 0.50);
    CHECK(kv[6] == 0.75);
}

TEST_CASE("a Bezier segment is the degenerate uniform case", "[nurbs][knots]") {
    const KnotVector kv = KnotVector::uniform_clamped(3, 4);

    CHECK(kv.size() == 8);
    CHECK(kv.num_spans() == 1);
    CHECK(kv.knots() == std::vector<double>{0, 0, 0, 0, 1, 1, 1, 1});
}

TEST_CASE("invalid knot vectors are rejected at construction", "[nurbs][knots]") {
    using Catch::Matchers::ContainsSubstring;

    SECTION("degree must be at least 1") {
        CHECK_THROWS_AS((KnotVector{0, {0, 1}}), std::invalid_argument);
    }

    SECTION("too few knots for the degree") {
        // p = 3 needs at least 2 * (p + 1) = 8 knots.
        CHECK_THROWS_WITH((KnotVector{3, {0, 0, 0, 0, 1, 1, 1}}), ContainsSubstring("at least 8"));
    }

    SECTION("knots must be non-decreasing") {
        CHECK_THROWS_WITH((KnotVector{1, {0, 0, 0.5, 0.2, 1, 1}}),
                          ContainsSubstring("non-decreasing"));
    }

    SECTION("the domain must not be degenerate") {
        CHECK_THROWS_WITH((KnotVector{1, {0, 0, 0, 0, 0, 0}}), ContainsSubstring("degenerate"));
    }

    SECTION("an interior knot may not repeat more than degree times") {
        // Multiplicity p + 1 in the interior disconnects the curve.
        CHECK_THROWS_WITH((KnotVector{2, {0, 0, 0, 1, 1, 1, 2, 2, 2}}),
                          ContainsSubstring("multiplicity"));
    }
}
