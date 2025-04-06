#pragma once
#include "common_types.h"
#include <functional>
#include <optional> // Use optional for potential failure

namespace NumericalSolvers {

    // Brent's method for finding root of f(x) = 0 in [a, b]
    std::optional<Real> brent(std::function<Real(Real)> func, Real a, Real b,
                              Real tol = 1e-9, int max_iter = 100);

    // Simple 3x3 matrix solver (Gaussian elimination or Cramer's rule)
    // Solves Ax = b for x
    // Returns true on success, false if matrix is singular
    bool solve_3x3(const std::array<std::array<Real, 3>, 3>& A, const std::array<Real, 3>& b, std::array<Real, 3>& x);

} // namespace NumericalSolvers