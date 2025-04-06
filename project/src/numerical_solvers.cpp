#include "numerical_solvers.h"
#include <cmath>
#include <algorithm> // For std::swap, std::abs
#include <iostream> // For warnings

namespace NumericalSolvers {

// Brent's method implementation (standard algorithm)
std::optional<Real> brent(std::function<Real(Real)> func, Real a, Real b, Real tol, int max_iter) {
    Real fa = func(a);
    Real fb = func(b);

    if (std::abs(fa) < tol) return a;
    if (std::abs(fb) < tol) return b;

    // Check if initial bounds are valid for Brent's method (root must be bracketed)
    // However, the entropic objective might not strictly bracket the root if alpha=1 is the minimum.
    // We rely on the caller (EntropicCollision) to handle cases where fa*fb >= 0.
    // if (fa * fb >= 0) {
    //     std::cerr << "Warning (Brent): Root not bracketed or function invalid at bounds. f(a)=" << fa << ", f(b)=" << fb << std::endl;
    //     // Return nullopt or handle based on context (e.g., check if |fa| or |fb| is small)
    //      return std::nullopt;
    // }
     if (fa * fb > 0 && std::abs(fa) > tol && std::abs(fb) > tol) {
         // If no sign change and neither bound is close to zero, Brent won't work well.
         // The caller should have handled this (e.g., defaulting alpha=1).
         // Returning nullopt indicates failure to the caller.
         return std::nullopt;
     }


    if (std::abs(fa) < std::abs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }

    Real c = a;
    Real fc = fa;
    bool mflag = true;
    Real d = 0.0, s = 0.0;

    for (int iter = 0; iter < max_iter; ++iter) {
        if (std::abs(fb) < tol || std::abs(b-a) < tol) { // Check using tolerance relative to b? std::abs(b-a) < tol*std::abs(b) + tol?
            return b;
        }

        if (fa != fc && fb != fc) {
            // Inverse quadratic interpolation
            s = (a * fb * fc / ((fa - fb) * (fa - fc))) +
                (b * fa * fc / ((fb - fa) * (fb - fc))) +
                (c * fa * fb / ((fc - fa) * (fc - fb)));
        } else {
            // Secant method
            s = b - fb * (b - a) / (fb - fa);
        }

        Real delta = tol; // Or some small value relative to b

        bool condition1 = (s < (3 * a + b) / 4.0 || s > b); // s not between (3a+b)/4 and b
        bool condition2 = mflag && (std::abs(s - b) >= std::abs(b - c) / 2.0);
        bool condition3 = !mflag && (std::abs(s - b) >= std::abs(c - d) / 2.0);
        bool condition4 = mflag && (std::abs(b - c) < delta);
        bool condition5 = !mflag && (std::abs(c - d) < delta);


        if (condition1 || condition2 || condition3 || condition4 || condition5) {
            // Bisection method
            s = (a + b) / 2.0;
            mflag = true;
        } else {
            mflag = false;
        }

        Real fs = func(s);
         if (std::abs(fs) < tol) return s;

        // Check for invalid function evaluation (e.g., from non-positive f in alpha objective)
        if (fs >= REAL_MAX) {
             // Function evaluation failed, likely means s is out of valid range.
             // Try bisection instead? Or just fail.
             // Let's try bisection towards 'a' if fs is invalid.
             s = (a + b) / 2.0;
             fs = func(s);
             if (fs >= REAL_MAX || std::abs(fs) < tol) {
                 // Still failing or converged at midpoint
                 if (std::abs(fs) < tol) return s;
                 // std::cerr << "Warning (Brent): Function evaluation failed repeatedly." << std::endl;
                 return std::nullopt; // Indicate failure
             }
             mflag = true; // Force bisection next time if needed
        }


        d = c; // d gets old c
        c = b; // c gets old b
        fc = fb;

        if (fa * fs < 0) { // Root is in [a, s]
            b = s;
            fb = fs;
        } else { // Root is in [s, b]
            a = s;
            fa = fs;
        }

        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
         if (iter == max_iter - 1) {
              // std::cerr << "Warning (Brent): Maximum iterations reached. Current estimate: " << b << std::endl;
         }
    }

    return std::nullopt; // Failed to converge within max_iter
}


// Simple 3x3 solver using Cramer's rule (or implement Gaussian elimination)
bool solve_3x3(const std::array<std::array<Real, 3>, 3>& A, const std::array<Real, 3>& b, std::array<Real, 3>& x) {
    Real detA = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
                A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

    if (std::abs(detA) < REAL_EPSILON * REAL_EPSILON) { // Check for singularity
        return false;
    }

    Real invDetA = 1.0 / detA;

    // Calculate determinant of A with column 0 replaced by b
    Real detA0 = b[0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
                 A[0][1] * (b[1] * A[2][2] - A[1][2] * b[2]) +
                 A[0][2] * (b[1] * A[2][1] - A[1][1] * b[2]);

    // Calculate determinant of A with column 1 replaced by b
    Real detA1 = A[0][0] * (b[1] * A[2][2] - A[1][2] * b[2]) -
                 b[0] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                 A[0][2] * (A[1][0] * b[2] - b[1] * A[2][0]);

    // Calculate determinant of A with column 2 replaced by b
    Real detA2 = A[0][0] * (A[1][1] * b[2] - b[1] * A[2][1]) -
                 A[0][1] * (A[1][0] * b[2] - b[1] * A[2][0]) +
                 b[0] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

    x[0] = detA0 * invDetA;
    x[1] = detA1 * invDetA;
    x[2] = detA2 * invDetA;

    return true;
}


} // namespace NumericalSolvers