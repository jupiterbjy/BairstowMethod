"""
Bairstow's Method for finding polynomial roots in iterative way.

Original Python 2+ recursive code by https://github.com/PowerUpMasters
Fork by https://github.com/jupiterbjy

This is more optimized, simplified and heavily modernized rewrite
intended for lightweight use cases and help better understanding of algorithm.
"""

import math
import cmath
from typing import Sequence, MutableSequence, List


__all__ = ["bairstow", "bairstow_complex"]


# --- Utilities ---


def _synthetic_division(
    coeffs: MutableSequence[float], u: float, v: float
) -> List[float]:
    """Perform synthetic division of polynomial via `x^2 + ux + v`

    Args:
        coeffs: Polynomial Coefficients (constant term first)
        u: first degree coefficient of divisor
        v: constant term of divisor

    Returns:
        Resulting polynomial coeffs in reversed order
    """

    n = len(coeffs)
    result: List[float] = [0] * n

    result[n - 1] = coeffs[n - 1]
    result[n - 2] = coeffs[n - 2] + u * result[n - 1]
    for i in range(n - 3, -1, -1):
        result[i] = coeffs[i] + u * result[i + 1] + v * result[i + 2]

    return result


# --- Logics ---


# noinspection DuplicatedCode
def bairstow_complex(coeffs: Sequence[float], max_iter: int = 50) -> List[complex]:
    """
    Bairstow's Method for finding polynomial roots.

    Args:
        coeffs: Coefficients of the polynomial (constant term last, natural order).
        max_iter: Maximum number of iterations to perform.

    Returns:
        List storing the found roots, none-sorted due to existence of complex roots.

    Raises:
        ZeroDivisionError: If calculation encountered division by zero.
    """

    roots: List[complex] = []

    u: float = 0
    v: float = 0

    # initial guess based on leading 3 coeffs as shown on wikipedia
    # if degree is less than 3 then skip, we don't even need these in that case
    if len(coeffs) >= 3:
        u = coeffs[1] / coeffs[0]
        v = coeffs[2] / coeffs[0]

    # flip coeffs to fit into this algorithm
    coeffs = coeffs[::-1]

    for _ in range(max_iter):
        deg = len(coeffs) - 1

        # perform hardwired calculations for deg < 3
        if deg < 1:
            return roots

        if deg == 1 and coeffs[1] != 0:
            roots.append(-coeffs[0] / coeffs[1])
            return roots

        if deg == 2:
            d_sqrt = cmath.sqrt(coeffs[1] ** 2 - 4 * coeffs[2] * coeffs[0])
            roots.append((-coeffs[1] - d_sqrt) / (2 * coeffs[2]))
            roots.append((-coeffs[1] + d_sqrt) / (2 * coeffs[2]))
            return roots

        # deg >= 3, perform Bairstow's method

        # do synthetic divisions
        b = _synthetic_division(coeffs, u, v)
        c = _synthetic_division(b, u, v)

        # calculate delta and update u & v accordingly
        d = c[2] * c[2] - c[3] * c[1]
        u += (c[2] * -b[1] - c[3] * -b[0]) / d
        v += (-c[1] * -b[1] + c[2] * -b[0]) / d

        # check for convergence or iteration limit, infinite loops do happen
        if abs(b[0]) > 1e-14 or abs(b[1]) > 1e-14:
            continue

        # if degree is still large, then extract roots of a quadratic factor & continue
        if deg >= 3:
            d_sqrt = cmath.sqrt(u**2 - 4 * (-v))
            roots.append((u - d_sqrt) / 2)
            roots.append((u + d_sqrt) / 2)

            coeffs = b[2:]

    # welp couldn't break within max_iter
    return roots


# noinspection DuplicatedCode
def bairstow(coeffs: Sequence[float], max_iter: int = 50) -> List[float]:
    """
    Bairstow's Method for finding polynomial roots, non-complex roots only.

    Args:
        coeffs: Coefficients of the polynomial (constant term last, natural order).
        max_iter: Maximum number of iterations to perform.

    Returns:
        List storing the found roots, sorted in ascending order.

    Raises:
        ZeroDivisionError: If calculation encountered division by zero.
    """

    roots: List[float] = []

    u: float = 0
    v: float = 0

    # initial guess based on leading 3 coeffs as shown on wikipedia
    # if degree is less than 3 then skip, we don't even need these in that case
    if len(coeffs) > 3:
        u = coeffs[1] / coeffs[0]
        v = coeffs[2] / coeffs[0]

    # flip coeffs to fit into this algorithm
    coeffs = coeffs[::-1]

    for _ in range(max_iter):
        deg = len(coeffs) - 1

        # perform hardwired calculations for deg < 3
        if deg < 1:
            return sorted(roots)

        if deg == 1 and coeffs[1] != 0:
            roots.append(-coeffs[0] / coeffs[1])
            return sorted(roots)

        if deg == 2:
            d = coeffs[1] ** 2 - 4 * coeffs[2] * coeffs[0]

            if d >= 0:
                d_sqrt = math.sqrt(d)
                roots.append((-coeffs[1] - d_sqrt) / (2 * coeffs[2]))
                roots.append((-coeffs[1] + d_sqrt) / (2 * coeffs[2]))
            return sorted(roots)

        # deg >= 3, perform Bairstow's method

        # do synthetic divisions
        b = _synthetic_division(coeffs, u, v)
        c = _synthetic_division(b, u, v)

        # update u & v
        d = c[2] * c[2] - c[3] * c[1]
        u += (c[2] * -b[1] - c[3] * -b[0]) / d
        v += (-c[1] * -b[1] + c[2] * -b[0]) / d

        # check for convergence or iteration limit, infinite loops do happen
        if abs(b[0]) > 1e-14 or abs(b[1]) > 1e-14:
            continue

        # if degree is still large, then extract roots of a quadratic factor & continue
        if deg >= 3:
            d = u**2 - 4 * (-v)

            if d >= 0:
                d_sqrt = math.sqrt(d)
                roots.append((u - d_sqrt) / 2)
                roots.append((u + d_sqrt) / 2)

            coeffs = b[2:]

    # welp couldn't break within max_iter
    return sorted(roots)


# --- Tests ---

if __name__ == "__main__":
    _polynomials = [
        [1, -9, 20, -12, 0],
        [6, 11, -33, -44, 11, 6],
        [1, 2, 1],
        [1, 3, 2],
        [1, 2, 3, 4],
        [1, 2, 3, 4, 5],
    ]

    for _polynomial in _polynomials:
        print("\ncoeffs:", _polynomial)

        _roots_dict = {
            "Complex": bairstow_complex(_polynomial),
            "Real": bairstow(_polynomial),
        }

        for _name, _roots in _roots_dict.items():
            print(f"  {_name} roots:")

            for _i, _root in enumerate(_roots):
                print(f"    R{_i} = {_root}")
