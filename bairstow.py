"""
Bairstow's Method for finding polynomial roots in iterative way.

Original Python 2+ recursive code by https://github.com/PowerUpMasters
Fork by https://github.com/jupiterbjy

This is more optimized, simplified and heavily modernized rewrite
intended for lightweight use cases and help better understanding of algorithm.
"""

import cmath
from typing import Sequence, MutableSequence, List


# --- Utilities ---


def _synthetic_division(
    coefficients: MutableSequence[float], u: float, v: float
) -> List[float]:
    """Perform synthetic division of polynomial via `x^2 + ux + v`

    Args:
        coefficients: Polynomial coefficients
        u: first degree coefficient
        v: constant term

    Returns:
        Resulting polynomial coefficients
    """

    n = len(coefficients)
    result: List[float] = [0] * n

    result[n - 1] = coefficients[n - 1]
    result[n - 2] = coefficients[n - 2] + u * result[n - 1]
    for i in range(n - 3, -1, -1):
        result[i] = coefficients[i] + u * result[i + 1] + v * result[i + 2]

    return result


def _filter_n_sort_complex(complexes: List[complex]) -> List[float]:
    """Filter out complex numbers and return sorted real numbers.

    Args:
        complexes: List of complex numbers

    Returns:
        ascending sorted list of real numbers
    """

    return sorted([c.real for c in complexes if not c.imag])


# --- Logics ---


def bairstow(coefficients: Sequence[float], max_iterations: int = 100) -> List[complex]:
    """
    Bairstow's Method for finding polynomial roots.

    Args:
        coefficients: Coefficients of the polynomial (constant term last, natural order).
        max_iterations: Maximum number of iterations to perform.

    Returns:
        List storing the found roots sorted in ascending order.

    Raises:
        ZeroDivisionError: If calculation encountered division by zero.
    """

    roots: List[complex] = []

    u: float = 0
    v: float = 0

    # initial guess based on leading 3 coefficients as shown on wikipedia
    # if degree is less than 3 then skip, we don't even need these in that case.
    if len(coefficients) >= 3:
        u = coefficients[1] / coefficients[0]
        v = coefficients[2] / coefficients[0]

    # flip coefficients to fit into this algorithm
    coefficients = coefficients[::-1]

    while True:
        max_iterations -= 1
        deg = len(coefficients) - 1

        # perform hardwired calculations for deg < 3
        if deg < 1:
            return roots

        if deg == 1 and coefficients[1] != 0:
            roots.append(-coefficients[0] / coefficients[1])
            return roots

        if deg == 2:
            d = coefficients[1] ** 2 - 4 * coefficients[2] * coefficients[0]
            roots.append((-coefficients[1] - cmath.sqrt(d)) / (2 * coefficients[2]))
            roots.append((-coefficients[1] + cmath.sqrt(d)) / (2 * coefficients[2]))
            return roots

        # deg >= 3, perform Bairstow's method

        # do synthetic divisions
        b = _synthetic_division(coefficients, u, v)
        c = _synthetic_division(b, u, v)

        # Update u & v using linear corrections
        denominator = c[2] * c[2] - c[3] * c[1]
        try:
            u += (c[2] * -b[1] - c[3] * -b[0]) / denominator
            v += (-c[1] * -b[1] + c[2] * -b[0]) / denominator

        except ZeroDivisionError as err:
            raise ValueError("demonstrator was zero") from err

        # Check for convergence
        if (abs(b[0]) > 1e-14 or abs(b[1]) > 1e-14) and max_iterations > 0:
            continue

        # if degree is still large then extract roots of quadratic factor & continue
        if deg >= 3:
            discriminant_sqrt = cmath.sqrt(u**2 - 4 * 1 * (-v))
            roots.append((u - discriminant_sqrt) / 2)
            roots.append((u + discriminant_sqrt) / 2)

            coefficients = b[2:]


# Test run
if __name__ == "__main__":
    _polynomials = [
        [1, -9, 20, -12, 0],
        [6, 11, -33, -44, 11, 6],
        [1, 2, 1],
        [1, 3, 2],
        [1, 2, 3, 4],
        # ^ this should trigger complex roots
    ]

    for _polynomial in _polynomials:
        print("\nCoefficients:", _polynomial)
        _roots = bairstow(_polynomial)

        for i, root in enumerate(_filter_n_sort_complex(_roots)):
            print(f"R{i} = {root}")
