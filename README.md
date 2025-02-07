# Bairstow's Method Python3/GDScript Implementation

Bairstow's Method for finding polynomial roots in iterative way.

Original Python 2+ recursive code by https://github.com/PowerUpMasters  
Fork by https://github.com/jupiterbjy

This is more optimized, simplified and heavily modernized rewrite
intended for lightweight use cases and help better understanding of algorithm.


## Examples

```python
from bairstow import bairstow, bairstow_complex

_polynomials = [
    [1, -9, 20, -12, 0],
    [6, 11, -33, -44, 11, 6],
    [1, 2, 1],
    [1, 3, 2],
    [1, 2, 3, 4],
    [1, 2, 3, 4, 5],
]

for _polynomial in _polynomials:
    print("\nCoefficients:", _polynomial)

    _roots_dict = {
        "Complex": bairstow_complex(_polynomial),
        "Real": bairstow(_polynomial),
    }

    for _name, _roots in _roots_dict.items():
        print(f"  {_name} roots:")

        for _i, _root in enumerate(_roots):
            print(f"    R{_i} = {_root}")
```

```text
Coefficients: [1, -9, 20, -12, 0]
  Complex roots:
    R0 = 0j
    R1 = (6+0j)
    R2 = (1+0j)
    R3 = (2+0j)
  Real roots:
    R0 = 0.0
    R1 = 1.0
    R2 = 2.0
    R3 = 6.0

Coefficients: [6, 11, -33, -44, 11, 6]
  Complex roots:
    R0 = (-0.28926555025561557+0j)
    R1 = (0.44051767714393747+0j)
  Real roots:
    R0 = -0.28926555025561557
    R1 = 0.44051767714393747

Coefficients: [1, 2, 1]
  Complex roots:
    R0 = (-1+0j)
    R1 = (-1+0j)
  Real roots:
    R0 = -1.0
    R1 = -1.0

Coefficients: [1, 3, 2]
  Complex roots:
    R0 = (-2+0j)
    R1 = (-1+0j)
  Real roots:
    R0 = -2.0
    R1 = -1.0

Coefficients: [1, 2, 3, 4]
  Complex roots:
    R0 = (-0.17468540428030588-1.5468688872313963j)
    R1 = (-0.17468540428030588+1.5468688872313963j)
    R2 = -1.6506291914393882
  Real roots:
    R0 = -1.6506291914393882

Coefficients: [1, 2, 3, 4, 5]
  Complex roots:
    R0 = (-1.287815479557648-0.85789675832849j)
    R1 = (-1.287815479557648+0.85789675832849j)
    R2 = (0.28781547955764863-1.4160930801719085j)
    R3 = (0.28781547955764863+1.4160930801719085j)
  Real roots:
```
