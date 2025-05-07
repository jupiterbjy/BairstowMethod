class_name Bairstow
## Math class with ballistics, polynominal solvers & utilities.
## Would really appreciate if someone rewrite this into GDNative C++ plugin.
##
## https://github.com/jupiterbjy/BairstowMethod


# --- References ---

static var GRAVITY_VEC: Vector3 = (
	ProjectSettings.get_setting("physics/3d/default_gravity_vector")
	* ProjectSettings.get_setting("physics/3d/default_gravity")
)


# --- Utilities ---

## Perform synthetic division. Pass inverted coeffs.
static func synthetic_division(
	coeffs: PackedFloat32Array, u: float, v: float
) -> PackedFloat32Array:

	var n := len(coeffs)
	var result: PackedFloat32Array
	result.resize(n)

	result[n - 1] = coeffs[n - 1]
	result[n - 2] = coeffs[n - 2] + u * result[n - 1]
	for i in range(n - 3, -1, -1):
		result[i] = coeffs[i] + u * result[i + 1] + v * result[i + 2]

	return result


## Bairstow polynomial solver, assuming degree >= 3
## Coefficients' order starts from constant term so make sure to reverse it.
static func bairstow(
	coeffs: PackedFloat32Array, u: float, v: float, max_iter: int = 30, eps: float = 1e-6
) -> Array[float]:
	# https://github.com/jupiterbjy/BairstowMethod/blob/master/bairstow.py

	var roots: Array[float]

	# used later in calc, just for keeping duplicated ref
	var d: float
	var d_sqrt: float
	var deg: int
	var b: PackedFloat32Array
	var c: PackedFloat32Array

	for _step: int in range(max_iter):
		#print("%s\n%s\n%s\n%s\n%s\n\n" % [coeffs, u, v, b, c])
		deg = len(coeffs) - 1

		if deg < 1:
			#print("Took %s cycle" % _step)
			roots.sort()
			return roots

		if deg == 1 and coeffs[1] != 0:
			roots.append(-coeffs[0] / coeffs[1])

			#print("Took %s cycle" % _step)
			roots.sort()
			return roots

		if deg == 2:
			d = coeffs[1] ** 2.0 - 4.0 * coeffs[2] * coeffs[0]
			if d >= 0.0:
				d_sqrt = sqrt(d)
				roots.append((-coeffs[1] - d_sqrt) / (2.0 * coeffs[2]))
				roots.append((-coeffs[1] + d_sqrt) / (2.0 * coeffs[2]))

			#print("Took %s cycle" % _step)
			roots.sort()
			return roots

		# deg >= 3, perform Bairstow's method

		# do synthetic divisions
		b = synthetic_division(coeffs, u, v)
		c = synthetic_division(b, u, v)

		# calculate delta and update u & v accordingly
		d = c[2] * c[2] - c[3] * c[1]
		u += (c[2] * -b[1] - c[3] * -b[0]) / d
		v += (-c[1] * -b[1] + c[2] * -b[0]) / d

		# check for convergence or iteration limit, infinite loops do happen
		if abs(b[0]) > eps or abs(b[1]) > eps:
			continue

		# if degree is still large then extract roots of quadratic factor & continue
		if deg >= 3:
			d = u ** 2.0 - 4 * (-v)
			if d >= 0.0:
				d_sqrt = sqrt(d)
				roots.append((u - d_sqrt) / 2.0)
				roots.append((u + d_sqrt) / 2.0)

			coeffs = b.slice(2)

	# welp couldn't break within max_iter
	#print("Took %s cycle" % max_iter)
	roots.sort()
	return roots


## Calculate ballistic arc. Returns all possible answers.
static func calc_ballistic_arc(
	effective_tgt_pos: Vector3,
	effective_tgt_vel: Vector3,
	proj_speed: float,
	gravity: Vector3 = GRAVITY_VEC,
	max_iter: int = 50,
	eps: float = 1e-8,
) -> PackedFloat32Array:
	# https://www.forrestthewoods.com/blog/solving_ballistic_trajectories/
	# https://playtechs.blogspot.com/2007/04/aiming-at-moving-target.html

	# prep reversed coefficents equation
	var coeffs: PackedFloat32Array = [
		effective_tgt_pos.length_squared(),
		2.0 * (effective_tgt_pos.dot(effective_tgt_vel)),
		(effective_tgt_pos.dot(gravity) + effective_tgt_vel.length_squared() - proj_speed ** 2),
		effective_tgt_vel.dot(gravity),
		0.25 * gravity.length_squared()
	]

	# if you're using this in game then you pretty much have more than 3 deg poly to solve
	#if len(coeffs) >= 3:
		#u = coeffs[1] / coeffs[0]
		#v = coeffs[2] / coeffs[0]

	# make a wikipedia guess
	var u: float = coeffs[-2] / coeffs[-1]
	var v: float = coeffs[-3] / coeffs[-1]

	# get all(actually maximum two) possible projectile travel time
	var all_ts := bairstow(coeffs, u, v, max_iter, eps)

	# if fails to converge within iteration limit give it last ditch effort with u=1 v=1
	if not all_ts:
		all_ts = bairstow(coeffs, 1.0, 1.0, max_iter, eps)

	# filter out all negs
	var roots: Array[float]
	for t: float in all_ts:
		if t >= 0.0:
			roots.append(t)

	return roots


## Calculate firing solution in global space in relative to shooter to hit target at given proj.
static func calc_aim_from_t(
	effective_tgt_pos: Vector3,
	effective_tgt_vel: Vector3,
	travel_time: float,
	gravity: Vector3 = GRAVITY_VEC,
) -> Vector3:

	return (
		effective_tgt_pos + (effective_tgt_vel * travel_time) - (0.5 * gravity * travel_time ** 2)
	)
