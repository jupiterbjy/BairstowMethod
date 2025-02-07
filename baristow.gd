class_name Baristow
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

## Perform synthetic division. Pass inverted coefficients.
static func synthetic_division(
	coefficients: PackedFloat32Array, u: float, v: float
) -> PackedFloat32Array:

	var n := len(coefficients)
	var result: PackedFloat32Array
	result.resize(n)

	result[n - 1] = coefficients[n - 1]
	result[n - 2] = coefficients[n - 2] + u * result[n - 1]
	for i in range(n - 3, -1, -1):
		result[i] = coefficients[i] + u * result[i + 1] + v * result[i + 2]

	return result


## Bairstow polynomial solver, assuming degree >= 3
## Coefficients order starts from constant term so make sure to reverse it.
static func bairstow(coefficients: PackedFloat32Array, max_iterations: int = 15) -> Array[float]:
	# https://github.com/jupiterbjy/BairstowMethod/blob/master/bairstow.py

	var roots: Array[float]

	var u: float = coefficients[-2] / coefficients[-1]
	var v: float = coefficients[-3] / coefficients[-1]

	# used later in calc, just for keeping duplicated ref
	var d: float
	var d_sqrt: float
	var deg: int
	var b: PackedFloat32Array
	var c: PackedFloat32Array

	#if len(coefficients) >= 3:
		#u = coefficients[1] / coefficients[0]
		#v = coefficients[2] / coefficients[0]

	for _step: int in range(max_iterations):
		#print("%s\n%s\n%s\n%s\n%s\n\n" % [coefficients, u, v, b, c])
		deg = len(coefficients) - 1

		if deg < 1:
			#print("Took %s cycle" % _step)
			roots.sort()
			return roots

		if deg == 1 and coefficients[1] != 0:
			roots.append(-coefficients[0] / coefficients[1])

			#print("Took %s cycle" % _step)
			roots.sort()
			return roots

		if deg == 2:
			d = coefficients[1] ** 2.0 - 4.0 * coefficients[2] * coefficients[0]
			if d >= 0.0:
				d_sqrt = sqrt(d)
				roots.append((-coefficients[1] - d_sqrt) / (2.0 * coefficients[2]))
				roots.append((-coefficients[1] + d_sqrt) / (2.0 * coefficients[2]))

			#print("Took %s cycle" % _step)
			roots.sort()
			return roots

		# deg >= 3, perform Bairstow's method

		# do synthetic divisions
		b = synthetic_division(coefficients, u, v)
		c = synthetic_division(b, u, v)

		# Update u & v
		d = c[2] * c[2] - c[3] * c[1]
		u += (c[2] * -b[1] - c[3] * -b[0]) / d
		v += (-c[1] * -b[1] + c[2] * -b[0]) / d

		# Check for convergence or iteration limit, infinite loop do happens.
		if abs(b[0]) > 1e-8 or abs(b[1]) > 1e-8:
			continue

		# if degree is still large then extract roots of quadratic factor & continue
		if deg >= 3:
			d = u ** 2.0 - 4 * (-v)
			if d >= 0.0:
				d_sqrt = sqrt(d)
				roots.append((u - d_sqrt) / 2.0)
				roots.append((u + d_sqrt) / 2.0)

			coefficients = b.slice(2)

	# welp couldn't break within max_iterations
	#print("Took %s cycle" % max_iterations)
	roots.sort()
	return roots


## Calculate ballistic arc. Returns all possible answers.
static func calc_ballistic_arc(
	shooter_pos: Vector3,
	shooter_vel: Vector3,
	tgt_pos: Vector3,
	tgt_vel: Vector3,
	proj_speed: float,
	gravity: Vector3 = GRAVITY_VEC,
) -> PackedFloat32Array:
	# https://playtechs.blogspot.com/2007/04/aiming-at-moving-target.html

	# calculate relative position
	var rel_tgt_pos := tgt_pos - shooter_pos
	var rel_tgt_vel := tgt_vel - shooter_vel

	# prep equation
	var coefficients_reversed: PackedFloat32Array = [
		rel_tgt_pos.length_squared(),
		2.0 * (rel_tgt_pos.dot(rel_tgt_vel)),
		(rel_tgt_pos.dot(gravity) + rel_tgt_vel.length_squared() - proj_speed ** 2),
		rel_tgt_vel.dot(gravity),
		0.25 * gravity.length_squared()
	]

	# get all(actually maximum two) possible projectile travel time
	var all_ts := bairstow(coefficients_reversed)

	# filter out all negs
	var roots: Array[float]
	for t: float in all_ts:
		if t >= 0.0:
			roots.append(t)

	return roots


## Calculate global-space firing solution to hit target at given proj. flight time
static func calc_aim_from_t(
	shooter_pos: Vector3,
	shooter_vel: Vector3,
	tgt_pos: Vector3,
	tgt_vel: Vector3,
	travel_time: float,
	gravity: Vector3 = GRAVITY_VEC,
) -> Vector3:
	# https://playtechs.blogspot.com/2007/04/aiming-at-moving-target.html

	# calculate relative position
	var rel_tgt_pos := tgt_pos - shooter_pos
	var rel_tgt_vel := tgt_vel - shooter_vel

	var aim: Vector3 = rel_tgt_pos + (rel_tgt_vel * travel_time) - (0.5 * gravity * travel_time ** 2)
	return shooter_pos + aim
