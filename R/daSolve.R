#' Rayleigh’s Method of Dimensional Analysis
#'
#' Performs Rayleigh’s method of dimensional analysis
#'
#' @param dv dependent variable
#' @param iv independent variables
#'
#' @return fundamental dimensions (MLT equations) and the solution
#'
#' @details
#'The study of the relationship between physical quantities with the help of dimensions and units of
#'measurement is termed as dimensional analysis. Dimensional analysis is essential because it keeps the
#'units the same, helping us perform mathematical calculations smoothly. Dimensional Analysis
#'(also called Factor-Label Method or the Unit Factor Method) is a problem-solving method that uses
#'the fact that any number or expression can be multiplied by one without changing its value.
#'
#'Rayleigh’s method of dimensional analysis is a conceptual tool used in physics, chemistry, and
#'engineering. This form of dimensional analysis expresses a functional relationship of some variables
#'in the form of an exponential equation. It was named after Lord Rayleigh.
#'
#'Here are the types and variables:
#'
#' |    Type     |      Variable      |   |
#' | ------------- |:-------------:| -----:|
#' | Geometric       |'length', 'area', 'volume', 'curvature', 'slope', 'angle', 'shape factore', 'diameter', 'distance'
#' | Kinematic      | 'time', 'linear velocity', 'angular velocity', 'velocity','frequency', 'linear acceleration', 'angular acceleration','gravitational acceleration', 'discharge per unit width', 'kinematic viscosity', 'circulation'|
#' | Dynamic      | 'mass', 'force', 'weight', 'density', 'specific weight', 'specific gravity', 'pressure','stress', 'shear stress', 'strain', 'dynamic viscosity', 'surface tension', 'modulus of elasticity', 'compressibility', 'impulse','momentum', 'work', 'energy', 'torque', 'power', 'weight rate of flow', 'viscosity'|
#'
#' @md
#'
#'
#' @examples
#'\donttest{
#'## Example 1:
#'daSolve(dv = "force",
#'        iv = c("mass", "velocity", "length"))
#'
#'## Example 2
#'daSolve(dv = "force",
#'        iv = c("velocity", "diameter",
#'               "density", "viscosity"))
#'}
#'
#' @export
#' @importFrom hash hash
#' @importFrom caracas def_sym as_sym solve_sys symbol
#'

daSolve <- function(dv, iv) {

  masterDict <-
    hash(
      "length" = c(0, 1, 0),
      "area" = c(0, 2, 0),
      "volume" = c(0, 3, 0),
      "curvature" = c(0,-1, 0),
      "slope" = c(0, 0, 0),
      "angle" = c(0, 0, 0),
      "shape factore" = c(0, 0, 0),

      "time" = c(0, 0, 1),
      "linear velocity" = c(0, 1,-1),
      "angular velocity" = c(0, 0,-1),
      "velocity" = c(0, 1,-1),
      "frequency" = c(0, 0,-1),
      "linear acceleration" = c(0, 1,-2),
      "angular acceleration" = c(0, 0,-2),
      "gravitational acceleration" = c(0, 1,-2),
      "discharge per unit width" = c(0, 2,-1),
      "kinematic viscosity" = c(0, 2,-1),
      "circulation" = c(0, 2,-1),
      "diameter" = c(0, 1, 0),

      "mass" = c(1, 0, 0),
      "force" = c(1, 1,-2),
      "weight" = c(1, 1,-2),
      "density" = c(1,-3, 0),
      "specific weight" = c(1,-2,-2),
      "specific gravity" = c(0, 0, 0),
      "pressure" = c(1,-1,-2),
      "stress" = c(1,-1,-2),
      "shear stress" = c(1,-1,-2),
      "strain" = c(0, 1, 0),
      "dynamic viscosity" = c(1,-1,-1),
      "surface tension" = c(1, 0,-2),
      "modulus of elasticity" = c(1,-1,-2),
      "compressibility" = c(-1, 1, 2),
      "impulse" = c(1, 1,-1),
      "momentum" = c(1, 1,-1),
      "work" = c(1, 2,-2),
      "energy" = c(1, 2,-2),
      "torque" = c(1, 2,-2),
      "power" = c(1, 2,-3),
      "weight rate of flow" = c(1, 1,-3),
      "viscosity" = c(1,-1,-1)
    )

  key = c(dv, iv)

  inDict = masterDict[key]

  nx = length(inDict) - 1

  symbs = def_sym(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z)

  mySymbols = symbs[1:nx]

  iniMatrixLft = as.list(inDict[key[!(key %in% key[1])]])[key[-1]]
  iniMatrixRft = as.list(inDict[key[1]])

  mx = iniMatrixLft[[1]][1] * symbol(as.character(mySymbols[[1]]))
  lx = iniMatrixLft[[1]][2] * symbol(as.character(mySymbols[[1]]))
  tx = iniMatrixLft[[1]][3] * symbol(as.character(mySymbols[[1]]))
  rc = 2

  for (rows in iniMatrixLft[-1]) {
    mx = mx + rows[1] * symbol(as.character(mySymbols[[rc]]))
    lx = lx + rows[2] * symbol(as.character(mySymbols[[rc]]))
    tx = tx + rows[3] * symbol(as.character(mySymbols[[rc]]))
    rc = rc + 1
  }

  fx = paste0(as.character(mx), " = ", iniMatrixRft[[1]][1])
  gx = paste0(as.character(lx), " = ", iniMatrixRft[[1]][2])
  hx = paste0(as.character(tx), " = ", iniMatrixRft[[1]][3])

  lhs <- cbind(mx, lx, tx)
  rhs <-
    t(as_sym(c(
      iniMatrixRft[[1]][1], iniMatrixRft[[1]][2], iniMatrixRft[[1]][3]
    )))

  soln <- solve_sys(lhs, rhs, symbs[1:3])

  return(list(
    M = fx,
    L = gx,
    T = hx,
    Soultion = soln
  ))

}
