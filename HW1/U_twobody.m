% Universal Variable two-body orbit propagator from Curtis 3.3 and 3.4
  function [R, V] = U_twobody(R0, V0, t, Mu)
%{
  This function uses kepler_U, stumpS, stumpC, f_an_g, fDot_and_gDot, and rv_from_r0v0 to solve the universal two-body problem.

  mu   - gravitational parameter (km^3/s^2)
  x    - the universal anomaly (km^0.5)
  dt   - time since x = 0 (s)
  t  - the time elapsed since R0 (s)
  ro   - radial position (km) when x = 0
  vro  - radial velocity (km/s) when x = 0
  a    - reciprocal of the semimajor axis (1/km)
  z    - auxiliary variable (z = a*x^2)
  C    - value of Stumpff function C(z)
  S    - value of Stumpff function S(z)
  n    - number of iterations for convergence
  nMax - maximum allowable number of iterations
%}

%Currently this function just wraps rv_from_r0v0 but in the future I might allow for more complex logic internally.
global mu = Mu;

[R, V] = rv_from_r0v0(R0, V0, t);


% Internal functions I would like to thank Howard D. Curtis for:
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function x = kepler_U(dt, ro, vro, a)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
  This function uses Newton's method to solve the universal
  Kepler equation for the universal anomaly.

  mu   - gravitational parameter (km^3/s^2)
  x    - the universal anomaly (km^0.5)
  dt   - time since x = 0 (s)
  ro   - radial position (km) when x = 0
  vro  - radial velocity (km/s) when x = 0
  a    - reciprocal of the semimajor axis (1/km)
  z    - auxiliary variable (z = a*x^2)
  C    - value of Stumpff function C(z)
  S    - value of Stumpff function S(z)
  n    - number of iterations for convergence
  nMax - maximum allowable number of iterations

  User M-functions required: stumpC, stumpS
%}
% ----------------------------------------------
%global mu

%...Set an error tolerance and a limit on the number of iterations:
error = 1.e-6;
nMax  = 100;

sqrmu = sqrt(mu);

%...Starting value for x:
x = sqrmu*abs(a)*dt;

%...Iterate on Equation 3.65 until until convergence occurs within
%...the error tolerance:
n     = 0;
ratio = 1;
while abs(ratio) > error && n <= nMax
    n     = n + 1;
    C     = stumpC(a*x^2);
    S     = stumpS(a*x^2);
    F     = (ro*vro/sqrmu)*(x^2)*C + (1 - a*ro)*(x^3)*S + (ro*x) - (sqrmu*dt);
    dFdx  = (ro*vro/sqrmu)*x*(1 - a*x^2*S) + (1 - a*ro)*(x^2)*C + ro;
    ratio = F/dFdx;
    x     = x - ratio;
end

%...Deliver a value for x, but report that nMax was reached:
if n > nMax
    %fprintf('\n **No. iterations of Kepler''s equation = %g', n);
    %fprintf('\n   F/dFdx                              = %g\n', F/dFdx);
end

endfunction
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function s = stumpS(z)
% ~~~~~~~~~~~~~~~~~~~~~~
%{
  This function evaluates the Stumpff function S(z) according
  to Equation 3.52.

  z - input argument
  s - value of S(z)

  User M-functions required: none
%}
% ----------------------------------------------

% MODIFIED TO MINIMIZE CALLS TO SQRT()

if z > 0
    sqrz = sqrt(z);
    s = real((sqrz - sin(sqrz))/(sqrz)^3);
elseif z < 0
    sqrz = sqrt(z);
    nsqrz = (sqrz + i);
    s = real((sinh(nsqrz) - nsqrz)/(sqrz)^3);
else
    s = 1/6;
end

endfunction
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function c = stumpC(z)
% ~~~~~~~~~~~~~~~~~~~~~~
%{
  This function evaluates the Stumpff function C(z) according
  to Equation 3.53.

  z - input argument
  c - value of C(z)

  User M-functions required: none
%}
% ----------------------------------------------


if z > 0
    c = real((1 - cos(sqrt(z)))/z);
elseif z < 0
    c = real((cosh(sqrt(-z)) - 1)/(-z));
else
    c = 1/2;
end

endfunction
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function [f, g] = f_and_g(x, t, ro, a)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
  This function calculates the Lagrange f and g coefficients.

  mu - the gravitational parameter (km^3/s^2)
  a  - reciprocal of the semimajor axis (1/km)
  ro - the radial position at time to (km)
  t  - the time elapsed since ro (s)
  x  - the universal anomaly after time t (km^0.5)
  f  - the Lagrange f coefficient (dimensionless)
  g  - the Lagrange g coefficient (s)

  User M-functions required:  stumpC, stumpS
%}
% ----------------------------------------------

%global mu

z = a*x^2;

%...Equation 3.69a:
f = 1 - x^2/ro*stumpC(z);

%...Equation 3.69b:
g = t - 1/sqrt(mu)*x^3*stumpS(z);

endfunction
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function [fdot, gdot] = fDot_and_gDot(x, r, ro, a)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
  This function calculates the time derivatives of the
  Lagrange f and g coefficients.

  mu    - the gravitational parameter (km^3/s^2)
  a     - reciprocal of the semimajor axis (1/km)
  ro    - the radial position at time to (km)
  t     - the time elapsed since initial state vector (s)
  r     - the radial position after time t (km)
  x     - the universal anomaly after time t (km^0.5)
  fdot  - time derivative of the Lagrange f coefficient (1/s)
  gdot  - time derivative of the Lagrange g coefficient (dimensionless)

  User M-functions required:  stumpC, stumpS
%}
% --------------------------------------------------

% global mu

z = a*x^2;

%...Equation 3.69c:
fdot = sqrt(mu)/r/ro*(z*stumpS(z) - 1)*x;

%...Equation 3.69d:
gdot = 1 - x^2/r*stumpC(z);

endfunction
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function [R,V] = rv_from_r0v0(R0, V0, t)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
  This function computes the state vector (R,V) from the
  initial state vector (R0,V0) and the elapsed time.

  mu - gravitational parameter (km^3/s^2)
  R0 - initial position vector (km)
  V0 - initial velocity vector (km/s)
  t  - elapsed time (s)
  R  - final position vector (km)
  V  - final velocity vector (km/s)

% User M-functions required: kepler_U, f_and_g, fDot_and_gDot
%}
% ----------------------------------------------

% global mu

%...Magnitudes of R0 and V0:
r0 = norm(R0);
v0 = norm(V0);

%...Initial radial velocity:
vr0 = dot(R0, V0)/r0;

%...Reciprocal of the semimajor axis (from the energy equation):
alpha = 2/r0 - v0^2/mu;

%...Compute the universal anomaly:
x = kepler_U(t, r0, vr0, alpha);

%...Compute the f and g functions:
[f, g] = f_and_g(x, t, r0, alpha);

%...Compute the final position vector:
R = f*R0 + g*V0;

%...Compute the magnitude of R:
r = norm(R);

%...Compute the derivatives of f and g:
[fdot, gdot] = fDot_and_gDot(x, r, r0, alpha);

%...Compute the final velocity:
V            = fdot*R0 + gdot*V0;

endfunction

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


endfunction

































