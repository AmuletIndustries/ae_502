clear all;
close all;

# Part 1
# Need to find a, e, i for a Molniya orbit where perigee is > 600km and period is 8 hours
#Molnia orbits require argument of perigee to be zero under J2, and argument of ascending node to be <<1 under J2

J2 = 0.00108;
Rearth = 6370;
mu0 = 3.986e5;

function rv = rvprop(t, t0, r0, v0, mu4)
  #Computes the r and v state vectors in the 2 body problem given mu, t0 (usually zero), time in seconds, r0 and v0 starting vectors in km and km/s
  #semimajor axis from vis-viva at t0
  a = norm(v0)^2;
  a = a/mu4;
  a = a - (2/norm(r0));
  a = -a;
  a = 1/a;

  #e from h
  h = norm(cross(r0, v0));
  e = h^2;
  e = e/mu4;
  e = e/a;
  e = e - 1;
  e = -e;
  e = sqrt(e);

  #E0 from r0
  E0 = acos(((-norm(r0)/a)+1)/e);
  #Need E from modified keplers equation
  #n(t-t0) + (E0 - esin(E0)) - (E - esin(E)) = 0
  g = 1;
  i = 0;
  n = sqrt(mu4/(a^3));
  E = E0;
  E0term = E0-e*sin(E0);

  while abs(g) > 1e-13
    #loops until tolerance is met, this is newton's methods
    g = n*(t-t0) + E0term - (E - e*sin(E));
    dgdE = -1 + e*cos(E);
    Eg = E - g/dgdE;

    E = Eg;
    i++;
  end
  #i = i;
  #Calculate F and G
  F = (1-(a/norm(r0))*(1-cos(E-E0)));
  G = ((t-t0) - sqrt(a^3 / mu4) * ((E-E0)-sin(E-E0)));

  #Calculate r(t)
  r = F*r0 + G*v0;

  #Calculate Fdot and Gdot
  Fdot = (-(sqrt(mu4*a)/(norm(r)*norm(r0))) * sin(E-E0));
  Gdot = (1-(a/norm(r)) * (1-cos(E-E0)));

  #Calculate v(t)
  v = Fdot*r0 + Gdot*v0;

  rv = [r, v];
end




disp('Part 1:')
fprintf('\n')

#Starting with period of 8 hours we can get a semi-major axis from earths Mu
a = ((sqrt(mu0) * 8 * 3600) / (2 * pi))^(2/3)

#From a we can use the 600km perigee limit and the earths radius, and the semi-major axis, to get an eccentricity
r_per = 600 + Rearth;
c = a - r_per;
e = c/a

#Now we need the mean motion, n, of the spacecraft
n = (2 * pi)/(8 * 3600);  #rad/sec

#Now that we have a and e, we can optimize the equations for argument of perigee and argument of ascending node to get a value of i that works.
#Simplifying the restriction of the argiment of perigee change to be zero, we can remove all the extraneous bits from it and solve directly for an inclination:
i = (1 + acos(0)) #rad

#Then we can plug that into the change of argument of ascending node and test for a low rate of change
asc_node_rate = ((-3/2) * n * J2 * ((Rearth/a)^2) * (cos(i))/(1-(e^2))^2)

disp(['The lowest drift rate is: Omega_dot = (',num2str(asc_node_rate),' rad/sec, which is far less than 1'])

fprintf('\n')
disp('Part 2:')
fprintf('\n')

#We can follow the same basic process to solve this one as the last one.

J2mars = 0.00196;
Rmars = 3390;
mumars = 4.282e4;

#Starting with period of 8 hours we can get a semi-major axis from earths Mu
amars = ((sqrt(mumars) * ((24 * 3600) + (39 * 60) + 35)) / (2 * pi))^(2/3)

#From a we can use the 600km perigee limit and the earths radius, and the semi-major axis, to get an eccentricity
r_perm = 400 + Rmars;
c = a - r_perm;
e = c/a

#Now we need the mean motion, n, of the spacecraft
n = (2 * pi)/((24 * 3600) + (39 * 60) + 35);  #rad/sec

#Now that we have a and e, we can optimize the equations for argument of perigee and argument of ascending node to get a value of i that works.
#Simplifying the restriction of the argiment of perigee change to be zero, we can remove all the extraneous bits from it and solve directly for an inclination:
i = (1 + acos(0)) #rad

#Then we can plug that into the change of argument of ascending node and test for a low rate of change
asc_node_rate = ((-3/2) * n * J2 * ((Rmars/a)^2) * (cos(i))/(1-(e^2))^2)

disp(['The lowest drift rate is: Omega_dot = (',num2str(asc_node_rate),' rad/sec, which is far less than 1'])


fprintf('\n')
disp('Part 3:')
fprintf('\n')

#Perturbed equations of motion
#I can use the variation of each orbital element to determine the rate of change, then plot that effect over time.
#First I need to get h from the given parameters.

#Then implement each variation equation in the book

#Then plot them all over time

#We can convert the orbital elements given into a state vector:
#1. Find E using keplers with M and e
e = 0.74;
M = deg2rad(10);
a = 26600;
i = 1.10654;
Omega = deg2rad(90);
omega = deg2rad(5);

#taken from slides
E0 = M % Guess E0
g = 1;
itr = 0;
while abs(g) > 1e-10
  g = E0 - e * sin(E0) - M;
  dgdE = 1 - e * cos(E0);
  E1 = E0 - g / dgdE;
  % Update
  E0 = E1;
  itr = itr + 1;
end

E = E0;

#2. Find true anomaly with E
f = 2 * atan2(sqrt(1 + e) * tan(E0/2), sqrt(1 - e));

#3. Find Angular Momentum mag from a, e, and mu
h = sqrt(mu0 * a * (1 - (e^2)));

#4. Theta = omega + f
theta = f + omega;

#5. r mag
rmag = ((h^2)/mu0)/(1 + e * cos(f));

#6. r vector
r = rmag*[(cos(theta)*cos(Omega) - cos(i)*sin(Omega)*sin(theta)), (cos(theta)*sin(Omega) + cos(i)*cos(Omega)*sin(theta)), (sin(i)*sin(theta))];

#7. v vector
v = [-(mu0/h)*(cos(Omega)*(sin(theta) + e * sin(omega)) + sin(Omega)*(cos(theta) + e * sin(omega))*cos(i)), -(mu0/h)*(sin(Omega)*(sin(theta) + e * sin(omega)) - cos(Omega)*(cos(theta) + e * sin(omega))*cos(i)), (mu0/h)*((cos(theta) + e * cos(omega))*sin(i))];

#Initial state vector
X0 = [r, v];
t_start = 0;
t_end = 100 * 24 * 3600; % 100 days in seconds
t_step = 100; % Define an appropriate time step in seconds
tspan = t_start:t_step:t_end;
%









function [dXdt] = j2_perturbed_eom(t, X, mu, J2, R_e)
    % Unpack the state vector
    r = X(1:3);
    v = X(4:6);

    % Constants
    r_norm = norm(r);
    r_hat = r / r_norm;

    % J2 Perturbation
    P = (3/2) * J2 * (R_e^2) * mu / (r_norm^4);
    F_J2 = P * [(1 - 5*(r(3)^2)/(r_norm^2))*r(1);
                (1 - 5*(r(3)^2)/(r_norm^2))*r(2);
                (3 - 5*(r(3)^2)/(r_norm^2))*r(3)];

    % Perturbed equations of motion
    a = -mu * r_hat + F_J2;
    dXdt = [v; a];
end



function xdot = eom_twobody(t,x,mue,J2,R)
  #need to compute pj2 here
  #like 90% sure that x1,2,3 is the current r vector
  r = norm([x(1), x(2), x(3)]);
  pJ2 = (3/2) * ((J2 * mue * R^2)/(r^4)) * [((x(1)/r)*(5*(x(3)/r)^2) - 1), ((x(2)/r)*(5*(x(3)/r)^2) - 1), ((x(3)/r)*(5*(x(3)/r)^2) - 3)];

xdot(1) = x(4);
xdot(2) = x(5);
xdot(3) = x(6);
xdot(4) = -mue * x(1) / (norm(x(1:3))^3) + pJ2(1);
xdot(5) = -mue * x(2) / (norm(x(1:3))^3) + pJ2(2);
xdot(6) = -mue * x(3) / (norm(x(1:3))^3) + pJ2(3);

if t == 0
  format long
  pJ2;
  format shortEng
endif
% Output a column vector
xdot = xdot';
end

# Two-Body + J2
options = odeset('AbsTol',1e-10,'RelTol',1e-13);
[t,X] = ode45(@eom_twobody, tspan, X0, options, mu0, J2, Rearth);
#[t,X] = ode45(@j2_perturbed_eom, tspan, X0, options, mu0, J2, Rearth);

#X(end,:);

#part d
#compute orbital elements at each time
r1_mov = zeros(length(X),3);
v1_mov = zeros(length(X),3);
a1_mov = zeros(length(X),1);
e1_mov = zeros(length(X),3);
e1_mov_mag = zeros(length(X),1);
h1_mov = zeros(length(X),3);
h1_mov_mag = zeros(length(X),1);
i1_mov = zeros(length(X),1);
n1_mov = zeros(length(X),3);
n1_mov_mag = zeros(length(X),1);
Omega1_mov = zeros(length(X),1);
omega1_mov = zeros(length(X),1);

itt1 = 1;

while itt1 <= length(X)
  r1_mov(itt1,:) = [X(itt1,1), X(itt1,2), X(itt1,3)];
  v1_mov(itt1,:) = [X(itt1,4), X(itt1,5), X(itt1,6)];
  #1: find a using vis viva
  a1_mov(itt1) = 1/((2/norm(r1_mov(itt1,:))) - ((norm(v1_mov(itt1,:))^2)/mu0));

  #2: find e
  e1_mov(itt1,:) = ((((norm(v1_mov(itt1,:))^2)/mu0) - (1/norm(r1_mov(itt1,:))))*r1_mov(itt1,:)) - ((1/mu0)*(dot(r1_mov(itt1,:), v1_mov(itt1,:)))*v1_mov(itt1,:));

  e1_mov_mag(itt1) = norm(e1_mov(itt1,:));

  #2.5: find h
  h1_mov(itt1,:) = cross(r1_mov(itt1,:),v1_mov(itt1,:));
  h1_mov_mag(itt1) = norm(h1_mov(itt1,:));

  #3: node vector
  K1 = [0, 0, 1];
  n1_mov(itt1,:) = cross(K1, h1_mov(itt1,:));
  n1_mov_mag(itt1) = norm(n1_mov(itt1,:));

  #Inclination
  i1_mov(itt1) = acos(dot((h1_mov(itt1,:)/h1_mov_mag(itt1)), K1));

  #4: RAAN
  I1 = [1, 0, 0];
  J1 = [0, 1, 0];

  Omega1_mov(itt1) = acos((dot(n1_mov(itt1,:), I1))/n1_mov_mag(itt1));

  if dot(n1_mov(itt1,:), J1) < 0
    Omega1_mov(itt1) = 2*pi - Omega1_mov(itt1);
  endif



  #5: argument of periapse
  omega1_mov(itt1) = acos((dot(n1_mov(itt1,:), e1_mov(itt1,:)))/(norm(n1_mov(itt1,:))*norm(e1_mov(itt1,:))));

  if dot(e1_mov(itt1,:), K1) < 0
    omega1_mov(itt1) = 2*pi - omega1_mov(itt1);
  endif

  itt1 = itt1 + 1;
end #end while loop

disp(['Semimajor Axis (a) ', num2str(a1_mov(end))])
disp(['Eccentricity (e) ', num2str(e1_mov_mag(end))])
disp(['Inclination (i) ', num2str(rad2deg(i1_mov(end)))])
disp(['Omega (W) ', num2str(rad2deg(Omega1_mov(end)))])
disp(['omega (w) ', num2str(rad2deg(omega1_mov(end)))])

#plot orbital elements
figure(1)
subplot(5,1,1);
#axis equal
#grid minor
#view(30, 30)
#xlim([0, 5200])
#ylim([20,90])
plot(a1_mov, 'r')
title('Plot of Semimajor Axis (a) for given Molnia orbit with J2 Perturbations over 100 days')
xlabel('Time (sec)')
ylabel('Semimajor Axis (km)')

subplot(5,1,2);
#axis equal
#grid minor
#view(30, 30)
#xlim([0, 5200])
#ylim([20,90])
plot(e1_mov_mag, 'b')
title('Plot of Eccentricity (e) for given Molnia orbit with J2 Perturbations over 100 days')
xlabel('Time (sec)')
ylabel('Eccentricity')

subplot(5,1,3);
#axis equal
#grid minor
#view(30, 30)
#xlim([0, 5200])
#ylim([20,90])
plot(i1_mov, 'g')
title('Plot of Inclination (i) for given Molnia orbit with J2 Perturbations over 100 days')
xlabel('Time (sec)')
ylabel('Inclination (rad)')

subplot(5,1,4);
#axis equal
#grid minor
#view(30, 30)
#xlim([0, 5200])
#ylim([20,90])
plot(Omega1_mov, 'c')
title('Plot of RAAN (Omega) for given Molnia orbit with J2 Perturbations over 100 days')
xlabel('Time (sec)')
ylabel('RAAN (rad)')

subplot(5,1,5);
#axis equal
#grid minor
#view(30, 30)
#xlim([0, 5200])
#ylim([20,90])
plot(Omega1_mov, 'm')
title('Plot of Argument of Periapse (omega) for given Molnia orbit with J2 Perturbations over 100 days')
xlabel('Time (sec)')
ylabel('Argument of Periapse (rad)')

hold off
%





































