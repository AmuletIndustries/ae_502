# SGP4 Orbital determination
clear all; clc
csv_filename = "AE502_HW4.csv";
mu = 398600.435507;
test = 0;

#Step 1
 #  Find observer position (WGS84)
obsv_pos = zeros(9, 5);    %MJD (1), Sidereal Time (1), R (3). Makes a 9 x 5 matrix which stores the times and positions of the observer at each observation
obsv_data = dlmread(csv_filename, ",");

obsv_pos(:,1) = obsv_data(2:10, 5);   #Store the MJD time from the CSV into position matrix

 function localSiderealTime = mjd2siderealTime(mjd)
  % Convert MJD time to sidereal time

  % Input:
  %   mjd - Modified Julian Date

  % Output:
  %   siderealTime - Sidereal time in hours

  % Initialize variables

  jd = mjd + 2400000.5;

  % Calculate sidereal time

 T0 = (jd - 2451545.0) / 36525.0;

 #Universal time of the given MJD, in fractions of a day which is what is needed here
  UT_frac = mod(mjd, 1);

  GSdT = 100.4606184 + (36000.77004*T0) + (0.000387933*(T0^2)) - ((2.583e-8) * (T0^3));
  while GSdT < 0
    GSdT = GSdT + 360;
  endwhile
  if GSdT > 360
    GSdT = mod(GSdT, 360);
  endif

  GSdT = GSdT + (360.98564724 * UT_frac);  #UTC -5 for Apr 2023 in CU

  SdT = GSdT + (-88.2434);   # Add angle east, so negative for CU

    while SdT < 0
    SdT = SdT + 360;
    endwhile
  if SdT > 360
    SdT = mod(SdT, 360);
  endif

  localSiderealTime = SdT;

end

#Solve for sidereal time
obsv_pos(1, 2) = mjd2siderealTime(obsv_pos(1, 1));
obsv_pos(2, 2) = mjd2siderealTime(obsv_pos(2, 1));
obsv_pos(3, 2) = mjd2siderealTime(obsv_pos(3, 1));
obsv_pos(4, 2) = mjd2siderealTime(obsv_pos(4, 1));
obsv_pos(5, 2) = mjd2siderealTime(obsv_pos(5, 1));
obsv_pos(6, 2) = mjd2siderealTime(obsv_pos(6, 1));
obsv_pos(7, 2) = mjd2siderealTime(obsv_pos(7, 1));
obsv_pos(8, 2) = mjd2siderealTime(obsv_pos(8, 1));
obsv_pos(9, 2) = mjd2siderealTime(obsv_pos(9, 1));

#Solve for inertial pos vector
f = 1/298.257223563; #flattening
for i = 1:9
  obsv_pos(i, 3) = (((6371.0087714 / (sqrt(1- (2*f - f^2) * sind(40.1164)^2))) + 0.233)* cosd(40.1164)*cosd(obsv_pos(i, 2)));    #I
  obsv_pos(i, 4) = (((6371.0087714 / (sqrt(1- (2*f - f^2) * sind(40.1164)^2))) + 0.233)* cosd(40.1164)*sind(obsv_pos(i, 2)));    #J
  obsv_pos(i, 5) = ((((6371.0087714 * (1-f)^2) / (sqrt(1- (2*f - f^2) * sind(40.1164)^2))) + 0.233)* sind(40.1164));    #K
end


#Step 2
  #Initial Orbit Determination
  #Everything here needs to be done in a function so that it can be called while rolling
function RV = gauss_method(ra1, ra2, ra3, dec1, dec2, dec3, stime1, stime2, stime3, R1, R2, R3, time12, time13)
  #Gauss' method for initial orbital determination, uses three input data values and generates RV for them.
  #Page 276 in Curtis, Algorithm 5.11
  mu = 398600.435507;

  #disp("rho_hats: \n");
  rho1_hat = [cosd(dec1)*cosd(ra1), cosd(dec1)*sind(ra1), sind(dec1)];
  rho2_hat = [cosd(dec2)*cosd(ra2), cosd(dec2)*sind(ra2), sind(dec2)];
  rho3_hat = [cosd(dec3)*cosd(ra3), cosd(dec3)*sind(ra3), sind(dec3)];

  #5.11 step 1
  #disp("Step 1: \n");
  #tau1 = (stime1 - stime2)*24*60*60
  #tau3 = (stime3 - stime2)*24*60*60
  tau1 = 0 - time12;
  tau3 = time13 - time12;
  tau = tau3 - tau1;

  #5.11 step 2
  #disp("Step 2: \n");
  p1 = cross(rho2_hat, rho3_hat);
  p2 = cross(rho1_hat, rho3_hat);
  p3 = cross(rho1_hat, rho2_hat);

#5.11 step 3
#disp("Step 3: \n");
D0 = dot(rho1_hat, p1);

#5.11 step 4
#disp("Step 4: \n");
D11 = dot(R1, p1);
D12 = dot(R1, p2);
D13 = dot(R1, p3);

D21 = dot(R2, p1);
D22 = dot(R2, p2);
D23 = dot(R2, p3);

D31 = dot(R3, p1);
D32 = dot(R3, p2);
D33 = dot(R3, p3);

#5.11 step 5
#disp("Step 5: \n");
#error starts to be visible here
A = (1/D0) * ((-D12 * tau3/tau) + D22 + (D32*tau1/tau));

B = (1/(6*D0)) * ((D12 * (tau3^2 - tau^2)*tau3/tau) + (D32 * (tau^2 - tau1) * tau1/tau));

#5.11 step 6
#disp("Step 6: \n")
E = dot(R2, rho2_hat);
R2_2 = dot(R2, R2);

#5.11 step 7
#disp("Step 7: \n");
a = -(A^2 + 2*A*E + R2_2);
b = -2*mu*B*(A+E);
c = -(mu^2)*B^2;

#5.11 step 8
#disp("Step 8: \n");
eqn5116 = [1, 0, a, 0, 0, b, 0, 0, c];
eqn5116_roots = roots(eqn5116);

#disp("The roots of equation 5.116 for this run of the Gauss solver are:\n")
#eqn5116_roots
ans_found = 0;
problem = 0;
for idx = 1:rows(eqn5116_roots)
  #print formatted stuff
  if (ans_found == 1) && isreal(eqn5116_roots(idx)) && (eqn5116_roots(idx) > 0)
    #This is also likely the answer, requires human intervention
    problem = 1;
  endif
  if (ans_found == 0) && isreal(eqn5116_roots(idx)) && (eqn5116_roots(idx) > 0)
    #This is likely the answer
    ans_found = 1;
    best_choice_idx = idx;
  endif
end

if problem
  for idx = 1:rows(eqn5116_roots)
    #print formatted stuff
    fprintf('\n Index=%g: %g', idx, eqn5116_roots(idx));
  end
  best_choice_idx = input("The automated root finder has detected more than one root. Please enter the index of the best choice root: ");
end

if problem == 0
  #No problem! Roots finder worked probably!
  #disp("The roots finder found this to be the best root of equation 5.116 for this run of the Gauss solver:")
  #fprintf('\n Index=%g: %g\n', best_choice_idx, eqn5116_roots(best_choice_idx));
end

r2 = eqn5116_roots(best_choice_idx);

#5.11 step 9
#disp("Step 9: \n");
  #eqn 5.113
  rho1 = (1/D0) * (((6 * ((D31 * tau1/tau3) + (D21 * tau/tau3)) * r2^3) + (mu*D31*(tau^2 - tau1^2) * tau1/tau)) / (6*(r2^3) + mu*(tau^2 - tau3^2)) - D11);

  #eqn 5.112a
  rho2 = A + (mu*B)/(r2^3);

  #eqn 5.114
  rho3 = (1/D0) * ((((6 * ((D13 * tau3/tau1) - (D23 * tau/tau1)) * r2^3) + (mu*D13*(tau^2 - tau3^2) * (tau3/tau1))) / (6*(r2^3) + mu*(tau^2 - tau1^2))) - D33);

#5.11 step 10
#disp("Step 10: \n");
#things are close but not like, that close
#eqn 5.86
br1 = R1 + rho1*rho1_hat;
br2 = R2 + rho2*rho2_hat;
br3 = R3 + rho3*rho3_hat;

#5.11 step 11
#disp("Step 11: \n");
f1 = 1-(0.5 * (mu/(r2^3)) * (tau1^2));
f3 = 1-(0.5 * (mu/(r2^3)) * (tau3^2));

g1 = tau1 - ((1/6) * (mu/(r2^3)) * (tau1^3));
g3 = tau3 - ((1/6) * (mu/(r2^3)) * (tau3^3));

#5.11 step 12
#disp("Step 12: \n");
v2 = (1/((f1*g3) - (f3*g1))) * ((-f3 * br1) + (f1 * br3));

#5.11 step 13
#N/A

RV = [br2(1), br2(2), br2(3), v2(1), v2(2), v2(3)];

#Basically this isn't super far off but it's not spot on either

endfunction

#Step 2.5
#Parameter generation
#Take RV and generate parameters from that

#elements = [e, a, i, Omega, omega, f], this is a 1x6 matrix
function elements = rv2elements(RV)
  mu = 398600.435507;

  #elements = [e, a, i, Omega, omega, f], this is a 1x6 matrix
  r1 = [RV(1), RV(2), RV(3)];
  v1 = [RV(4), RV(5), RV(6)];

  #1: find a using vis viva
  a1 = 1/((2/norm(r1)) - ((norm(v1)^2)/mu));

  #2: find e
  e1 = ((((norm(v1)^2)/mu) - (1/norm(r1)))*r1) - ((1/mu)*(dot(r1, v1))*v1);

  norm(e1);

  #2.5: find h
  h1 = cross(r1,v1);

  #3: node vector
  K1 = [0, 0, 1];

  i1 = acos(dot((h1/norm(h1)), K1));

  n1 = cross(K1, h1);
  norm(n1);

  #4: RAAN
  I1 = [1, 0, 0];
  J1 = [0, 1, 0];

  Omega1 = acos((dot(n1, I1))/norm(n1));
  if dot(n1, J1) < 0
    Omega1 = 2*pi - Omega1;
  endif

  #5: argument of periapse
  omega1 = acos((dot(n1, e1))/(norm(n1)*norm(e1)));
  if dot(e1, K1) < 0
    omega1 = 2*pi - omega1;
  endif

  #6: true anomaly
  f1 = acos((dot(r1, e1))/(norm(r1)*norm(e1)));

  if dot(r1, v1) < 0
    f1 = 2*pi - f1;
  endif

  elements = [norm(e1), a1, rad2deg(i1), rad2deg(Omega1), rad2deg(omega1), rad2deg(f1)];

endfunction


#Step 3
  #Differential Correction
  #Could probably do this myself by taking the set of RV from Step 2 and performing an average of their resulting parameters
#Test of functions
if test == 1
  time12 = 118.10;
  time13 = 237.58;
  test_RV = gauss_method(43.537, 54.420, 64.318, -8.7833, -12.074, -15.105, 44.506, 45.000, 45.499, [3489.8, 3430.2, 4078.5], [3460.1, 3460.1, 4078.5], [3429.9, 3490.1, 4078.5], time12, time13);
  test_elements = rv2elements(test_RV);
endif


#Loop through different combinations of elements
  #7 loops through to get each shingled block of tracking data
RV = [0,0,0,0,0,0];
if test != 1
for itt = 1:7
  #Get the current data
  #Get the current position
  #Use the gauss_method function
  time12 = (obsv_data(itt+2, 5) - obsv_data(itt+1, 5)) * 24 * 60 * 60;    #Seconds between observation 1 and 2
  time13 = (obsv_data(itt+3, 5) - obsv_data(itt+1, 5)) * 24 * 60 * 60;    #Seconds between observation 1 and 3;

##  if itt == 1
##    RV = gauss_method(obsv_data(itt+1,1), obsv_data(itt+2,1), obsv_data(itt+3,1), obsv_data(itt+1,2), obsv_data(itt+2,2), obsv_data(itt+3,2), obsv_pos(itt, 2), obsv_pos(itt+1, 2), obsv_pos(itt+2, 2), obsv_pos(itt, 3:5), obsv_pos(itt+1, 3:5), obsv_pos(itt+2, 3:5), time12, time13);
##  elseif (itt == 4 || itt == 7)
##    RV = vertcat(RV, gauss_method(obsv_data(itt+1,1), obsv_data(itt+2,1), obsv_data(itt+3,1), obsv_data(itt+1,2), obsv_data(itt+2,2), obsv_data(itt+3,2), obsv_pos(itt, 2), obsv_pos(itt+1, 2), obsv_pos(itt+2, 2), obsv_pos(itt, 3:5), obsv_pos(itt+1, 3:5), obsv_pos(itt+2, 3:5), time12, time13));
##  else #(itt != 1 || itt != 4 || itt != 7)
##    RV = vertcat(RV, zeros(1, 6));
##  endif
  #fprintf('\n Done with gauss_method, itt: %g', itt);

  #pass that to the rv2elements function
  if itt == 1
    #elements = rv2elements(RV);
    elements = Example_5_11(obsv_data(itt+1,1), obsv_data(itt+2,1), obsv_data(itt+3,1), obsv_data(itt+1,2), obsv_data(itt+2,2), obsv_data(itt+3,2), obsv_pos(itt, 2), obsv_pos(itt+1, 2), obsv_pos(itt+2, 2), time12, time13);
  elseif (itt == 4 || itt == 7)
    #elements = vertcat(elements,  rv2elements(RV(itt, :)));
    elements = vertcat(elements,  Example_5_11(obsv_data(itt+1,1), obsv_data(itt+2,1), obsv_data(itt+3,1), obsv_data(itt+1,2), obsv_data(itt+2,2), obsv_data(itt+3,2), obsv_pos(itt, 2), obsv_pos(itt+1, 2), obsv_pos(itt+2, 2), time12, time13));
  else
    #elements = vertcat(elements,  zeros(1, 6));
  endif
  #fprintf('\n Done with rv2elements, itt: %g\n', itt);
end
endif

#Average the orbital elements
avg_results = mean([elements(1,:); elements(2,:); elements(3,:)]);

















