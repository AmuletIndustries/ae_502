clear all;
close all;

pkg load parallel

parcellfun_set_nproc(0);

% Part 1
%Uses my U_twobody function to solve for the position and velocity given a time elapsed from the position and velocity epoch, and mu
%Uses curtis lambert.m function to solve lamberts problem

from_file = 0;    #Bool that tells program to pull from file instead of calculating display vectors itself
to_file3 = 0;      #Bool that tells program to save the display vector data to a file
to_file4 = 0;      #Bool that tells program to save the display vector data to a file
do_lambert3 = 0;   %bool to control lambert loop
do_lambert4 = 0;   %bool to control lambert loop
lambert_to_file3 = 0;  %bool for controlling saving of lambert solutions to file.
lambert_to_file4 = 0;  %bool for controlling saving of lambert solutions to file.
lambert_from_file = 0;%bool to control pulling lambert solutions from file.
file_name = 'object_positions.mat';
lambert_file_name = 'lambert_solutions.mat';
lambert_file_name2 = 'lamber_solutions.mat';

%control flags to produce the displays for q3 and q4
make_display_1 = 0;
make_display_2 = 0;

%1I State Vectors
r1Iinit = [3.515868886595499e-2, -3.162046390773074, 4.493983111703389];
v1Iinit = [-2.317577766980901e-3, 9.843360903693031e-3, -1.541856855538041e-2];

%Earth State Vectors
rEinit = [-1.796136509111975e-1, 9.667949206859814e-1, -3.668681017942158e-5];
vEinit = [-1.720038360888334e-2, -3.211186197806460e-3, 7.927736735960840e-7];

%2I State Vectors
r2Iinit = [7.249472033259724, 14.61063037906177, 14.24274452216359];
v2Iinit = [-8.241709369476881e-3, -1.156219024581502e-2, -1.317135977481448e-2];

Musun = 1.32712440018e20;

%AU to km
r1Iinit = r1Iinit * 149597870.691;
rEinit = rEinit * 149597870.691;
r2Iinit = r2Iinit * 149597870.691;

%au/day to km/sec
v1Iinit = (v1Iinit * 149597870.691)/(24*3600)
v2Iinit = (v2Iinit * 149597870.691)/(24*3600)
vEinit = (vEinit * 149597870.691)/(24*3600);

% Need to generate two pork chop plots (date from earth departure) from:
%  Jan 2017 to Dec 2017
% vs arrival dates in range:
%  Aug 2017 to Jan 2019
% Where total DeltaV is color coded
%One plot is for fly-by mission (<20km/s total), other plot is for rendezvous (<50 km/s)

q3data = zeros(364, 548, 2);    %dday X aday X DVflyby, DVrendez
q4data = zeros(1307, 1856, 2);    %dday X aday X DVflyby, DVrendez

object_positions = zeros(1307, 1856, 18);    %dday X aday X rE, vE, rI1, vI1, r2I, 2I

if to_file3 == 1
  %load object_positions so we don't lose data
  load(file_name);

  % if we're saving the position data to the file, we need to compute it now.
  for dday = 1:364    %Days from Jan 1 2017 to Dec 31 2017
    %Loop through each day of arrival
    fprintf('\n dday: %g', dday);
    for aday = 212:760  %Days from Aug 1 2017 to Jan 31 2019 dated such that 1 is Jan 1 2017
      % Escape impossible arrival before departure
      if aday >= (dday + 60)
        %Determine location of Earth at departure date
        if dday == 1
          %Special case where initial date is correct and doesn't need to be computed
          rE = rEinit;
          vE = vEinit;
        else
          %Use U_twobody to determine Earth vectors at departure date
          tsec = ((dday - 1)*24*3600);    %Seconds between epoch and departure date
          [rE, vE] = U_twobody(rEinit, vEinit, tsec, Musun);
          %fprintf('\n rE: %g, vE: %g', rE, vE);
        endif
        %Determine location of Target at arrival date
        tsec = (aday*24*3600);    %Seconds between epoch and departure date
        %[r1I, v1I] = pararrayfun(15, @(r, v, t, Mu) U_twobody(r, v, t, Mu), r1Iinit, v1Iinit, tsec, Musun, "Vectorized", true, "ChunksPerProc", 1);
        [r1I, v1I] = U_twobody(r1Iinit, v1Iinit, tsec, Musun);
        %fprintf('\n r1I: %g, v1I: %g', r1I, v1I);

        %fprintf('\n dday: %g, aday: %g', dday, aday);
        %store data to appropriate location
        object_positions(dday, (aday-211), 1) = rE(1);
        object_positions(dday, (aday-211), 2) = rE(2);
        object_positions(dday, (aday-211), 3) = rE(3);
        object_positions(dday, (aday-211), 4) = vE(1);
        object_positions(dday, (aday-211), 5) = vE(2);
        object_positions(dday, (aday-211), 6) = vE(3);

        object_positions(dday, (aday-211), 7) = r1I(1);
        object_positions(dday, (aday-211), 8) = r1I(2);
        object_positions(dday, (aday-211), 9) = r1I(3);
        object_positions(dday, (aday-211), 10) = v1I(1);
        object_positions(dday, (aday-211), 11) = v1I(2);
        object_positions(dday, (aday-211), 12) = v1I(3);
      endif
    end
  end
  %Save the object_positions variable to the file
  save(file_name, 'object_positions');
endif

if to_file4 == 1
  %load object_positions so we don't lose data
  load(file_name);

    % if we're saving the position data to the file, we need to compute it now.
  for dday = 1:1307    %Days from Jan 1 2017 to Dec 31 2017
    %Loop through each day of arrival
    fprintf('\n dday: %g', dday);
    for aday = 881:1856  %Days from Aug 1 2017 to Jan 31 2019 dated such that 1 is Jan 1 2017
      % Escape impossible arrival before departure
      if aday >= (dday + 60)
        %Determine location of Earth at departure date
        if dday == 1
          %Special case where initial date is correct and doesn't need to be computed
          rE = rEinit;
          vE = vEinit;
        else
          %Use U_twobody to determine Earth vectors at departure date
          tsec = ((dday - 1)*24*3600);    %Seconds between epoch and departure date
          [rE, vE] = U_twobody(rEinit, vEinit, tsec, Musun);
          %fprintf('\n rE: %g, vE: %g', rE, vE);
        endif

        %Determine location of Target at arrival date
        tsec = (aday*24*3600);    %Seconds between epoch and departure date
        %[r1I, v1I] = pararrayfun(15, @(r, v, t, Mu) U_twobody(r, v, t, Mu), r1Iinit, v1Iinit, tsec, Musun, "Vectorized", true, "ChunksPerProc", 1);
        [r2I, v2I] = U_twobody(r2Iinit, v2Iinit, tsec, Musun);
        %fprintf('\n r1I: %g, v1I: %g', r1I, v1I);

        %fprintf('\n dday: %g, aday: %g', dday, aday);
        %store data to appropriate location
        object_positions(dday, (aday-880), 1) = rE(1);
        object_positions(dday, (aday-880), 2) = rE(2);
        object_positions(dday, (aday-880), 3) = rE(3);
        object_positions(dday, (aday-880), 4) = vE(1);
        object_positions(dday, (aday-880), 5) = vE(2);
        object_positions(dday, (aday-880), 6) = vE(3);

        object_positions(dday, (aday-880), 13) = r2I(1);
        object_positions(dday, (aday-880), 14) = r2I(2);
        object_positions(dday, (aday-880), 15) = r2I(3);
        object_positions(dday, (aday-880), 16) = v2I(1);
        object_positions(dday, (aday-880), 17) = v2I(2);
        object_positions(dday, (aday-880), 18) = v2I(3);
      endif
    end
  end
  %Save the object_positions variable to the file
  save(file_name, 'object_positions');
endif


if from_file == 1
  load(file_name);
end


% Loop through each day of departure
if do_lambert3 == 1
  for dday = 1:364    %Days from Jan 1 2017 to Dec 31 2017
    %Loop through each day of arrival
    fprintf('\n dday: %g\n', dday);
    for aday = 212:760  %Days from Aug 1 2017 to Jan 31 2019 dated such that 1 is Jan 1 2017
      % Escape impossible arrival before departure
      if aday >= (dday + 60)
        %Feed all this to Lambert solver
        tofsec = ((aday - dday)*24*3600);    %Time of flight in seconds
        tofday = (aday - dday);    %Time of flight in days

        rE(1) = object_positions(dday, (aday - 211), 1);
        rE(2) = object_positions(dday, (aday - 211), 2);
        rE(3) = object_positions(dday, (aday - 211), 3);
        vE(1) = object_positions(dday, (aday - 211), 4);
        vE(2) = object_positions(dday, (aday - 211), 5);
        vE(3) = object_positions(dday, (aday - 211), 6);

        r1I(1) = object_positions(dday, (aday - 211), 7);
        r1I(2) = object_positions(dday, (aday - 211), 8);
        r1I(3) = object_positions(dday, (aday - 211), 9);
        v1I(1) = object_positions(dday, (aday - 211), 10);
        v1I(2) = object_positions(dday, (aday - 211), 11);
        v1I(3) = object_positions(dday, (aday - 211), 12);

        external_distances = [0,0];
        exitflag = 0;
        [lambV1, lambV2, extremal_distances, exitflag] = robust_lambert_solver(rE, r1I, tofday, 1, Musun);
        if exitflag == -1
          i = 2;
          while exitflag == -1 && i < 4
            [lambV1, lambV2, extremal_distances, exitflag] = robust_lambert_solver(rE, r1I, tofday, i, Musun);
            i = i + 1;
          end
          if i >= 4
            %Means unsolveable
          endif
        endif

        %[lambV1, lambV2] = lambert(rE, r1I, tofsec, 'retro');     %Prograde orbit
        %fprintf('\n dday: %g, aday: %g', dday, aday);

        %Subtract Planetary velocities from lambert velocities
        DVflyby = abs(norm(lambV1 - vE));
        DVrendez = abs(norm(lambV1 - vE)) + abs(norm(v1I - lambV2));

        %Put data into matrix
        q3data(dday, (aday-211), 1) = DVflyby;
        q3data(dday, (aday-211), 2) = DVrendez;
      endif
    end
  end
endif

if lambert_to_file3
  save(lambert_file_name, 'q3data');
endif

if do_lambert4 == 1
  for dday = 1:1307    %Days from Jan 1 2017 to July 31 2020
    %Loop through each day of arrival
    fprintf('\n dday: %g\n', dday);
    for aday = 881:1856  %Days from June 1 2019 to Jan 31 2022 dated such that 1 is Jan 1 2017             975
      % Escape impossible arrival before departure
      if aday >= (dday + 60)
        %Feed all this to Lambert solver
        tofsec = ((aday - dday)*24*3600);    %Time of flight in seconds
        tofday = (aday - dday);    %Time of flight in days

        rE(1) = object_positions(dday, (aday - 880), 1);
        rE(2) = object_positions(dday, (aday - 880), 2);
        rE(3) = object_positions(dday, (aday - 880), 3);
        vE(1) = object_positions(dday, (aday - 880), 4);
        vE(2) = object_positions(dday, (aday - 880), 5);
        vE(3) = object_positions(dday, (aday - 880), 6);

        r2I(1) = object_positions(dday, (aday - 880), 13);
        r2I(2) = object_positions(dday, (aday - 880), 14);
        r2I(3) = object_positions(dday, (aday - 880), 15);
        v2I(1) = object_positions(dday, (aday - 880), 16);
        v2I(2) = object_positions(dday, (aday - 880), 17);
        v2I(3) = object_positions(dday, (aday - 880), 18);

        [lamb2V1, lamb2V2, extremal_distances, exitflag] = robust_lambert_solver(rE, r2I, tofday, 1, Musun);
        if exitflag == -1
          i = 2;
          while exitflag == -1 && i < 4
            [lamb2V1, lamb2V2, extremal_distances, exitflag] = robust_lambert_solver(rE, r2I, tofday, i, Musun);
            i = i + 1;
          end
          if i >= 4
            %Means unsolveable

          endif
        endif
        %[lambV1, lambV2] = lambert(rE, r1I, tofsec, 'retro');     %Prograde orbit
        %fprintf('\n dday: %g, aday: %g', dday, aday);

        %Subtract Planetary velocities from lambert velocities
        DVflyby2 = abs(norm(lamb2V1 - vE));
        DVrendez2 = abs(norm(lamb2V1 - vE)) + abs(norm(v2I - lamb2V2));

        %Put data into matrix
        q4data(dday, (aday-880), 1) = DVflyby2;
        q4data(dday, (aday-880), 2) = DVrendez2;
      endif
    end
  end
endif

if lambert_to_file4
  save(lambert_file_name, 'q4data');
endif

if lambert_from_file
  if make_display_1
    load(lambert_file_name2);
  endif

  if make_display_2
    load(lambert_file_name);
  endif
endif

%Display plots
if make_display_1
  %Display first set
  %flyby

  %figure;
  %q3_image = uint8(255*cat(3,));
  subplot(1, 2, 1);
  %  surf(q3data(:,:,1), 'EdgeColor', 'None')
   % view(2)

    image_mx = q3data(:,:,1);
    imagesc(image_mx), axis equal tight, colorbar
    ylabel("Departure date, days since Jan 1 2017");
    xlabel("Arrival date, days since Aug 1 2017");
    title("DeltaV Required for Flyby from Earth to ’Oumouamoua");
    colormap(jet);

  subplot(1, 2, 2);
  %  surf(q3data(:,:,2), 'EdgeColor', 'None')
   % view(2)

    image_mx = q3data(:,:,2);
    imagesc(image_mx), axis equal tight, colorbar
    ylabel("Departure date, days since Jan 1 2017");
    xlabel("Arrival date, days since Aug 1 2017");
    title("DeltaV Required for Rendezvous from Earth to ’Oumouamoua");
    %colormap(jet);
endif

if make_display_2
  %Display second set
  subplot(1, 2, 1);
  %  surf(q3data(:,:,1), 'EdgeColor', 'None')
   % view(2)

    image_mx = q4data(:,:,1);
    imagesc(image_mx), axis equal tight, colorbar
    ylabel("Departure date, days since Jan 1 2017");
    xlabel("Arrival date, days since June 1 2019");
    title("DeltaV Required for Flyby from Earth to Borisov");
    colormap(jet);

  subplot(1, 2, 2);
  %  surf(q3data(:,:,2), 'EdgeColor', 'None')
   % view(2)

    image_mx = q4data(:,:,2);
    imagesc(image_mx), axis equal tight, colorbar
    ylabel("Departure date, days since Jan 1 2017");
    xlabel("Arrival date, days since June 1 2019");
    title("DeltaV Required for Rendezvous from Earth to Borisov");
    %colormap(jet);
endif















