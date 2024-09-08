clear;clc;

% PARAMETERS THAT ARE THE SAME FOR ALL CASES
V1 = 1.00 + 0.0i; % Voltage of slack bus bar (pu)
base_power = 100; % Base power (MVA) delivered to loads

% Load power at busbars
SD2 = -0.75 - 0.50i;
SD3 = -0.50 - 0.25i;
SD4 = -0.25 - 0.75i;
SD5 = -0.25 - 0.25i;
SD6 = -0.75 - 0.25i;

% Loads (in pu) connected to busbars 2–6
Z = zeros(6,6);
Z(1,2) = 0.010 + 0.02i;
Z(1,3) = 0.020 + 0.03i;
Z(2,3) = 0.010 + 0.04i;
Z(3,4) = 0.020 + 0.03i;
Z(4,5) = 0.010 + 0.03i;
Z(1,6) = 0.020 + 0.04i;
Z(5,6) = 0.020 + 0.03i;

for i = 6:-1:1
    for j = 6:-1:1
        if i == j
           Z(i,j) = 0;
        else
           Z(i,j) = Z(j,i);
        end
    end
end
Z;

% Bus admittance matrix (Y Bus)
Y = zeros(6,6);
for i = 1:6
    for j = 1:6
        if i == j
           for k = 1:6
               if Z(i,k) == 0
                   % neglect this one
               else
                   Y(i,j) = Y(i,j) + 1/Z(i,k);
               end
           end
        else
           if Z(i,j) == 0
                % neglect this one
           else
               Y(i,j) = -1/Z(i,j);
           end
        end
    end
end
Y

% Create a matrix to summarise all the cases (used later for line loss table)
cases = zeros(8,6);

% Summarise locations of the new generator and capacitor bank for each case
% (used later for generated power at busbars)
SG = [0, 0, 0, 0, 0;
     (0.5 + 0.5i), 0, 0, 0, 0;
     0, (0.5 + 0.5i), 0, 0, 0;
     0, 0, (0.5 + 0.5i), 0, 0;
     0, 0, 0, (0.5 + 0.5i), 0;
     0, 0, 0, 0, (0.5 + 0.5i)];

for c = 1:6
    % Generated power at busbars (using summary created earlier)
    SG2 = SG(c,1);
    SG3 = SG(c,2);
    SG4 = SG(c,3);
    SG5 = SG(c,4);
    SG6 = SG(c,5);

    % Total power at busbars
    S2 = SG2 + SD2;
    S3 = SG3 + SD3;
    S4 = SG4 + SD4;
    S5 = SG5 + SD5;
    S6 = SG6 + SD6;

    % Voltage at each busbar (using Gauss-Seidle iteration method)
    % start with non-zero values
    V2(1) = 0.5 + 0i;
    V3(1) = 0.5 + 0i;
    V4(1) = 0.5 + 0i;
    V5(1) = 0.5 + 0i;
    V6(1) = 0.5 + 0i;

    for j=1:1:99
        V2(j+1) = (1/Y(2,2))*(conj(S2/V2(j)) - Y(2,1)*V1 - Y(2,3)*V3(j) - Y(2,4)*V4(j) - Y(2,5)*V5(j) - Y(2,6)*V6(j));
        V3(j+1) = (1/Y(3,3))*(conj(S3/V3(j)) - Y(3,1)*V1 - Y(3,2)*V2(j) - Y(3,4)*V4(j) - Y(3,5)*V5(j) - Y(3,6)*V6(j));
        V4(j+1) = (1/Y(4,4))*(conj(S4/V4(j)) - Y(4,1)*V1 - Y(4,2)*V2(j) - Y(4,3)*V3(j) - Y(4,5)*V5(j) - Y(4,6)*V6(j));
        V5(j+1) = (1/Y(5,5))*(conj(S5/V5(j)) - Y(5,1)*V1 - Y(5,2)*V2(j) - Y(5,3)*V3(j) - Y(5,4)*V4(j) - Y(5,6)*V6(j));
        V6(j+1) = (1/Y(6,6))*(conj(S6/V6(j)) - Y(6,1)*V1 - Y(6,2)*V2(j) - Y(6,3)*V3(j) - Y(6,4)*V4(j) - Y(6,5)*V5(j));
    end

    % Use these final values
    V2(100);
    V3(100);
    V4(100);
    V5(100);
    V6(100);

    % Power generated at busbar 1 (S1 = V1 x I1*)
    % (use the final voltage values Vj(100))
    S1 = V1*(Y(1,1)*V1 + Y(1,2)*V2(100) + Y(1,3)*V3(100) + Y(1,4)*V4(100) + Y(1,5)*V5(100) + Y(1,6)*V6(100))

    % Current flow (pu)
    % Note: add negative sign to admittance as it is negative in the Y BUS
    I12 = -Y(1,2)*(V1 - V2(100));
    I21 = -I12;
    I13 = -Y(1,3)*(V1 - V3(100));
    I31 = -I13;
    I16 = -Y(1,6)*(V1 - V6(100));
    I61 = -I16;
    I23 = -Y(2,3)*(V2(100) - V3(100));
    I32 = -I23;
    I34 = -Y(3,4)*(V3(100) - V4(100));
    I43 = -I34;
    I45 = -Y(4,5)*(V4(100) - V5(100));
    I54 = -I45;
    I56 = -Y(5,6)*(V5(100) - V6(100));
    I65 = -I56;

    % Power flow (pu)
    S12 = V1*conj(I12);
    S21 = V2(100)*conj(I21);
    S13 = V1*conj(I13);
    S31 = V3(100)*conj(I31);
    S16 = V1*conj(I16);
    S61 = V6(100)*conj(I61);
    S23 = V2(100)*conj(I23);
    S32 = V3(100)*conj(I32);
    S34 = V3(100)*conj(I34);
    S43 = V4(100)*conj(I43);
    S45 = V4(100)*conj(I45);
    S54 = V5(100)*conj(I54);
    S56 = V5(100)*conj(I56);
    S65 = V6(100)*conj(I65);

    % Line losses (pu)
    SL12_pu = S12 + S21;
    SL13_pu = S13 + S31;
    SL16_pu = S16 + S61;
    SL23_pu = S23 + S32;
    SL34_pu = S34 + S43;
    SL45_pu = S45 + S54;
    SL56_pu = S56 + S65;

    % Line losses (MW & MVAR)
    SL12 = SL12_pu*base_power;
    SL13 = SL13_pu*base_power;
    SL16 = SL16_pu*base_power;
    SL23 = SL23_pu*base_power;
    SL34 = SL34_pu*base_power;
    SL45 = SL45_pu*base_power;
    SL56 = SL56_pu*base_power;

    % Total Line loss (pu)
    SL_total_pu = SL12_pu + SL13_pu + SL16_pu + SL23_pu + SL34_pu + SL45_pu + SL56_pu;

    % Total Line loss (MW & MVAR)
    SL_total = SL12 + SL13 + SL16 + SL23 + SL34 + SL45 + SL56;

    % Summarise the line losses for this case
    cases(1,c) = SL12;
    cases(2,c) = SL13;
    cases(3,c) = SL16;
    cases(4,c) = SL23;
    cases(5,c) = SL34;
    cases(6,c) = SL45;
    cases(7,c) = SL56;
    cases(8,c) = SL_total;
end

% TABULATE line losses for each case
T = array2table(cases,...
    'VariableNames',{'non-existent' 'at busbar 2' 'at busbar 3' 'at busbar 4' 'at busbar 5' 'at busbar 6'}, ...
    'RowNames',{'SL12';'SL13';'SL16';'SL23';'SL34';'SL45';'SL56';'SL_total'});
disp(T)

