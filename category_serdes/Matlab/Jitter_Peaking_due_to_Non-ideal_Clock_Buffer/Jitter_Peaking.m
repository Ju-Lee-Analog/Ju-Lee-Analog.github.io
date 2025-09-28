clc; clear; close all;

color_order = [
    0.0, 0.0, 0.0;   % Black
    0.9, 0.6, 0.0;   % Orange
    0.6, 0.1, 0.2;   % Red
    0.0, 0.6, 0.5;   % Green
    0.8, 0.4, 0.0;   % Brown
    0.95, 0.9, 0.25; % Yellow
    0.0, 0.45, 0.7;  % Blue
];
set(groot, 'defaultAxesColorOrder', color_order, ...
    'defaultFigureColor', 'white', ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultTextFontName', 'Times New Roman', ...
    'defaultLegendFontName', 'Times New Roman', ...
    'defaultColorbarFontName', 'Times New Roman');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_cycles = 10;                                       % Run n_cycles Square Wave                                                                                         
n_points = 100*n_cycles;                             % n_points in 1 cycle

A0 = 1;                                              % Clock Amplitude
f1 = 14e9;                                           % Clock Frequency
t = linspace(0, 1/f1*n_cycles, n_points);            % Time Scale = n_cycles cycles

CK = A0 * square (2*pi*f1*t);                        % Generate Square Wave Clock

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inject_cycle = 3;                                    % Jitter Injection at Timing = #inject_cycle Period
jitter_time = inject_cycle/f1;                       % Jitter Injection 'Time'
jitter_amp = 0.1/f1;                                 % Jitter Amplitude = jitter_amp (sec)
[~, idx] = min(abs(t-jitter_time));                  % idx = points at injection
[~, idx2] = min(abs(t-jitter_time-jitter_amp));      % idx2 = points after injection
t_skew = t;
t_skew(idx:end) = t_skew(idx:end) - jitter_amp;      % Shift CK from idx to idx2

CK_SKEW = square(2*pi*f1*t_skew);                   
CK_SKEW(idx2:end) = CK(idx2:end);                    % Generate Jittered CK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flpf = 3e9;                                         % LPF Corner = flpf (Hz)
omega_lpf = 2*pi*flpf;                              % LPF Corner = omega_lpf (rad/s)
LPF = tf(omega_lpf, [1, omega_lpf]);                % LPF Transfer Function

[CK_LPF, ~] = lsim(LPF, CK, t);                     % CK_LPF = LPF(CK)
[CK_SKEW_LPF, ~] = lsim(LPF, CK_SKEW, t);           % CK_SKEW_LPF = LPF(CK_SKEW)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CK_EDGE = find(diff(sign(CK_LPF))~=0);                              % Find Edge (Rising/Falling) of CK_LPF --> Discrete Signal, Length = # of Edges
CK_SKEW_EDGE = find(diff(sign(CK_SKEW_LPF))~=0);                    % Find Edge (Rising/Falling) of CK_SKEW_LPF --> Discrete Signal
PHASE_DIFF = (CK_SKEW_EDGE-CK_EDGE).*n_cycles./n_points.*2*pi;      % PHASE_DIFF = Phase Shifted (rad) --> Discrete Signal
EDGE = t(CK_EDGE);                                                  % EDGE = Ideal Time of CK EDGE --> Discrete Signal

PHASE_DIFF_CT = zeros(n_points, 1);                                 % Fill zeros to PHASE_DIFF --> CT Signal, Length = n_cycles * n_points

A1 = 0;             % A1 = First Peak Value of PHASE_DIFF, A1 > 0, unit: rad
A2 = 0;             % A2 = Second Peak Value of PHASE_DIFF, A2 > 0, unit: rad
T1 = 0;             % T1 = Time of A1
T2 = 0;             % T2 = Time of A2
IDX1 = 0;           % TDX1 = Discrete Time of A1
IDX2 = 0;           % TDX2 = Discrete Time of A2
n_A1 = 0;           % n_A1 = Continuse "Point" of A1, n_A1 ~= T1 / (1/f1 * n_cycles/n_points)
n_A2 = 0;           % n_A2 = Continuse "Point" of A2, n_A2 ~= T2 / (1/f1 * n_cycles/n_points)

cnt = 0;
j = 1;
for i = 1:n_points
    if (j <= n_cycles*2)
        if (t(i) == EDGE(j))
            PHASE_DIFF_CT(i) = PHASE_DIFF(j);
            if (cnt == 0 && PHASE_DIFF(j)>0)
                A1 = PHASE_DIFF(j);
                T1 = EDGE(j);
                n_A1 = i;
                IDX1 = j;
                cnt = cnt + 1;
            elseif (cnt == 1 && PHASE_DIFF(j)>0)
                A2 = PHASE_DIFF(j);
                T2 = EDGE(j);
                n_A2 = i;
                IDX2 = j;
                cnt = cnt + 1;
            end
            j = j + 1;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = linspace(1, n_cycles*2, n_cycles*2);

input_idx = inject_cycle*n_points/n_cycles;                       % input_idx = Continuse "Point" of Jitter Injection
output_idx = n_A1;                                                % output_idx = n_A1 = Continuse "Point" of A1, n_A1 ~= T1 / (1/f1 * n_cycles/n_points)
in = zeros(n_points, 1);                                          % Generate Continuse Time Input Impulse Signal, in
in(input_idx) = jitter_amp*2*pi*f1;                               % in only has value at input_idx, value = jitter_amp * 2 * pi * f1 (rad)

delta = log(A1/A2);                                               %%% Suppose Jitter Transfer is 2nd-order System %%%
zeta = delta/sqrt(4*pi^2+delta^2);
wn = 2*pi/(T2-T1);

x = wn/sqrt(1-zeta^2)*exp(-1*zeta*wn*t).*sin(wn*sqrt(1-zeta^2)*t);      % x = Impulse Response of a 2nd-order System
x = x./max(x).*A1;                                                      % Normalized x to the amplitude of PHASE_DIFF_CT
x_fit = [zeros(1, n_A1-(n_A2-n_A1)/4), x];                              % x_fit: Shift timing offset of x to fit PHASE_DIFF_CT
x_fit = x_fit(1:length(x));

data = iddata(x_fit', in, 1/f1*n_cycles/n_points);                      % data: input = in, output = x_fit 
JP = tfest(data, 2, 1);                                                 % Find a Transfer Function JP which fits data
                                                                        % Using 1 zero, 2 pole System

[mag, phase, w] = bode(JP);                                             % Export JP Bode Plot Values
mag = squeeze(mag);
phase = squeeze(phase);
w = squeeze(w);
f_Hz = w / (2 * pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Ideal CK and Jittered CK
figure(1);
plot(t, CK, 'LineWidth', 1.5);
hold on;
plot(t, CK_SKEW, 'LineWidth', 1.5);
legend('CK_{IN}', 'CK_{SKEW}');
title('Ideal Input Clock and Jittered Clock')
xlabel('Time (sec)');
ylabel('Voltage (V)');
grid on;

% Plot Ideal CK and Jittered CK after LPF
figure(2);
plot(t, CK_LPF, 'LineWidth', 1.5);
hold on;
plot(t, CK_SKEW_LPF, 'LineWidth', 1.5);
legend('CK_{LPF}', 'CK_{SKEW LPF}');
title('Output Clock and Jittered Clock (Non-Clamped)')
xlabel('Time (sec)');
ylabel('Voltage (V)');
grid on;

% Plot Discrete Time Phase Error Due to Limited Bandwidth
figure(3);
stem(EDGE, PHASE_DIFF, 'LineWidth', 1.5);
legend('Phase Noise');
title('Jitter-induced Phase Noise')
xlabel('Time (sec)');
ylabel('Phase Error (rad)');
grid on;

% Plot CT Jitter Injection, CT Fitting Phase Error and Real DT Phase Error
figure(4)
plot(t, x_fit, 'LineWidth', 1.5);
hold on;
stem(EDGE, PHASE_DIFF, 'LineWidth', 1.5);
hold on;
stem(t, in, 'LineWidth', 1.5);
legend('Phase Error Fit', 'Phase Error', 'Jitter Injection');
title("Jitter Transfer Injection Response");
xlabel("Time (sec)");
ylabel("Phase Error (rad)");
grid on;

% Plot Fitting Correlation between Transfer Function and Fitting Phase Error
figure(5)
compare(data, JP);
title("Jitter Transfer Function Fitting Results");
xlabel("Time (sec)");
ylabel("Phase Error (rad)");

% Plot Impulse Response of the Transfer Function
figure(6)
impulse(JP);
title("Jitter Transfer Impulse Response (Fitting)");
xlabel("Time (sec)");
ylabel("Phase Error (rad)");
grid on;

% Plot Bode of Transfer Function
figure(7);
h = bodeplot(JP);
setoptions(h, 'FreqUnits', 'Hz');
title("Jitter Transfer Bode Plot");
grid on;


