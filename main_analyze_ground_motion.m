% This script performs a basic analysis of an example ground motion, to
% demonstrate how the accompanying functions work. 

% Jack Baker
% October 19th, 2006

filename = 'example_ground_motion.acc'; % an acceleration time history

% the following block of code reads the specified file and integrates it to
% determine the velocity time history
file_in = fopen(filename, 'r');
for i=1:6
    ans = fgetl(file_in);
end
[ans ans] = fscanf(file_in, '%5c', 1);
[ans ans] = fscanf(file_in, '%5c', 1);
[record_dt ans] = fscanf(file_in, '%f', 1);
ans = fgetl(file_in);
[A, np] = fscanf(file_in, '%f');
fclose(file_in);
% integrate to get velocity time history in units of cm/s
signal = cumsum(A)' .* record_dt .* 981;


% perform the wavelet analysis
pulse_data = fn_analyze_record(signal, record_dt);

% output the classification results to the screen
pulse_indicator = pulse_data.pulse_indicator
late = pulse_data.late
PGV = pulse_data.PGV
is_pulse = pulse_data.is_pulse
Tp = pulse_data.Tp


% plot the velocity time history, extracted pulse, and residual
fn_plot_pulse(pulse_data, signal, record_dt);

