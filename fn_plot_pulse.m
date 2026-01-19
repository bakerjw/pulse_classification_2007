function [] = fn_plot_pulse(pulse_data, signal, record_dt)

% plot time history, extracted pulse, and residual
% Jack Baker
% August 22, 2006

% create the extracted pulse
pulse_th = fn_build_pulse(pulse_data);

% create the residual time history 
resid_th = signal - pulse_th;

% plot the record, pulse and residual
np = length(signal);
time = record_dt:record_dt:record_dt*np;
figure
subplot(3,1,1)
plot(time, signal, '-k')
set(gca, 'xticklabel', [])
yl= [-max(abs(get(gca, 'ylim'))) max(abs(get(gca, 'ylim')))];
set(gca, 'ylim', yl)

subplot(3,1,2)
plot(time, pulse_th, '-b')
hy = ylabel('Velocity [cm/s]');
set(gca, 'xticklabel', [])
set(gca, 'ylim', yl)

subplot(3,1,3)
plot(time, resid_th , '-k')
hx = xlabel('Time [s]');
set(gca, 'ylim', yl)