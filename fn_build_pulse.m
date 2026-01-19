function [pulse_th] = fn_build_pulse(pulse_data)

% reassemble the extracted velocity pulse, based on the wavelet
% coefficients and locations provided in the data structure created by the
% function fn_analyze_record

% Jack Baker
% August 22, 2006

% INPUT
% pulse_data    Structure containing relevant information about the
%               extracted pulse. See the function fn_analyze_record

% OUTPUT
% pulse_th      Velocity time history of the extracted pulse

% get relevant data from the pulse_data structure
rec_dt = pulse_data.dt;
scale = pulse_data.pulse_scale;
n = length(pulse_data.rows); % reduce this if coarser pulse recreation is desired
max_t = pulse_data.num_pts*rec_dt;
time = rec_dt:rec_dt:max_t;


% reconstruct pulse
wname = pulse_data.wavelet_name;
iter = 4; % how detailed is the wavelet deconstruction
wtype = wavemngr('type',wname);
switch wtype
    case 1 , [nul,psi,xval] = wavefun(wname,iter);
    case 2 , [nul,psi,nul,nul,xval] = wavefun(wname,iter);
    case 3 , [nul,psi,xval] = wavefun(wname,iter);
    case 4 , [psi,xval] = wavefun(wname,iter);
    case 5 , [psi,xval] = wavefun(wname,iter);
end
pulse_th = zeros(size(time));
for i=1:n
    basis = xval.*(scale);
    basis = basis + (pulse_data.rows(i) - median(basis));
    basis = basis.*rec_dt;
    y_vals = psi * pulse_data.coefs(i) / sqrt(scale);
    delta = basis(2) - basis(1);

    num_pads = ceil((max(time) - max(basis))/delta);
    final_basis = [0:delta:min(basis-0.00001) basis (max(basis)+delta:delta:max(basis)+delta*num_pads)];
    final_yvals = [zeros(size(0:delta:min(basis-0.00001))) y_vals zeros(1,num_pads)];

    pulse_th = pulse_th + interp1(final_basis,final_yvals,time);
end


% Occasionally the interpolation function above 
% thinks it is extrapolating due to numerical round-off 
% error, and inserts a 'nan' value at the end of the signal. 
% Replace the nan with a zero. 
pulse_th(find(isnan(pulse_th))) = 0; 

