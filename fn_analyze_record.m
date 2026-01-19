function [pulse_data, Tp, signals] = fn_analyze_record(signal, dt)

% by Jack Baker
% January 12, 2006
% revised August 22, 2006

% INPUT
% signal        Vector of the ground motion's velocity time history, in
%               units of cm/s^2
% dt            Time step of the velocity time history

% OUTPUT
% pulse_data    Structure containing relevant information about the
%               extracted pulse. See the bottom of this function for
%               information about what is included
% Tp            Period of the extracted pulse
% signals       Structure containing the pulse time history and residual.
%               This is currently commented out to speed up the function.


% defined parameters
num_coefs = 10; % number of coefficients to extract
wname = 'db4'; % use Daubechies Wavelet of fourth order for the transform
iter = 4; % level of detail of the wavelet deconstruction

% first analysis, to find Tp and location
[coefs, pulse_scale, rows, Tp] = fn_extract_one_wavelet(signal, dt);

% extract residual
    time = dt:dt:length(signal)*dt;
    wtype = wavemngr('type',wname);
    switch wtype
      case 1 , [nul,psi,xval] = wavefun(wname,iter);
      case 2 , [nul,psi,nul,nul,xval] = wavefun(wname,iter);
      case 3 , [nul,psi,xval] = wavefun(wname,iter);
      case 4 , [psi,xval] = wavefun(wname,iter);
      case 5 , [psi,xval] = wavefun(wname,iter);
    end
    pulse_th = zeros(size(signal));
    % compute pulse signal and residual
    basis = xval.*(pulse_scale);
    basis = basis + (rows(1) - median(basis));
    basis = basis.*dt;
    y_vals = psi * coefs(1) / sqrt(pulse_scale);
    delta = basis(2) - basis(1);
    num_pads = ceil((max(time) - max(basis))/delta);
    final_basis = [0:delta:min(basis-0.00001) basis (max(basis)+delta:delta:max(basis)+delta*num_pads)];
    final_yvals = [zeros(size(0:delta:min(basis-0.00001))) y_vals zeros(1,num_pads)];
    pulse_th = pulse_th + interp1(final_basis,final_yvals,time);
    pulse_th(find(isnan(pulse_th))) = 0; %remove nan values due to round-off error leading interp to think it is extrapolating

    first_pulse_th = pulse_th;
    resid_th = signal - pulse_th;
    first_resid_th = resid_th;

% get additional coefficients in the neighborhood of the pulse
for i = 2:num_coefs
    [coefs(1,i), ans, rows(1,i), ans] = fn_extract_one_wavelet(resid_th, dt, pulse_scale, rows(1));
    
    % compute the residual
    basis = xval.*(pulse_scale);
    basis = basis + (rows(1,i) - median(basis));
    basis = basis.*dt;
    y_vals = psi * coefs(1,i) / sqrt(pulse_scale);
    num_pads = ceil((max(time) - max(basis))/delta);
    final_basis = [0:delta:min(basis-0.00001) basis (max(basis)+delta:delta:max(basis)+delta*num_pads)];
    final_yvals = [zeros(size(0:delta:min(basis-0.00001))) y_vals zeros(1,num_pads)];
    pulse_th = pulse_th + interp1(final_basis,final_yvals,time);
    pulse_th(find(isnan(pulse_th))) = 0; %remove nan values due to round-off error leading interp to think it is extrapolating
    resid_th = signal - pulse_th;
end


% criterion to eliminate late-arriving pulses
    % energy fraction threshold
    orig_threshold = 0.2;
    pulse_threshold = 0.1;

    % times at which the critical energy fraction thresholds are reached
    cum_energy_orig = cumsum(signal.^2) ./ sum(signal.^2);
    orig_t = time(min(find(cum_energy_orig>orig_threshold)));
    cum_energy_pulse = cumsum(pulse_th.^2) ./ sum(pulse_th.^2);
    pulse_t = time(min(find(cum_energy_pulse>pulse_threshold)));

    delta = (orig_t - pulse_t)/(dt*length(signal)); %fraction of record duration between orig_t and pulse_t
    late = orig_t < pulse_t;



% sum of the coefficients of the DWT of the original signal and the final
% residual. This sum is also equal to the energy of the ground motion.
    pulse_data.dwt_orig = sum(sum(abs(dwt(signal,wname))));
    pulse_data.dwt_resid = sum(sum(abs(dwt(resid_th,wname))));
    pulse_data.dwt_squared_orig = sum(sum( dwt(signal,wname).^2  ));
    pulse_data.dwt_squared_resid = sum(sum( dwt(resid_th,wname).^2  ));

% compute the "Pulse Indicator" based on the PGV ratio and energy ratio
    PGV = max(abs(signal));
    PGV_resid = max(abs(resid_th));
    pulse_indicator = (1./(1+exp(-(23.297  -14.579*(PGV_resid/PGV) -20.477*(pulse_data.dwt_squared_resid/pulse_data.dwt_squared_orig) ))));


% identify whether the record contains a significant pulse
    is_pulse = (pulse_indicator>0.85) & (PGV>30) & ~late;


% build a structure to hold the decomposition data
    pulse_data.dt = dt; % time step of the ground motion recording
    pulse_data.num_pts = length(signal); % number of points in the recording
    pulse_data.Tp = Tp; % pulse period
    pulse_data.wavelet_name = wname; % name of the wavelet used for analysis
    pulse_data.pulse_scale = pulse_scale; % scale at which the largest wavelet was found
    pulse_data.rows = rows; % rows where wavelets were extracted
    pulse_data.coefs = coefs; % coefficients for the extracted wavelets
    pulse_data.PGV = PGV; % Peak ground velocity of the original ground motion
    pulse_data.PGV_resid = PGV_resid; % Peak ground velocity of the residual ground motion
    pulse_data.delta_energy_t= delta; % if this is negative, then the pulse is late-arriving
    pulse_data.late = late; % is the pulse late arriving
    pulse_data.pulse_indicator = pulse_indicator; % the score from the logistic regression used for classification
    pulse_data.is_pulse = is_pulse; % classification of the record

% % build a structure to hold the reconstructed signals
% % this step is somewhat slow, so leave it commented out unless it is needed
% signals.first_pulse_th = first_pulse_th;
% signals.pulse_th = pulse_th;
% signals.first_resid_th = first_resid_th;
% signals.resid_th = resid_th;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract one additional wavelet
function [coef, pulse_scale, row, Tp] = fn_extract_one_wavelet(signal, dt, pulse_scale, pulse_row)

    % defined parameters
    Tp_min = 0.25; % minimum period to consider for pulses
    Tp_max = 15; % maximum period to consider for pulses
    row_range = 10/25; % range of row locations in which to look for supplementary wavelets
    wname = 'db4'; % use Daubechies Wavelet of fourth order for the transform

    % set the scales at which to perform the CWT
    if nargin < 3
        % set an initial coarse range of scales to inspect
        num_scales = 50; % optimal would be sqrt(2* (scale_max-scale_min)
        scale_min = floor(Tp_min/1.4/dt); % the factor of 1.4 was back-calculated from the scal2frq function below
        scale_step = ceil((Tp_max/1.4/dt - scale_min)/num_scales);
        scale_max = scale_min + num_scales*scale_step;
        scales = scale_min:scale_step:scale_max;
    else
        % scale and range of rows to consider are already known
        scales = pulse_scale;
        row_indices = [max(1, pulse_row - ceil(scales*row_range)) : min(length(signal), pulse_row + ceil(scales*row_range))];
    end

    % peform the Continuous Wavelet Transform
    cwt_coefs = cwt(signal,scales, wname);

    % find max coefficient
    if nargin < 3 % critial scale has not yet been identified
        % find the neighborhood of the maximum coefficient based on the coarse scales
        z=max(max(abs(cwt_coefs)));
        col = find(z==max(abs(cwt_coefs)'));
        row = find(z==(abs(cwt_coefs(col,:))));

        % revise the scales to a finer grid and repeat CWT to find the exact period of the maximum coeficient
        scales = [scales(max(1, col-1)): scales(min(num_scales+1, col+1))];
        cwt_coefs = cwt(signal,scales, wname);
        z=max(max(abs(cwt_coefs)));
        col = find(z==max(abs(cwt_coefs)'));
        row = find(z==(abs(cwt_coefs(col,:))));
        pulse_scale = scales(col);
    else
        z=max(max(abs(cwt_coefs(1,row_indices))));
        col = 1;
        row = find(z==(abs(cwt_coefs(1,:))));
    end

    % the value of the coeficient associated with the extracted pulse
    coef = cwt_coefs(col, row);

    % pulse period (Tp)
    Tp = 1./scal2frq(pulse_scale,wname,dt);







