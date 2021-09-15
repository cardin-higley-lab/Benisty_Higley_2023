%% main script for extracting diffusion embedding of correlation traces and modeling behavior

% simulation parameters
winhopSec=0.1;
Kfolds = 10;
winLen = 3;
J = 20;

% load imaging and behavior data
load('inputdata.mat');
imaging_data1.time = imaging_data.time;
imaging_data1.data = imaging_data.data;
imaging_data1.fsample = imaging_data.fsample;
% extract instantaneous correlation matrices 
tC = sliding_window_corr_reg(imaging_data1, winLen, winhopSec);
[X, t_win] = sliding_window_mean(imaging_data1, winLen, winhopSec);


L = min([length(behavior_data.face_time), length(behavior_data.facemap)]);
t_win=round(t_win);

face_resampled = interp1(behavior_data.face_time(1:L), double(behavior_data.facemap(1:L)), imaging_data1.time(t_win));
pupil_resampled = interp1(behavior_data.face_time(1:L), double(behavior_data.pupil(1:L)), imaging_data1.time(t_win));
wheel_resampled = interp1(behavior_data.wheeltime, double(behavior_data.wheel_speed), imaging_data1.time(t_win));


behavior_traces{1} = pupil_resampled;
behavior_traces{2} = face_resampled;
behavior_traces{3} = wheel_resampled;



[recon_signals_corr, R2_corr] = predict_behavior_from_corr(t_win, Kfolds, tC, imaging_data1.time, ...
    behavior_traces, J);


[recon_signals_activity, R2_activity] = predict_behavior_from_activity(t_win, Kfolds, X, imaging_data1.time, ...
    behavior_traces, J);
    
ttls = {'pupil','facemap','wheel'};
figure;
for i = 1:3
subplot(3,1,i); plot(imaging_data1.time(t_win), behavior_traces{i});
hold all;
plot(imaging_data1.time(t_win), recon_signals_activity{i});
plot(imaging_data1.time(t_win), recon_signals_corr{i});
ylabel(ttls{i});axis tight;
end
xlabel('Time [sec]');
legend('Behavior','\Phi_a', '\Phi_c');




