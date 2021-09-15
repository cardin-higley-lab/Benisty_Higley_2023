function [tC, t_win] = sliding_window_corr_reg(imaging_data, winsizeSec, winhopSec)

Np = size(imaging_data.data,1);
T = floor(winsizeSec*imaging_data.fsample);


% add tau as the noise level
T = winsizeSec*imaging_data.fsample;
s=svd(imaging_data.data(:,1:round(T)));
th = quantile(s,0.3);
tau1 =median(s(s<=th));
[tC, t_win] = sliding_window_corr(imaging_data, winsizeSec, winhopSec, 2*tau1);
end

function [Wmat, t_win] = sliding_window_corr(imaging_data, winsizeSec, winhopSec, tau)

if ~exist('tau','var')
    tau=0;
end


fsampleimaging = imaging_data.fsample;
dFoF_parcells = imaging_data.data;
imaging_time = imaging_data.time;
winsize = round(fsampleimaging*winsizeSec);
winhop = round(fsampleimaging*winhopSec);


winst = 1:winhop:length(imaging_time);
winen = winst+winsize-1;
winst=winst(winen<size(dFoF_parcells,2));
winen=winen(winen<size(dFoF_parcells,2));
t_win = round((winst+winen)/2);

Np = size(dFoF_parcells,1);
Wmat = zeros(Np, Np, length(winst));
for i=1:length(winst)
    W = corr(dFoF_parcells(:, winst(i):winen(i))')+eye(Np)*tau;
    Wmat(:,:,i) = (W+W')/2;
end
end