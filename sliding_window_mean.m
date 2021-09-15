function [X, t_win] = sliding_window_mean(imaging_data, winsizeSec, winhopSec)





fsampleimaging = imaging_data.fsample;
dFoF_parcells = imaging_data.data;
imaging_time = imaging_data.time;
winsize = round(fsampleimaging*winsizeSec);
winhop = round(fsampleimaging*winhopSec);


winst = 1:winhop:length(imaging_time);
winen = winst+winsize;
winst=winst(winen<size(dFoF_parcells,2));
winen=winen(winen<size(dFoF_parcells,2));
t_win = (winst+winen)/2;
X = zeros(size(dFoF_parcells,1),  length(winst));
for i=1:length(winst)
    %disp(i/length(winst))
    W = corr(dFoF_parcells(:, winst(i):winen(i))');
    X(:,i) = nanmean(dFoF_parcells(:, winst(i):winen(i)),2);
end
