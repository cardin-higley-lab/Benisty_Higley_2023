function [recon_signals, R2] = predict_behavior_from_activity(t_win, Kfolds, X, t_imaging, ...
    behavior_traces, J)

Len = t_imaging(t_win(end))-t_imaging(t_win(1));
chunksize = floor(Len/Kfolds);
Na = length(behavior_traces);

recon_signals = cell(Na,1);
for n=1:Na
    recon_signals{n} = nan(size(behavior_traces{n}));
end

for chunk_i = 1:Kfolds
    disp(chunk_i);
    teststart = t_imaging(1) + (chunk_i - 1)*chunksize;
    testend = t_imaging(1)+chunk_i*chunksize;
    
    testinds = findClosestDouble(teststart,t_imaging(t_win)):findClosestDouble(testend,t_imaging(t_win));
    traininds = setdiff(1:length(t_win), testinds);
    
    arousal_resampled_tr = cell(Na,1);
    for n=1:Na
        arousal_resampled_tr{n} = behavior_traces{n}(traininds);
    end
    
    %% train
    
    
    Xtr = X(:, traininds);
    par.knn=5000;
    [initAll, sig] = CalcInitAff2D( Xtr, par );
    dParams.maxInd = min(size(initAll,1),J+1);
    [~, Lambda, Psi_tr] = calcDiffusionMap(initAll,dParams);
    
    diffmap_tr=Psi_tr(:,2:end).';
    
    
    
    
    mdlglobal_phi=cell(1,Na);
    
    for n=1:Na
        mdlglobal_phi{n} = fitlm(diffmap_tr(1:J,:)', arousal_resampled_tr{n});
    end
    
    %% test
    
    isnan_te = squeeze(isnan(X(1,testinds)));
    Xte=X(:,  testinds(~isnan_te));
    Psi_te = nystrom_extension(Psi_tr, Xtr, Xte, sig, Lambda)';
    
    
    diffmap_te=Psi_te(:,2:end).';
    
    for n=1:Na
        
        recon_signals{n}(testinds(~isnan_te)) = predict(mdlglobal_phi{n}, diffmap_te(1:J,:)');
        
    end
    
    
    
    
end



R2 = max(getstats(recon_signals, behavior_traces), 0);






end
