function [recon_signals, R2] = predict_behavior_from_corr(t_win, Kfolds, tC, t_imaging, ...
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
    behavior_resampled_tr = cell(Na,1);
    for n=1:Na
        behavior_resampled_tr{n} = behavior_traces{n}(traininds);
    end
    
    %% train
    
    mRiemannianMean = RiemannianMean(tC(:,:, traininds));
    mX_tr = proj_R1(mRiemannianMean^(-1/2),tC(:,:, traininds));
    
    
    par.knn=max(round(size(mX_tr,2)*0.01),20);
    [ initAll, sig ] = CalcInitAff2D( mX_tr, par );
    dParams.maxInd = min(size(initAll,1),1+J);
    [~, Lambda, Psi_tr] = calcDiffusionMap(initAll,dParams);
    diffmap_tr=Psi_tr(:,2:end).';
    
    
    mdlglobal_phi=cell(Na, 1);
    
    
    for n=1:Na
        mdlglobal_phi{n} = fitlm(diffmap_tr(1:J,:)', behavior_resampled_tr{n});
    end
    %% test
    isnan_te = squeeze(isnan(tC(1,1,testinds)));
    mX_te = proj_R1(mRiemannianMean^(-1/2),tC(:, :, testinds(~isnan_te)));
    Psi_te = nystrom_extension(Psi_tr, mX_tr, mX_te, sig, Lambda)';
    diffmap_te=Psi_te(:,2:end).';
    for n=1:Na
        
        recon_signals{n}(testinds(~isnan_te)) = predict(mdlglobal_phi{n}, diffmap_te(1:J,:)');
        
    end
    
    
end
R2 = max(getstats(recon_signals, behavior_traces), 0);




