function R2 = getstats(reconY, behavior_traces)

for n=1:length(reconY)    
    notnans = ~isnan(reconY{n});    
    R2(n) = getRsq(reconY{n}(notnans), behavior_traces{n}(notnans));
end
end
function rsquare = getRsq(newy, y)

%y=nanzscore(y);
%newy=nanzscore(newy);
newy=newy(:);
y=y(:);
L=min(length(y), length(newy));
y=y(1:L);
newy=newy(1:L);
notnans = ~isnan(y) & ~isnan(newy);
SSE = mean((newy(notnans)-y(notnans)).^2);
TSS = var(y(notnans));
rsquare = 1 - SSE/TSS;
end