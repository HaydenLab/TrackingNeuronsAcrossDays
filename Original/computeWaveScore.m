function score = computeWaveScore(wmean)
score = cell(numel(wmean)-1,1);

% Interpolate up x10
for iid=1:length(wmean)
    for unit1=1:length(wmean{iid})
        maxlength = length(wmean{iid}{unit1}(:));
        wmean{iid}{unit1} = interp1(1:maxlength,wmean{iid}{unit1}(:),1:.1:maxlength); % interp1(1:32,wmean{iid}{unit1}(:),1:.1:32);
    end
end

for iid=1:length(wmean)-1
    score{iid} = nan(length(wmean{iid}),length(wmean{iid+1}));
    for unit1=1:length(wmean{iid})
        for unit2=1:length(wmean{iid+1})
            w1 = wmean{iid}{unit1};
            w2 = wmean{iid+1}{unit2};
            w1 = w1 / norm(w1);
            w2 = w2 / norm(w2);
            score{iid}(unit1,unit2) = atanh(max(conv(w1,w2(end:-1:1),'same')));
            %score{iid}(unit1,unit2) = atanh(max(xcorr(wmean{iid}{unit1},wmean{iid+1}{unit2},ceil(numel(wmean{iid}{unit1})/4),'coeff')));
        end
    end
end