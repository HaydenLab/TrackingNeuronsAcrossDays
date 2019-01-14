function base = computeDailyBaserates(spiketimes)
base = cell(size(spiketimes));

for day=1:length(spiketimes)
    base{day} = nan(1,length(spiketimes{day}));
    for iic=1:length(spiketimes{day})
        base{day}(iic) = log(length(spiketimes{day}{iic})./range(spiketimes{day}{iic}));
    end
end