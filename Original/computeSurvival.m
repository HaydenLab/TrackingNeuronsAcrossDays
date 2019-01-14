function [survival, score] = computeSurvival(channel, unit, corrscore, wavescore, autoscore, basescore, survival)

sameChan = cell(size(corrscore));
sameUnit = cell(size(corrscore));
for iid=1:length(corrscore)
    sameChan{iid} = bsxfun(@eq,channel{iid}(:),channel{iid+1}(:)');
    sameUnit{iid} = sameChan{iid} & bsxfun(@eq,unit{iid}(:),unit{iid+1}(:)');
end

data = [unroll(corrscore) unroll(wavescore) unroll(autoscore) unroll(basescore)];
good = ~isnan(sum(data));
fprintf('Using %d classifiers, %d others are nan\n',sum(good),sum(~good));
% Remove missing classifiers
data = data(:,good);

% Approximates expectation-maximization with partly specified labels
% The E-step is making a hard assignment which is not exactly right but in
% the context of this type of data where most labels are specified it makes
% no difference
unrollSameChan = unroll(sameChan); % removed only same channel
C = unroll(survival);
for i=1:10
    [C,err,P] = classify(data, data, C,'quadratic');
    % Only same-channel is possible
    C = C & unrollSameChan;
end

% Identify a threshold for a 1% FP rate
negative = data(~unroll(sameChan),:);
[Cneg,err,Pneg] = classify(negative,data,C,'quadratic');
negative = Pneg(:,2);
threshold = quantile(negative,.95);
C = P(:,2)>threshold;
C = C & unrollSameChan;

survival = reroll(corrscore, C);
score = reroll(corrscore, P(:,2));

% Occasionally survival will indicate something impossible, like the same 
% unit becomes two different units or vice versa.  We need to fix all those 
% cases.
for iid=1:length(survival)
    for iic=channel{iid}(:)'
        left = channel{iid}==iic;
        right = channel{iid+1}==iic;
        survival{iid}(left,right) = takeBest(score{iid}(left,right), threshold);
    end
end
end

function survival = takeBest(similarity, thresh)
% Posterior must be > thresh
[pre,post] = find(similarity > thresh);
pre = unique(pre); post = unique(post);

similarity(similarity < thresh) = nan;
sim = similarity(pre,post);

survival = zeros(size(similarity));
survival(pre,post) = eye(numel(pre),numel(post));
survival = survival > 0;
if(numel(pre)>numel(post))
    P = perms(pre);
    score = nan(size(P,1),1);
    for i=1:size(P,1)
        score(i) = sum(sum(sim(survival(P(i,:),post))));
    end
    P = P(find(score==max(score),1),:);
    survival(pre,post) = survival(P,post);
else
    P = perms(post);
    score = nan(size(P,1),1);
    for i=1:size(P,1)
        score(i) = sum(sum(sim(survival(pre,P(i,:)))));
    end
    P = P(find(score==max(score),1),:);
    survival(pre,post) = survival(pre,P);
end
end
