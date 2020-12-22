function [AUC, tpr, fpr, p, AUC_shuffle] = calc_roc(targets, outputs, n_shuffle)
%CALC_ROC  Compute ROC curve and the area under the curve (AUC)
%   CALC_ROC(targets, outputs, n_shuffle) is based on Fawcett, 2003;
%   Algorithm 4 and additionally computes shuffled AUC values and p-value
%   (permutation test)
%   It returns:
%   AUC          area under the ROC curve ("choice probability")
%   tpr          true positive rate
%   fpr          false positive rate
%   p            p value
%   AUC_shuffle  AUC values for shuffled targets (used to establish
%                significance bounds)


% usage and self test
if nargin == 0
    
    % --- first example
    
    x1 = 15 + 2.5 * randn(1,40);
    x2 = 20 + 2.0 * randn(1,50);
    
    x = [x1 x2];                            % e.g. vector containing spike count of each trial
    r = [zeros(size(x1)) ones(size(x2))];   % e.g. vector of the monkeys choice in each trial (0, 1)
    
    [CP, tpr, fpr, p] = calc_roc(r, x, 1000);
    fprintf('example 1:\n');
    fprintf('CP: %2.5f (p=%2.5f)\n', CP, p);
    
    figure
    xmin = min([x1 x2]);
    xmax = max([x1 x2]);
    x = linspace(xmin,xmax,20);
    a1 = hist(x1,x);
    a2 = hist(x2,x);
    
    subplot(121)
    bar(x,[a1; a2]',1.5);
    colormap(gray);
    set(gca,'xlim',[xmin xmax])
    axis square
    legend('x1','x2')
    legend(gca,'boxoff')
    
    subplot(122)
    plot(fpr,tpr,'x');
    xlabel('false positive')
    ylabel('true positive')
    hold on;
    plot([0,1],[0,1],'k');
    axis square
    axis([0 1 0 1])
    
    % --- second example
    
    mean_x1 = 15;
    mean_x2 = mean_x1 * 1.1;
    std_x1 = sqrt(mean_x1);
    std_x2 = sqrt(mean_x2);
    x1 = randn(1,1e6) * std_x1 + mean_x1;
    x2 = randn(1,1e6) * std_x2 + mean_x2;
    
    x = [x1 x2];                            % e.g. vector containing spike count of each trial
    r = [zeros(size(x1)) ones(size(x2))];   % e.g. vector of the monkeys choice in each trial (0, 1)
    
    CP = calc_roc(r, x);
    fprintf('\nexample 2:\n');
    fprintf('area under ROC curve (numerical): %2.5f\n', CP);
    fprintf('area under ROC curve (theoretical): %2.5f\n', ...
        0.5 * erfc((mean_x1 - mean_x2)/(sqrt(2)*sqrt(std_x1^2+std_x2^2))));
    
    return
end


% ------------------------------------------------------------------------

TT = targets(:)';

if nargin < 3
    n_shuffle = 0;
end

targets = targets(:)';
outputs = outputs(:)';

% --- sort outputs
[L sort_idx] = sort(outputs,'descend');
targets = targets(sort_idx);

% ---- normalize to [0 1]
L = L - L(end);     % L - min(L)
if L(1) > 0         % do not normalize if all outputs are equal (-> division by zero)
    L = L ./ L(1);  % L ./ max(L)
end

% --- calculate true and false positive rate
TP = [0 cumsum(targets == 1)];
FP = [0 cumsum(targets == 0)];

idx = [true, diff(L) < 0, true];        % add first and last index as well (generates the points (0,0) and (1,1) on the ROC curve

fpr = FP(idx);
tpr = TP(idx);

fpr = fpr / FP(end);      % divide by total number of pos / neg targets
tpr = tpr / TP(end);

% calculate area under the curve
% AUC = trapz(fpr,tpr);
AUC = (tpr(1:end-1) + tpr(2:end)) * diff(fpr).' / 2;        % more efficient thatn trapz


if nargout >= 4
    
    % --- calculate shuffled ROCs
    AUC_shuffle = zeros(1,n_shuffle);
    for i=1:n_shuffle
        targets_ = TT(randperm(numel(TT)));             % shuffle targets
        targets_ = targets_(sort_idx);                  % sort according to sort_idx
        % sorting is important for doing trial-based shuffling (where the
        % same random numbers are using for each time bin but the sort
        % order is different) => otherwise we'll get the same shuffled
        % values in each time bin!
        % sorting does not matter if the random number generator is not
        % reset always before calling calc_roc.m
        
        TP_ = [0 cumsum(targets_ == 1)];
        FP_ = [0 cumsum(targets_ == 0)];
        fpr_ = FP_(idx) / FP_(end);
        tpr_ = TP_(idx) / TP_(end);
        
        AUC_shuffle(i) = (tpr_(1:end-1) + tpr_(2:end)) * diff(fpr_).' / 2;        % more efficient thatn trapz
    end
    
    % --- calculate p value corresponding to a two-tailed test (the null
    % hyptothesis is that s = 0.5)
    AUC_left = 0.5 - abs(AUC-0.5);
    AUC_right = 0.5 + abs(AUC-0.5);
    p = sum(AUC_shuffle <= AUC_left) / n_shuffle  + sum(AUC_shuffle >= AUC_right) / n_shuffle;
    
    if (AUC_left == AUC_right)
        % subtract the shuffled trials that are counted twice if the two
        % thresholds are equal (pathological case where s_left == s_right ==
        % s == s_shuffle)
        p = p - sum(AUC_shuffle == AUC_right)/n_shuffle;
    end
    
end


