
function [out, pvals] = a_peaks (resdata, opt_plot)

% out = a_peaks (resdata)
%
%   resdata -> Structure array with a fly's data in each row.
%
%   opt_plot -> (Optional) If set to 1, the function will produce a plot
%               with the results. It is set to 0 by default.
%
% This function will calculate and add to the original structure array 
% fields containing calculations of the peaks of the main "velocity" 
% parameters in three particular time periods. 
%
% It will also perform statistical tests to compare the different peaks and
% can return (optionally) the p-values obtained. The matrix returned will
% have a row for each comparison (see below) and a column for each
% parameter.


out = resdata; 
if nargin < 2, ploplo = 0; elseif nargin==2, ploplo = opt_plot; end

if ploplo == 1, figure('position',[264 102 111 809]); end


for fly = 1:length(out) % Iterates through each fly in the structure
    
    % Resets variables to make sure it substitutes any pre-existing version
    out(fly).ppeaks = [];
    out(fly).vpeaks = [];
    out(fly).vypeaks = [];
    out(fly).angvpeaks = [];
    out(fly).curvpeaks = [];
    out(fly).tpeaks = [];
    
    % ---------------------------------------------------------------------
    % Time windows are defined
    
    stind = find(out(fly).stimulus(:,1) ~= 0); % Detects stimulus time    
    
    if isempty(stind) % If there is no stimulus (i.e. blank trial) arbitrary
                      % time windows are defined.
        prea=1; preb=3500;
        dura=1; durb=3500;
        posa=1; posb=3500;
    else % If there is a stimulus
        
        prea = 1; % All time before odor
        preb = stind(1)-1; 

        dura = stind(1) + 100; % 2 to 3 seconds after odor onset
        durb = stind(1) + 100 + 150;
        
        posa = stind(end) + 50; % 1 to 3 seconds after odor offset
        posb = stind(end) + 50 + 100;
    end

    % ---------------------------------------------------------------------
    
    
    % Iterates through each trial and collects the peak measurements
    for trial = 1:size(out(fly).v,2)
     
        % PROBABILITY OF MOVEMENT peaks
        out(fly).ppeaks(trial,1) = nanmean(out(fly).pmovew(prea:preb,trial));
        out(fly).ppeaks(trial,2) = nanmean(out(fly).pmovew(dura:durb,trial));
        out(fly).ppeaks(trial,3) = nanmean(out(fly).pmovew(posa:posb,trial));
        
        % UPWIND VELOCITY peaks
        out(fly).vypeaks(trial,1) = nanmean(out(fly).vymovew(prea:preb,trial));
        out(fly).vypeaks(trial,2) = nanmean(out(fly).vymovew(dura:durb,trial));
        out(fly).vypeaks(trial,3) = nanmean(out(fly).vymovew(posa:posb,trial));
        
        % GROUND SPEED peaks
        out(fly).vpeaks(trial,1) = nanmean(out(fly).vmovew(prea:preb,trial));
        out(fly).vpeaks(trial,2) = nanmean(out(fly).vmovew(dura:durb,trial));
        out(fly).vpeaks(trial,3) = nanmean(out(fly).vmovew(posa:posb,trial));
        
        % ANGULAR VELOCITY peaks
        out(fly).angvpeaks(trial,1) = nanmean(out(fly).angvturnw(prea:preb,trial));
        out(fly).angvpeaks(trial,2) = nanmean(out(fly).angvturnw(dura:durb,trial));
        out(fly).angvpeaks(trial,3) = nanmean(out(fly).angvturnw(posa:posb,trial));
        
        % CURVATURE peaks
        out(fly).curvpeaks(trial,1) = nanmean(out(fly).curvaturew(prea:preb,trial));
        out(fly).curvpeaks(trial,2) = nanmean(out(fly).curvaturew(dura:durb,trial));
        out(fly).curvpeaks(trial,3) = nanmean(out(fly).curvaturew(posa:posb,trial));
%         out(fly).curvpeaks(trial,3) = nanmean(out(fly).curvature(posa:posb,trial)) - out(fly).curvpeaks(trial,1)for;
        
        % PROBABILITY OF TURNING peaks
        out(fly).tpeaks(trial,1) = nanmean(out(fly).pturnw(prea:preb,trial));
        out(fly).tpeaks(trial,2) = nanmean(out(fly).pturnw(dura:durb,trial));
        out(fly).tpeaks(trial,3) = nanmean(out(fly).pturnw(posa:posb,trial));
    end
    
    
    % Saves fly averages
    ppeaks(fly,:) = nanmean(out(fly).ppeaks);
    vpeaks(fly,:) = nanmean(out(fly).vpeaks);
    vypeaks(fly,:) = nanmean(out(fly).vypeaks);
    tpeaks(fly,:) = nanmean(out(fly).tpeaks);
    angvpeaks(fly,:) = nanmean(out(fly).angvpeaks);
    curvpeaks(fly,:) = nanmean(out(fly).curvpeaks);
    
    
    % Plots if it was requested
        if ploplo == 1
            
            subplot(6,1,1);
            hold on;
            plot(nanmean(out(fly).ppeaks),'color',[.7 .7 .7]);
            
            subplot(6,1,2);
            hold on;
            plot(nanmean(out(fly).vypeaks),'color',[.7 .7 .7]);
            
            subplot(6,1,3);
            hold on;
            plot(nanmean(out(fly).vpeaks),'color',[.7 .7 .7]);
            
            subplot(6,1,4);
            hold on;
            plot(nanmean(out(fly).angvpeaks),'color',[.7 .7 .7]);
            
            subplot(6,1,5);
            hold on;
            plot(nanmean(out(fly).curvpeaks),'color',[.7 .7 .7]);
            
            subplot(6,1,6);
            hold on;
            plot(nanmean(out(fly).tpeaks),'color',[.7 .7 .7]);
        end     
end


% IT PERFORMS STATISTICAL ANALYSES comparing all averages from each fly
dats = {'ppeaks','vypeaks','vpeaks','angvpeaks','curvpeaks','tpeaks'};
titles = {'Probability of movement','Upwind velocity','Ground speed',...
    'Angular velocity','Curvature','Turn probability'};
units = {'Probability','mm/s','mm/s','deg/s','deg/mm','probability'};
alpha = 0.0167; % alpha level, adjusted according to three multiple comparisons
                % by the Bonferroni method (or adjust to whatever strikes your fancy)
 

for i = 2:6 % Iterates through each parameter and plots and does stats
    
        subplot(6,1,i);
        box off;
        dat = eval(dats{i}); %vertcat(out.(dats{i}));
        if ~isempty(dat == inf)
            dat(dat==inf) = nan; % Substitutes Inf values
        end
        
%         hold on
        plot(nanmean(dat),'k','linewidth',2);
        title(titles{i});
        % Paired Wilcoxon signed rank test
        pvals(1,i) = signrank(dat(:,1),dat(:,2)); % Before VS During 
        pvals(2,i) = signrank(dat(:,1),dat(:,3)); % Before VS After
        pvals(3,i) = signrank(dat(:,2),dat(:,3)); % During VS After
        % Some figure adjustment
        ylabel(units{i})
        set(gca,'XTick',[1, 2, 3])
        set(gca,'XTickLabel',{'Pre-odor';'Odor';'Post-odor'})
        xtickangle(30)
        axis tight
        xlim([0.7 3.3]);
        % Plots significance
        yl = ylim;
        if  pvals(1,i) <= alpha
            plot(1.5,yl(2)-((yl(2)-yl(1))/5),'.g','markersize',10); % Before VS During significance is plotted as a GREEN dot
        end
        if pvals(2,i) <= alpha
            plot(2,yl(1)+((yl(2)-yl(1))/6),'.b','markersize',10); % Before VS After significance is plotted as a BLUE dot
        end
        if pvals(3,i) <= alpha
            plot(2.5,yl(2)-((yl(2)-yl(1))/5),'.g','markersize',10); % During VS After significance is plotted as a GREEN dot
        end  
        
end
  
% If the opt_figure is requested as a 2, additionally it produces an
% additional plot for figures
if ploplo < 1
    close;
end    
% if ploplo == 2
%     figurepeaks (out);
% end
        
        
        
        
        
        
        