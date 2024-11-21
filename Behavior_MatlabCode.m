%% Behavior analysis
% Written by Rita Gil and Mafalda Valente
% Modified last in 15.12.2023

% For every batch calculate the percent of pulsating light reports across sessions
% Join batches
% Fit of the average response and bootstrapp with replacement for the std calculation

% Correlation with fMRI steady-state data

%% Loading data

for batch=1:2
    if batch==1
        load('DataBatch1.mat')
        val=[-1,1,1,1,1,1];%Code for the rewarded port (-1 continuous port; 1 pulsating port).
        %Order: 2Hz, 8Hz, 15Hz, 20Hz, 25Hz and continuous light
        
    else
        load('DataBatch2.mat')
        
        val=[-1,1,1,1,1,1,1,1,1,1,1,1]; %Code for the rewarded port (-1 continuous port; 1 pulsating port).
        %Order: 1Hz ,2Hz ,8Hz ,12.5Hz ,15Hz ,16Hz ,18Hz ,20Hz ,25Hz ,40Hz ,50Hz and continuous light
    end
    
    out6 = out6(out6.training_level==4,:);
    animals = unique(out6.animal);
    cond = unique(out6.LED_index);
    session = unique(out6.session);
    
    percContsess = nan(length(animals),length(session),length(cond));
    
    for an=1:length(animals)
        
        for a = 1:length(session)
            
            for c=1:length(cond)
                percContsess(an-1,a,c)=sum(out6.success==val(c) & out6.LED_index==cond(c) & out6.session==session(a) & out6.animal==animals(an))*100 ...
                    /sum(out6.success~=0 & out6.session==session(a) & out6.LED_index==cond(c) & out6.animal==animals(an));
                
            end
        end
    end
    
    if batch==1
        
        behaviorDatasess=percContsess(:,:,[2 3 4 5 6 1]); %organizing frequencies in an ascending order
        
        behaviorDatasessRESH=reshape(behaviorDatasess, [length(animals)*length(session),length(cond)]);
        behaviorDatasessRESH(behaviorDatasessRESH==0)=nan;
        %Save data
        save('BehaviorBatch1.mat','behaviorDatasess','behaviorDatasessRESH')
    else
        
        behaviorDatasess2=percContsess(:,:,[2 3 4 5 6 7 8 9 10 11 12 1]); %organizing frequencies in an ascending order
        behaviorDatasessRESH2=reshape(behaviorDatasess2, [length(animals)*length(session),length(cond)]);
        behaviorDatasessRESH2(behaviorDatasessRESH2==0)=nan;
        %Save data
        save('BehaviorSess2.mat','behaviorDatasess2','behaviorDatasessRESH2')
    end
end

%% Join both batches

totanimal=7;%total number of animals in both batches
totcond=12; %total number of tested frequencies
totsess=size(behaviorDatasessRESH,1)+size(behaviorDatasessRESH2,1); %total number of sessions

finalsessionBatch=nan(totsess,totcond); %Vector with individual session information per condition
finalsessionBatch(1:length(behaviorDatasessRESH),[2,3,5,8,9,12])= behaviorDatasessRESH;
finalsessionBatch(length(behaviorDatasessRESH2)+1:totsess,:)=behaviorDatasessRESH2;

totavean=size(behaviorDatasess,1)+size(behaviorDatasess2,1);

finalsessionBatchan=nan(totanimal,totcond); %Vector with individual animal information per condition
finalsessionBatchan(1:3,[2,3,5,8,9,12])= squeeze(nanmean(behaviorDatasess,2));
finalsessionBatchan(4:7,:)=squeeze(nanmean(behaviorDatasess2,2));

save('BehaviorSessFINAL.mat','finalsessionBatch','finalsessionBatchan')
%% Fit of the average response and bootstrapp for the std calculation

xfreq=[1,2,8,12.5,15,16,18,20,25,40,50]; %frequency vector

%fitting a psychometric curve to the data (used in Wichmann, Felix A., and N. Jeremy Hill., Perception & psychophysics, 2001)
clear xData yData ymeans thresh1 d
ymeans=nanmean(finalsessionBatchan(:,1:end-1),1)/100; %use averaged animal information
[xData, yData]=prepareCurveData(xfreq,ymeans);
q1=fittype(@(a,b,c,d,x) d + (c-d)./(1+exp(-(2*a)*(x-b))),'independent',{'x'},'dependent',{'y'});
[f,~]=fit(xData,squeeze(yData),q1,'Lower', [0 xfreq(1) 0 0], 'Upper', [Inf xfreq(end) 1 1],'startpoint',[1,0,0.05, 0.95]);
thresh=f.b-(1/(2*f.a))*log(((f.c-f.d)/(0.50-f.d))-1); %Calculation of threshold of 50% choices to pulsating port ('chance level')


%Plotting fitted curve
close all
figure(1)
plot(xfreq,ymeans,'-o')
hold on
plot(f)
ylim([0 1])
fig1=figure(1);
fig1.Renderer='Painters';
saveas(gcf,'Fit','epsc')

%Bootstrapp with replacement for standard deviation calculation
it=50; %Number of iterations

d=finalsessionBatch(:,1:end-1); %use individual session information
bootstat  = bootstrp(it,@nanmean,d);

clear xData yData ymeans1 threshBOOT d
for a=1:it
    a
    ymeans1(a,:)=bootstat(a,:)/100;
    [xData, yData(a,:)]=prepareCurveData(xfreq,ymeans1(a,:));
    q1=fittype(@(a,b,c,d,x) d + (c-d)./(1+exp(-(2*a)*(x-b))),'independent',{'x'},'dependent',{'y'});
    [f,~]=fit(xData,squeeze(yData(a,:))',q1,'Lower', [0 xfreq(1) 0 0], 'Upper', [Inf xfreq(end) 1 1],'startpoint',[1,0,0.05, 0.95]);
    threshBOOT(a)=f.b-(1/(2*f.a))*log(((f.c-f.d)/(0.5-f.d))-1);
end

finalthreshhigh=mean(thresh);
finalthreshhighSTDBOOT=std(threshBOOT,0,2);

%% Final behavior plots

close all

%plot averaged and std
xfreqfinal=[1,2,8,12.5,15,16,18,20,25,40,50,100]; %frequencies tested in batch 2
xfreqfinal1=[2,8,15,20,25,100]; %frequencies tested in batch 1

%individual animals
close all
figure(1)
plot(xfreqfinal1,squeeze(nanmean(behaviorDatasess,2)),'o')
hold on
plot(xfreqfinal,squeeze(nanmean(behaviorDatasess2,2)))
plot(xfreqfinal,nanmean(finalsessionBatchan,1),'Linewidth',1.5)
ylim([0 100])
fig1=figure(1);
fig1.Renderer='Painters';
saveas(gcf,'Finalplot1','epsc')

%final fitted curve
close all
figure
plot(xfreqfinal,ymeans,'-o')
hold on
plot(f)
ylim([0 1])
fig1=figure(1);
fig1.Renderer='Painters';
saveas(gcf,'Finalplot2','epsc')


%overlap both plots and save plots

%% Correlation with fMRI (run this code for each ROI)

load('SteadyState_fMRIData.mat') %load fMRI data for each specific ROI
testedfreq=6; %correlated frequencies

BehaviorData=finalsessionBatch(:,[1,2,5,8,9,12]); %Behavior data using only frequencies that were used in fMRI
fMRIData=nan(length(BehsAreas),6);
fMRIData(1:length(fMRIData),:)=ROIX_fMRIsteadystateData(:,[2:7]); %we do not inlcude the 0.25Hz condition

%mean and std for fMRI and behavior data
MEANfMRI=[nanmean(fMRIData,1)];
MEANBEHAVIOR=[nanmean(BehaviorData,1)];

STDfMRI=[nanstd(fMRIData,0,1)];
STDBEHAVIOR=[nanstd(BehaviorData,0,1)];

%Plot the x and Y error bars
close all
figure('units','normalized','outerposition',[0 0 1 1])

for i=1:testedfreq
    errorbarxy(MEANfMRI(i), MEANBEHAVIOR(i), STDfMRI(i), STDBEHAVIOR(i))
    hold on
    pause
end

ylim([0 100]) %0 to 100% of flicker reports
xlim([a b]) %define x axis based on the error bars (equal for each ROI)
fig1=figure(1);
fig1.Renderer='Painters';
saveas(gcf,'ROIX_stdXY','epsc')

%Pearson correlation
close all
figure('units','normalized','outerposition',[0 0 1 1])
[RHO,PVAL,H] = corrplot([MEANfMRI',MEANBEHAVIOR'],'Type','Pearson');%manually modify x and y axis to correspond to the first plot with x and y error bars
[R,P,RL,RU] = corrcoef([MEANfMRI',MEANBEHAVIOR']);

fig1=figure(1);
fig1.Renderer='Painters';
saveas(gcf,'ROIX_corr','epsc')

%Resulting image will be an overlap of 'ROIX_corrCloud' and 'ROIX_stdXY'
