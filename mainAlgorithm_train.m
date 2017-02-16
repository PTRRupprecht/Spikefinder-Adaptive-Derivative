% example matlab script for loading spikefinder data
%
% for more info see https://github.com/codeneuro/spikefinder

clear Kurtouis
for j = 1:10
    dataset = num2str(j);
    calcium_train = csvread([dataset '.train.calcium.csv']);
    spike_train = csvread([dataset '.train.spikes.csv']);
    
    for n = 1:size(calcium_train,2)
        L_trace = calcium_train(:,n);
        indizes = find(~isnan(L_trace) & ((L_trace~=0 | circshift(L_trace,1)~=0)) );
        L_trace = L_trace(indizes(2:end));
        
        L_trace = (L_trace-median(L_trace))/std(L_trace);
        Kurtouis(n,j) = kurtosis(L_trace);
    end
end

optimal_delay = [28 22 29 32 49 9 7 9 17 8]; % optimal delay measured for the training dataset
optimal_delay = [24    26    26    34    37     4     1     1    19     1];
Kurtouis(Kurtouis==0) = NaN;
figure(839), hold on;
for j = 1:10
    plot(nanmean(Kurtouis(:,j)),optimal_delay(j)*2+2,'.');
    text(nanmean(Kurtouis(:,j))+0.023,optimal_delay(j)*2+2,num2str(j));
end

[xData, yData] = prepareCurveData( nanmean(Kurtouis), optimal_delay );
% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( ft );
opts.SmoothingParam = 0.0136;
[fitresult, gof] = fit( xData, yData, ft, opts );

figure(839), hold on;
for j = 1:10
    plot(nanmean(Kurtouis(:,j)),optimal_delay(j),'.');
    text(nanmean(Kurtouis(:,j))+0.023,optimal_delay(j),num2str(j));
end
hold on; plot(fitresult)

clear simDataset
%% load dataset
for j = 1:10
    dataset = num2str(j);
    calcium_train = csvread([dataset '.train.calcium.csv']);
    spike_train = csvread([dataset '.train.spikes.csv']);

    clear simNeuron final_prediction
    for n = 1:size(calcium_train,2)
        L_trace = calcium_train(:,n);
        final_prediction(:,n) = L_trace;
        indizes = find(~isnan(L_trace) & ((L_trace~=0 | circshift(L_trace,1)~=0)) );
        indizesM = find(isnan(L_trace) | ((L_trace==0 & circshift(L_trace,1)==0)) );
        L_trace = L_trace(indizes(2:end));
        L_trace = (L_trace-median(L_trace))/std(L_trace);
        S_trace = spike_train(indizes(2:end),n);
        
        delay = fitresult(nanmean(Kurtouis(:,j)));
        
        clear prediction
        CL = L_trace;
        for jj = 1:(12*2+1)
            delay1 = round(delay-(12+1)+jj);
            pprediction = CL - circshift(CL,delay1);%
            pprediction = circshift(pprediction,-round(delay1/2));
            pprediction = smooth(pprediction,5);
            pprediction( pprediction < 0*std(pprediction) ) = 0;
            prediction(:,jj) = pprediction/mean(pprediction);
        end
        if delay > 15
            prediction = nanmean(prediction')';
        else
            prediction = nanmedian(prediction')';
        end
        final_prediction(indizes(2:end),n) = prediction;
        final_prediction(indizesM,n) = NaN;
        simNeuron(n) = corr(conv(S_trace,[1 1 1 1],'same'),conv(prediction,[1 1 1 1],'same'));
    end

    spike_pred = final_prediction;     % fill spike_pred with your own predictions
    csvwrite([dataset '.train.pred.csv'], spike_pred);
    
    simDataset(j) = mean(simNeuron);
end

mean(simNeuron-red)


STM = [0.46 0.44 0.44 0.44 0.28 0.52 0.49 0.40 0.41 0.55];
OOP = [0.35 0.35 0.38 0.32 0.21 0.60 0.64 0.53 0.34 0.67];
Suite2p = [0.45 0.45 0.47 0.48 0.48 0.57 0.63 0.59 0.45 0.72];

mean(simDataset-STM)
mean(simDataset-OOP)
mean(simDataset-Suite2p)

figure(9), imagesc(SSDx)
% end
% end
figure(5), imagesc(SS,[0.2 0.6]); 


figure, imagesc(SSD)






% saving predictions
spike_pred = calcium_train;     % fill spike_pred with your own predictions
csvwrite([dataset '.train.pred.csv'], spike_pred);
