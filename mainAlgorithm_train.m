% For more info on the competition, see https://github.com/codeneuro/spikefinder
% Code written by Peter Rupprecht (2017), ptrrupprecht.wordpress.com

% kurtosis for each dataset can predict the time lag that is lateron used
% for prediction of spike timing
clear Kurtouis
for j = 1%:10
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

% The "optimal delay" is the delay that has been chosen in a parameter test
% for the time lag; the kurtosis is used as a proxy to predict this
% "optimal delay", as it was measured for the 10 training datasets, for yet
% unknown datasets
% optimal_delay = [28 22 29 32 49 9 7 9 17 8];
optimal_delay = [24    26    26    34    37     4     1     1    19     1];
Kurtouis(Kurtouis==0) = NaN;
figure(839), hold on;
for j = 1:10
    plot(nanmean(Kurtouis(:,j)),optimal_delay(j)*2+2,'.');
    text(nanmean(Kurtouis(:,j))+0.023,optimal_delay(j)*2+2,num2str(j));
end
% Fit optimal_delay(kurtosis) -- since no model underlies this observation,
% I simply used a spline for fitting
[xData, yData] = prepareCurveData( nanmean(Kurtouis), optimal_delay );
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

mean(simNeuron)
