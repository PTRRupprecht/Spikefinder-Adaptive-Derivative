
% For more detailed explanations, see the comments on 'mainAlgorith_train.m'

%% Choosing the best parameters for predictions

load('Kurtosis.mat');

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

%% main part of the program
clear simDataset
for j = 4%1:10
    % load dataset
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
