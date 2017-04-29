% For more info on the competition, see https://github.com/codeneuro/spikefinder
% Code written by Peter Rupprecht (2017), ptrrupprecht.wordpress.com


%% preliminary part for choosing the best parameter space (mainly to estimate
% the estimated delay of a calcium transient upon an action potential event


% The "optimal delay" is the delay that has been chosen in a parameter test
% for the time lag; the kurtosis is used as a proxy to predict this
% "optimal delay", as it was measured for the 10 training datasets, for yet
% unknown datasets; "optimal window" is roughly the temporal window over which the
% derivative is averaged
% optimal_window = [29    27    32    45    42    15    17    22    18    15];
% optimal_delay = [24    26    26    34    37     4     1     1    19     1];
optimal_delay = [24    24    20    52    40    10     2     2     18     2];
optimal_window = [27    18    15    45    18    12    12    18    19    15];

% mean kurtosis of the different datasets
Parameter = [ 4.8549    9.3126    6.9904    3.6841    4.6713   26.9578   54.8355   33.0293    7.9486   18.1246];

% Fit optimal_delay(kurtosis) -- since no model underlies this observation,
% I simply used a spline for fitting
[xData, yData] = prepareCurveData( Parameter, optimal_delay );
ft = fittype( 'smoothingspline' );
opts = fitoptions( ft );
opts.SmoothingParam = 0.0136;
[fitresult, gof] = fit( xData, yData, ft, opts );
% check that the fit is neither overfitting nor way off
figure(839),
for j = 1:10
    plot(Parameter(j),optimal_delay(j),'.');hold on;
    text(Parameter(j)+0.023,optimal_delay(j),num2str(j));
end
hold on; plot(fitresult); hold off;

[xData, yData] = prepareCurveData( Parameter, optimal_window );
% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( ft );
% opts.SmoothingParam = 0.9036;
opts.SmoothingParam = 0.0046;
[fitresult2, gof] = fit( xData, yData, ft, opts );

figure(859),
for j = 1:10
    plot(Parameter(j),optimal_window(j),'.');hold on;
    text(Parameter(j)+0.023,optimal_window(j),num2str(j));
end
hold on; plot(fitresult2); hold off;



%% Main part of the program
% For each dataset, 'simNeuron' will contain the predictive/correlative
% values for each neuron
% Those will be combined (averaged/mean) in the vector 'simDataset',
% containing one value for each dataset

% the j-loop goes through all datasets

% the n-loop goes through each neuron of the chosen dataset

% the jj-loop serves to average across a finite time-window around the
% optimal delay point, which in turn is set by the fit based on the
% kurtosis values, as described above

clear simDataset final_prediction
for j = 1:10
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
        
        delay = fitresult(nanmean(Parameter(j)));
        window = fitresult2(nanmean(Parameter(j)));
        
        clear prediction
        CL = L_trace;
        for jj = 1:(window*2+1)
            delay1 = round(delay-(window+1)+jj);
            pprediction = CL - circshift(CL,delay1);%
            pprediction = circshift(pprediction,-round(delay1/2));
            pprediction = smooth(pprediction,6);
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
    
    simDataset(j) = median(simNeuron);
end
