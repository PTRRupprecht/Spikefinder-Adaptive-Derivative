% For more info on the competition, see https://github.com/codeneuro/spikefinder
% Code written by Peter Rupprecht (2017), ptrrupprecht.wordpress.com

Parameter = [5.7304    8.7100    6.8567    2.9975    6.9790]; %Kurtosis of test data set

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
for j = 1:5
    % load dataset
    dataset = num2str(j);
    calcium_train = csvread([dataset '.test.calcium.csv']);

    clear simNeuron final_prediction
    for n = 1:size(calcium_train,2)
        L_trace = calcium_train(:,n);
        final_prediction(:,n) = L_trace;
        indizes = find(~isnan(L_trace) & ((L_trace~=0 | circshift(L_trace,1)~=0)) );
        indizesM = find(isnan(L_trace) | ((L_trace==0 & circshift(L_trace,1)==0)) );
        L_trace = L_trace(indizes(2:end));
        L_trace = (L_trace-median(L_trace))/std(L_trace);
        
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
    end

    spike_pred = final_prediction;     % fill spike_pred with your own predictions
    csvwrite([dataset '.test.spikes.csv'], spike_pred);
end



