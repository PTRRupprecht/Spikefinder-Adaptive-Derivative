% example matlab script for loading spikefinder data
% for more info see https://github.com/codeneuro/spikefinder
% modified and enhanced by Peter Rupprecht (2017), ptrrupprecht.wordpress.com

% load dataset; normally, there are 1:10 datasets; here, only one dataset is included in the repository
dataset = '1';
calcium_train = csvread([dataset '.train.calcium.csv']);
spike_train = csvread([dataset '.train.spikes.csv']);

% plot example neuron
% figure;
% n = 5; % index of neuron
% t = (0:length(calcium_train(:,n))-1)/100;     % 100Hz sampling rate, each bin 10 ms
% plot(t,zscore(calcium_train(:,n)),'k'); hold on
% plot(t,spike_train(:,5)-2,'r')
% xlim([t(1) 400]); ylim([-2.5 7]); xlabel('Time (s)'); ylabel('Fluorescence / Spike rate')
% legend('calcium','spikes')


%% load dataset and calculate filter / calcium response function using a cross-correlation function between spikes and calcium traces
counter = 1; clear filterXX CCC
for j = 1 % datasets : normally, there are 1:10 datasets; here, only one dataset is included in the repository
    dataset = num2str(j);
    calcium_train = csvread([dataset '.train.calcium.csv']);
    spike_train = csvread([dataset '.train.spikes.csv']);

    for n = 1:size(calcium_train,2)
        L_trace = calcium_train(:,n);
        indizes = find(~isnan(L_trace) & ((L_trace~=0 | circshift(L_trace,1)~=0)) );
        L_trace = L_trace(indizes);
        L_trace = (L_trace-median(L_trace))/std(L_trace);
        S_trace = spike_train(indizes,n);

        timet = (( 1:((numel(L_trace)-1)*2+1) ) - (numel(L_trace)))/100;

        ffx = xcorr(L_trace,S_trace,'unbiased');

        filterX = ffx((numel(L_trace)-20):(numel(L_trace)+250)); % window of 2 sec
        filterX(21) = (filterX(20)+filterX(22))/2;
        filterXX(counter,:) = (filterX - mean(filterX))/std(filterX);
        filterXXX = nanmean(filterXX);
        
        prediction = conv(S_trace,filterX,'same');
        prediction = circshift(prediction,115);
        prediction = (prediction-median(prediction))/std(prediction);

        CCC(j,n) = corr(prediction,L_trace);
        
        F = zeros(size(L_trace));
        for k = 1:(numel(L_trace)-numel(filterX)+1)
            F(k) = corr(filterX,L_trace(k:(k+numel(filterX)-1)));
        end
        
        counter = counter + 1;
    end
end
% plot the filters
figure, imagesc(([1:250]-20)/100,[],filterXX);
xlabel('time (sec)');
ylabel('neuron index');

% This plot makes more sense once j = 1:10 (all datasets included)
% it shows the scatter of correlation coefficients (indicating quality of
% the measurements) for neurons in one dataset
figure(243)
for j = 1:size(CCC,1)
    CC = CCC(j,:);
    CC = CC(CC~=0);
    plot(j*ones(size(CC))+ (rand(size(CC))-0.5)*0.5, CC,'.'); hold on
end
CCC(CCC==0) = NaN;
figure(244); boxplot(CCC')

