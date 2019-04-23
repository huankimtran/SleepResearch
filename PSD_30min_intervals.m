% Code to get power in 30-min intervals.
% In addition to typical frequency bands (SWO, Delta, Theta, Alpha, Beta,
% Gamma) also want to do calculations for high delta (2.75 - 5 Hz) and low
% alpha (6.25 - 10 Hz)

%% Parameters

subs=1; % MAX 5804
m_per=20; % MAX number of periods
EDFPATH='shhs1-%d.edf';
STAGEPATH='shhs1-%d-staging.csv';
FILESAVEPATH='30minPSDdata.mat';
BASE=200000;

%% Variable declarations
periods=NaN(subs,1); % # of 30 min periods each person has

% Matrix for each power band. Matrix position (i,j) corresponds to
% (subject, period)

% Typical frequency bands (raw and normalized)
swo_raw=NaN(subs,m_per); % (0-1 Hz)
delta_raw=NaN(subs,m_per); % (1-4 Hz)
theta_raw=NaN(subs,m_per); % (4-8 Hz)
alpha_raw=NaN(subs,m_per); % (8-13 Hz)
beta_raw=NaN(subs,m_per); % (13-30 Hz)
gamma_raw=NaN(subs,m_per); % (30-50 Hz)

swo_norm=NaN(subs,m_per);
delta_norm=NaN(subs,m_per);
theta_norm=NaN(subs,m_per);
alpha_norm=NaN(subs,m_per);
beta_norm=NaN(subs,m_per);
gamma_norm=NaN(subs,m_per);

% Couple of special ones here, to check our what that paper had to say
high_delta_raw=NaN(subs,m_per); % (2.75 - 5 Hz)
low_alpha_raw=NaN(subs,m_per);

high_delta_norm=NaN(subs,m_per); 
low_alpha_norm=NaN(subs,m_per); % (6.25 - 10 Hz)

REM_raw=NaN(subs,m_per); % All REM
NONREM_raw=NaN(subs,m_per); %NonREM

REM_norm=NaN(subs,m_per); 
NONREM_norm=NaN(subs,m_per);

% Total power
total_power=NaN(subs,m_per);
%% Calculations
for i=1:subs
    fprintf('Processing no. %d of %d\n',i,subs);
    
    % Check if PSG(.edf) and staging(.csv) files exist
%     a=exist(strcat('/project/bsheth/sleepdata/nsrr/shhs/polysomnography/edfs/shhs1/shhs1-', int2str(i+200000), '.edf'), 'file');
%     b=exist(strcat('/project/bsheth/sleepdata/nsrr/shhs/polysomnography/annotations-staging/shhs1/shhs1-', int2str(i+200000), '-staging.csv'), 'file');
%    a=exist(strcat('G:\nsrr\shhs\polysomnography\edfs\shhs1\shhs1-', int2str(i+200000), '.edf'),'file');
%    b=exist(strcat('G:\nsrr\shhs\polysomnography\annotations-staging\shhs1\shhs1-',int2str(i+200000),'-staging.csv'),'file');
    a=exist(sprintf(EDFPATH,i+BASE),'file');
    b=exist(sprintf(STAGEPATH,i+BASE),'file');
    
    % If either does not exist - go to next loop iteration (next subject)
    if(a==0||b==0)
        fprintf('Empty file with subject %d',i);
        continue;
    end
    
    % Load in edf.... hdr has information, record has actual data
%     [hdr,record]=edfread(strcat('/project/bsheth/sleepdata/nsrr/shhs/polysomnography/edfs/shhs1/shhs1-', int2str(i+200000), '.edf'));
    [hdr,record]=edfread(sprintf(EDFPATH,i+BASE));
    
    % Find EEG channel using channel labels in header
    eeg=record(find(ismember(hdr.label,'EEG')==1),:);
    % Pick out sampling rate
    sr=hdr.samples(find(ismember(hdr.label,'EEG')==1));
    
    % Load in staging data
%     stages=csvread(strcat('/project/bsheth/sleepdata/nsrr/shhs/polysomnography/annotations-staging/shhs1/shhs1-',int2str(i+200000),'-staging.csv'),1,1);
    stages=csvread(sprintf(STAGEPATH,i+BASE),1,1);
    
    % Figure out how many FULL 30-min periods the person has.
    start=find(stages~=0,1,'first'); % Gives index of first sleep epoch
    finish=find(stages~=0,1,'last'); % Gives index of last sleep epoch
    % Filter through
    for j=1:m_per
        if ((j-1)*30*2)+start>finish || (j*60)+start>finish
            % IF checks both ends of 30-min period. TRUE if hits last
            % period
            
            % Note: 60 epochs in 30-min period
            periods(i)=j-1; % Record # of 30 min periods
            break; % Jump out of for loop, no need to keep going
        end
    end
    
    %Calculate power for each period
    for j=1:periods(i)
        nfft=30*sr/2;
        epochs=60; % 60 epochs per 30-min period
        
        pxx=zeros(floor(nfft/2)+1,1); % Accumulative power vector
        
        % Loop for power calculations. Doing power over each epoch
        % individually (to avoid unnecessary error), storing sum into pxx.
        % Since going in continuous 30-min periods, really don't need to go
        % epoch by epoch, but no harm in doing so.
        for k=1:epochs
            m=j*60+k;
            data=eeg((m-1)*30*sr+1:m*30*sr+1);
            %temppwr is power spectral density, f is frequency
            [temppwr,f]=pwelch(data,kaiser(nfft),floor(nfft/2),nfft,sr,'power');
            pxx=pxx+temppwr;
        end
        
        % Take average of power in 30-min interval by dividing by # of
        % epochs summed over.
        pxx=pxx./epochs;
        
        zero=find(f==0);
        zero=zero+1; % Take out DC bias level
        fifty=find(f==50);
        
        total_power(i,j)=sum(pxx(zero:fifty,1)); % Will use this for normalized power calculations
        
        % SWO band (0-1 Hz)
        low=find(f==0);
        low=low+1; % Again, ignore DC level
        high=find(f==1);
        swo_raw(i,j)=sum(pxx(low:high,1));
        swo_norm(i,j)=swo_raw(i,j)/total_power(i,j);
        
        % Delta band (1-4 Hz)
        low=high;
        high=find(f==4);
        delta_raw(i,j)=sum(pxx(low:high,1));
        delta_norm(i,j)=delta_raw(i,j)/total_power(i,j);
        
        % Theta band (4-8 Hz)
        low=high;
        high=find(f==8);
        theta_raw(i,j)=sum(pxx(low:high,1));
        theta_norm(i,j)=theta_raw(i,j)/total_power(i,j);
        
        % Alpha band (8-13 Hz)
        low=high;
        high=find(f==13);
        alpha_raw(i,j)=sum(pxx(low:high,1));
        alpha_norm(i,j)=alpha_raw(i,j)/total_power(i,j);
        
        % Beta band (13-30 Hz)
        low=high;
        high=find(f==30);
        beta_raw(i,j)=sum(pxx(low:high,1));
        beta_norm(i,j)=beta_raw(i,j)/total_power(i,j);
        
        % Gamma band (30-50 Hz)
        low=high;
        high=find(f==50);
        gamma_raw(i,j)=sum(pxx(low:high,1));
        gamma_norm(i,j)=gamma_raw(i,j)/total_power(i,j);
        
        % High delta (2.75-5 Hz)
%         low=find(f==2.75); % 2.75 will not be in f, so need to find closest
        [~,low]=min(abs(f-2.75)); % Will find closest value in f to 2.75
        high=find(f==5);
        high_delta_raw(i,j)=sum(pxx(low:high,1));
        high_delta_norm(i,j)=high_delta_raw(i,j)/total_power(i,j);
        
        % Low alpha (6.25-10 Hz)
%         low=find(f==6.25); % Similar situation to f==2.75 above
        [~,low]=min(abs(f-6.25)); % Will find closest value in f to 6.25
        high=find(f==10);
        low_alpha_raw(i,j)=sum(pxx(low:high,1));
        low_alpha_norm(i,j)=low_alpha_raw(i,j)/total_power(i,j);
    end
end    
save(FILESAVEPATH,'periods','swo_raw','swo_norm','delta_raw','delta_norm',...
    'theta_raw','theta_norm','alpha_raw','alpha_norm','beta_raw',...
    'beta_norm','gamma_raw','gamma_norm','high_delta_raw',...
    'high_delta_norm','low_alpha_raw','low_alpha_norm','total_power');
fprintf('DONE');

