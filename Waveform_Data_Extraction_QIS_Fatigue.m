%%%%%% This code is designed to find waveform files for valid signals and measure their information entropy %%%%%%%%%%%%% 

clc;
clear all;

load('AE_Data_Extraction_Fatigue_04-Mar-2021')
lib = load('QIS_Test_lib','Test');
lib = lib.Test;

%% Get the list of waveform files in different waveform folders

for i= 1:size(AE_data,2) % Extract all the information needed
    
    for j=1:size(lib,2) % Find the index of it
        if contains(lib(j).name, AE_data(i).test_name) == 1
            ind_fldr = j;
        end
    end
    
    AE_loc = [lib(ind_fldr).AE_loc,'\'];
    
    d = dir(AE_loc);
    dfolders = d([d(:).isdir]);
    dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
    
    for k=1:size(dfolders,1)
        
        wave_fldr_loc = [AE_loc,dfolders(k).name,'\'];
        
        t_name = struct2table(dir(wave_fldr_loc)).name; % Find all the file names in the address
        t_name(ismember(t_name,'.')) = []; t_name(ismember(t_name,'..')) = [];
        
        Ch1_waves(k,:) = length(t_name(contains(t_name,'_1_')));
        Ch2_waves(k,:) = length(t_name(contains(t_name,'_2_')));
        
    end
    
    AE_wave(i).test_name = AE_data(i).test_name;
    
    AE_wave(i).address = AE_loc;
    
    name = cell2mat(extractfield(dfolders,'name').');
    AE_wave(i).folder_name = name;
    
    % This field will be used later to correct address when reading waveforms
    % Each row corresponds to the foldar name (waveform1, ...), column 1
    % corresponds to the number of Ch1 waveforms in the folder and column 2
    % corresponds to the number of Ch2 waveforms in the folder
    AE_wave(i).wavecount = [cumsum(Ch1_waves),cumsum(Ch2_waves)]; 
    
    AE_wave(i).sum = cumsum(AE_data(i).sum); % This is the cumulative sum of waveforms for each .DTA file
    
end

clear i j k ind AE_loc d dfolders wave_fldr_loc t_name Ch1_waves Ch2_waves name ind_fldr

%% Find correct waveforms and record their location

%%%%%%%%%%%%%%%%%%% Only waveforms recorded in Channel 2 will be used %%%%%%%%%%%%%%%%%%%%%%%%%
% The reason behind using Ch2 data is the duration for AE signals. When
% compared to Ch1, Ch2 signals show lower duration which mean more of the
% signal (or all of it in more than 95% of the cases) is captured by the
% waveform. The difference between duration and waveforms in Ch1 and Ch2 is
% due to the higher suseptibility of Ch1 to record noise becuase it is
% closer to actuator. The paper on AE Entropy in metals clearly show that
% the difference in the waveforms do not make a huge difference in entropy
% measurement but the other AE characteristics differ substantially.
% However, it should be remembered that in the case of that paper, the bin
% width was fixed. 

% depending on how fast is waveform finding for Ch2, it is then decided to
% either look at Ch1 or not. 
tic
for i = 1:size(AE_data,2)
    disp(['Test ',AE_data(i).test_name,' --> Start ...'])
    
    for n = 1:size(AE_wave(i).folder_name,1)
        wave_source(n,:) = [AE_wave(i).address,AE_wave(i).folder_name(n,:),'\'];
    end
    
    y = ones(size(AE_data(i).Ch2_dtf,1),1);
    
    AEwsum = AE_wave(i).sum;
    AEwcount = AE_wave(i).wavecount;
    wave_num = [];
    
    for j = 1:size(AE_data(i).Ch2_dtf,1) % 1:size(AE_data(i).Ch2_dtf,1)
        
        if rem(j,5000) == 0
            fprintf('\t i = %d - Progress %.1f percent\n',i, j./size(AE_data(i).Ch2_dtf,1)*100)
        end
        
        init_wnum = AE_data(i).Ch2_dtf(j,18);

        [Wave_filename,WaveAddress,wloc] = Wave_fname_loc(init_wnum,AEwsum,AEwcount,wave_source);

        wave_time = importdata(WaveAddress,' ',10);
        wave_time = wave_time.data;

        wave_num(j,:) = wloc;

        if j>1 && y(j)~=y(j-1)
            y(j) = y(j-1)-1;
        end
        
        if wave_time - AE_data(i).Ch2_dtf(j,2) < -2e-5
            while wave_time - AE_data(i).Ch2_dtf(j,2) < -2e-5
                
                wnum = AE_data(i).Ch2_dtf(j,18)+y(j);

                [Wave_filename,WaveAddress,wloc] = Wave_fname_loc(wnum,AEwsum,AEwcount,wave_source);
                
                wave_time = importdata(WaveAddress,' ',10);
                wave_time = wave_time.data;
                
                y(j) = y(j)+1;

            end
            wave_num(j,:) = wloc;
        end

        p=1;
        if wave_time - AE_data(i).Ch2_dtf(j,2) > 2e-5
            while wave_time - AE_data(i).Ch2_dtf(j,2) > 2e-5
                
                wnum = AE_data(i).Ch2_dtf(j,18)+y(j)-p;
                
                [Wave_filename,WaveAddress,wloc] = Wave_fname_loc(wnum,AEwsum,AEwcount,wave_source);
                
                wave_time = importdata(WaveAddress,' ',10);
                wave_time = wave_time.data;
                
                p=p+1;
                
            end
            y(j) = y(j)-(p-1);
            wave_num(j,:) = wloc;
        end

        if abs(wave_time - AE_data(i).Ch2_dtf(j,2)) > 2e-5
            
            fprintf('No waveform found for i=%d and j=%d \n',i,j)
            wave_num(j,:) = NaN;
            
        end
        
    end
    
    % Using the three # in the wave_num the correct address for each
    % waveform can be rebuilt
    AE_wave(i).wave_loc_col_header(:,1) = {'Folder #'};
    AE_wave(i).wave_loc_col_header(:,2) = {'Test #'};
    AE_wave(i).wave_loc_col_header(:,3) = {'Waveform #'};

    AE_wave(i).Ch2_wave_loc = wave_num;

    disp(['Test ',AE_data(i).test_name,' --> Done!'])
    clear wave_num wave_time

end

clear AEwcount AEwsum i j n p w_order wave_num wave_source wave_time wloc wnum y init_wnum WaveAddress Wave_filename ind_sum

toc 

%% Waveform verification, HLT selection, Binwidth measurement and Entropy calculation

disp('-------------------------')
disp('1- Double check waveforms are selected correctly')
disp('2- Find Hit Lockout Time')
disp('3- Find optimum binwidth')
disp('4- Calculate entropy based on optimum BinWidth')
disp(' ')

skew = 'G:\My Drive\Research\Fatigue on composites (with Dr Modarres)\TDA Project\Phase 2 Option\Tests\5RAA07_102218\AE\Results\Bin Width Study\Skewnes Coefficient.csv';
kurt = 'G:\My Drive\Research\Fatigue on composites (with Dr Modarres)\TDA Project\Phase 2 Option\Tests\5RAA07_102218\AE\Results\Bin Width Study\Kutosis Coefficient.csv';
sk_coeff = csvread(skew);
kur_coeff = csvread(kurt);

for i= 1:size(AE_data,2)
    
    disp(['Test ',AE_data(i).test_name,' --> Start ...'])
    
    wave_loc = AE_wave(i).Ch2_wave_loc;
    wave_loc_nan = length(wave_loc(isnan(wave_loc(:,1)))); % number of hits with no waveforms 
    
    fprintf('\tThis test has %d hits with no waveform! \n', wave_loc_nan)
    
    entropy = zeros(size(AE_wave(i).Ch2_wave_loc,1),1); % Generate variable for entropy
    
    for k = 1:size(AE_wave(i).Ch2_wave_loc,1)
        
        % Exclude hits with no waveforms (NaN values)
        if ~isnan(wave_loc(k,1))
            fldr_loc = wave_loc(k,1); t_loc = wave_loc(k,2); w_loc = wave_loc(k,3);
        else
            continue
        end
        
        % Generate the address for a waveform
        if t_loc < 10
            WaveAddress = [AE_wave(i).address,'Waveform',num2str(fldr_loc),'\Test_0',num2str(t_loc),'_2_',num2str(w_loc),'.txt'];
        else
            WaveAddress = [AE_wave(i).address,'Waveform',num2str(fldr_loc),'\Test_',num2str(t_loc),'_2_',num2str(w_loc),'.txt'];
        end
        
        % Import data
        wave_time = importdata(WaveAddress,' ',10);
        waveform = importdata(WaveAddress,' ',12);
        
        % Report progress
        if rem(k,5000) == 0
            fprintf('\t i = %d - Progress %.1f percent\n',i, k./size(AE_wave(i).Ch2_wave_loc,1)*100)
        end
        
        % Extract waveform information 
        waveform = waveform.data;
        wave_time = wave_time.data;
        
        % Double check waveforms are selected correctly
        if wave_time ~= AE_data(i).Ch2_dtf(k,2)
            disp('Error ! Wrong waveform !!')
            disp(['i = ',num2str(i),' k = ',num2str(k)])
            break
        end
        
        % Find Hit Lockout Time
        fl_waveform = flip(waveform);
        HLT_ind(k,:) = length(fl_waveform) - find(fl_waveform ~= 0,1) + 1;
        
        
        % Find optimum binwidth
        waveform = waveform(1:HLT_ind(k,1));
        
        % Find standard deviation, skewness, kurtosis for the waveform after HLT removal
        sd = std(waveform);
        sk = skewness(waveform);
        kur = kurtosis(waveform);
        n = length(waveform);

    	% Calcualte initial bin width before applying correction coefficients  
        bin = 3.49*sd*(n.^-0.33);
        
        % Find correction factors
        for j=1:length(sk_coeff)-1
           if abs(sk) < sk_coeff(1,1)
               wave_sk_coeff = sk_coeff(1,2);
               continue
           elseif (sk_coeff(j,1) < abs(sk) & abs(sk) < sk_coeff(j+1,1))
               wave_sk_coeff = sk_coeff(j+1,2);
               continue
           elseif sk_coeff(length(sk_coeff),1) < abs(sk)
               wave_sk_coeff = sk_coeff(length(sk_coeff),2);
               continue
           end
        end
        
        for m=1:length(kur_coeff)-1
           if abs(kur) < kur_coeff(1,1)
               wave_kur_coeff = kur_coeff(1,2);
               continue
           elseif (kur_coeff(m,1) < abs(kur) & abs(kur) < kur_coeff(m+1,1))
               wave_kur_coeff = kur_coeff(m+1,2);
               continue
           elseif kur_coeff(length(kur_coeff),1) < abs(kur)
               wave_kur_coeff = kur_coeff(length(kur_coeff),2);
               continue
           end
        end
        
        % Modify the binwodth according to correction factors
        bin_mod(k,1) = bin.*wave_sk_coeff.*wave_kur_coeff;
        
        
        %%%%% Measure entropy based on the optimum bin width
        
        BinWidth = bin_mod(k,1);
        
        [prob,edge] = histcounts(waveform,'Binwidth',BinWidth,'Normalization','probability');
                
        ent_temp = zeros(length(prob),3);
        ent_temp(:,1) = prob.';
        
        for n = 1:size(ent_temp,1) %length(ent_temp)
            if ent_temp(n,1) == 0
                ent_temp(n,2) = 0;
            else
                ent_temp(n,2) = log(ent_temp(n,1));
            end
            ent_temp(n,3) = -ent_temp(n,1).*ent_temp(n,2);
        end
    
%       AE entropy is calculated by sumation of all entries of column 3
        entropy(k) = sum(ent_temp(:,3)); % Keep Entropy value in column 1
        
    end
    
    AE_wave(i).wave_col_header(:,1) = {'HLT_index'};
    AE_wave(i).Ch2_wave(:,1) = HLT_ind;
    
    AE_wave(i).wave_col_header(:,2) = {'Optimum_BinWidth'};
    AE_wave(i).Ch2_wave(:,2) = bin_mod;
    
    AE_wave(i).wave_col_header(:,3) = {'Entropy_SOBW [nats]'};
    AE_wave(i).Ch2_wave(:,3) = entropy;
    
    disp(['Test ',num2str(i),' --> Done!'])
end

clear sd sk kur n bin m j skew sk_coeff kurt kur_coeff wave_kur_coeff wave_sk_coeff i k wave_loc wave_loc_nan WaveAddress
clear ent_temp j prob edge entropy wave_time waveform fl_waveform HLT_ind bin_mod t_loc fldr_loc w_loc BinWidth
       
%% save

name = ['Waveform_data_QIS_',date]; % To save the result of all tests

save(name,'-v7.3')
