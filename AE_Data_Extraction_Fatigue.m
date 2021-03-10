%%%%%%% This code is designed to extract valid signals from all the available signals and prepare them for waveform analysis %%%%%%%%%%%%

clc;
clear all;

lib = load('QIS_Test_lib','Test');
lib = lib.Test;

%% Select the tests
% Select the test names
test_name = ['QIS_08';'QIS_09';'QIS_10';'QIS_16';...
             'QIS_17';'QIS_20';'QIS_21';'QIS_22';
             'QIS_23';'QIS_24';'QIS_25';'QIS_26';...
             'QIS_27';'QIS_28']; 
o=1;
for j =1:size(test_name,1)
    fprintf('Test name is %s \n', test_name(j,:))
    for i=1:size(lib,2) % Find the index of it
        if contains(lib(i).name, test_name(j,:)) == 1
            t_ind(o,:) = i;

            fprintf('\t Test name is %s and index is %d \n', lib(i).name, o)

            o = o+1;
        end
    end
end
clear i o 

%% Load AE data

for i= 1:length(t_ind) % Extract all the information needed
    AE_loc = [lib(t_ind(i)).AE_loc,'\'];
    
    AE_data(i).test_name = lib(t_ind(i)).name;
    
    t_name = struct2table(dir(AE_loc)).name; % Find all the file names in the address
    
    t_name(ismember(t_name,'.')) = []; t_name(ismember(t_name,'..')) = []; % Remove the . and .. in the test name list
    
    txt_files = t_name(contains(t_name,'.TXT'));
    test_fname = txt_files(contains(txt_files,'Test'));
    para_fname = txt_files(contains(txt_files,'Para'));
    sum_fname = txt_files(contains(txt_files,'Sum'));
    
    for j=1:length(test_fname)
        th = importdata([AE_loc,char(test_fname(j))],' ',8);
        temp_hit = th.data;
        
        tp = importdata([AE_loc,char(para_fname(j))],' ',8);
        temp_para = tp.data;
        
        ts = importdata([AE_loc,char(sum_fname(j))],' ',4);
        temp_sum = ts.data(:,2).';
        
        if j==1
            t_hit = temp_hit;
            t_para = temp_para;
            t_sum = temp_sum;
        else
            t_hit = [t_hit;temp_hit];
            t_para = [t_para;temp_para];
            t_sum = [t_sum;temp_sum];
        end
    end
    
    AE_data(i).hit = t_hit;
    AE_data(i).hit_colheaders = th.colheaders;
    
    AE_data(i).para = t_para;
    AE_data(i).para_colheaders = tp.colheaders;
    
    AE_data(i).sum = t_sum; % First column shows hits in Channel 1, second columns Channel 2
    AE_data(i).sum_colheaders = ts.colheaders;

    AE_data(i).DeltaT = lib(t_ind(i)).DeltaT;
    
end

clear i j th temp_hit tp temp_para t_hit t_para txt_files test_fname sum_fname para_fname t_ind AE_loc t_name ts temp_sum t_sum

%% Find initial and end index number of loading blocks

for i= 1:size(AE_data,2)
    
    block_inds = find(AE_data(i).para(:,3)<0.5); % Initial tensile loadings in blocks
    block_init = block_inds([0;diff(block_inds)]>1); % Find the borders of loading blocks
    
    figure()
    plot(AE_data(i).para(:,2),AE_data(i).para(:,3),'DisplayName','Para_raw')
    hold on
    scatter(AE_data(i).para(block_inds,2),AE_data(i).para(block_inds,3),'DisplayName','block_inds')
    hold on
    scatter(AE_data(i).para(block_init,2),AE_data(i).para(block_init,3),'x','DisplayName','block_init')
    hold on
    legend('Interpreter','none')

    % The LENGTH OF THIS VARIABLE is equal to number of loading blocks and
    % the first and second column hold the index of start and end of the
    % block
    block_inits = [[1;block_init(1:end-1)],block_init]; 

    for j=1:size(block_inits,1)
        
        block_data = AE_data(i).para(block_inits(j,1):block_inits(j,2),:); % Extract data of the block
        
        plot(block_data(:,2),block_data(:,3),'DisplayName',['Block_data ',num2str(j)])
        hold on
        
        block_avg_load = mean([max(block_data(:,3)),0.1*max(block_data(:,3))]); % Mean stress in a block
            
        first_ind = find(block_data(:,3)>block_avg_load*1.03,1); % First fatigue cycle in a block
        last_ind = find(block_data(:,3)>block_avg_load*1.03,1,'last'); % Last fatigue cycle in a block
        
        block_inits(j,3:4) = [first_ind,last_ind]; % Keep the begining and end index of valid data in a block in columns 3 and 4
        
        scatter(block_data(block_inits(j,3):block_inits(j,4),2),block_data(block_inits(j,3):block_inits(j,4),3),'DisplayName',['Valid data block ',num2str(j)])
        hold on
        
    end
    
    hold off
    
    AE_data(i).block_inits = block_inits;
    AE_data(i).block_inits_colheaders = {'block_init_ind','block_end_ind','block_data_init_ind','block_data_end_ind'};
        
end

clear block_avg_load block_data block_inds block_init block_inits cycle filter first_ind i j last_ind t_cycle t_filter filter

%% Extract cycle number from Para file

for i= 1:size(AE_data,2)
    
    filter = [];
    cycle = [];

    for j=1:size(AE_data(i).block_inits,1) %length(block_init)
        
        block_data = AE_data(i).para(AE_data(i).block_inits(j,1):AE_data(i).block_inits(j,2),:); % Read block index values from the structure filed
        
        first_ind = AE_data(i).block_inits(j,3); % Find the initial index of valid data in the block from the field in structure
        last_ind = AE_data(i).block_inits(j,4); % Find the final index of valid data in the block from the field in structure
        
        t_filter = false(length(block_data),1); % A filter to extract useful data of the Para file
        t_filter(first_ind:last_ind) = 1;
        
        t_cycle = zeros(length(block_data),1); % Record cycle numbers 
        t_cycle(1:first_ind-1) = 0;
        t_cycle(first_ind:last_ind) = fix((block_data(first_ind:last_ind,2) - block_data(first_ind,2)).*10); % Assign cycle numbers based on f=10Hz
        t_cycle(last_ind+1:end) = 0; 
        
        filter = [filter;t_filter];

        if isempty(cycle)
            t_cycle(t_filter) = t_cycle(t_filter);
        else
            t_cycle(t_filter) = t_cycle(t_filter) + max(cycle);
        end
        cycle = [cycle;t_cycle];
        
    end
    
    if length(AE_data(i).para) > length(cycle)
        cycle(AE_data(i).block_inits(j,2)+1:length(AE_data(i).para)) = 0;
        filter(AE_data(i).block_inits(j,2)+1:length(AE_data(i).para)) = 0;
    end
    
    filter = logical(filter);
    
    AE_data(i).para_colheaders(1,5) = {'Cycle'};
    AE_data(i).para(:,5) = cycle;
    AE_data(i).para_filter = filter;
        
end

clear block_avg_load block_data block_inds block_init block_inits cycle filter first_ind i j last_ind t_cycle t_filter filter

%% Plot for cycle number adding to Para and Para filtering

i=14; % Row number in AE_data

figure()
plot(AE_data(i).para(:,2),AE_data(i).para(:,3),'DisplayName','Para_raw')
hold on
plot(AE_data(i).para(AE_data(i).para_filter,2),AE_data(i).para(AE_data(i).para_filter,3),'DisplayName','Filtered_data')
hold off

legend('Interpreter','none')
xlabel('Time [s]')
ylabel('Load [kN]')

figure()
plot(AE_data(i).para(:,5),AE_data(i).para(:,3),'DisplayName','Filtered_data')

legend('Interpreter','none')
xlabel('Cycle')
ylabel('Load [kN]')

clear i

%% Check if number of hits recorded matches number of waveforms

ch1cnt_hit = sum(AE_data.hit(:,5)==1); ch1cnt_sum = sum(AE_data.sum(:,1));
ch2cnt_hit = sum(AE_data.hit(:,5)==2); ch2cnt_sum = sum(AE_data.sum(:,2));

fprintf('\n------- \n CH1: \n\t Number of hit data available for Ch1 is %d number of waveforms available for the channel is %d \n\t Difference is %d \n',...
        ch1cnt_hit, ch1cnt_sum, abs(ch1cnt_hit - ch1cnt_sum))
fprintf('\n------- \n CH2: \n\t Number of hit data available for Ch1 is %d number of waveforms available for the channel is %d \n\t Difference is %d \n',...
        ch2cnt_hit, ch2cnt_sum, abs(ch2cnt_hit - ch2cnt_sum))

clear ch1cnt_hit ch1cnt_sum ch2cnt_hit ch2cnt_sum
%% Separate channels and add order
    
for i= 14 %1:size(AE_data,2)

    AE_data(i).hit_colheaders(:,18) = {'Order'};
    
    AE_data(i).Ch1_data = AE_data(i).hit(AE_data(i).hit(:,5)==1,:);
    for k=1:length(AE_data(i).Ch1_data) %add order to column 18
        AE_data(i).Ch1_data(k,18) = k;
    end
    
    AE_data(i).Ch2_data = AE_data(i).hit(AE_data(i).hit(:,5)==2,:);
    for k=1:length(AE_data(i).Ch2_data) %add order to column 18
        AE_data(i).Ch2_data(k,18) = k;
    end
end

clear i k 

%% Apply delta T filter
% We do this loading block by loading block.

for i= 1:size(AE_data,2)
    
    deltaT = AE_data(i).DeltaT; % Delta T value recorded at the begining of test
    
    Ch1_dtf = []; % DeltaT filtered data in Ch1
    Ch2_dtf = []; % DeltaT filtered data in Ch2
    
    for j=1:size(AE_data(i).block_inits,1)
        
        block_data = AE_data(i).para(AE_data(i).block_inits(j,1):AE_data(i).block_inits(j,2),:); % Read block index values from the structure filed
        
        first_ind = AE_data(i).block_inits(j,3); % Find the initial index of valid data in the block from the field in structure
        last_ind = AE_data(i).block_inits(j,4); % Find the final index of valid data in the block from the field in structure
        
        init_time = block_data(first_ind,2); % time value for the block start
        end_time = block_data(last_ind,2); % time value for the block end
        
        block_hit_1 = AE_data(i).Ch1_data(AE_data(i).Ch1_data(:,2)>=init_time,:); % Ch1 hit data of the block
        block_hit_1 = block_hit_1(block_hit_1(:,2) <= end_time,:); % Ch1 hit data of the block
        
        block_hit_2 = AE_data(i).Ch2_data(AE_data(i).Ch2_data(:,2)>=init_time,:); % Ch1 hit data of the block
        block_hit_2 = block_hit_2(block_hit_2(:,2) <= end_time,:); % Ch1 hit data of the block
        
        % deltaT filter in Channel 1
        ch1kd = [];
        for k=1:size(block_hit_2,1)
            temp_ind = abs(block_hit_1(:,2) - block_hit_2(k,2))<deltaT;
            ch1kd = [ch1kd;block_hit_1(temp_ind,:)];
        end
        
        % deltaT filter in Channel 2
        ch2kd = [];
        for k=1:size(ch1kd,1)
            temp_ind = abs(block_hit_2(:,2) - ch1kd(k,2))<deltaT;
            ch2kd = [ch2kd;block_hit_2(temp_ind,:)];
        end
        
        Ch1_dtf = [Ch1_dtf;ch1kd]; % Keep all the Ch1 dtf data 
        Ch2_dtf = [Ch2_dtf;ch2kd]; % Keep all the Ch2 dtf data 
        
        fprintf('\n Progress ... %.1f percent \n',j*100/size(AE_data(i).block_inits,1))
        
    end
    
    fprintf('----------- \nTest %s done!\n',AE_data(i).test_name)
    
    AE_data(i).Ch1_dtf = Ch1_dtf; % Add the dtf data to structure
    AE_data(i).Ch2_dtf = Ch2_dtf; % Add the dtf data to structure
    
end

% Check to make sure dtf data are the same size
fprintf('\n ------ \n\t # of Ch1 data points: %d \n\t # of Ch2 data points: %d\n',length(Ch1_dtf), length(Ch2_dtf))

clear i j k ch1kd ch2kd Ch1_dtf Ch2_dtf block_data first_ind last_ind 
clear init_time end_time block_hit_1 block_hit_2 temp_ind deltaT

%% Add cycle number to dtf data

% Find cycle numbers for 1 channel and add the same to the other. Since
% delta T filtereing is already performed, this is correct and saves time.

for i= 1:size(AE_data,2)
    
    cycle = [];
    
    for j=1:size(AE_data(i).block_inits,1)
        % block_data includes para info of the block
        block_data = AE_data(i).para(AE_data(i).block_inits(j,1):AE_data(i).block_inits(j,2),:); % Read block index values from the structure filed
        
        first_ind = AE_data(i).block_inits(j,3); % Find the initial index of valid data in the block from the field in structure
        last_ind = AE_data(i).block_inits(j,4); % Find the final index of valid data in the block from the field in structure
        
        init_time = block_data(first_ind,2); % time value for the block start
        end_time = block_data(last_ind,2); % time value for the block end
        
        block_hit_1 = AE_data(i).Ch1_dtf(AE_data(i).Ch1_dtf(:,2)>=init_time,:); % Ch1 hit data of the block
        block_hit_1 = block_hit_1(block_hit_1(:,2) <= end_time,:); % Ch1 hit data of the block
        
        % Find cycle number based on minimum difference between hit time
        % and para time. Add cycle based on para cycle value in column 5. 
        t_cycle = [];
        for k=1:size(block_hit_1,1)
           temp_ind = abs(block_hit_1(k,2) - block_data(:,2)) == min(abs(block_hit_1(k,2) - block_data(:,2)),1);
           
           if sum(temp_ind)>1 % This is to prevent more than 1 cycle number for a hit
               temp_ind = find(temp_ind==1,1,'first');
           end
           
           t_cycle = [t_cycle;block_data(temp_ind,5)];
        end
        
        cycle = [cycle;t_cycle]; % Keep the cycle numbers
        
        fprintf('\n Progress ... %.1f percent \n',j*100/size(AE_data(i).block_inits,1))
        
    end
    
    fprintf('----------- \nTest %s done!\n',AE_data(i).test_name)
    
    AE_data(i).hit_colheaders(1,19) = {'Cycle'};
    
    AE_data(i).Ch1_dtf(:,19) = cycle; % Add cycle numbers to dtf data
    AE_data(i).Ch2_dtf(:,19) = cycle; % Add cycle numbers to dtf data
    
end

clear i j k t_cycle cycle block_data first_ind last_ind 
clear init_time end_time block_hit_1 temp_ind

%% Now remove the hit data with cycle number = 0
% When a cycle number 0 is assigned to a hit, it means the hit was recorded
% in tensile loading not fatigue loading

for i= 1:size(AE_data,2)
    
   AE_data(i).Ch1_dtf(AE_data(i).Ch1_dtf(:,19) == 0,:) = [];
   AE_data(i).Ch2_dtf(AE_data(i).Ch2_dtf(:,19) == 0,:) = [];
   
end

clear i

%% Now remove hits that have load values less than 1 [kN]

for i= 1:size(AE_data,2)
    
   AE_data(i).Ch1_dtf(AE_data(i).Ch1_dtf(:,3) < 1,:) = [];
   AE_data(i).Ch2_dtf(AE_data(i).Ch2_dtf(:,3) < 1,:) = [];
   
end

clear i

%% Statistics of valid hits 

for i= 1:size(AE_data,2)
    fprintf('----------------\n Test %s \n',AE_data(i).test_name(1:6))
    
    ch1 = length(AE_data(i).Ch1_dtf)/length(AE_data(i).Ch1_data);
    ch2 = length(AE_data(i).Ch2_dtf)/length(AE_data(i).Ch2_data);
    
    fprintf('\tRatio of the valid data in Ch1 is %.1f percent \n\tRatio of the valid data in Ch2 is %.1f percent\n ', ch1*100, ch2*100)
end

clear i ch1 ch2
%% Plot

x_axis = 2;
y_axis = 3;

for i= 1:size(AE_data,2)
    
    figure()
    yyaxis left
    plot(AE_data(i).para(AE_data(i).para_filter,2),AE_data(i).para(AE_data(i).para_filter,3),'DisplayName','Loading profile')
    hold on
    
    ylim([0 inf])
    ylabel('Load [kN]')
    
    yyaxis right
    scatter(AE_data(i).Ch1_dtf(:,x_axis),AE_data(i).Ch1_dtf(:,y_axis),'DisplayName',[AE_data(i).test_name(1:6), ' Ch1'])
    hold on
    scatter(AE_data(i).Ch2_dtf(:,x_axis),AE_data(i).Ch2_dtf(:,y_axis),'x','DisplayName',[AE_data(i).test_name(1:6), ' Ch2'])
    hold on
    
    ylim([0 inf])
    
    xlabel(AE_data(i).hit_colheaders(x_axis))
    ylabel(AE_data(i).hit_colheaders(y_axis))
    
    legend('location','nw','Interpreter','None')
    
end

clear x_axis y_axis i 

%% Save

% name = 'AE_Data_Extraction_Fatigue';
% name = [name,'_',date];

save(name,'-v7.3')
