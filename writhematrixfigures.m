%% adolescent 

patient = [9 14 21]; 
vertnames = [arrayfun(@(i) ['T' num2str(i)], 1:12, 'UniformOutput', 0), ...
    arrayfun(@(i) ['L' num2str(i)], 1:5, 'UniformOutput', 0)];
segmentnames = arrayfun(@(i) [vertnames{i-1} '-' vertnames{i}], 2:length(vertnames), ...
    'UniformOutput', 0);

figure; 

for I = 1:length(patient)
    subplot(1,3,I); 
    load('spines_XYZ.mat')
    idx = patient(I);

    xyz = spinesXYZ{idx, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    
    cm = [x; y; z]';
    writheM; 
    
    M = M(2:end, 2:end); 
    heatmap(segmentnames, segmentnames, M);
    title(['Patient ' num2str(idx)]);
    xlabel('Segment 1'); ylabel('Segment 2'); 
end

%% adult 

patient = [6, 28, 1, 3, 12, 15, 22]; NN = length(patient);
Range = [3:7; 1:5; 1:5; 3:7; 3:7; 3:7; 3:7]; 
vertnames = arrayfun(@(i) ['L' num2str(i)], 1:5, 'UniformOutput', 0);
segmentnames = arrayfun(@(i) [vertnames{i-1} '-' vertnames{i}], 2:length(vertnames), ...
    'UniformOutput', 0);

%figure; 

for I = 1:length(patient)
    %subplot(1,NN,i); 
    figure('Position', [100, 100, 400, 300]);
    idx = patient(I);
    load(['points_pat' num2str(idx) '.mat']);

    cm = cm(Range(I,:),:);
    writheM; 
    
    M = M(2:end, 2:end); 
    heatmap(segmentnames, segmentnames, M);
    title(['Patient ' num2str(I)]);
    xlabel('Segment 1'); ylabel('Segment 2'); 
end