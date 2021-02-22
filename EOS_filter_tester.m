clear;
[fn, fp] = uigetfile('*.mat', 'Select Scan Directory File'); 
load([fp, '\', fn]); 

patients_avail = [];
for p = patient_list 
    fnp = [base_fp, num2str(p), img_fp];
    if exist([fnp, 'patient',num2str(p),' EOSoutline data.mat'], 'file')
        patients_avail = [patients_avail, p];
        
        %{
        load([fnp, 'patient',num2str(p),' EOSoutline data.mat']);
        p
        size(splSclSmp)
        %}
    end
end

% Sag = x; Cor = y

filtfilt3d = @(filtobj, sig3d) [filtfilt(filtobj, sig3d(:,1)), filtfilt(filtobj, sig3d(:,2)), sig3d(:,3)];
plot3d = @(sig3d) plot3(sig3d(:,1), sig3d(:,2), sig3d(:,3), 'LineWidth', 1);

figure('Position', [50, 100, 1400, 800]);
ntot = length(patients_avail);
ncol = 5; nrow = ceil(ntot/ncol);
for i = 1:ntot
    subplot(nrow, ncol, i)
    
    p = patients_avail(i);
    fnp = [base_fp, num2str(p), img_fp, 'patient',num2str(p),' EOSoutline data.mat'];
    load(fnp);
    
    plot3(splSclSmp(:,1), splSclSmp(:,2), splSclSmp(:,3), ':');
    grid on; hold on;
    plot3(splSclSmp([1,end],1), splSclSmp([1,end],2), splSclSmp([1,end],3), 'ob', 'LineWidth', 2);
    plot3(splSclSmp(splSclBnd,1), splSclSmp(splSclBnd,2), splSclSmp(splSclBnd,3), 'or', 'LineWidth', 2);
    %wr = getWrithe(splSclSmp); title(num2str(wr));
    
    fcutoff = 1/100; %1/mm
    
    fs = 1/mean(diff(splSclSmp(:,3))); %1/mm
    filt_fir = designfilt('lowpassfir', 'CutoffFrequency',fcutoff, 'SampleRate',fs, ...
        'FilterOrder',100);
    filt_iir = designfilt('lowpassiir', 'HalfPowerFrequency',fcutoff, 'SampleRate',fs, ...
        'FilterOrder',5, 'DesignMethod','butter');
    filt_chb = designfilt('lowpassiir', 'PassbandFrequency',fcutoff, 'SampleRate',fs, ...
        'StopbandFrequency',1/80, 'DesignMethod','cheby2', ...
        'StopbandAttenuation',80, 'PassbandRipple',1);
    
    splFIR = filtfilt3d(filt_fir, splSclSmp); 
    splIIR = filtfilt3d(filt_iir, splSclSmp);
    splCHB = filtfilt3d(filt_chb, splSclSmp);
    plot3d(splFIR); plot3d(splIIR); plot3d(splCHB);
    
    splFIR = splFIR(splSclBnd(2):splSclBnd(1),:);
    splIIR = splIIR(splSclBnd(2):splSclBnd(1),:);
    splCHB = splCHB(splSclBnd(2):splSclBnd(1),:);
    
    xlim([min(splSclSmp(:,1)), max(splSclSmp(:,1))]);
    ylim([min(splSclSmp(:,2)), max(splSclSmp(:,2))]);
    view(2);
    
    splfilt = splCHB;
    %Wr = 0;
    Wr = getWrithe(splfilt); %wr2 = levittWrithe(splf);
    title([num2str(patients_avail(i)), '|', num2str(Wr)]);
    
    fnp2 = [base_fp, num2str(p), img_fp, 'patient',num2str(p),' filtered data.mat'];
    save(fnp2, 'splfilt', 'Wr');
    
    %legend('orig','FIR','butter','cheby');
    
    %{
    [~,corFilt,~,corFS,corPoly,corZ] = ...
        piecepolyfilter(fcutoff, splCorObj.p.breaks, splCorObj.p.coefs, splSclRng);
    [~,sagFilt,~,sagFS,sagPoly,sagZ] = ...
        piecepolyfilter(fcutoff, splSagObj.p.breaks, splSagObj.p.coefs, splSclRng);
    size(corZ)
    size(sagZ)
    plot3(sagFilt, corFilt, sagZ);
    plot3(sagFS, corFS, sagZ);
    %}
    
    i/ntot
end