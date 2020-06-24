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


figure('Position', [50, 100, 1400, 800]);
ntot = length(patients_avail);
ncol = 4; nrow = ceil(ntot/ncol);
for i = 1:ntot
    subplot(nrow, ncol, i)
    
    p = patients_avail(i);
    fnp = [base_fp, num2str(p), img_fp, 'patient',num2str(p),' EOSoutline data.mat'];
    load(fnp);
    
    plot3(splSclSmp(:,1), splSclSmp(:,2), splSclSmp(:,3));
    grid on;
    wr = getWrithe(splSclSmp); title(num2str(wr));
end