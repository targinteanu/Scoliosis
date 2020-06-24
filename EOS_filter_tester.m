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

figure('Position', [50, 100, 1400, 800]);
ntot = length(patients_avail);
ncol = 4; nrow = ceil(ntot/ncol);
for i = 1:ntot
    subplot(nrow, ncol, i)
    
    p = patients_avail(i);
    fnp = [base_fp, num2str(p), img_fp, 'patient',num2str(p),' EOSoutline data.mat'];
    load(fnp);
    
    plot3(splSclSmp(:,1), splSclSmp(:,2), splSclSmp(:,3));
    grid on; hold on;
    %wr = getWrithe(splSclSmp); title(num2str(wr));
    
    fcutoff = 1/40; %1/mm
    
    [~,corFilt,~,corFS,corPoly,corZ] = ...
        piecepolyfilter(fcutoff, splCorObj.p.breaks, splCorObj.p.coefs, splSclRng);
    [~,sagFilt,~,sagFS,sagPoly,sagZ] = ...
        piecepolyfilter(fcutoff, splSagObj.p.breaks, splSagObj.p.coefs, splSclRng);
    size(corZ)
    size(sagZ)
    plot3(sagFilt, corFilt, sagZ);
    plot3(sagFS, corFS, sagZ);
end