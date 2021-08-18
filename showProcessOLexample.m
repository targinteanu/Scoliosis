%function [imgFiltered, splineObj, splineSample, splineObjBound, splineSampleBound] ...
%    = processOL(img, outln, boundX, boundY, zscl, xscl)

p = 42; 
imgtype = 'sag';

fp = ['C:\Users\Toren\Desktop\scoliosis\EOS patient ',num2str(p),'\selected DICOM\'];
img = dicomread([fp,imgtype]); ifo = dicominfo([fp,imgtype]);
img = double(img) / double(2^ifo.BitsStored - 1);
load([fp,'patient',num2str(p),' EOSoutline data']);
    zscl = ifo.PixelSpacing(1); 
    xscl = ifo.PixelSpacing(2);

if strcmp(imgtype,'sag')
    outln = SagOL;
end
if strcmp(imgtype,'cor')
    outln = CorOL;
end
        
% img is a double img [0, 1]; outln is a BW outline

% filtering 
im = imreconstruct(double(~outln), img); 
imgFiltered = img - im; 
imgFiltered = adapthisteq(imgFiltered);
imgFiltered = imgFiltered .* outln;
figure; imshow(imgFiltered, []);

%%
load([fp,'patient',num2str(p),' filtered data']);

sig = splSclSmp;
fs = 1/mean(diff(sig(:,3)));

filtobj = designfilt('lowpassiir', 'SampleRate',fs, ...
    'DesignMethod','butter', ...
    'PassbandFrequency',.004, 'StopbandFrequency',.005, ...
    'PassbandRipple',1, 'StopbandAttenuation',5);

figure; 
[wunit, w, pwr, pwr_dB] = plotFilterAndSignal(sig(:,1:2), fs);
pwr = pwr/(1e12);
stem(w, pwr(:,1), ':ob', 'LineWidth', 1.25); 
grid on; hold on; 
stem(w, pwr(:,2), '--sr', 'LineWidth', 1.25);
%xlabel(['Frequency, ' num2str(wunit) '/mm']);
xlim([0, 100]); 
xlabel('Frequency (1/m)'); ylabel('Power (m^4)'); title('Spine Magnitude Power Spectrum');
legend('Coronal Spine', 'Sagittal Spine');
set(gca, 'YScale', 'log');
fvtool(filtobj); 
%xlabel(['Frequency, ' num2str(wunit) '/mm']);
xlabel('Frequency (1/m)'); title('Filter Magnitude Response (dB)'); 
xlim([0, 100]); ylim([-125, 10]);

lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.25;
end

function [wunit, w, pwr, pwr_dB] = plotFilterAndSignal(x, fs)
X = abs(fftshift(fft(x)));
w = linspace(-fs/2, fs/2, size(X,1));

wmax = w(end-1); 
wordr = log10(wmax);
wunitlog = 3*floor(wordr/3); wunit = 10^wunitlog;

w = w/wunit;
pwr = X.^2; pwr_dB = 10*log10(pwr./sum(pwr));

end