[fn, fp] = uigetfile('*.mat', 'Select Scan Directory File'); 
load([fp, '\', fn]);

for p = patient_list
    fnp = [base_fp, num2str(p), img_fp];
    %{
    src_fp = [fnp, 'patient',num2str(p),' EOSoutline data.mat'];
    if exist(src_fp, 'file')
        dst_fp = ['C:\Users\Toren\Desktop\newfolder\','patient',num2str(p),' EOSoutline data.mat'];
        copyfile(src_fp, dst_fp);
    end
    %}
    %%{
    src_fp = [fnp, 'patient',num2str(p),' filtered data.mat'];
    if exist(src_fp, 'file')
        dst_fp = ['C:\Users\Toren\Desktop\newfolder\','patient',num2str(p),' filtered data.mat'];
        copyfile(src_fp, dst_fp);
    end
    %}
end