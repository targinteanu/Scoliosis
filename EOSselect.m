function EOSselect

%% setting up the window and background
width = 1280; windowheight = 750; margin = 50; menuwidth = 100; 
ax = []; btn = []; txt = [];

[fn, fp] = uigetfile('*.mat', 'Select Scan Directory File'); 
dta = load([fp, '\', fn]); 
base_fp = dta.base_fp; img_fp = dta.img_fp; patient_list = dta.patient_list;
current_patient = 1;

% create the GUI but make it invisible for now
window = figure('Visible','off','Position',[360,500,width,windowheight]);

%% objects 
modselect = uicontrol('Style','popupmenu','Units','Pixels', ...
                'Position',[width/2 - menuwidth/2, windowheight - margin, menuwidth, margin/2],...
                'Callback',@setupPatient);
            
btnFwr = uicontrol('Style', 'pushbutton', 'String', 'Next >',...
        'Position', [width*.75 - menuwidth/2, windowheight - margin, menuwidth, margin/2],...
        'Callback', @funFwr);
btnPrv = uicontrol('Style', 'pushbutton', 'String', '< Prev',...
        'Position', [width*.25 - menuwidth/2, windowheight - margin, menuwidth, margin/2],...
        'Callback', @funPrv);
    
titletxt = uicontrol('Style','text',...
        'Position',[margin, windowheight - margin, menuwidth, margin/2]);
    
handleNewPatient(current_patient);
    
%% display the window
% Move the GUI to the center of the screen.
   movegui(window,'center')
   % Make the GUI visible.
   window.Visible='on';

%% functions 
    function handleNewPatient(patnum)
        current_patient = patnum;
        titletxt.String = ['Patient ' num2str(patient_list(current_patient))];
        getDicomTypes();
        if exist([base_fp, num2str(patient_list(current_patient)), img_fp], 'dir')
            showPatient();
        end
    end
    
    function getDicomTypes()
        dirfp = [base_fp, num2str(patient_list(current_patient))];
        [~,modselect.String] = parseDicomdir(fullfile(dirfp, 'DICOMDIR'));
    end

    function clearprevdata()
        delete(ax);
        delete(btn);
        delete(txt);
        ax = []; btn = []; txt = [];
    end
    
    function showPatient()
        % configure window for current patient
        imgs = dir([base_fp, num2str(patient_list(current_patient)), img_fp, '*.dcm']);
        N = length(imgs);
        imwidth = (width-(N+1)*margin)/N;
        clearprevdata();
        for i = 1:N
            ax(i) = axes('Units','Pixels','Position', ...
                [(i-1)*imwidth + i*margin, 2*margin, imwidth, windowheight-4*margin]);
            imgpath = fullfile(imgs(i).folder, imgs(i).name);
            img = dicomread(imgpath);
            imshow(img, []);
            btn(i) = uicontrol('Style','popupmenu','Units','Pixels', ...
                'Position',[(i-1)*imwidth + i*margin, margin, imwidth, margin/2],...
                'Callback',@assignImg, 'UserData',imgpath, ...
                'String', {'none', 'cor', 'sag'});
            txt(i) = uicontrol('Style','text',...
                'Position',[(i-1)*imwidth + i*margin, windowheight-2*margin, imwidth, margin],...
                'String',[num2str(size(img,1)), ' x ', num2str(size(img,2)), newline, imgs(i).name]);
        end
    end

%% callbacks 
    function funFwr(source, event)
        handleNewPatient( min(length(patient_list), current_patient+1) );
    end
    function funPrv(source, event)
        handleNewPatient( max(1, current_patient-1) );
    end

    function setupPatient(source,event)
        modality = source.String{source.Value};
        dirfp = [base_fp, num2str(patient_list(current_patient))];
        copySeriesFromDicomdir(dirfp, modality, [dirfp, img_fp], false);
        showPatient();
    end

    function assignImg(source,event)
        % copy file to new name if appropriate 
        val = source.String{source.Value};
        if ~strcmp(val, 'none')
            copyfile(source.UserData, [base_fp, num2str(patient_list(current_patient)), img_fp, val]);
        end
    end

end