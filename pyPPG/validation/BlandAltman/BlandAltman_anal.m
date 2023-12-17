clear all
close all

% Define folder
input_folder='2023_12_16_22_4';

% Set analysis: 
%   detection vs reference (1), 
%   reference_1 vs reference_2 (0)
isdet=0;

% Run Bland Altman 
run_BA(input_folder,isdet)


function run_BA(input_folder,isdet)
    if isdet
        output_folder=strcat(input_folder,'_DET-REF12');
    else
        output_folder=strcat(input_folder,'_REF1-REF2');
    end
    
    % Get fiducial data
    [fps1,fps2]=get_fps(input_folder, isdet);
    names=fps1.Properties.VariableNames;
    
    % % plot results 
    plot_BA(fps1,fps2,isdet,names,output_folder);
    
    % % crop BA plots
    crop_BA(names,output_folder);
    
    % merge BA plots
    merge_BA(names,output_folder);
end

function merge_BA(names,output_folder)
   
    imageFiles = cell(1, length(names));
    for i = 1:length(names)
        tmp_fp=names(i);
            
        % Full file path
        Filename=strcat(tmp_fp,'.jpeg');
        croppedFolder=strcat('figures\',output_folder,'\cropped');
        FilePath = string(fullfile(croppedFolder, Filename));

        imageFiles{i} = FilePath;
    end
    
    % Read all the images
    images = cell(1, length(names));
    for i = 1:length(names)
        imData = imread(imageFiles{i});

        % Define the region to crop (in the format [xmin, ymin, width, height])
        crop_region = [720, 180, 4700, 3900];
        cropped_image=imcrop(imData, crop_region);
        imshow(cropped_image);

        images{i} = cropped_image;

    end
    
    % Create a montage with specified parameters
    montage(images, 'Size', [5, 3], 'ThumbnailSize', [NaN NaN], 'BorderSize', [0, 0]);
    
    % Save the composite image with high resolution
    outputFile = ['figures\',output_folder,'\merged_BA.jpeg'];
    print(gcf, outputFile, '-dpng', '-r4400'); % Adjust resolution as needed

end

function [fps1,fps2]=get_fps(output_folder,isdet)
    % Load data
    if isdet
        load(strcat('..\results\',output_folder,'\MG-pyPPG.mat'))
        load(strcat('..\results\',output_folder,'\PC-pyPPG.mat'))
    else 
        load(strcat('..\results\',output_folder,'\MG-PC.mat'))
    end
    
    % define annotated fiducial points
    mg_fps=struct2table(MG_fps);
    mg_fps.('dp')=[];
    
    pc_fps=struct2table(PC_fps);
    pc_fps.('dp')=[];

    if isdet
        tmp_fps=(mg_fps{:,:} + pc_fps{:,:}) / 2;
        ref_fps = array2table(tmp_fps, 'VariableNames', mg_fps.Properties.VariableNames);
        
        ref_fps_on=ref_fps.on;
        tmp_fps=(ref_fps{:,:} - ref_fps.on)+1;
        ref_fps = array2table(tmp_fps, 'VariableNames', mg_fps.Properties.VariableNames);
        
        % define detected fiducial points
        det_fps=struct2table(pyPPG_fps);
        det_fps.('dp')=[];
        
        tmp_fps=(det_fps{:,:} - ref_fps_on)+1;
        det_fps = array2table(tmp_fps, 'VariableNames', det_fps.Properties.VariableNames);
        
        fps1=det_fps;
        fps2=ref_fps;
    else
        fps1=mg_fps;
        fps2=pc_fps;
    end
end

function crop_BA(names,folder_name)
    % Check if the folder exists; if not, create it
    output_folder=strcat('figures\',folder_name,'\cropped')
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    for i=1:length(names)
        tmp_fp=names(i);
    
        % Full file path
        Filename=strcat(tmp_fp,'.jpeg');
        originFolder=strcat('figures\',folder_name,'\origin');
        FilePath = string(fullfile(originFolder, Filename));
    
        % Read the image file using imread
        imData = imread(FilePath);
    
        % Define the region to crop (in the format [xmin, ymin, width, height])
        crop_region = [1200, 200, 900, 720];
    
        % Crop the specified region
        cropped_image = imcrop(imData, crop_region);
        imshow(cropped_image)

        % Calculate the position dynamically based on the size of the cropped image
        imageSize = size(cropped_image);
        titlePosition = [imageSize(2)/2, imageSize(1)*0.025]; % Middle of the width, 2.5% from the top

        % Add the title using the text function
        text(titlePosition(1), titlePosition(2), tmp_fp, Color='black', FontSize=20,FontWeight='bold');     

        % Save the figure with high quality and added bold text
        outputFilePath=string(strcat(output_folder,'\',tmp_fp,'.jpeg'));
        print(gcf, outputFilePath, '-dpng', '-r1000', '-opengl');

    end 
end

function plot_BA(fps1,fps2,isdet,names,folder_name)
    for i = 1:length(names)
        FID_f=i;
        states = names(FID_f);
        
        %% Example 1
        % data1
        data1(1:height(fps1),:,1:length(FID_f)) =  NaN;
        data1(1:height(fps1),:,1:length(FID_f)) =  table2array(fps1(:,FID_f));
        
        % data2
        data2(1:height(fps2),:,1:length(FID_f)) =  NaN;
        data2(1:height(fps2),:,1:length(FID_f)) =  table2array(fps2(:,FID_f));
        
        % BA plot paramters
        tit = strcat(states); % figure title
        gnames = {states}; % names of groups in data {dimension 1 and 2}
        
        if isdet
            label = {'Reference','pyPPG','ms'}; % Names of data sets
        else
            label = {'Annotator 1','Annotator 2','ms'}; % Names of data sets
        end
        corrinfo = {'n','RMSE','r2','eq'};  % stats to display of correlation scatter plot
        BAinfo = {'LOA(%)','n','RMSE'};     % stats to display on Bland-ALtman plot
        limits = 'auto';                    %'tight';%'tight';%'auto'; % how to set the axes limits
        colors = 'rbgmcky';                 % character codes
        colors = colors(1:length(states));
        
        % Generate figure with symbols
        [cr, fig, statsStruct] = BlandAltman_func(data1, data2,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors,'markerSize',4,'showFitCI',' on');
        
        % Check if the folder exists; if not, create it
        output_folder=strcat('figures\',folder_name,'\origin')
        if ~exist(output_folder, 'dir')
            mkdir(output_folder);
        end
        
        % Full file path
        outputFilename=strcat(states,'.jpeg')
        outputFilePath = string(fullfile(output_folder,outputFilename));

        % Save figure
        saveas(gcf, outputFilePath);

        % Close the figure
        close(gcf);

    end
end
