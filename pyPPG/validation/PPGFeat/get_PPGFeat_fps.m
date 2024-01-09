% The MIT License (MIT)
% 
% Copyright (c) 2023 Marton A. GODA, PhD.
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function get_PPGFeat_fps(results_date)

    % Load test referece data
    input_folder=strcat('..',filesep,'results',filesep,results_date,filesep,'MG_PC',filesep);
    load(strcat(input_folder,'MG_PC.mat'))
    input=strcat('..',filesep,'PPG-BP_annot',filesep,'PPG-BP_ref1.mat');
    load(input)
    
    % Create struct for fiducial points 
    fieldNames=fieldnames(MG_fps);
    PPGFeat_fps = struct();
    for i = 1:length(fieldNames)
        PPGFeat_fps.(fieldNames{i}) = [];
        MG_PPGFeat_diff.(fieldNames{i}) = [];
        PC_PPGFeat_diff.(fieldNames{i}) = [];
    end
    
    %%
    % File name of pulse wave start points
    filename='start_sig.csv';
    
    % Use xlsread to read the data into a matrix
    win_start = table2array(readtable(filename));
    
    %%
    % File name of fiducial points
    filename='ID_min1_min2.xlsx';
    
    % Specify the sheet name and range if needed
    sheet = 'Sheet1';  % Change to the actual sheet name
    range = 'B1:B219';  % Change to the actual range
    
    % Use xlsread to read the data into a matrix
    sig_start = xlsread(filename, sheet, range);
    
    %%
    % File name of fiducial points
    filename='PPG_features.xlsx';
    
    % Specify the sheet name and range if needed
    sheet = 'Sheet1';  
    range = 'P1:AD219'; 
    
    % Use xlsread to read the data into a matrix
    PPG_fps = xlsread(filename, sheet, range);
    
    %%
    % Get absolute location
    PPG_init=PPG_fps(:,1);
    PPG_fps_new=PPG_fps+win_start+sig_start-PPG_init;
    
    % Detection of other pulse wave
    tmp_inds1=find(PPG_fps_new(:,2)>[MG_fps.off]');
    pw_len=[MG_fps.off]'-[MG_fps.on]';
    for ind=tmp_inds1
        PPG_fps_new(ind,:)=PPG_fps(ind,:)+win_start(ind)+sig_start(ind)-PPG_init(ind)-pw_len(ind);
    end
    
    tmp_inds2=find(PPG_fps_new(:,2)<[MG_fps.on]');
    pw_len=[MG_fps.off]'-[MG_fps.on]';
    for ind=tmp_inds2
        PPG_fps_new(ind,2:end)=PPG_fps(ind,2:end)+win_start(ind)+sig_start(ind)-PPG_init(ind)+pw_len(ind);
    end
    
    % Update fiducial struct
    PPGFeat_fps.on=PPG_fps_new(:,1);
    PPGFeat_fps.sp=PPG_fps_new(:,2);
    PPGFeat_fps.dn=PPG_fps_new(:,3);
    PPGFeat_fps.dp=PPG_fps_new(:,4);
    PPGFeat_fps.off=PPG_fps_new(:,5);
    PPGFeat_fps.u=PPG_fps_new(:,6);
    PPGFeat_fps.v=PPG_fps_new(:,8);
    PPGFeat_fps.w=PPG_fps_new(:,9);
    
    for i = 9:length(fieldNames)-2
        PPGFeat_fps.(fieldNames{i}) = PPG_fps_new(:,i+1);
    end
    
    PPGFeat_fps.p1=nan(219,1);
    PPGFeat_fps.p2=nan(219,1);
    
    PPGFeat_fps=struct2table(PPGFeat_fps);
    PPGFeat_fps_copy=PPGFeat_fps;
    
    %%
    for tmp_name=["MG","PC"]
        % Calculate the difference
        ref_fps=eval([char(tmp_name),'_fps']);
        PPGFeat_fps=PPGFeat_fps_copy;
        diff_fps=table2array(struct2table(ref_fps)) - table2array(PPGFeat_fps);
        
        %%
        tmp_fname=strcat(char(tmp_name),'_PPGFeat');
        output_folder=char(strcat('..',filesep,'results',filesep,results_date,filesep,'PPGFeat',filesep,tmp_fname,filesep));
        if exist(output_folder, 'dir') ~= 7
            abs_loc=pwd;
            rel_loc=strcat(filesep,'PPGFeat');
            tmp_i=strfind(abs_loc,rel_loc);
            output_folder=strcat(abs_loc(1:tmp_i),output_folder(4:end));
            mkdir(output_folder);
        end
        
        % Save data
        PPGFeat_fps=table2struct(PPGFeat_fps);
    
        for i = 1:length(fieldNames)
            eval([char(tmp_name),'_PPGFeat_diff.(fieldNames{i})=diff_fps(:,i);']);
        end
    
        eval([char(tmp_name),'_PPGFeat_diff=table2struct(struct2table(',char(tmp_name),'_PPGFeat_diff','));'])
    
        mat_name=strcat(output_folder,tmp_fname,'.mat');
        ref_name=strcat(char(tmp_name),'_fps');
        diff_name=strcat(char(tmp_name),'_PPGFeat_diff');
        save(mat_name,ref_name, 'PPGFeat_fps',diff_name);
    end
    
    disp(strcat('Output folder: ..',filesep,'results',filesep,results_date,filesep,'PPGFeat'));
    disp('PPGFeat has been finished!')
end
