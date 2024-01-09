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


function get_PA_fps(results_date)
    input_folder=strcat('..',filesep,'results',filesep,results_date,filesep,'MG_PC',filesep);
    load(strcat(input_folder,'MG_PC.mat'))
    input=strcat('..',filesep,'PPG-BP_annot',filesep,'PPG-BP_ref1.mat');
    load(input)
    
    fieldNames=fieldnames(MG_fps);
    
    PulseAnal_fps = struct();
    for i = 1:length(fieldNames)
        PulseAnal_fps.(fieldNames{i}) = [];
        MG_PulseAnal_diff.(fieldNames{i}) = [];
        PC_PulseAnal_diff.(fieldNames{i}) = [];
    end
    
    for tmp_name=["MG","PC"]
        ref_fps=eval([char(tmp_name),'_fps']);
    
        for i=1:length(ppg_data)
            %% default setup
            fs=ppg_data(i).fs;
            [ppg,vpg,apg,jpg]=get_signals(ppg_data(i).sig, fs);
            S.v = ppg; 
            S.fs = fs;
            S.d1 = vpg;
            S.d2 = apg;
            S.d3 = jpg;
            options.do_plot = 0;
            options.exclude_low_quality_data = 0;
            options.do_filter = 0;
            options.do_quality = 0;
            
            pulses.onsets=[ref_fps(i).on,ref_fps(i).off]+1;
            pulses.peaks=ref_fps(i).sp+1;
            
            %% fiducial points by PulseAnalyse2
            try
                [pw_inds, fid_pts, pulses, sigs] = PulseAnalyse2(S,options,pulses);
                fid_pts.ind.sp=fid_pts.ind.s;
                fid_pts.ind.on=fid_pts.ind.os;
                fid_pts.ind.off=NaN;
                fid_pts.ind.dn=fid_pts.ind.dic;
                fid_pts.ind.dp=NaN;
                fid_pts.ind.u=fid_pts.ind.ms;
                fid_pts.ind.v=NaN;
                fid_pts.ind.w=NaN;
                fid_pts.ind.p1=fid_pts.ind.p1pk;
                fid_pts.ind.p2=fid_pts.ind.p2pk;
        
                for j = 1:length(fieldNames)
                    if isempty(fid_pts.ind.(fieldNames{j}))
                        PulseAnal_fps(i).(fieldNames{j}) = NaN;
                    else
                        PulseAnal_fps(i).(fieldNames{j}) = fid_pts.ind.(fieldNames{j});
                    end
                end
            catch
                for j = 1:length(fieldNames)
                    PulseAnal_fps(i).(fieldNames{j}) = NaN;
                end
            end
        end
    
        diff_fps=eval([char(tmp_name),'_PulseAnal_diff']);
        for j = 1:length(fieldNames)
            tmp_ref_fps=[ref_fps.(fieldNames{j})];
            numElements=length(tmp_ref_fps);
            [diff_fps(1:numElements).(fieldNames{j})] = deal([]);
            for k = 1:numElements
                if isempty(PulseAnal_fps(k).(fieldNames{j})) || isempty(tmp_ref_fps(k))
                    diff_fps(k).(fieldNames{j})=NaN;
                else
                    tmp_diff=tmp_ref_fps(k)-PulseAnal_fps(k).(fieldNames{j});
                    diff_fps(k).(fieldNames{j}) = tmp_diff;
                end
            end
        end
        eval([char(tmp_name),'_PulseAnal_diff=diff_fps;']);
        
        tmp_fname=strcat(char(tmp_name),'_PulseAnal');
        
        output_folder=char(strcat('..',filesep,'results',filesep,results_date,filesep,'PulseAnal',filesep,tmp_fname,filesep));
        if exist(output_folder, 'dir') ~= 7
            abs_loc=pwd;
            rel_loc=strcat(filesep,'PulseAnal');
            tmp_i=strfind(abs_loc,rel_loc);
            output_folder=strcat(abs_loc(1:tmp_i),output_folder(4:end));
            mkdir(output_folder);
        end
    
        mat_name=strcat(output_folder,tmp_fname,'.mat');
        ref_name=strcat(char(tmp_name),'_fps');
        diff_name=strcat(char(tmp_name),'_PulseAnal_diff');
        save(mat_name,ref_name, 'PulseAnal_fps',diff_name);
    
        disp(strcat('Output folder: ..',filesep,'results',filesep,results_date,filesep,'PulseAnal'));
        disp('PulseAnalysis has been finished!')
    end

end


function [ppg,vpg,apg,jpg]=get_signals(signal,Fs)
    %% Moving average filter 50 ms and 10 ms window
    Fn = Fs/2;
    fL=0.5;
    fH=12;
    order=4;
    [A,B,C,D] = cheby2(order,20,[fL,fH]/Fn);
    [filter_SOS,g] = ss2sos(A,B,C,D);
    filt_sig_cb2=filtfilt(filter_SOS,g,signal);
        
    win=round(Fs*0.05);
    B = 1/win*ones(win,1);
    filt_sig_MAFcb2=filtfilt(B,1,filt_sig_cb2);
    ppg=filt_sig_MAFcb2;
    %%

    med_max = median(findpeaks(ppg,'MinPeakDistance',Fs/2));
    med_min = median(-findpeaks(-ppg,'MinPeakDistance',Fs/2));
    norm_factor=(med_max-med_min);

    vpg = diff(ppg);
    win=round(Fs*0.01);
    B = 1/win*ones(win,1);
    
    vpg = filtfilt(B,1,vpg);
    vpg = vpg/(max(vpg)-min(vpg))*norm_factor;
    vpg = vpg+mean(ppg);

    apg = diff(vpg);
    apg = filtfilt(B,1,apg);
    apg = apg/(max(apg)-min(apg))*norm_factor;
    apg = apg+mean(ppg);

    jpg = diff(apg);
    jpg = filtfilt(B,1,jpg);
    jpg = jpg/(max(jpg)-min(jpg))*norm_factor;
    jpg = jpg+mean(ppg);
end

