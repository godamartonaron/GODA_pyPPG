function [pw_inds, fid_pts, pulses, sigs] = PulseAnalyse2(S, options, pulses)
% PULSEANALYSE  Extracts pulse wave (PW) indices from pulse waves.
%   
%  Inputs:
%
%    S        -  a pulsatile signal, consisting of ueither a single pulse, or a
%                 signal containing several pulses. S should be a structure, containing
%                 a vector of amplitudes, S.v, and the sampling frequency (in Hz), S.fs,
%                 and optionally the subject's height in metres, S.ht.
%    options  -  (optional) a structure of options which determine the settings used for the analysis:
%                    options.do_plot                    - (1 or 0, default value of 1) A logical indicating whether or not to make plots
%                    options.exclude_low_quality_data   - (default of 1) A logical indicating whether or not to exclude low quality pulse waves from the calculation of median values of pulse wave indices
%                    options.plot_third_deriv           - (default of 0) A logical indicating whether or not to plot the third derivative
%                    options.close_figures              - (default of 1) A logical indicating whether or not to close all Matlab figures before running the analysis
%                    options.do_filter                  - (default of 1) A logical indicating whether or not to filter pulse waves prior to analysis
%                    options.do_quality                 - (default of 1) A logical indicating whether or not to assess the quality of pulse waves
%                    options.save_folder                - (default of '') A string containing the path of a directory in which to save plots
%                    options.save_file                  - (default of '') A string containing the filename under which to save plots (without an extension)
%                    options.beat_detector              - (default of 'IMS') The type of beat detector to use when analysing a prolonged recording of multiple pulse waves. Definitions: IMS - [incremental merge segmentation algorithm](https://doi.org/10.1109/EMBC.2012.6346628) (implemented by M.A.F Pimentel as part of the [RRest](https://github.com/peterhcharlton/RRest) toolbox).
%                    options.verbose                    - (default of 0) A logical indicating whether or not to provide text outputs for any warnings
%                    options.plot_areas                 - (default of 0) A logical indicating whether or not to plot the systolic and diastolic areas on the pulse wave
%                    options.plot_pw_only               - (default of 0) A logical indicating whether or not to plot only the pulse wave and not its derivatives
%                    options.normalise_pw               - (default of 1) A logical indicating whether or not to normalise the pulse wave to occupy a range of 0 to 1
%                    options.calc_pw_inds               - (default of 1) A logical indicating whether or not to identify fiducial points and calculate PW indices
%                    options.downsample                 - (default of 0) A logical indicating whether or not to downsample the PPG signal prior to analysis
%                    options.downsample_freq            - (default of []) The downsampling frequency (only used if options.downsample is true).
%
%  Outputs:
%
%    pw_inds  -  a structure containing the calculated pulse wave indices. 
%                   For instance, pw_inds.AI.v is the value of the
%                   augmentation index (AI). If the input S contains
%                   several pulses, then  pw_inds.AI.v is the median AI
%                   value, and pw_inds.AI.raw is a vector of raw values for
%                   each pulse.
%    fid_pts  -  a structure containing:
%                   fid_pts.ind: the indices of the fiducial points
%                   fid_pts.amp: the ampltiudes of the fiducial points
%                   fid_pts.amp_norm: the amplitudes of the fiducial points (when the pulse wave is normalised to occupy a range of 0 to 1) 
%                   fid_pts.t: the times of the fiducial points
%                   For instance, fid_pts.ind.dic provides the indices of the
%                   dicrotic notch(es).
%    pulses   -  a structure containing information on the pulses:
%                   pulses.peaks    - the indices of systolic peaks
%                   pulses.onsets   - the indices of pulse onsets (i.e. beginning of systole)
%                   pulses.quality  - a logical indicating the signal quality of each pulse (true indicates high quality)
%    sigs     -  the following signals are provided:
%                   sigs.orig - the original signal
%                   sigs.filt - a filtered signal (if asked for)
%                   sigs.first_d - first derivative
%                   sigs.second_d - second derivative
%                   sigs.third_d - third derivative
%
%  Exemplary usage:
%
%    [pw_inds, fid_pts, pulses, sigs] = PulseAnalyse            performs a demo pulse wave analysis (of a photoplethysmogram, PPG, pulse wave)
%    [pw_inds, fid_pts, pulses, sigs] = PulseAnalyse([],'pressure_example')     performs a demo pressure pulse wave analysis
%    [pw_inds, fid_pts, pulses, sigs] = PulseAnalyse(S)         analyses the pulsatile signal specified in S
%    [pw_inds, fid_pts, pulses, sigs] = PulseAnalyse(S, options)                uses options to specify the settings for the analysis
%
%  For further information please see the accompanying manual: https://github.com/peterhcharlton/pulse-analyse/wiki/Examples
%
%   Licence:
%       Available under the GNU public license - please see the accompanying
%       file named "LICENSE"
%
%   Sources:
%       This script contains items either copied or modified from the RRest
%       toolbox which is covered by the GNU public licence (<a href="http://github.com/peterhcharlton/RRest/">link</a>).
%
% Author: Peter H. Charlton, University of Cambridge, April 2021.
% Version: 1.3beta, and is still in development.

%% Setup
% - check to see whether a signal has been provided
if nargin == 2 && strcmp(options, 'pressure_example'), fprintf('\n - Running basic pressure example'), S = basic_example_S('pressure'); clear options, options.normalise_pw = 0; options.do_filter = 0; end
if nargin == 0 || isempty(S) || ~sum(strcmp(fieldnames(S),'v')) || isempty(S.v), fprintf('\n - No input signal provided, so running basic PPG example'), S = basic_example_S('ppg'); end
% - create options structure if it hasn't been provided
if nargin < 2, options = struct; end
% - insert nan for height if it's not provided
if ~sum(strcmp(fieldnames(S), 'ht')), S.ht = nan; end
% - convert height into metres if it's in cm
if S.ht > 20, S.ht = S.ht/100; end
% - make sure S.v is a column vector
S.v = S.v(:);
% - setup universal parameters
up = setup_up(S, options);

%% Downsample
S = do_downsample(S, up);

% %% Beat Detection
% pulses = beat_detection(S, up);

%% Quality Assessment
pulses.quality = assess_signal_quality(S, pulses, up);

%% Filter Signal (if required)
[sigs, pulses] = filter_signal(S, pulses, up);

%% Calibrate (if required)
sigs = calib_signal(sigs, pulses, up);

%% Generate central pressure (if required)
[sigs, pulses] = generate_central_pressure(sigs, pulses, up);

%% Calculate Average Pulse Wave (if required)
[sigs, pulses, up] = calculate_average_pw(sigs, pulses, up);

%% Calculate Derivatives
sigs = calculate_derivs(sigs, up, S);

%% Identify Fiducial Points, and Calculate Fiducial Point Timings and Amplitudes
fid_pts = identify_fiducial_point_indices(sigs, pulses, up);
%%% MG 13/06/2023
up.analysis.no_of_pulses = 'multiple';
p = beat_detection(S, up);
diff_os=pulses.peaks(1)-p.onsets;
tempi_os=find(diff_os>0);
[v_os,i_os]=min(diff_os(tempi_os));
fid_pts.ind.os=p.onsets(i_os);
%%%


%% Calculate Pulse Wave Indices
pw_inds = calculate_pw_inds(fid_pts, sigs, pulses, S.ht, up);

%% Plot Results
plot_results(sigs, pulses, fid_pts, pw_inds, up);

%% Make pre-processing plot
make_pre_processing_plot(sigs, pulses, up);

end

function options = setup_options(options)

if isempty(fieldnames(options)) || ~sum(strcmp(fieldnames(options), 'do_plot'))
    options.do_plot = 1;
end

if ~strcmp(fieldnames(options), 'exclude_low_quality_data')
    options.exclude_low_quality_data = 1;
end
if ~strcmp(fieldnames(options), 'plot_third_deriv')
    options.plot_third_deriv = 0;
end
if ~strcmp(fieldnames(options), 'close_figures')
    options.close_figures = 1;
end
if ~strcmp(fieldnames(options), 'do_filter')
    options.do_filter = 1;
end
if ~strcmp(fieldnames(options), 'do_quality')
    options.do_quality = 1;
end
if ~strcmp(fieldnames(options), 'save_folder')
    options.save_folder = '';
end
if ~strcmp(fieldnames(options), 'save_file')
    options.save_file = '';
end
if ~strcmp(fieldnames(options), 'beat_detector')
    options.beat_detector = 'IMS';
    %options.beat_detector = 'AMPD'; %%%%%%%%%%%%%% CHANGE ^%%%%%%%%%%%
end
if ~strcmp(fieldnames(options), 'verbose')
    options.verbose = false;
end
if ~strcmp(fieldnames(options), 'plot_areas')
    options.plot_areas = false;
end
if ~strcmp(fieldnames(options), 'plot_pw_only')
    options.plot_pw_only = false;
end
if ~strcmp(fieldnames(options), 'normalise_pw')
    if sum(strcmp(fieldnames(options), 'calibrate_pw')) & options.calibrate_pw
        options.normalise_pw = false;
    else
        options.normalise_pw = true;
    end
end
if ~strcmp(fieldnames(options), 'calibrate_pw')
    options.calibrate_pw = false;
end
if ~strcmp(fieldnames(options), 'calc_average_pw')
    options.calc_average_pw = false;
end
if ~strcmp(fieldnames(options), 'tran_func')
    options.tran_func = 0;
end
if ~strcmp(fieldnames(options), 'sig_type')
    options.sig_type = '';
end
if ~strcmp(fieldnames(options), 'calc_pw_inds')
    options.calc_pw_inds = 1;
end
if ~strcmp(fieldnames(options), 'plot_f_pt')
    options.plot_f_pt = 0;
end
if ~strcmp(fieldnames(options), 'downsample')
    options.downsample = 0;
end
if ~strcmp(fieldnames(options), 'downsample_freq')
    options.downsample_freq = [];
end
if ~strcmp(fieldnames(options), 'do_beat_filter')
    options.do_beat_filter = 1;
end

% insert filename if needed
if ~isempty(options.save_folder) && isempty(options.save_file)
    options.save_file = 'PulseAnalyse_outputs';
end

end

function S = basic_example_S(signal_type)

if strcmp(signal_type, 'ppg')
    % Data from the PPGDiary Pilot Study (a single pulse wave acquired using an infrared PPG sensor on the thumb)
    S.v = [-14880;-14508;-13390;-11384;-8573;-5236;-1745;1569;4502;6979;8992;10536;11600;12196;12396;12328;12123;11860;11547;11160;10700;10218;9789;9458;9211;8990;8742;8465;8206;8016;7894;7778;7573;7223;6747;6239;5814;5554;5478;5556;5744;6010;6339;6720;7137;7561;7961;8309;8580;8760;8840;8825;8732;8580;8384;8151;7884;7595;7306;7042;6817;6629;6460;6286;6093;5874;5630;5364;5077;4771;4453;4132;3822;3534;3273;3041;2828;2622;2410;2187;1967;1771;1617;1501;1391;1243;1029;760;486;269;140;79;28;-69;-231;-425;-599;-716;-784;-847;-949;-1104;-1291;-1475;-1639;-1800;-1993;-2239;-2524;-2798;-3011;-3150;-3256;-3395;-3613;-3898;-4189;-4417;-4562;-4658;-4771;-4942;-5167;-5410;-5645;-5865;-6077;-6270;-6414;-6491;-6521;-6568;-6694;-6911;-7167;-7385;-7519;-7588;-7653;-7778;-7982;-8245;-8523;-8777;-8992;-9169;-9323;-9475;-9648;-9855;-10093;-10344;-10577;-10768;-10908;-11012;-11116;-11267;-11502;-11815;-12151;-12427;-12585;-12644;-12698;-12845;-13082;-13241];
    S.v = [-368;-355;-303;-206;-75;69;202;310;388;439;468;478;471;451;421;385;346;305;266;228;196;170;150;136;127;119;109;96;78;61;47;40;38;39;39;39;44;54;70;86;98;104;104;101;97;91;82;70;55;38;22;7;-5;-13;-19;-27;-38;-51;-63;-70;-73;-74;-79;-92;-108;-120;-126;-124;-123;-128;-142;-160;-176;-184;-187;-188;-191;-198;-207;-215;-221;-227;-235;-248;-265;-282;-296;-303;-304;-306;-316;-334;-352];
    S.v = repmat(S.v, [10,1]);
    S.v = movmean(S.v, 7);
    S.fs = 100; % in Hz
elseif strcmp(signal_type, 'pressure')
    % A simulated carotid pressure pulse wave, from the pulse wave database described at: https://peterhcharlton.github.io/pwdb/'
    S.v = [74.83;74.86;75.02;75.35;75.89;76.6;77.44;78.36;79.31;80.27;81.2;82.12;82.99;83.83;84.62;85.36;86.05;86.71;87.32;87.92;88.49;89.06;89.61;90.16;90.71;91.26;91.8;92.34;92.88;93.41;93.95;94.48;95.02;95.56;96.11;96.67;97.22;97.77;98.32;98.86;99.39;99.9;100.4;100.9;101.3;101.8;102.2;102.5;102.9;103.2;103.5;103.7;103.9;104.1;104.2;104.3;104.4;104.4;104.4;104.4;104.4;104.3;104.3;104.2;104.1;103.9;103.8;103.7;103.5;103.4;103.2;103;102.9;102.7;102.5;102.3;102.1;102;101.8;101.6;101.4;101.2;101.1;100.9;100.7;100.5;100.4;100.2;100.1;99.94;99.81;99.69;99.59;99.49;99.4;99.32;99.25;99.19;99.15;99.11;99.08;99.06;99.05;99.05;99.06;99.08;99.11;99.14;99.18;99.23;99.28;99.34;99.4;99.46;99.53;99.6;99.67;99.74;99.81;99.88;99.95;100;100.1;100.1;100.2;100.2;100.2;100.2;100.2;100.2;100.2;100.1;100.1;100;99.9;99.78;99.63;99.46;99.27;99.04;98.79;98.52;98.23;97.93;97.61;97.29;96.97;96.65;96.35;96.11;95.96;95.92;95.99;96.14;96.29;96.42;96.49;96.5;96.46;96.37;96.25;96.11;95.95;95.77;95.59;95.41;95.23;95.06;94.89;94.72;94.55;94.39;94.24;94.09;93.94;93.81;93.68;93.57;93.47;93.37;93.3;93.23;93.17;93.12;93.07;93.03;92.99;92.96;92.92;92.88;92.84;92.79;92.75;92.7;92.64;92.58;92.52;92.46;92.39;92.32;92.25;92.17;92.1;92.02;91.94;91.85;91.77;91.69;91.61;91.53;91.45;91.36;91.28;91.2;91.12;91.04;90.96;90.88;90.8;90.72;90.63;90.54;90.45;90.36;90.27;90.17;90.08;89.98;89.88;89.78;89.69;89.59;89.49;89.39;89.3;89.2;89.11;89.02;88.93;88.84;88.75;88.67;88.58;88.5;88.42;88.34;88.26;88.19;88.12;88.04;87.97;87.9;87.83;87.77;87.7;87.64;87.57;87.51;87.45;87.38;87.32;87.26;87.2;87.14;87.08;87.02;86.96;86.9;86.84;86.78;86.72;86.66;86.6;86.53;86.47;86.41;86.34;86.28;86.21;86.15;86.08;86.01;85.94;85.87;85.8;85.73;85.66;85.58;85.51;85.43;85.36;85.28;85.2;85.12;85.03;84.95;84.87;84.78;84.69;84.61;84.52;84.43;84.34;84.25;84.15;84.06;83.96;83.87;83.77;83.68;83.58;83.48;83.39;83.29;83.19;83.09;82.99;82.9;82.8;82.7;82.6;82.51;82.41;82.31;82.21;82.12;82.02;81.92;81.83;81.73;81.64;81.54;81.45;81.35;81.26;81.16;81.07;80.98;80.88;80.79;80.7;80.61;80.52;80.42;80.33;80.24;80.15;80.06;79.97;79.88;79.79;79.71;79.62;79.53;79.44;79.35;79.27;79.18;79.09;79.01;78.92;78.84;78.75;78.66;78.58;78.49;78.41;78.33;78.24;78.16;78.07;77.99;77.91;77.82;77.74;77.66;77.57;77.49;77.41;77.33;77.24;77.16;77.08;77;76.92;76.84;76.76;76.68;76.6;76.52;76.44;76.36;76.28;76.2;76.12;76.04;75.96;75.88;75.81;75.73;75.65;75.57;75.5;75.42;75.35;75.27;75.19;75.12;75.04;74.97;74.9];
    S.fs = 500; % in Hz
end

end

function S = do_downsample(S, up)

% Skip if not required
if ~up.options.downsample
    return
end

%% Downsample

% Work out whether to downsample or interpolate
if rem(S.fs, up.options.downsample_freq) == 0
    % Downsample
    ds_factor = S.fs/up.options.downsample_freq;
    S.v = downsample(S.v, ds_factor);
    S.fs = up.options.downsample_freq;
else
    use_resampling = 1;
    if use_resampling
        % Resample
        S.v(isnan(S.v)) = median(S.v(~isnan(S.v)));
        S.v = resample(S.v, up.options.downsample_freq, S.fs);
    else
        % Interpolate
        t = 1:length(S.v);
        t_new = 1:(S.fs/up.options.downsample_freq):t(end);
        S.v = interp1(t, S.v, t_new);
    end
    S.fs = up.options.downsample_freq;
end

end

function pulses = beat_detection(S, up)

if sum(strcmp(fieldnames(up.options), 'do_beats')) && ~up.options.do_beats
    [pulses.peaks, pulses.onsets] = deal([]);
    return
end


if strcmp(up.analysis.no_of_pulses, 'multiple')
    
    if up.options.do_beat_filter
        %% Filter to remove high frequencies
        filt_characteristics = up.filtering.beat_detection.elim_high_freqs;
        s_filt = elim_vhfs3(S, filt_characteristics);
        
        %% Filter to remove low frequencies
        filt_characteristics = up.filtering.beat_detection.elim_low_freqs;
        s_filt = elim_vlfs(s_filt, filt_characteristics, up);
    else
        s_filt = S;
    end
    
    %% Detect beats
    pulses = detect_beats(s_filt, up);
    
else

    %% Locate onsets and peak in single pulse wave
    pulses.onsets = [1;length(S.v)];
    [~, pulses.peaks] = max(S.v);
    desired_ht = mean(S.v([pulses.onsets(1), pulses.peaks(1)]));
    [~, pulses.mid_amps] = min(abs(S.v(pulses.onsets(1):pulses.peaks(1)) - desired_ht));
    
end

% Check

% plot(s_filt.v), hold on, pts = pulses.onsets; plot(pts, s_filt.v(pts), 'ok'), pts = pulses.peaks; plot(pts, s_filt.v(pts), 'or'), pts = pulses.mid_amps; plot(pts, s_filt.v(pts), '*b')


end

function up = setup_up(S, options)

% setup options (using defaults if any of the options aren't specified)
up.options = setup_options(options);

if up.options.close_figures
    close all
end

%% Analysis settings

% Filter characteristics: Eliminate VLFs (below resp freqs): For 4bpm cutoff
up.paramSet.elim_vlf.Fpass = 0.157;  % in Hz
up.paramSet.elim_vlf.Fstop = 0.02;   % in Hz     (0.157 and 0.02 provide a - 3dB cutoff of 0.0665 Hz)
up.paramSet.elim_vlf.Dpass = 0.05;
up.paramSet.elim_vlf.Dstop = 0.01;

% Filter characteristics: Eliminate VHFs (above frequency content of signals)
up.paramSet.elim_vhf.Fpass = 38.5;  % in HZ
up.paramSet.elim_vhf.Fstop = 33.12;  % in HZ   (33.12 and 38.5 provide a -3 dB cutoff of 35 Hz)
% up.paramSet.elim_vhf.Fpass = 20;  % in HZ
% up.paramSet.elim_vhf.Fstop = 15;  % in HZ
up.paramSet.elim_vhf.Dpass = 0.05;
up.paramSet.elim_vhf.Dstop = 0.01;

% No of times to repeat a single pulse to perform VHF filtering
up.paramSet.no_pulse_repeats = 5;

%% Filter Characteristics

if sum(strcmp(fieldnames(up.options),'beat_detection_filtering'))
    
    up.filtering.beat_detection = up.options.beat_detection_filtering;
    
else
    % - Beat detection: eliminating high freqs
    up.filtering.beat_detection.elim_high_freqs.Fpass = 38.5;  % in HZ
    up.filtering.beat_detection.elim_high_freqs.Fstop = 33.12;  % in HZ   (33.12 and 38.5 provide a -3 dB cutoff of 35 Hz)
    up.filtering.beat_detection.elim_high_freqs.Dpass = 0.05;
    up.filtering.beat_detection.elim_high_freqs.Dstop = 0.01;
    
    % - Beat detection: eliminating low freqs
    up.filtering.beat_detection.elim_low_freqs.Fpass = 0.3;  % in Hz
    up.filtering.beat_detection.elim_low_freqs.Fstop = 0.1;  % in Hz
    up.filtering.beat_detection.elim_low_freqs.Dpass = 0.05;
    up.filtering.beat_detection.elim_low_freqs.Dstop = 0.01;
    
end


if sum(strcmp(fieldnames(up.options),'fiducial_point_filtering'))
    
    up.filtering.fiducial_points = up.options.fiducial_point_filtering;
    
else
    
    % - Fiducial Points: eliminating high freqs
    up.filtering.fiducial_points.elim_high_freqs.Fpass = 20;  % in HZ
    up.filtering.fiducial_points.elim_high_freqs.Fstop = 15;  % in HZ
    up.filtering.fiducial_points.elim_high_freqs.Dpass = 0.05;
    up.filtering.fiducial_points.elim_high_freqs.Dstop = 0.01;
    
    % - Fiducial Points: eliminating low freqs
    up.filtering.fiducial_points.elim_low_freqs.Fpass = 0.3;  % in Hz
    up.filtering.fiducial_points.elim_low_freqs.Fstop = 0.1;  % in Hz
    up.filtering.fiducial_points.elim_low_freqs.Fpass = 0.8;  % in Hz
    up.filtering.fiducial_points.elim_low_freqs.Fstop = 0.6;  % in Hz
    up.filtering.fiducial_points.elim_low_freqs.Dpass = 0.05;
    up.filtering.fiducial_points.elim_low_freqs.Dstop = 0.01;
    
end

% - Derivative: number of points to use in Savitzky-Golay filter
up.filtering.derivatives.s_g_filt_len_no_filtering = 5;
up.filtering.derivatives.s_g_filt_len_no_filtered = 9;

% - Pulse onsets: tolerance in which to search for minima
up.paramSet.onset_tol = 0.05; % duration of search region either side of current minimum (in secs)

%% Number of pulses

% Threshold signal duration to distinguish between a single pulse or multiple pulses: 
up.analysis.max_duration_of_single_pulse = 2.5;   % in secs

% Determine how many pulses there are:
up.analysis.no_of_pulses = determine_no_of_pulses(S, up);

end

function no_of_pulses = determine_no_of_pulses(S, up)
% DETERMINE_NO_OF_PULSES  Determines whether this is a single pulse wave,
% of a pulsatile signal containing multiple pulse waves.

signal_duration = (length(S.v)-1)/S.fs;

% If duration of signal is greater than a threshold then assume this is
% multiple pulse waves:
if signal_duration > up.analysis.max_duration_of_single_pulse
    no_of_pulses = 'multiple';
else
    no_of_pulses = 'single';
end

end

function s_filt = elim_vlfs(s, filt_characteristics, up)
%% Filter pre-processed signal to remove frequencies below resp
% Adapted from RRest

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Make filter
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord([filt_characteristics.Fstop filt_characteristics.Fpass]/(s.fs/2), [1 0], [filt_characteristics.Dstop filt_characteristics.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% % Gives a -3 dB cutoff at ? Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 4.435e-3;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(s.fs/2);

if length(s.v) > (length(AMfilter.numerator)-1)*3
    % - Tukey Window to avoid edge effects
    win_edge_durn = 0.5; % in secs
    prop_win = win_edge_durn*2/((length(s.v)-1)/s.fs);
    tw = tukeywin(length(s.v),prop_win); 
    s_filt.v = filtfilt(AMfilter.numerator, 1, s.v.*tw);
    s_filt.v = s.v-s_filt.v;
else
    s_filt.v = s.v;
end
s_filt.fs = s.fs;
end

function pulses = detect_beats(s, up)

%% Detect peaks in pulse signal
switch up.options.beat_detector
    case 'AMPD'
        pulses.peaks = ampd_peak_detector(s.v,s.fs,up);
    case 'IMS'
        [pulses.peaks, pulses.onsets, ~] = adaptPulseSegment(s.v,s.fs);
        for wave_no = 1 : length(pulses.peaks)-1
            [~, temp] = max(s.v(pulses.onsets(wave_no):pulses.onsets(wave_no+1)));
            pulses.peaks(wave_no) = temp+pulses.onsets(wave_no)-1;
        end
    case 'ABD'
        pulses.peaks = abd_algorithm(s);
    case 'MSPTD'
        [pulses.peaks, pulses.onsets] = bishop_peak_detector(s);
    case 'Pulses'
        [pulses.peaks, ~, peaks_filt, ~] = PPG_pulses_detector(s.v, s.fs);
        pulses.peaks = pulses.peaks(:);
        % note that I haven't used peaks_filt
    case 'qppg'
        [pulses.onsets] = qppg(s.v,s.fs);
        pulses.onsets = pulses.onsets(:);
    case 'HeartPy'
        pulses.peaks = heartpy_pc(s);
    case 'COppg'
        [pulses.peaks, pulses.onsets] = COppg_peak_detector(s.v, s.fs);
    case 'SPAR3'
        % uses external file
        N = 3;
        [beat_indices_all,cycle_times_all] = SPARcycleTimesPPG(s.v,s.fs,N);
        pulses.peaks = nan(length(beat_indices_all)-1,1);
        for beat_no = 1 : length(beat_indices_all)-1
            curr_ind = beat_indices_all(beat_no);
            [~, temp] = max(s.v(curr_ind:beat_indices_all(beat_no+1)));
            pulses.peaks(beat_no) = curr_ind+temp-1;
        end
    case 'SPAR7'
        % uses external file
        N = 7;
        [beat_indices_all,cycle_times_all] = SPARcycleTimesPPG(s.v,s.fs,N);        
        pulses.peaks = nan(length(beat_indices_all)-1,1);
        for beat_no = 1 : length(beat_indices_all)-1
            curr_ind = beat_indices_all(beat_no);
            [~, temp] = max(s.v(curr_ind:beat_indices_all(beat_no+1)));
            pulses.peaks(beat_no) = curr_ind+temp-1;
        end
    case 'ERMA'
        pulses.peaks = erma(s);
    case 'PWD'
        [pulses.onsets, pulses.peaks] = pwd(s);
    case 'PDA'
        pulses.peaks = pda(s);
    case 'WFD'
        pulses.peaks = wfd(s);
end

% tidy up peaks and onsets
if sum(strcmp(fieldnames(pulses),'onsets')) && sum(strcmp(fieldnames(pulses),'peaks'))
    
    % If there are two consecutive peaks, then insert a trough between them
    peak_log = [true(size(pulses.peaks)); false(size(pulses.onsets))];
    [els,order] = sort([pulses.peaks;pulses.onsets]);
    peak_log = peak_log(order);
    bad_els = find(diff(peak_log)==0 & peak_log(1:end-1)); % repeated peaks
    if ~isempty(bad_els)  % if there is a repeated peak
        new_troughs = nan(length(bad_els),1);
        for bad_el_no = 1 : length(bad_els)   % cycle through each repeated peak
            curr_pks = [els(bad_els(bad_el_no)),els(bad_els(bad_el_no)+1)];
            [~,temp]=min(s.v(curr_pks(1):curr_pks(2)));
            new_troughs(bad_el_no,1) = curr_pks(1) -1 + temp;
        end
        pulses.onsets = sort([pulses.onsets; new_troughs]);
    end
    clear new_troughs temp curr_pks bad_el_no bad_els peak_log order els
    
    % If there are two consecutive troughs, then insert a peak between them
    trough_log = [false(size(pulses.peaks)); true(size(pulses.onsets))];
    [els,order] = sort([pulses.peaks;pulses.onsets]);
    trough_log = trough_log(order);
    bad_els = find(diff(trough_log)==0 & trough_log(1:end-1)); % repeated troughs
    if ~isempty(bad_els)  % if there is a repeated peak
        new_peaks = nan(length(bad_els),1);
        for bad_el_no = 1 : length(bad_els)   % cycle through each repeated peak
            curr_trs = [els(bad_els(bad_el_no)),els(bad_els(bad_el_no)+1)];
            [~,temp]=max(s.v(curr_trs(1):curr_trs(2)));
            new_peaks(bad_el_no,1) = curr_trs(1) -1 + temp;
        end
        pulses.peaks = sort([pulses.peaks; new_peaks]);
    end
    clear new_peaks temp curr_trs bad_el_no bad_els trough_log order els
    
    % Make sure that the first onset is before the first peak, and the last peak is after the last onset
    if ~isempty(pulses.onsets) && ~isempty(pulses.peaks)
        if pulses.onsets(1)>pulses.peaks(1)
            pulses.peaks(1) = [];
        end
        if pulses.peaks(end)<pulses.onsets(end)
            pulses.onsets(end) = [];
        end
    end
end

% Add onsets if required
if ~sum(strcmp(fieldnames(pulses),'onsets'))
    pulses.onsets = nan(length(pulses.peaks)-1,1);
    for wave_no = 1 : length(pulses.peaks)-1
        [~, temp] = min(s.v(pulses.peaks(wave_no):pulses.peaks(wave_no+1)));
        pulses.onsets(wave_no) = temp + pulses.peaks(wave_no) - 1;
    end
    pulses.peaks = pulses.peaks(2:end); % so that there are the same number of peaks and onsets
end

% Add peaks if required
if ~sum(strcmp(fieldnames(pulses),'peaks'))
    pulses.peaks = nan(length(pulses.onsets)-1,1);
    for wave_no = 1 : length(pulses.onsets)-1
        [~, temp] = max(s.v(pulses.onsets(wave_no):pulses.onsets(wave_no+1)));
        pulses.peaks(wave_no,1) = temp + pulses.onsets(wave_no) - 1;
    end
    pulses.onsets = pulses.onsets(1:end-1); % so that there are the same number of peaks and onsets
end

% Add mid-amplitude point if required
if ~sum(strcmp(fieldnames(pulses),'mid_amps'))
    pulses.mid_amps = nan(length(pulses.onsets),1);
    for wave_no = 1 : length(pulses.onsets)
        desired_ht = mean(s.v([pulses.onsets(wave_no), pulses.peaks(wave_no)]));
        [~, temp] = min(abs(s.v(pulses.onsets(wave_no):pulses.peaks(wave_no)) - desired_ht));
        pulses.mid_amps(wave_no,1) = temp + pulses.onsets(wave_no) - 1;
    end
end

end

function [overall_peaks, overall_onsets] = COppg_peak_detector(ppg, fs)
%COPPG_PEAK_DETECTOR detects PPG pulse peaks
%
% This has been adapted by Peter Charlton from code supplied by Christina
% Orphanidou.
%
% Purpose: To detect PPG peaks.
%
% Original source: https://github.com/peterhcharlton/RRest/blob/master/RRest_v3.0/Algorithms/extract_resp_sig/feat_based_extraction/COr_peak_detector/co_ppg_peak_detector.m
% Modifications from original by PC:
%  - Removed pre-processing
%  - Hard-coded constants
%
% With grateful thanks to Christina Orphanidou and Alexander Darrell.
%

%% Constants
% Hard-coded by PC
up.paramSet.CO_peak_det.upctl = 0.9;
up.paramSet.CO_peak_det.lpctl = 0.1;

%% Pre-processing
% removed by PC (see original for details)
ppgfilt.v = ppg;
ppgfilt.t = [0:length(ppgfilt.v)-1]./fs;

%% Setup
% Segmentation into 10s windows
win_length = 10;
overlap = 3;
win_starts = ppgfilt.t(1):win_length:ppgfilt.t(end);
win_ends = win_starts + win_length + overlap;
win_ends(end) = win_ends(end) - overlap;
% setup variables
[overall_peaks, overall_onsets] = deal([]);

%% Perform peak detection on each 10s window
for win_no = 1 : length(win_starts)
    
    rel_data = struct;
    
    rel_data.i = find(ppgfilt.t >= win_starts(win_no) & ppgfilt.t <= win_ends(win_no));
    rel_data.t = ppgfilt.t(rel_data.i);
    rel_data.v = ppgfilt.v(rel_data.i);
    
    %% Calculate thresholds
    thresh1=quantile(rel_data.v,up.paramSet.CO_peak_det.upctl);
    thresh2=quantile(rel_data.v,up.paramSet.CO_peak_det.lpctl);
    thresh3=thresh2+0.3*(thresh1-thresh2);
    thresh4=thresh2+0.7*(thresh1-thresh2);
    
    %% Find all peaks
    % identify peaks
    diffs_on_left_of_pt = diff(rel_data.v); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt>0);
    diffs_on_right_of_pt = diff(rel_data.v); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); level_on_right_of_pt = logical(diffs_on_right_of_pt==0); diffs_on_right_of_pt = logical(diffs_on_right_of_pt<0);
    diffs_on_second_right_of_pt = [1; diff(rel_data.v(2:end))]; diffs_on_second_right_of_pt = diffs_on_second_right_of_pt(2:end); diffs_on_second_right_of_pt = logical(diffs_on_second_right_of_pt<0);
    peaks.i = find((diffs_on_left_of_pt & diffs_on_right_of_pt) | ...
        (diffs_on_left_of_pt & level_on_right_of_pt & diffs_on_second_right_of_pt))+1;
    peaks.i = peaks.i';
    % Take maximum value in window if no peaks were found
    if isempty(peaks.i)
        peaks.i = find(max(rel_data.v) == rel_data.v);
    end
    % Extract signal values at peaks
    peaks.t = rel_data.t(peaks.i);
    peaks.v = rel_data.v(peaks.i);
    
    %% Classify peaks
    % Identify relevant peaks according to amplitude
    upperdiff = abs(peaks.v-thresh1);
    middlehighdiff = abs(peaks.v-thresh4);
    middlelowdiff = abs(peaks.v-thresh3);
    lowerdiff = abs(peaks.v-thresh2);
    upper_pks = find(upperdiff<middlehighdiff & upperdiff<middlelowdiff & upperdiff<lowerdiff);
    PPG_PKS.i = peaks.i(upper_pks);
    PPG_PKS.v = peaks.v(upper_pks);
    PPG_PKS.t = peaks.t(upper_pks);
    % eliminate peaks which are too close together in time
    good_els = find(diff(PPG_PKS.i)>=fs/3)+1;
    PPG_PKS.i = PPG_PKS.i(good_els);
    PPG_PKS.v = PPG_PKS.v(good_els);
    PPG_PKS.t = PPG_PKS.t(good_els);
    
    %% Find troughs
    
    PPG_TRS.i = nan(length(PPG_PKS.t)-1,1);
    for s = 1 : (length(PPG_PKS.t)-1)
        start_el = PPG_PKS.i(s);
        [~, additional_el] = min(rel_data.v(PPG_PKS.i(s):PPG_PKS.i(s+1)));
        PPG_TRS.i(s) = start_el - 1 + additional_el;
    end
    PPG_TRS.v = rel_data.v(PPG_TRS.i);
    PPG_TRS.t = rel_data.t(PPG_TRS.i);
    
    temp_pk_is = rel_data.i(PPG_PKS.i(:));
    temp_on_is = rel_data.i(PPG_TRS.i(:));
    
    overall_peaks = [overall_peaks; temp_pk_is(:)]; overall_peaks = unique(overall_peaks);
    overall_onsets = [overall_onsets; temp_on_is(:)]; overall_onsets = unique(overall_onsets);
    
    clear upper lower middlehigh middlelow temp_pk_is temp_on_is PPG_TRS PPG_PKS additional_el s good_els peaks lower_pks upper_pks lowerdiff middlelowdiff middlehighdiff upperdiff ind
    
end

end

function p = abd_algorithm(data)
% ABD_ALGORITHM  identifies pulse peaks in a pulsatile signal. It is an
% adaptation of the algorithm described in:
%
%    Aboy, M. et al., An automatic beat detection algorithm for pressure
%    signals. IEEE Trans. Biomed. Eng. 2005, 52, 1662?1670,
%      http://doi.org/10.1109/TBME.2005.855725
%   
%  Inputs:
%
%    data     -  a pulsatile signal formatted as a structure, containing
%                 a vector of amplitudes, data.v, and the sampling
%                 frequency (in Hz), data.fs.
%    options  -  (optional) a structure of options, which may contain any of:
%                    options.                    - a logical (true or false)
%
%  Outputs:
%
%    p        -  a vector of pulse indices
%
%  Exemplary usage:
%
%    p = abd_algorithm(data)                              extracts pulse indices from the pulsatile signal, data.
%    p = abd_algorithm(data, options)                     uses options to specify the analysis.
%    [p, S_filt] = abd_algorithm(___)                     also outputs the filtered signals.
%
%  For further information please see the accompanying manual.
%
%  This script contains items either copied or modified from the RRest
%  toolbox which is covered by the GNU public licence (<a href="http://github.com/peterhcharlton/RRest/">link</a>).
%
% Peter H. Charlton, King's College London, August 2017

% Changes from original algorithm description:
%  1) PSD calculation method may not be exactly the same
%  2) Not conducted on windows of 10 s
%  3) Band-pass filtering may not produce exactly the right cut-offs
%  4) Wasn't sure what the lower and upper HR values were, so used 30 and 200 bpm
%  5) Changed the proportion of the lower HR value at which to draw the lower cut-off (from 50% to 80%)
%  6) Changed the percentile threshold for identifying peaks in the derivative from 90% to 75%
%  7) Haven't implemented harmonic PSD
%  8) HR estimation only considers HRs within the range of plausible HRs

% load data
if nargin<1
    load('/Users/petercharlton/Downloads/ppg_data.mat');
    data.v = -1*data.v;
end

% inputs
x = data.v;  % signal
fs = data.fs; % sampling freq
up = setup_up_abd_algorithm; % settings
w = fs*10; % window length (number of samples)
win_starts = 1:round(0.8*w):length(x);
win_starts(win_starts>=length(x)-w+1) = [];
win_starts = [win_starts, length(x)+1-w];

% before pre-processing
px = DetectMaxima(x,0);  % detect all maxima
if isempty(px)
    p = [];
    return
end

% detect peaks in windows
all_p4 = [];
all_hr = nan(length(win_starts)-1,1);
for win_no = 1 : length(win_starts)-1
    curr_els = win_starts(win_no):win_starts(win_no)+w-1;
    curr_x = x(curr_els);
    y1 = Bandpass(curr_x, fs, 0.9*up.fl_hz, 3*up.fh_hz); % Filter no.1
    hr = EstimateHeartRate(y1, fs, up); % Estimate HR from weakly filtered signal
    all_hr(win_no) = hr;
    y2 = Bandpass(curr_x, fs, 0.9*up.fl_hz, 2.5*hr/60); % Filter no.2
    y2_deriv = EstimateDeriv(y2); % Estimate derivative from highly filtered signal
    p2 = DetectMaxima(y2_deriv,up.deriv_threshold); % Detect maxima in derivative
    % plot(x), hold on, plot(p2, x(p2), 'or')
    % plot(y2_deriv), hold on, thresh = prctile(y2_deriv, up.deriv_threshold); plot([0,length(y2_deriv)], thresh*[1,1])
    y3 = Bandpass(curr_x, fs, 0.9*up.fl_hz, 10*hr/60);
    p3 = DetectMaxima(y3,60); % Detect maxima in moderately filtered signal
    % plot(x), hold on, plot(p3, x(p3), 'or')
    p4 = find_pulse_peaks(p2,p3);
%     plot(curr_x), hold on, plot(p4, curr_x(p4), 'or')
    all_p4 = [all_p4;win_starts(win_no)+p4-1];
end

all_p4 = unique(all_p4);

if ~isempty(all_p4)
    [p, fn] = IBICorrect(all_p4, px, median(all_hr), fs, up);
    p = unique(p);
else
    p = all_p4;
end
% plot(x), hold on, plot(p, x(p), 'or')

% plot(y1), hold on, plot(p4, y1(p4), 'or')


end

function up = setup_up_abd_algorithm

% plausible HR limits
up.fl = 30; % lower bound for HR
up.fh = 200; % upper bound for HR
up.fl_hz = up.fl/60;
up.fh_hz = up.fh/60;

% Thresholds
up.deriv_threshold = 75;            % originally 90
up.upper_hr_thresh_prop = 2.25;     % originally 1.75
up.lower_hr_thresh_prop = 0.5;     % originally 0.75

% Other parameters
up.win_size = 10; % in secs

end

function mt = DetectMaxima(sig,percentile)

% Table VI pseudocode

tr = prctile(sig, percentile);
ld = length(sig);

m = 1+find(sig(3:end) < sig(2:end-1) & sig(1:end-2) < sig(2:end-1));
mt = m(sig(m)>tr);

end

function bpf_sig = Bandpass(sig, fs, lower_cutoff, upper_cutoff)

% Filter characteristics: Eliminate VLFs (below resp freqs): For 4bpm cutoff
up.paramSet.elim_vlf.Fpass = 1.3*lower_cutoff;  % in Hz
up.paramSet.elim_vlf.Fstop = 0.8*lower_cutoff;   % in Hz
up.paramSet.elim_vlf.Dpass = 0.05;
up.paramSet.elim_vlf.Dstop = 0.01;

% Filter characteristics: Eliminate VHFs (above frequency content of signals)
up.paramSet.elim_vhf.Fpass = 1.2*upper_cutoff;  % in HZ
up.paramSet.elim_vhf.Fstop = 0.8*upper_cutoff;  % in HZ
up.paramSet.elim_vhf.Dpass = 0.05;
up.paramSet.elim_vhf.Dstop = 0.03;

% perform BPF
s.v = sig;
s.fs = fs;
s_evlf = elim_vlfs_abd(s, up);
s_filt = elim_vhfs(s_evlf, up);
bpf_sig = s_filt.v;

end

function s_filt = elim_vlfs_abd(s, up)
%% Filter pre-processed signal to remove frequencies below resp
% Adapted from RRest

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Make filter
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord([up.paramSet.elim_vlf.Fstop up.paramSet.elim_vlf.Fpass]/(s.fs/2), [1 0], [up.paramSet.elim_vlf.Dstop up.paramSet.elim_vlf.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% % Gives a -3 dB cutoff at ? Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.0266;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(fs/2);

try
    s_filt.v = filtfilt(AMfilter.numerator, 1, s.v);
    s_filt.v = s.v-s_filt.v;
catch
    s_filt.v = s.v;
end
s_filt.fs = s.fs;
end

function s_filt = elim_vhfs(s, up)
%% Filter signal to remove VHFs
% Adapted from RRest

s_filt.fs = s.fs;

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Check to see if sampling freq is at least twice the freq of interest
if (up.paramSet.elim_vhf.Fpass/(s.fs/2)) >= 1
    % then the fs is too low to perform this filtering
    s_filt.v = s.v;
    return
end

%% Create filter
% parameters for the low-pass filter to be used
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord([up.paramSet.elim_vhf.Fstop up.paramSet.elim_vhf.Fpass]/(s.fs/2), [1 0], [up.paramSet.elim_vhf.Dstop up.paramSet.elim_vhf.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% % Gives a -3 dB cutoff at cutoff_freq Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.3355;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(s.fs/2);

%% Remove VHFs
s_dt=detrend(s.v);
s_filt.v = filtfilt(AMfilter.numerator, 1, s_dt);
end

function hr = EstimateHeartRate(sig, fs, up)

% Estimate PSD
blackman_window = blackman(length(sig), 'periodic');
[pxx,f] = periodogram(sig, blackman_window,length(sig), fs);
% [ph, fh] = harmonicPSD(pxx,f);
ph = pxx; fh = f;

% Extract HR
rel_els = fh >= up.fl_hz & fh <= up.fh_hz;
rel_p = ph(rel_els);
rel_f = fh(rel_els);
[~,max_el] = max(rel_p);
hr = rel_f(max_el)*60;

end

function [ph, fh] = harmonicPSD(pxx,f)

% settings
n = 11;
alpha = 2;

ph = nan(length(f),1);
for freq_no = 1 : length(f)
    temp = 0;
    for k = 1 : n
        if freq_no*k > length(f)
            harmonic_p = 0;
        else
            harmonic_p = pxx(freq_no*k);
        end
        temp = temp + min([alpha*pxx(freq_no), harmonic_p]);
    end
    ph(freq_no) = temp;
end

fh = f;

end

function deriv = EstimateDeriv(sig)

% Savitzky Golay
deriv_no = 1;
win_size = 5;
deriv = savitzky_golay_abd(sig, deriv_no, win_size);

end

function deriv = savitzky_golay_abd(sig, deriv_no, win_size)

%% assign coefficients
% From: https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter#Tables_of_selected_convolution_coefficients
% which are calculated from: A., Gorry (1990). "General least-squares smoothing and differentiation by the convolution (Savitzky?Golay) method". Analytical Chemistry. 62 (6): 570?3. doi:10.1021/ac00205a007.

switch deriv_no
    case 0
        % - smoothing
        switch win_size
            case 5
                coeffs = [-3, 12, 17, 12, -3];
                norm_factor = 35;
            case 7
                coeffs = [-2, 3, 6, 7, 6, 3, -2];
                norm_factor = 21;
            case 9
                coeffs = [-21, 14, 39, 54, 59, 54, 39, 14, -21];
                norm_factor = 231;
            otherwise
                error('Can''t do this window size')
        end
    case 1
        % - first derivative
        switch win_size
            case 5
                coeffs = -2:2;
                norm_factor = 10;
            case 7
                coeffs = -3:3;
                norm_factor = 28;
            case 9
                coeffs = -4:4;
                norm_factor = 60;
            otherwise
                error('Can''t do this window size')
        end
        
    case 2
        % - second derivative
        switch win_size
            case 5
                coeffs = [2,-1,-2,-1,2];
                norm_factor = 7;
            case 7
                coeffs = [5,0,-3,-4,-3,0,5];
                norm_factor = 42;
            case 9
                coeffs = [28,7,-8,-17,-20,-17,-8,7,28];
                norm_factor = 462;
            otherwise
                error('Can''t do this window size')
        end
        
    case 3
        % - third derivative
        switch win_size
            case 5
                coeffs = [-1,2,0,-2,1];
                norm_factor = 2;
            case 7
                coeffs = [-1,1,1,0,-1,-1,1];
                norm_factor = 6;
            case 9
                coeffs = [-14,7,13,9,0,-9,-13,-7,14];
                norm_factor = 198;
            otherwise
                error('Can''t do this window size')
        end
        
    case 4
        % - fourth derivative
        switch win_size
            case 7
                coeffs = [3,-7,1,6,1,-7,3];
                norm_factor = 11;
            case 9 
                coeffs = [14,-21,-11,9,18,9,-11,-21,14];
                norm_factor = 143;
            otherwise
                error('Can''t do this window size')
        end
        
    otherwise
        error('Can''t do this order of derivative')        
end

if rem(deriv_no, 2) == 1
    coeffs = -1*coeffs;
end

A = [1,0];
filtered_sig = filter(coeffs, A, sig);
s=length(sig);
half_win_size = floor(win_size*0.5);
deriv=[filtered_sig(win_size)*ones(half_win_size,1);filtered_sig(win_size:s);filtered_sig(s)*ones(half_win_size,1)];
deriv = deriv/norm_factor;

end

function p4 = find_pulse_peaks(p2,p3)
p4 = nan(length(p2),1);
for k = 1 : length(p2)
    rel_el = find(p3>p2(k),1);
    if ~isempty(rel_el)
        p4(k) = p3(rel_el);
    end
end
p4 = p4(~isnan(p4));
end

function [pc, fn] = IBICorrect(p, m, hr, fs, up)

% Correct peaks' location error due to pre-processing
pc = nan(length(p),1);
for k = 1 : length(p)
    [~,rel_el] = min(abs(m-p(k)));
    pc1(k,1) = m(rel_el);    
end

% Correct false positives
% - identify FPs
d = diff(pc1)/fs;  % interbeat intervals in secs
fp = find_reduced_IBIs(d, median(hr), up);
% - remove FPs
pc2 = pc1(~fp);

% Correct false negatives
d = diff(pc2)/fs;  % interbeat intervals in secs
fn = find_prolonged_IBIs(d, median(hr), up);

pc = pc1;

end

function fn = find_prolonged_IBIs(IBIs, med_hr, up)

IBI_thresh = up.upper_hr_thresh_prop*60/med_hr;
fn = IBIs > IBI_thresh;

end

function fp = find_reduced_IBIs(IBIs, med_hr, up)

IBI_thresh = up.lower_hr_thresh_prop*60/med_hr;
fp = IBIs < IBI_thresh;

end

function [onsets, peaks] = pwd(s)

% The Pulse Wave Delineator described in:
%    BN Li, MC Dong, MI Vai (2010) On an automatic delineator for arterial blood pressure waveforms. Biomedical Signal Processing and Control 5(1) 76-81.
%
% Aavailable at: https://www.mathworks.com/matlabcentral/fileexchange/29484-pulse-waveform-delineator

[onsets,peaks,~] = delineator(s.v,s.fs);
onsets = onsets(:);
peaks = peaks(:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This section is 'Pulse Waveform Delineator' code
%
% from: https://www.mathworks.com/matlabcentral/fileexchange/29484-pulse-waveform-delineator
%
% Copyright (c) 2010, LI Bing Nan
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [onsetp,peakp,dicron] = delineator(abpsig,abpfreq)
% This program is intended to delineate the fiducial points of pulse waveforms
% Inputs:
%   abpsig: input as original pulse wave signals;
%   abpfreq: input as the sampling frequency;
% Outputs:
%   onsetp: output fiducial points as the beginning of each beat;
%   peakp: output fiducial points as systolic peaks;
%   dicron: output fiducial points as dicrotic notches;
% Its delineation is based on the self-adaptation in pulse waveforms, but
% not in the differentials.
% Reference:
%   BN Li, MC Dong & MI Vai (2010) 
%   On an automatic delineator for arterial blood pressure waveforms
%   Biomedical Signal Processing and Control 5(1) 76-81.
% LI Bing Nan @ University of Macau, Feb 2007
%   Revision 2.0.5, Apr 2009
%Initialization
peakIndex=0;
onsetIndex=0;
dicroIndex=0;
stepWin=2*abpfreq;
closeWin=floor(0.1*abpfreq);    %invalide for pulse beat > 200BPM
sigLen=length(abpsig);
peakp=[];
onsetp=[];
dicron=[];
%lowpass filter at first
coh=25;                     %cutoff frequency is 25Hz
coh=coh*2/abpfreq;
od=3;                       %3rd order bessel filter
[B,A]=besself(od,coh);
abpsig=filter(B,A,abpsig);
abpsig=10*abpsig;
abpsig=smooth(abpsig);
%Compute differentials
ttp=diff(abpsig);
diff1(2:sigLen)=ttp;
diff1(1)=diff1(2);
diff1=100*diff1;
clear ttp;
diff1=smooth(diff1);
if sigLen>12*abpfreq
    tk=10;
elseif sigLen>7*abpfreq
    tk=5;
elseif sigLen>4*abpfreq
    tk=2;
else
    tk=1;
end
%Seek avaerage threshold in original signal
if tk>1             %self-learning threshold with interval sampling
    tatom=floor(sigLen/(tk+2));
    for ji=1:tk       %search the slopes of abp waveforms
        sigIndex=ji*tatom;
        tempIndex=sigIndex+abpfreq;
        [tempMin,jk,tempMax,jl]=seeklocales(abpsig,sigIndex,tempIndex);
        tempTH(ji)=tempMax-tempMin;
    end
    abpMaxTH=mean(tempTH);
else
    [tempMin,jk,tempMax,jl]=seeklocales(abpsig,closeWin,sigLen);
    abpMaxTH=tempMax-tempMin;
end
clear j*;
clear t*;
abpMaxLT=0.4*abpMaxTH;
%Seek pulse beats by MinMax method
% diffIndex=1;
diffIndex=closeWin;             %Avoid filter distortion
while diffIndex<sigLen
    tempMin=abpsig(diffIndex);   %Initialization
    tempMax=abpsig(diffIndex);
    tempIndex=diffIndex;
    tpeakp=diffIndex;        %Avoid initial error
    tonsetp=diffIndex;      %Avoid initial error
    while tempIndex<sigLen
        %If no pulses within 2s, then adjust threshold and retry
        if (tempIndex-diffIndex)>stepWin
%             tempIndex=diffIndex-closeWin;
            tempIndex=diffIndex;
            abpMaxTH=0.6*abpMaxTH;
            if abpMaxTH<=abpMaxLT
                abpMaxTH=2.5*abpMaxLT;
            end
            break;
        end
        if (diff1(tempIndex-1)*diff1(tempIndex+1))<=0  %Candidate fiducial points
            if (tempIndex+5)<=sigLen
                jk=tempIndex+5;
            else
                jk=sigLen;
            end
            if (tempIndex-5)>=1
                jj=tempIndex-5;
            else
                jj=1;
            end
            %Artifacts of oversaturated or signal loss?
            if (jk-tempIndex)>=5
                for ttk=tempIndex:jk
                    if diff1(ttk)~=0
                        break;
                    end
                end
                if ttk==jk
                    break;          %Confirm artifacts
                end
            end
            if diff1(jj)<0          %Candidate onset
                if diff1(jk)>0
                    [tempMini,tmin,ta,tb]=seeklocales(abpsig,jj,jk);
                    if abs(tmin-tempIndex)<=2
                        tempMin=tempMini;
                        tonsetp=tmin;
                    end
                end
            elseif diff1(jj)>0      %Candidate peak
                if diff1(jk)<0
                    [tc,td,tempMaxi,tmax]=seeklocales(abpsig,jj,jk);
                    if abs(tmax-tempIndex)<=2
                        tempMax=tempMaxi;
                        tpeakp=tmax;
                    end
                end
            end
            if ((tempMax-tempMin)>0.4*abpMaxTH)   %evaluation
                if ((tempMax-tempMin)<2*abpMaxTH)
                    if tpeakp>tonsetp
                        %If more zero-crossing points, further refine!
                        ttempMin=abpsig(tonsetp);
                        ttonsetp=tonsetp;
                        for ttk=tpeakp:-1:(tonsetp+1)
                            if abpsig(ttk)<ttempMin
                                ttempMin=abpsig(ttk);
                                ttonsetp=ttk;
                            end
                        end
                        tempMin=ttempMin;
                        tonsetp=ttonsetp;
                            
                        if peakIndex>0
                            %If pulse period less than eyeclose, then artifact
                            if (tonsetp-peakp(peakIndex))<(3*closeWin)
                                %too many fiducial points, then reset
                                tempIndex=diffIndex;                                
                                abpMaxTH=2.5*abpMaxLT;
                                break;
                            end
                            
                            %If pulse period bigger than 2s, then artifact
                            if (tpeakp-peakp(peakIndex))>stepWin
                                peakIndex=peakIndex-1;
                                onsetIndex=onsetIndex-1;
                                if dicroIndex>0
                                    dicroIndex=dicroIndex-1;
                                end
                            end
                            if peakIndex>0
                                %new pulse beat
                                peakIndex=peakIndex+1;
                                peakp(peakIndex)=tpeakp;
                                onsetIndex=onsetIndex+1;
                                onsetp(onsetIndex)=tonsetp;
                                tf=onsetp(peakIndex)-onsetp(peakIndex-1);
                                to=floor(abpfreq./20);   %50ms
                                tff=floor(0.1*tf);
                                if tff<to
                                    to=tff;
                                end
                                to=peakp(peakIndex-1)+to;
                                te=floor(abpfreq./2);   %500ms
                                tff=floor(0.5*tf);
                                if tff<te
                                    te=tff;
                                end
                                te=peakp(peakIndex-1)+te;
                                tff=seekdicrotic(diff1(to:te));
                                if tff==0
                                    tff=te-peakp(peakIndex-1);
                                    tff=floor(tff/3);
                                end
                                dicroIndex=dicroIndex+1;
                                dicron(dicroIndex)=to+tff;
                                tempIndex=tempIndex+closeWin;
                                break;
                            end
                        end
                        
                        if  peakIndex==0   %new pulse beat
                            peakIndex=peakIndex+1;
                            peakp(peakIndex)=tpeakp;
                            onsetIndex=onsetIndex+1;
                            onsetp(onsetIndex)=tonsetp;
                            tempIndex=tempIndex+closeWin;
                            break;
                        end
                    end
                end
            end
        end
        tempIndex=tempIndex+1;      %step forward
    end
%     diffIndex=tempIndex+closeWin;    %for a new beat
    diffIndex=tempIndex+1;
end
if isempty(peakp),return;end
%Compensate the offsets of lowpass filter
sigLen=length(peakp);
for diffIndex=1:sigLen          %avoid edge effect
    tempp(diffIndex)=peakp(diffIndex)-od;
end
ttk=tempp(1);
if ttk<=0
    tempp(1)=1;
end 
clear peakp;
peakp=tempp;
clear tempp;
sigLen=length(onsetp);
for diffIndex=1:sigLen
    tempp(diffIndex)=onsetp(diffIndex)-od;
end
ttk=tempp(1);
if ttk<=0
    tempp(1)=1;
end 
clear onsetp;
onsetp=tempp;
clear tempp;
if isempty(dicron),return;end
sigLen=length(dicron);
for diffIndex=1:sigLen
    if dicron(diffIndex)~=0
        tempp(diffIndex)=dicron(diffIndex)-od;
    else
        tempp(diffIndex)=0;
    end
end
clear dicron;
dicron=tempp;
clear tempp;

end

function [mini,minip,maxi,maxip]=seeklocales(tempsig,tempbegin,tempend)
tempMin=tempsig(tempbegin);
tempMax=tempsig(tempbegin);
minip=tempbegin;
maxip=tempbegin;
for j=tempbegin:tempend
    if tempsig(j)>tempMax
        tempMax=tempsig(j);
        maxip=j;
    elseif tempsig(j)<tempMin
        tempMin=tempsig(j);
        minip=j;
    end
end
mini=tempMin;
maxi=tempMax;

end

function [dicron]=seekdicrotic(tempdiff)
izcMin=0;
izcMax=0;
itemp=3;
tempLen=length(tempdiff)-3;
dicron=0;
tempdiff=smooth(tempdiff);
while itemp<=tempLen
    if (tempdiff(itemp)*tempdiff(itemp+1))<=0
        if tempdiff(itemp-2)<0
            if tempdiff(itemp+2)>=0
                izcMin=izcMin+1;
                tzcMin(izcMin)=itemp;
            end
        end
%         if tempdiff(itemp-2)>0
%             if tempdiff(itemp+2)<=0
%                 izcMax=izcMax+1;
%                 tzcMax(izcMax)=itemp;
%             end
%         end
    end
    itemp=itemp+1;
end
if izcMin==0     %big inflection
    itemp=3;
    tempMin=tempdiff(itemp);
    itempMin=itemp;
    
    while itemp<tempLen
        if tempdiff(itemp)<tempMin
            tempMin=tempdiff(itemp);
            itempMin=itemp;
        end
        itemp=itemp+1;
    end
    itemp=itempMin+1;
    while itemp<tempLen
        if tempdiff(itemp+1)<=tempdiff(itemp-1)
            dicron=itemp;
            return;
        end
        itemp=itemp+1;
    end
elseif izcMin==1
    dicron=tzcMin(izcMin);
    return;
else
    itemp=tzcMin(1);
    tempMax=tempdiff(itemp);
    itempMax=itemp;
    
    while itemp<tempLen
        if tempdiff(itemp)>tempMax
            tempMax=tempdiff(itemp);
            itempMax=itemp;
        end
        itemp=itemp+1;
    end
    for itemp=izcMin:-1:1
        if tzcMin(itemp)<itempMax
            dicron=tzcMin(itemp);
            return;
        end
    end
end

end

function [diap]=seekdiap(tempabp)
diap=0;
[tt,ti]=max(tempabp);
if ti==0
    diap=floor(length(tempabp)./2);
else
    diap=ti;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function peaks = wfd(s)

v = 0; % set to 1 to enable visualiation
peaks = wavelet(s.v,s.fs,v);

end

function ibis = wavelet(ppg,fs,v)
%% Detects inter-beat intervals using a Wavelet approach
% Citation: Conn NJ, Borkholder DA (2013) Wavelet based photoplethysmogram foot delineation for heart rate 
% variability applications. IEEE Signal Processing in Medicine and Biology Symposium (SPMB) 2013:1-5. 
% doi: 10.1109/SPMB.2013.6736782.
%
% Inputs:   PPG signal, ppg
%           sampling rate, fs
%           visualization option, v (1: plot, 0: don't plot)
% Outputs:  position of the starting points of inter-beat intervals, in number of samples, ibis
% 
% Developed by: Elisa Mejía Mejía
%               City, University of London
% Version:      1.0 -   May, 2021

    %% Beat segmentation
    b = fir1(round(fs/10),[0.5 8]/(fs/2),'bandpass');   % Designs the filter
    ppgf = filtfilt(b,1,ppg);                           % Filters the signal   
    
    if v == 1
        color = [   0, 0.4470, 0.7410; ...
                    0.8500, 0.3250, 0.0980; ...
                    0.9290, 0.6940, 0.1250; ...
                    0.4940, 0.1840, 0.5560; ...
                    0.4660, 0.6740, 0.1880; ...
                    0.3010, 0.7450, 0.9330; ...
                    0.6350, 0.0780, 0.1840; ...
                    0.75, 0, 0.75; ...
                    0.5, 0.5, 0.5];
                
        figure('position',[240 140 1200 625],'name','Wavelet');
        ax(1) = subplot(2,1,1);
        plot(ppg,'k','linewidth',1.2);
        hold on;
        xlim([1 length(ppg)]);
        title('Peak detection');
        xlabel('Samples');
        ylabel('PPG');
        grid on;
        box off;
        
        ax(2) = subplot(2,1,2);
        yyaxis left;
        plot(ppgf,'-k','linewidth',1.2);
        hold on;
        xlim([1 length(ppg)]);
        title('Wavelet');
        xlabel('Samples');
        ylabel('PPG');
        grid on;
        box off;  
        
        linkaxes(ax,'x'); 
    end
    
    fsi = 250;                                          % Defines the interpolation frequency
    ti = (0:length(ppgf) - 1)/fs;                       % Generates the original time array
    tf = ti(1):1/fsi:ti(end);                           % Generates the time array for the spline
    ppgss = spline(ti,ppgf,tf);                         % Subsamples the PPG signal to 250 Hz
    
    [~,q] = wavelet_filters(length(ppgss),fsi,0);       % Generates the Wavelet filters at 250 Hz
    w = wavelet_decomposition(ppgss,q,0);               % Applies the Wavelet decomposition to the 
                                                        % subsampled signal
    ppgw = w{5};                                        % Takes the fifth Wavelet scale
    ppgw = spline(tf,ppgw,ti);                          % Restores sampling rate
    
    if v == 1
        ax(2) = subplot(2,1,2);
        yyaxis left;
        plot(ppgw,'-','color',color(2,:),'linewidth',1.2);
        xlim([1 length(ppgw)]);
        
        linkaxes(ax,'x');  
    end
    
    ppgrs = ppgw;                                       % Creates a copy of the signal
    ppgrs(ppgrs < 0) = 0;                               % Rectifies the signal
    ppgrs = ppgrs.^2;                                   % Squares the rectified signal
    
    if v == 1
        ax(2) = subplot(2,1,2);
        yyaxis left;
        plot(ppgrs,'-','color',color(3,:),'linewidth',1.2);
        xlim([1 length(ppgrs)]); 
        
        linkaxes(ax,'x');  
    end
    
    b = fir1(round(length(ppgrs)/10),1/(fs/2),'low');   % Designs a Hamming FIR filter
    th = filtfilt(b,1,ppgrs);                           % Applies the Hamming FIR filter
    
    if v == 1
        ax(2) = subplot(2,1,2);
        yyaxis left;
        plot(th,'--','color',color(4,:),'linewidth',1.2);
        xlim([1 length(th)]);
        
        linkaxes(ax,'x');  
    end
    
    zc = [0 diff(ppgrs >= th)];                         % Finds the points where the Wavelet signal crosses the
                                                        % threshold
    zc_up = find(zc == 1);                              % Finds the left boundaries of the areas of interest
    zc_down = find(zc == -1);                           % Finds the right boundaries of the areas of interest
    if zc_up(1) > zc_down(1)                            % Verifies if the first left boundary is higher than
                                                        % the first right boundary
        zc_down(1) = [];                                % Removes the first right boundary
    end
    if zc_up(end) > zc_down(end)                        % Verifies if the last left boundary is higher than
                                                        % the last right boundary
        zc_up(end) = [];                                % Removes the last left boundary
    end
    aoi = [zc_up' zc_down'];                            % Determines the areas of interest
        
    %% Determination of PPG foot location
    b = firls(51,[0 0.9],[0 0.9*pi],'differentiator');  % Designs the differentiator FIR filter
    ppg1d = filtfilt(b,1,ppgf);                         % Obtaines the 1st derivative
    ppg2d = filtfilt(b,1,ppg1d);                        % Obtaines the 2nd derivative
    ppg3d = filtfilt(b,1,ppg2d);                        % Obtaines the 3rd derivative
    
    if v == 1
        ax(2) = subplot(2,1,2);
        yyaxis right;
        plot(ppg1d/max(ppg1d),':','color',color(7,:),'linewidth',1.2);
        plot(ppg2d/max(ppg2d),':','color',color(8,:),'linewidth',1.2);
        plot(ppg3d/max(ppg3d),'-','color',color(1,:),'linewidth',1.2);
        xlim([1 length(ppg3d)]);
        
        linkaxes(ax,'x');   
    end
    
    ons = zeros(size(aoi,1),1);                         % Initializes variable
    for i = 1:size(aoi,1)                               % Iterates through the areas of interest
        aux_ppg3d = ppg3d(aoi(i,1):aoi(i,2));           % Takes the portion of the 3rd derivative in the i-th AOI
        j = 2;                                          % Initializes counter
        zc = 0;                                         % Initializes variable
        
        while j <= length(aux_ppg3d)                    % Iterates through the portion of the signal
            if aux_ppg3d(j - 1) <= 0 && aux_ppg3d(j) > 0
                                                        % Verifies if the point corresponds to a zero-crossing
                zc = j;                                 % Updates variables
                j = length(aux_ppg3d);                  % Updates counter
            end
            j = j + 1;                                  % Updates counter
        end
        
        if zc ~= 0
            zc = zc + aoi(i,1) - 1;                     % Adds offset
            if ppg2d(zc) >= 0 && ppg1d(zc) > 0          % Verifies if the 1st and 2nd derivatives are positive at the point
                ons(i) = zc;
            else
                [~,ind] = max(ppg2d(aoi(i,1):aoi(i,2)));% Finds the location of the maximum of the second 
                                                        % derivative in the i-th AOI
                ons(i) = ind + aoi(i,1) - 1;            % Adds offset
            end
        else
            [~,ind] = max(ppg2d(aoi(i,1):aoi(i,2)));    % Finds the location of the maximum of the second 
                                                        % derivative in the i-th AOI
            ons(i) = ind + aoi(i,1) - 1;                % Adds offset
        end
    end
    
    if v == 1
        ax(1) = subplot(2,1,1);
        plot(ons,ppg(ons),'ok','MarkerFaceColor',color(1,:));
        
        ax(2) = subplot(2,1,2);
        yyaxis left;
        plot(ons,ppgf(ons),'ok','MarkerFaceColor',color(9,:));
        xlim([1 length(ppg3d)]);
        
        linkaxes(ax,'x');   
    end
    
    if v == 1
        ax(1) = subplot(2,1,1);
        miny = min(ppg);
        if miny < 0
            miny = 1.2*miny;
        else
            miny = 0.8*miny;
        end
        maxy = max(ppg);
        if maxy < 0
            maxy = 0.8*maxy;
        else
            maxy = 1.2*maxy;
        end
        yline = linspace(miny,maxy,1000);
        plot(aoi(1,1).*ones(size(yline)),yline,'--','color',color(5,:),'linewidth',1.2);
        plot(aoi(1,2).*ones(size(yline)),yline,'--','color',color(6,:),'linewidth',1.2);        
        plot(aoi(2:end,1).*ones(size(yline)),yline,'--','color',color(5,:),'linewidth',1.2);
        plot(aoi(2:end,2).*ones(size(yline)),yline,'--','color',color(6,:),'linewidth',1.2);
        legend('PPG','Onset','AOI: Left boundaries','AOI: Right boundaries','location','northeastoutside');
        
        ax(2) = subplot(2,1,2);
        yyaxis right;
        miny = min([min(ppgf) min(ppgw) min(ppgrs)]);
        if miny < 0
            miny = 1.2*miny;
        else
            miny = 0.8*miny;
        end
        maxy = max([max(ppgf) max(ppgw) max(ppgrs)]);
        if maxy < 0
            maxy = 0.8*maxy;
        else
            maxy = 1.2*maxy;
        end
        yline = linspace(miny,maxy,1000);
        plot(aoi(1,1).*ones(size(yline)),yline,'--','color',color(5,:),'linewidth',1.2);
        plot(aoi(1,2).*ones(size(yline)),yline,'--','color',color(6,:),'linewidth',1.2);
        plot(aoi(2:end,1).*ones(size(yline)),yline,'--','color',color(5,:),'linewidth',1.2);
        plot(aoi(2:end,2).*ones(size(yline)),yline,'--','color',color(6,:),'linewidth',1.2);
        legend( 'Filtered PPG','Wavelet scale','Rectified and squared','Threshold','Onset', ...
                '1D PPG','2D PPG','3D PPG','AOI: Left boundaries','AOI: Right boundaries',...
                'location','northeastoutside');
        
        linkaxes(ax,'x');  
    end
    ibis = ons;                                         % Sets the output
end

function h = waveleth(w)
%% Constructs the LPF required for the wavelet-based delineator at a sampling frequency of 250 Hz
% Inputs:   array with the frequency points, in radians, that will be used to construct the filter. 
%           Must be between 0 and 2pi
% Outputs:  array with the LPF coefficients

    aux1 = exp((1j)*w/2);
    aux2 = cos(w/2).^3;
    h = aux1.*aux2;

end

function g = waveletg(w)
%% Constructs the HPF required for the wavelet-based delineator at a sampling frequency of 250 Hz
% Inputs:   array with the frequency points, in radians, that will be used to construct the filter. 
%           Must be between 0 and 2pi
% Outputs:  array with the HPF coefficients

    aux1 = (4*(1j));
    aux2 = exp((1j)*w/2);
    aux3 = sin(w/2);
    g = (aux1*aux2).*aux3;

end

function [q, q250] = wavelet_filters(n,fs,v)
%% Creates the filters required to make the wavelet decomposition using the algorithme-a-trous
% Inputs:   number of samples of the signal that will be decomposed, n
%           sampling frequency of the signal that will be decomposed, fs
%           visualization option, v
% Outputs:  cell object with the five filters required to make the wavelet decomposition, q

    %% Generation of filters
    m = 250*n/fs;
    w = 0:2*pi/m:2*pi - 2*pi/m;
    
    q = cell(5,1);
    q{1} = waveletg(w);
    if v == 1
        figure;
        subplot(1,2,1);
        f = (0:length(q{1}) - 1)*(250/length(q{1}));
        plot(f,abs(q{1}));
        hold on;
        grid on; box off;
%         xlim([0 floor(250/2)]);
    end
    for k = 2:5
        g = waveleth(w);
        h = 1;
        while h < k - 1
            g = g.*waveleth((2^h).*w);
            h = h + 1;
        end
        g = waveletg(((2^(k - 1))*w)).*g;
        q{k} = g;
        if v == 1
            subplot(1,2,1);
            f = (0:length(q{k}) - 1)*(250/length(q{k}));
            plot(f,abs(q{k}));
            hold on;
            grid on; box off;
%           xlim([0 floor(250/2)]);
        end
    end
    q250 = q;
    
    
    %% Resampling to fs
    for k = 1:length(q)
        aux = real(ifft(q{k}));
        len = length(aux);
        xi = (0:length(aux) - 1)/250;
        xf = 0:1/fs:len/250; xf(end) = [];
        aux = spline(xi,aux,xf);
        q{k} = fft(aux); 
    end
    
    if v == 1
        subplot(1,2,2);
        hold on;
        for k = 1:length(q)
            f = (0:length(q{k}) - 1)*(fs/length(q{k}));
            plot(f,abs(q{k}));
        end
        grid on; box off;
%         xlim([0 fs/2]);
    end
end

function w = wavelet_decomposition(x,q,v)
%% Performs the Wavelet decomposition of a signal using the algorithme-a-trous
% Inputs:   signal to be decomposed, x
%           cell object containing the filters that will decompose the signal, q
%           visualization option, v
% Outputs:  cell object containing the Wavelet decomposition of the signal
%           at scales 2^1, 2^2, ..., 2^5, w

    %% Application of the Wavelet filters
    if v == 1
        figure('position',[250   90  980  680]);
        ax(1) = subplot(6,1,1);
        plot(x);
        grid on; box off;
        title('Original signal');
        linkaxes(ax,'x');
    end
    
    w = cell(size(q));
    for i = 1:length(q)
        aux = fft(x);
        aux_q = q{i};
        if size(aux,1) > size(aux,2)
            aux_q = aux_q';
        end
        aux = aux.*aux_q;
        aux = ifft(aux);
        w{i} = real(aux);
        
        if v == 1
            ax(i + 1) = subplot(6,1,i + 1);
            plot(w{i});
            grid on; box off;
            title(strcat(['Scale ' num2str(i)]));
            linkaxes(ax,'x');
        end
    end
    
end

function peaks = pda(s)

v = 0; % set to 1 to enable visualiation
peaks = upslopes(s.v,v);

end

function ibis = upslopes(ppg,v)
%% Detects inter-beat intervals using the upslopes
% Citation: Argüello Prada EJ, Serna Maldonado RD. (2018 ) A novel and low-complexity peak detection algorithm for 
% heart rate estimation from low-amplitude photoplethysmographic (PPG) signals. J Med Eng Technol 42(8):569-577. 
% doi: 10.1080/03091902.2019.1572237. Epub 2019 Mar 28. PMID: 30920315.
%
% Inputs:   PPG signal, ppg
%           visualization option, v (1: plot, 0: don't plot)
% Outputs:  position of the starting points of inter-beat intervals, in number of samples, ibis
% 
% Developed by: Elisa Mejía Mejía
%               City, University of London
% Version:      1.0 -   May, 2021

    %% Detection of peaks 
    if v == 1
        color = [   0, 0.4470, 0.7410; ...
                    0.8500, 0.3250, 0.0980; ...
                    0.9290, 0.6940, 0.1250; ...
                    0.4940, 0.1840, 0.5560; ...
                    0.4660, 0.6740, 0.1880; ...
                    0.3010, 0.7450, 0.9330; ...
                    0.6350, 0.0780, 0.1840; ...
                    0.75, 0, 0.75];
                
        figure('position',[240 140 1200 625],'name','Upslopes');
        ax(1) = subplot(2,1,1);
        plot(ppg,'k','linewidth',1.2);
        hold on;
        xlim([1 length(ppg)]);
        title('Peak detection');
        xlabel('Samples');
        ylabel('PPG');
        grid on;
        box off;
        
        ax(2) = subplot(2,1,2);
        plot(ppg,'k','linewidth',1.2);
        hold on;
        xlim([1 length(ppg)]);
        title('Upslopes');
        xlabel('Samples');
        ylabel('PPG');
        grid on;
        box off;  
        linkaxes(ax,'x');   
    end
    
    th = 6;                                     % Initializes threshold
    pks = [];                                   % Initializes variable  
    pos_peak = [];                              % Initializes variable
    pos_peak_b = 0;                             % Initializes variable
    n_pos_peak = 0;                             % Initializes counter
    n_up = 0;                                   % Initializes counter
    for i = 2:length(ppg)                       % Iterates through signal
        if ppg(i) > ppg(i - 1)                  % Verifies if it is the upslope
            n_up = n_up + 1;                    % Adds one to the coutner
            
            if v == 1
                ax(2) = subplot(2,1,2);
                plot(i,ppg(i),'.','color',color(5,:));
                linkaxes(ax,'x');
            end
        else
            if n_up >= th                       % Checks if the number of data in the upslope is higher than the
                                                % threshold
                pos_peak = [pos_peak; i];       % Adds the position to the possible peaks 
                pos_peak_b = 1;                 % Sets the value as 1 for the detection of a new possible peak
                n_pos_peak = n_pos_peak + 1;    % Refreshes counter
                n_up_pre = n_up;                % Stores the previous value of n_up
                
                if v == 1
                    ax(2) = subplot(2,1,2);
                    plot(pos_peak(n_pos_peak),ppg(pos_peak(n_pos_peak)),'ok','MarkerFaceColor',color(1,:));
                    linkaxes(ax,'x');
                end
            else
                if pos_peak_b == 1              % Verifies if a peak has been found
                    if ppg(i - 1) > ppg(pos_peak(n_pos_peak))
                                                % Verifies if the previous sample is higher than
                                                % the one selected as peak
                        pos_peak(n_pos_peak) = i - 1;
                                                % Updates the value of the position of the possible
                                                % peak
                        if v == 1
                            ax(2) = subplot(2,1,2);
                            plot(pos_peak(n_pos_peak),ppg(pos_peak(n_pos_peak)),'ok','MarkerFaceColor',color(2,:));
                            linkaxes(ax,'x');
                        end
                    else
                        pks = [pks; pos_peak(n_pos_peak)];
                                                % Refreshes array of peaks
                        if v == 1
                            ax(2) = subplot(2,1,2);
                            plot(pks(end),ppg(pks(end)),'ok','MarkerFaceColor',color(3,:));
                            linkaxes(ax,'x');
                        end
                    end
                    th = 0.6*n_up_pre;          % Updates value of threshold
                    pos_peak_b = 0;             % Updates value of variable
                end
            end
            n_up = 0;                           % Resets counter
        end
    end
    
    if v == 1
        ax(1) = subplot(2,1,1);
        plot(pks,ppg(pks),'ok','MarkerFaceColor',color(1,:));
        linkaxes(ax,'x');
    end
    
    ibis = pks;                                 % Sets the output
    
end

function peaks = erma(s)

%% Detects inter-beat intervals using D2Max
% Citation: Elgendi M, Norton I, Brearley M, Abbott D, Schuurmans D (2013) Systolic Peak Detection in Acceleration
% Photoplethysmograms Measured from Emergency Responders in Tropical Conditions. PLoS ONE 8(10): e76585.
% doi:10.1371/journal.pone.0076585
%
% Inputs:   PPG signal, ppg
%           sampling rate, fs
%           visualization option, v (1: plot, 0: don't plot)
% Outputs:  position of the starting points of inter-beat intervals, in number of samples, ibis
%
% Developed by: Elisa Mejía Mejía
%               City, University of London
%
% Adapted by: Peter H. Charlton
%
% Version:      1.1 -   Nov, 2021

%% setup 
% (inserted by PHC)
ppg = s.v;
fs = s.fs;
clear s
v = 0; % set to 1 to enable visualiation

%% Bandpass filter
% if length(ppg) < 4098                                   % Removed by PC
%     ppg = [ppg; zeros(4098 - length(ppg) + 1,1)];
% end
[b,a] = butter(2,[0.5, 8]/(fs/2),'bandpass');   % Designs the bandpass, butterworth filter
s = filtfilt(b,a,ppg);                          % Applies the zero-phase filter
z = s;                                          % Creates a new signal based on the filtered signal
z(z < 0) = 0;                                   % Clips the signal

if v == 1
    color = [   0, 0.4470, 0.7410; ...
        0.8500, 0.3250, 0.0980; ...
        0.9290, 0.6940, 0.1250; ...
        0.4940, 0.1840, 0.5560; ...
        0.4660, 0.6740, 0.1880; ...
        0.3010, 0.7450, 0.9330; ...
        0.6350, 0.0780, 0.1840; ...
        0.75, 0, 0.75];
    
    figure('position',[240 140 1200 625],'name','D2Max');
    ax(1) = subplot(2,1,1);
    plot(ppg,'k','linewidth',1.2);
    hold on;
    xlim([1 length(ppg)]);
    title('Peak detection');
    xlabel('Samples');
    ylabel('PPG');
    grid on;
    box off;
    
    ax(2) = subplot(2,1,2);
    plot(s,'k','linewidth',1.2);
    hold on;
    plot(z,'color',color(1,:),'linewidth',1.2);
    xlim([1 length(s)]);
    title('D2Max');
    xlabel('Samples');
    ylabel('Filtered PPG');
    grid on;
    box off;
    linkaxes(ax,'x');
end

%% Squaring
y = z.^2;                                       % Squares the filtered, clipped signal

if v == 1
    ax(2) = subplot(2,1,2);
    plot(y,'color',color(2,:),'linewidth',1.2);
    linkaxes(ax,'x');
end

%% Generating blocks of interest
w1 = (111e-3)*fs;                               % Determines the order of a first moving average
w1 = 2*floor(w1/2) + 1;                         % Rounds the order to the nearest odd integer
b = (1/w1)*ones(w1,1);                          % Designs the first moving average filter
ma_peak = filtfilt(b,1,y);                      % Applies the first moving average filter

w2 = (667e-3)*fs;                               % Determines the order of a second moving average
w2 = 2*floor(w2/2) + 1;                         % Rounds the order to the nearest odd integer
b = (1/w2)*ones(w2,1);                          % Designs the second moving average filter
ma_beat = filtfilt(b,1,y);                      % Applies the second moving average filter

if v == 1
    ax(2) = subplot(2,1,2);
    plot(ma_peak,'--','color',color(3,:),'linewidth',1.2);
    plot(ma_beat,'--','color',color(4,:),'linewidth',1.2);
    linkaxes(ax,'x');
end

%% Thresholding
alpha = 0.02*mean(y);                           % Determines offset level
th1 = ma_beat + alpha;                          % Determines a first threshold based on crossing of ma_beat
boi = ma_peak > th1;                            % Determines blocks of interest according to threshold 1

if v == 1
    ax(2) = subplot(2,1,2);
    plot(boi,'color',color(5,:),'linewidth',1.2);
    linkaxes(ax,'x');
end

th2 = w1;                                       % Determines a second threshold based on length of the
% blocks of interest
pos_blocks_init = find(diff(boi) > 0);          % Determines the location of the starting points for blocks
pos_blocks_init = pos_blocks_init + 1;          % Adds one for the delay
pos_blocks_end = find(diff(boi) < 0);           % Determines the location of the ending points for blocks
pos_blocks_end = pos_blocks_end + 1;            % Adds one for the delay
if pos_blocks_init(1) > pos_blocks_end(1)       % Verifies if the first block does not have a starting point
    pos_blocks_init = [1, pos_blocks_init];     % Adds a initial point for the first block
end
if pos_blocks_init(end) > pos_blocks_end(end)   % Verifies if the last block does not have an ending point
    pos_blocks_end = [pos_blocks_end, length(y)];
    % Adds an ending point for the last block
end

len_blocks = zeros(size(pos_blocks_init));      % Initializes variable
ibis = zeros(size(pos_blocks_init));            % Initializes variable
for i = 1:length(pos_blocks_init)               % Iterates through the blocks
    ind = find(pos_blocks_end > pos_blocks_init(i),1);
    % Determines the first ending position occurring after
    % the i-th initial position
    len_blocks(i) = pos_blocks_end(ind) - pos_blocks_init(i);
    % Measures the size of the block
    if len_blocks(i) >= th2                     % Verifies if the block is long enough
        [~,max_block] = max(ppg(pos_blocks_init(i):pos_blocks_end(ind)));
        % Detects the maximum in the signal within the block
        ibis(i) = max_block + pos_blocks_init(i) - 1;
        % Adds the offset and stores the results
        
    end
end
ind = find(len_blocks < th2);                   % Finds the blocks that are not long enough
if ~isempty(ind)                                % Checks if at least one block needs to be deleted
    for i = 1:length(ind)                       % Iterates through the blocks to be deleted
        boi(pos_blocks_init(i):pos_blocks_end(i)) = 0;
    end
end
ibis(ibis == 0) = [];                           % Removes the zeros from the array

if v == 1
    ax(1) = subplot(2,1,1);
    plot(ibis,ppg(ibis),'ok','MarkerFaceColor',color(1,:));
    ax(2) = subplot(2,1,2);
    plot(boi,'--','color',color(6,:),'linewidth',1.2);
    plot(ibis,s(ibis),'ok','MarkerFaceColor',color(7,:));
    linkaxes(ax,'x');
end

%% End matter
% (inserted by PHC)
peaks = ibis;

end

%%%%%%% start of heartpy_pc

function peaks = heartpy_pc(s)

% A Matlab implementation of HeartPy based on the description in:
%
%  van Gent P. et al., Analysing noisy driver physiology real-time using off-the-shelf sensors: Heart rate analysis software from the taking the fast lane project. J. Open Res. Softw. 2019, 7, doi:10.5334/jors.241.
%
% and v.1.2.5 of HeartPy downloaded on 7-June-2021 from: https://github.com/paulvangentcom/heartrate_analysis_python 

%% Setup
hrdata = s.v;
sample_rate = s.fs;
iterations = 2;

%% Pre-processing
% Peak enhancement
hrdata = enhance_peaks(hrdata);

%% Peak detection
% Process data
working_data = process(hrdata, sample_rate);

%% Output peaks
% Only output those that are kept
peaks = working_data.peaklist(working_data.binary_peaklist);

end

function working_data = process(hrdata, sample_rate)

% adapted from the 'process' function within 'heartpy.py'

% constants
bpmmin=40;
bpmmax=180;
windowsize=0.75;
reject_segmentwise=false;

working_data.hr = hrdata;
working_data.sample_rate = sample_rate;

% Calculate moving average
rol_mean = rolling_mean(hrdata, windowsize, sample_rate);

% Identify peaks
working_data = fit_peaks(hrdata, rol_mean, sample_rate, bpmmin, bpmmax);

% Calculate peak-to-peak intervals
working_data = calc_rr(working_data.peaklist, sample_rate, working_data);

% Check peaks
working_data = check_peaks(working_data.RR_list, working_data.peaklist, working_data.ybeat, reject_segmentwise, working_data);

% reminaing functions appear to be analysing peak-to-peak intervals, rather than identifying peaks, so not implemented here.

end

function working_data = fit_peaks(hrdata, rol_mean, sample_rate, bpmmin, bpmmax)
% from the 'fit_peaks' function within 'peakdetection.py'

% moving average values to test
ma_perc_list = [5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 150, 200, 300];

[rrsd.rrsd, rrsd.bpm, rrsd.ma_perc] = deal([]);
[valid_ma.rrsd, valid_ma.ma_perc] = deal([]);

update_dict = true;

for ma_val_no = 1 : length(ma_perc_list)
    ma_perc = ma_perc_list(ma_val_no);

    working_data = detect_peaks(hrdata, rol_mean, ma_perc, sample_rate, update_dict);
    bpm = ((length(working_data.peaklist)/(length(hrdata)/sample_rate))*60);
    
    rrsd.rrsd(end+1,1) = working_data.rrsd;
    rrsd.bpm(end+1,1) = bpm;
    rrsd.ma_perc(end+1,1) = ma_perc;

end

for ma_val_no = 1 : length(ma_perc_list)
    if rrsd.rrsd(ma_val_no) > 0.1 && rrsd.bpm(ma_val_no) >= bpmmin && rrsd.bpm(ma_val_no) <= bpmmax
        valid_ma.rrsd(end+1,1) = rrsd.rrsd(ma_val_no);
	    valid_ma.ma_perc(end+1,1) = rrsd.ma_perc(ma_val_no);
    end
end

if ~isempty(valid_ma.rrsd)
    [~, min_rrsd_el] = min(valid_ma.rrsd);
    working_data.best = valid_ma.ma_perc(min_rrsd_el);
    working_data = detect_peaks(hrdata, rol_mean, working_data.best, sample_rate, update_dict, working_data);
else
    warning('\n----------------\nCould not determine best fit for given signal. Please check the source signal.\n Probable causes:\n- detected heart rate falls outside of bpmmin<->bpmmax constraints\n- no detectable heart rate present in signal\n- very noisy signal (consider filtering and scaling)\nIf you\''re sure the signal contains heart rate data, consider filtering and/or scaling first.\n----------------\n');
end

end

function working_data = check_peaks(rr_arr, peaklist, ybeat, reject_segmentwise, working_data)
% excludes outlying peaks
%  an implementation of 'check_peaks' from 'peakdetection.py' in HeartPy.

% define RR range as mean +/- 30%, with a minimum of 300
mean_rr = mean(rr_arr);
thirty_perc = 0.3 * mean_rr;
if thirty_perc <= 300
    upper_threshold = mean_rr + 300;
    lower_threshold = mean_rr - 300;
else
    upper_threshold = mean_rr + thirty_perc;
    lower_threshold = mean_rr - thirty_perc;
end

% identify peaks to exclude based on RR interval
rem_idx = find((rr_arr <= lower_threshold) | (rr_arr >= upper_threshold));

working_data.removed_beats = peaklist(rem_idx);
working_data.removed_beats_y = ybeat(rem_idx);
working_data.binary_peaklist = true(size(peaklist));
working_data.binary_peaklist(rem_idx) = false;

if reject_segmentwise
    working_data = check_binary_quality(peaklist, working_data.binary_peaklist, working_data);
end

working_data = update_rr(working_data);

end

function working_data = check_binary_quality(peaklist, binary_peaklist, working_data)

% Constant
max_rejects = 3;

% The whole segment is rejected as it contains more than the specified 3 rejections per 10 beats.
idx = 0;
working_data.rejected_segments = [];
for i = 1: floor(length(binary_peaklist) / 10)
    if sum(binary_peaklist(idx:idx + 10)==0) > maxrejects
        binary_peaklist(idx:idx + 10) = 0;
        if idx + 10 < length(peaklist)
            working_data.rejected_segments(end+1,:) = [peaklist(idx), peaklist(idx + 10)];
        else
            working_data.rejected_segments(end+1,:) = [peaklist(idx), peaklist(end)];
        end
    end
    idx = idx+10;
end

end

function working_data = update_rr(working_data)
% Function that updates RR differences and RR squared differences based on corrected RR list
%  an implementation of 'update_rr' from 'analysis.py' in HeartPy

% - NB: this hasn't been completely implemented - it doesn't contain 'mask' variables. (and needs checking)

rr_source = working_data.RR_list;
b_peaklist = working_data.binary_peaklist;

rel_rrs = (b_peaklist(1:end-1) + b_peaklist(2:end)) == 2;
rr_list = rr_source(rel_rrs);
rr_diff = diff(rr_list);
rr_sqdiff = rr_diff.^2;

working_data.RR_list_cor = rr_list;
working_data.RR_diff = rr_diff;
working_data.RR_sqdiff = rr_sqdiff;

end

function working_data = detect_peaks(hrdata, rol_mean, ma_perc, sample_rate, update_dict, working_data)

rmean = rol_mean;

rol_mean = rmean + ((rmean ./ 100) .* ma_perc);
mn = mean(rmean ./ 100) .* ma_perc;
rol_mean = rmean + mn;

peaksx = find(hrdata > rol_mean);
peaksy = hrdata(peaksx);

peakedges = [1; find(diff(peaksx) > 1); length(peaksx)];
peaklist = [];

for i = 1: length(peakedges)
    try
        y_values = peaksy(peakedges(i):peakedges(i+1));
        [~, max_el] = max(y_values);
        peaklist = [peaklist; peaksx( peakedges(i) + max_el - 1)];
    catch
        % do nothing
    end
end

working_data.peaklist = peaklist;

if update_dict
    
    working_data.ybeat = hrdata(peaklist);
    working_data.rolling_mean = rol_mean;
    working_data = calc_rr(working_data.peaklist, sample_rate, working_data);
    if ~isempty(working_data.RR_list)
        working_data.rrsd = std(working_data.RR_list);
    else
        working_data.rrsd = inf;
    end
    
end

end

function working_data = calc_rr(peaklist, sample_rate, working_data)
% calculate peak-to-peak intervals
%  an implementation of 'calc_rr' from 'analysis.py' in HeartPy

% delete first peak if within first 150ms (signal might start mid-beat after peak)
if length(peaklist) > 0
    if peaklist(1) <= ((sample_rate / 1000.0) * 150)
        peaklist(1) = [];
        working_data.peaklist = peaklist;
        working_data.ybeat(1) = [];
    end
end

rr_list = (diff(peaklist) / sample_rate) * 1000.0;
rr_indices = [peaklist(1:end-1), peaklist(2:end)];
rr_diff = diff(rr_list);
rr_sqdiff = rr_diff.^2;
working_data.RR_list = rr_list;
working_data.RR_indices = rr_indices;
working_data.RR_diff = rr_diff;
working_data.RR_sqdiff = rr_sqdiff;

end

function rol_mean = rolling_mean(hrdata, windowsize, sample_rate)
% from the 'rolling_mean' function within 'datautils.py'


% calculate rolling mean
rol_mean = mean(sliding_window(hrdata, round(windowsize*sample_rate)), 1);


% need to fill 1/2 windowsize gap at the start and end
n_missvals = round((length(hrdata) - length(rol_mean))/2);
missvals_a = ones(1, n_missvals)*rol_mean(1);
missvals_b = ones(1, n_missvals)*rol_mean(end);
rol_mean = [missvals_a, rol_mean, missvals_b];
 
% only to catch length errors that sometimes unexplicably occur.
% Generally not executed, excluded from testing and coverage
%     if len(rol_mean) != len(data): % pragma: no cover
%         lendiff = len(rol_mean) - len(data)
%         if lendiff < 0:
%             rol_mean = np.append(rol_mean, 0)
%         else:
%             rol_mean = rol_mean[:-1]
%         end
%     end

rol_mean = rol_mean(:);

end

function window_mat = sliding_window(hrdata, win_samps)
% an implementation of a Python function

% assemble sliding window matrix
win_starts = win_samps:(length(hrdata)-win_samps+1);
window_mat = nan(win_samps,length(win_starts));
for s = 1 : length(win_starts)
    window_mat(:,s) = hrdata(win_starts(s):win_starts(s)+win_samps-1);
end

end

function hrdata = enhance_peaks(hrdata)
% Peak enhancement
% Implementation of the 'enhance_peaks' function from 'preprocessing.py'

% constants
iterations = 2;

hrdata = scale_data(hrdata);

for i = 1 : iterations
    hrdata = hrdata.^2;
    hrdata = scale_data(hrdata);
end

end

function data = scale_data(data)

% thresholds
lower = 0;
upper = 1024;

% calculate range
rng = max(data) - min(data);

% find minimum
minimum = min(data);

% normalise
data = (upper - lower) .* ((data - minimum) ./ rng) + lower;

end

%%%%%%%%%% end of heartpy_pc

function [peaks, troughs] = bishop_peak_detector(s)

data = s.v;

% adjusted slightly by PC from: https://link.springer.com/chapter/10.1007/978-3-319-65798-1_39

% ----------
% Physiology Feature Extraction Toolkit
% Dr Steven Bishop, 2015-16
% Division of Anaesthesia, University of Cambridge, UK
% Email: sbishop {AT} doctors.org.uk
% ----------
% PEAK_TROUGH_FINDER
%
% Based upon the algorithm by (with updates and optimisations):
% Scholkmann F, Boss J, Wolk M. An Efficient Algorithm for Automatic Peak
% Detection in Noisy Periodic and Quasi-Periodic Signals. Algorithms 2012
% (5), p588-603; doi:10.3390/a5040588
% ----------
% [peaks,troughs,maximagram,minimagram] = PEAK_TROUGH_FINDER(data, {max-interval})
% data: input data as vector
% sampling_frequency (optional): sampling frequency of input
% Returns: vectors [peaks, troughs, maximagram, minimagram] containing
% indices of the peaks and troughs and the maxima/minima scalograms

N = length(data);
if nargin == 2
    L = ceil(vargin{1}/2) - 1;
else
    L = ceil(N/2) - 1;
end
%Detrend the data
meanval = nanmean(data);
data(isnan(data)) = meanval;
data = detrend(data, 'linear');
Mx = false(N, L);   % Changed from 'zeros' by PC
Mn = false(N, L);   % Changed from 'zeros' by PC
%Produce the local maxima scalogram
use_orig_version = 1;
if use_orig_version
    % Original version
    for j=1:L
        k = j;
        for i=k+2:N-k+1
            if data(i-1) > data(i-k-1) && data(i-1) > data(i+k-1)
                Mx(i-1,j) = true;
            end
            if data(i-1) < data(i-k-1) && data(i-1) < data(i+k-1)
                Mn(i-1,j) = true;
            end
        end
    end
else
    % PC's version - to be completed
    for j=1:L
        % Maxima
        
        for i=j+2:N-j+1
            if data(i-1) > data(i-j-1) && data(i-1) > data(i+j-1)
                Mx(i-1,j) = true;
            end
        end
        
        % Minima
        for i=j+2:N-j+1
            if data(i-1) < data(i-j-1) && data(i-1) < data(i+j-1)
                Mn(i-1,j) = true;
            end
        end
    end
    
end
maximagram = Mx;
minimagram = Mn;
%Form Y the column-wise count of where Mx is 0, a scale-dependent distribution of
%local maxima. Find d, the scale with the most maxima (== most number
%of zeros in row). Redimension Mx to contain only the first d scales
Y = sum(Mx==true, 1);
[~, d] = max(Y);
Mx = Mx(:,1:d);
%Form Y the column-wise count of where Mn is 0, a scale-dependent distribution of
%local minima. Find d, the scale with the most minima (== most number
%of zeros in row). Redimension Mn to contain only the first d scales
Y = sum(Mn==true, 1);
[~, d] = max(Y);
Mn = Mn(:,1:d);
%Form Zx and Zn the row-rise counts of Mx and Mn's non-zero elements.
%Any row with a zero count contains entirely zeros, thus indicating
%the presence of a peak or trough
Zx = sum(Mx==false, 2);
Zn = sum(Mn==false, 2);
%Find all the zeros in Zx and Zn. The indices of the zero counts
%correspond to the position of peaks and troughs respectively
peaks = find(~Zx);
troughs = find(~Zn);

end

function peaks = ampd_peak_detector(sig,fs,up)

% based on the description in:
%?F. Scholkmann, J. Boss, and M. Wolf, ?An efficient algorithm for automatic peak detection in noisy periodic and quasi-periodic signals,? Algorithms, vol. 5, no. 4, pp. 588?603, 2012. https://doi.org/10.3390/a5040588

% split into overlapping 6 s windows
win_len = 6; % in secs
overlap = 0.2; % proportion of overlap between consecutive windows
no_samps_in_win = win_len*fs;
if length(sig) <= no_samps_in_win
    win_starts = 1;
    win_ends = length(sig);
else
    win_offset = round(no_samps_in_win*(1-overlap));
    win_starts = 1:win_offset:length(sig)-no_samps_in_win;
    win_ends = win_starts + no_samps_in_win;
    if win_ends(end) < length(sig)
        win_starts(end+1) = length(sig) - no_samps_in_win;
        win_ends(end+1) = length(sig);
    end % this ensures that the windows include the entire signal duration
    
end

% Set up downsampling if the sampling frequency is particularly high
min_fs = 2*up.filtering.beat_detection.elim_high_freqs.Fpass;
if fs > min_fs
    do_ds = 1;
    ds_factor = floor(fs/min_fs);
    ds_fs = fs/floor(fs/min_fs);
else
    do_ds = 0;
end
do_ds = 0;

% detect peaks in each window
peaks = [];
for win_no = 1 : length(win_starts)
    win_sig = sig(win_starts(win_no):win_ends(win_no));
    % downsample signal
    if do_ds
        rel_sig = downsample(win_sig, ds_factor);
    else
        rel_sig = win_sig;
    end
    win_peaks = detect_peaks_using_ampd(rel_sig);
    if do_ds
        win_peaks = win_peaks*ds_factor;
        % find highest point within tolerance either side of detected peaks
        for pk_no = 1 : length(win_peaks)
            curr_peak = win_peaks(pk_no);
            tol_start = curr_peak - ds_factor;
            tol_end = curr_peak + ds_factor;
            [~, temp] = max(win_sig(tol_start:tol_end));
            win_peaks(pk_no) = curr_peak - ds_factor + temp;
            clear temp curr_peak tol_start tol_end
            
        end
    end
    win_peaks = win_peaks + win_starts(win_no) -1;
    peaks = [peaks, win_peaks];  
end
peaks = unique(peaks(:));

% % check
% plot(sig), hold on, plot(peaks, sig(peaks), '*r')
% close all

end

function p = detect_peaks_using_ampd(sig)

do_orig = 0; % whether to do it the way described in the article, or the alternative way

%% Calculate Local Maxima Scalogram (LMS)

% - detrend signal
dt = detrend(sig);

% - initialise matrix
N = length(sig);
L = ceil(N/2)-1;
alpha = 1;
if do_orig
    % - original way
    m = rand(L,N)+alpha;
else
    % - alternative way (only works when using alternative way to find peaks)
    m = ones(L,N)+alpha;
end

% - Detect local maxima
for k = 1:L
        % - original way
        for i = (k+2) : (N-k+1)
            if (sig(i-1) > sig(i-k-1)) && (sig(i-1) > sig(i+k-1))
                m(k,i) = 0;
            end
        end
%         % - rewritten (not used as this takes longer)
%         i = (k+2) : (N-k+1);
%         rel_els = (sig(i-1) > sig(i-k-1)) & (sig(i-1) > sig(i+k-1));
%         m(k,i(rel_els)) = 0;
end

%% Row-wise summation
gamma = sum(m,2);
[~,lambda] = min(gamma);

%% Reshape LMS matrix
% - original interpretation
% k = 1:L;
% mr = m(k<=lambda,:);
% - re-written
mr = m(1:floor(lambda),:);

%% Identify peaks
if do_orig
    % % - original way
    sigma = std(mr);
    p = find(sigma==0);
else
    % - alternative way (knowing that all elements of mr are >= 0)
    p = find(sum(mr)==0);
end

%% (added) Slight correction
% peaks always seem to be one index too high
p = p-1;

% % check
% plot(sig), hold on, plot(p, sig(p), '*r')
% close all

end

function s_filt = elim_vhfs3(s, filt_characteristics)
%% Filter signal to remove VHFs
% Adapted from RRest

s_filt.fs = s.fs;

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Check to see if sampling freq is at least twice the freq of interest
if (filt_characteristics.Fpass/(s.fs/2)) >= 1
    % then the fs is too low to perform this filtering
    s_filt.v = s.v;
    return
end

%% Create filter
% parameters for the low-pass filter to be used
AMfilter = create_lpf(filt_characteristics, s);

%% Re-make filter if it requires too many samples for this signal
% check to see if it requires too many samples
req_sig_length = 3*(length(AMfilter.numerator)-1);
no_attempts = 0;
while no_attempts<4 && length(s.v)<=req_sig_length
    % change Fpass (i.e. the high frequency end of the filter band)
    filt_characteristics.Fpass = filt_characteristics.Fpass+0.75*(filt_characteristics.Fpass-filt_characteristics.Fstop);
    % re-make filter
    AMfilter = create_lpf(filt_characteristics, s);
    % update criterion
    req_sig_length = 3*(length(AMfilter.numerator)-1);
    % increment number of attempts
    no_attempts = no_attempts+1;
end
if length(s.v)<=req_sig_length
    fprintf('\n - Couldn''t perform high frequency filtering')
end

%% Check frequency response
% Gives a -3 dB cutoff at cutoff_freq Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.0512;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(s.fs/2);

%% Remove VHFs
if length(s.v)>req_sig_length
    s_filt.v = filtfilt(AMfilter.numerator, 1, s.v);
else
    s_filt.v = s.v;
end

end

function AMfilter = create_lpf(filt_characteristics, s)

[N,Wn,BETA,TYPE] = kaiserord([filt_characteristics.Fstop filt_characteristics.Fpass]/(s.fs/2), [1 0], [filt_characteristics.Dstop filt_characteristics.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), 'scale');
AMfilter = dfilt.dffir(b);

end

function [peaks,onsets,clipp] = adaptPulseSegment(y,Fs,annot)
%ADAPTPULSESEGMENT perform adaptive pulse segmentation and artifact detection 
%in ppg signals
%   [peaks,onsets,artif] = adaptPulseSegment(y,annot)
%
% Inputs:
%       y      vector, ppg signal [Lx1] or [1xL], in which L = length(signal)
%       Fs      scalar, sampling frequency in Hz
%       annot   vector, timestamps (in samples) with location of the peaks
%
% Outputs:
%       peaks   vector, locations of peaks (in samples)
%       onsets  vector, locations of onsets of the beats (in samples)
%       artif   vector, locations of peaks classified as artefacts
% 
% References:
%       Karlen et al., Adaptive Pulse Segmentation and Artifact Detection in 
%       Photoplethysmography for Mobile Applications, 34th Internat. Conf. IEEE-EMBS 2012
%       
% Written by Marco A. Pimentel

doOptimise = 1;
doPlot = 0;
if nargin < 3
    % no annotations are provided, therefore, no optimisation will take
    % place
    doOptimise = 0; 
end


% The algorithm in the paper is applied to signals sampled at 125 Hz...
% We do not resample our signal
%ys = resample(y,125,Fs);
%Fs = 125;
% if Fs ~= 300
%     ys = resample(y,300,Fs);
%     Fs = 300;
% else
    ys = y;
% end

% The paper is not clear about the selection of the range of search for m
% We define the range of m to be between [0.005 - 0.100] secs (5ms to 100ms)
% We define "m" in terms of samples
opt.bounds = 0.005:0.005:0.100;
opt.m = unique(ceil(opt.bounds*Fs));

opt.perf = zeros(length(opt.m),4); % store results of performance for each m

if doOptimise
    % Perform optimisation
    for i = 1 : length(opt.m)
        % Determine peaks and beat onsets
        [linez,linezSig] = pulseSegment(ys,Fs,opt.m(i));
        % Calculate performance of the peak detection
        opt.perf(i,:) = evalPerf(annot,linez(:,2));
    end
    
else
    % Do not perform optimization; fix m
    opt.m = 10;
    [peaks,artifs,clipp,linez,linezSig] = pulseSegment(ys,Fs,opt.m);
end

if doPlot
    colData = {'g','y','r'};
    figure; 
    h(1) = subplot(211);
    plot(ys); hold on;
    for i = 1 : size(linez,1)
        %if linezSig(i) > -1
        plot(linez(i,:),ys(linez(i,:)),'-x','Color',colData{linezSig(i)+2});
        %end
    end
    
    h(2) = subplot(212);
    plot(ys,'g'); hold on;
    for i = 1 : size(peaks,1)
        plot(peaks(i,:),ys(peaks(i,:)),'-xr');
    end
    if ~isempty(artifs)
    for i = 1 : size(artifs,1)
        plot(artifs(i,:),ys(artifs(i,:)),'--^b');
    end
    end
    if ~isempty(clipp)
    for i = 1 : size(clipp,1)
        plot(clipp(i,:),ys(clipp(i,:)),'-om');
    end
    end
    linkaxes(h,'x');
    
end

% Correct for the downsmapling performed during the peak detection
try
    onsets = peaks(:,1);
catch
    a = 1;
end
peaks  = peaks(:,2);
for i = 1 : size(peaks,1)
    [~,ind]  = min(ys(max([1 onsets(i)-opt.m]):min([length(ys) onsets(i)+opt.m])));
    onsets(i) = max([1 onsets(i)-opt.m]) + ind(1) - 1;
    [~,ind]  = max(ys(max([1 peaks(i)-opt.m]):min([length(ys) peaks(i)+opt.m])));
    peaks(i) = max([1 peaks(i)-opt.m]) + median(ind) - 1;
end

% Correct minimum value of onset of the beat
for i = 2 : length(onsets)
    [~,ind]   = min(ys(peaks(i-1):peaks(i)));
    onsets(i) = peaks(i-1) + ind - 1;
end

end

function [peaks,artifs,clipp,linez,linezSig] = pulseSegment(ys,Fs,m)
% Perform pulse segmentation in ys given m
% Inputs:
%       ys      vector, with ppg signal
%       m       scalar, with the length of each line (read paper for details)
% 
% Outputs:
%       linez      2-column vector, with location of beat onsets and peaks
%       linezSig   vector, with label of the slope of each line segment
%                  1 - positive slope; -1 - negative slope; 0 - constant
% 

% split signal in different segments
nseg = floor(length(ys)/m);    % number of segments
% nseg = round(length(ys)/m);    % number of segments
% intialize loop variables
seg = 1;    % segment counter
z = 1;      % line counter 
segInLine = 1;  % line controler
linez = zeros(nseg,2); linez(1,:) = [1,m];
% slope of segment/line
a = zeros(nseg,1); a(1) = slope(ys,linez(1,:));
% "classify" line segment according to the slope
linezSig = zeros(nseg,1); linezSig(1) = sign(a(1));
% Start loop over segments
z = z + 1; seg = seg + 1;
for i = 2 : nseg    % loop over segments
    linez(z,:) = [(seg-1)*m+1 seg*m];
    try
        a(z) = slope(ys,linez(z,:));
    catch
        a = 1;
    end
    linezSig(z) = sign(a(z));
    if sign(a(z)) == sign(a(z-1))
        linez(z-1,:) = [(seg-1-segInLine)*m+1 seg*m];
        seg = seg + 1;
        segInLine = segInLine + 1;
    else
        z = z + 1;
        seg = seg + 1;
        segInLine = 1;
    end
end

% remove extra spaces created in output variables
linezSig(sum(linez,2)==0,:) = [];
linez(sum(linez,2)==0,:) = [];

% Apply adaptive threshold algorithm
% For this algorithm to work, we need to first find a valide line segment 
% in order to intialize the thresholds! In order to this, we define a flag
% to control the intialization in the main loop
FOUND_L1 = 0;

% The algorithm includes the definition of 4 adaptation parameters
% We define the following adaptation parameters
% a =     | a_fast_low    a_fast_high |
%         | a_slow_low    a_slow_high |
% 
a = [0.5 1.6; ...
     0.6 2.0];
 
% Define fixed thresholds described in the paper
ThT  = 0.03 * Fs;    % Duration of the line
ThIB = 0.24 * Fs;    % Interbeat invertal (240 ms) 

% Define parameters used in the main loop
alpha = zeros(size(linez,1),1);
for i = 1 : size(linez,1)
    alpha(i) = slope(ys,linez(i,:));   % slopes of line segments
end
theta = diff(ys(linez),[],2);
durat = diff(linez,[],2);       % duration of line segments (in samples)

% remove lines that do not have the necessary duration
linez(durat<ThT,:) = [];
theta(durat<ThT,:) = [];
alpha(durat<ThT,:) = [];
horiz = horizontalLine(ys,linez,Fs);

FLAG = 0;
artifs = []; clipp = [];
% Select window for detect firs peaks!
wind = theta(theta>0);
try 
    wind = wind(1:10);
catch
    wind = wind;
end
ThAlow  = prctile(wind,95)*0.6;
ThAhigh = prctile(wind,95)*1.8;
peaks = [];
for z = 1 : size(linez,1)-1   % loop over line segments
    if FOUND_L1
        if alpha(z) > 0 && ... % slope must be positive
                alpha(z-1) ~= 0 && ...  % peaks before or after clipping are artefactual
                alpha(z+1) ~= 0
            if theta(z) >= ThAlow && theta(z) <= ThAhigh && ...
                    linez(z,2) >= peaks(end,2) + ThIB
                ThAlow  = (ThAlow + theta(z)*a(2,1))/2;
                ThAhigh = theta(z) * a(2,2);
                FLAG = 0;
                currTheta = [currTheta; theta(z)];
                peaks = [peaks; linez(z,:)];
            else
                if FLAG > 0
                    ThAlow  = (ThAlow + min(currTheta(max([1 end-4]):end))*a(1,1))/2;
                    ThAhigh = max(currTheta(max([1 end-4]):end)) * a(1,2);
                    %ThAlow  = (ThAlow + theta(z)*a(1,1))/2;
                    %ThAhigh = theta(z) * a(1,2);
                end
                FLAG = FLAG + 1;
                artifs = [artifs; linez(z,:)];
            end
        elseif theta(z) > 0 && ... 
                ((theta(z-1) ~= 0 || horiz(z-1) ~= 0) && ...
                (theta(z+1) ~= 0 || horiz(z+1) ~= 0))
            if theta(z) >= ThAlow && theta(z) <= ThAhigh && ...
                    linez(z,2) >= peaks(end,2) + ThIB
                ThAlow  = (ThAlow + theta(z)*a(2,1))/2;
                ThAhigh = theta(z) * a(2,2);
                FLAG = 0;
                currTheta = [currTheta; theta(z)];
                peaks = [peaks; linez(z,:)];
            else
                if FLAG > 0
                    %ThAlow  = (ThAlow + currTheta*a(1,1))/2;
                    %ThAhigh = currTheta * a(1,2);
                    ThAlow  = (ThAlow + min(currTheta(max([1 end-4]):end))*a(1,1))/2;
                    ThAhigh = max(currTheta(max([1 end-4]):end)) * a(1,2);
                    %ThAlow  = (ThAlow + theta(z)*a(1,1))/2;
                    %ThAhigh = theta(z) * a(1,2);
                end
                FLAG = FLAG + 1;
                artifs = [artifs; linez(z,:)];
            end
        elseif theta(z) == 0 && horiz(z) == 0
            artifs  = [artifs; linez(z,:)];
            clipp   = [clipp; linez(z,:)];
        end 
    else
        if alpha(z) > 0 && durat(z) >= ThT && ...
                theta(z) >= ThAlow && theta(z) <= ThAhigh 
            FOUND_L1 = 1;
            ThAlow  = theta(z)*0.5;
            ThAhigh = theta(z)*2.0;
            peaks = linez(z,:);    % loaction of onsets and peaks
            currTheta = theta(z);
        end
    end
end

end

function out = horizontalLine(ys,linez,Fs)
% Get horizontal lines from signal given linez
out = zeros(size(linez,1),1);
for i = 1 : size(linez,1)
    out(i) = median(abs(diff(ys(linez(i,1):linez(i,2)))));
    % check duration of the peaks
    if out(i) == 0 && diff(linez(i,:)) <= 0.200*Fs
        out(i) = 0.1;
    end
end

end

function out = slope(ys,interv)
start = interv(1); stop = interv(2);
out = sum(diff(ys([start:stop])))/(stop-start);
%out = median(gradient(ys(start:stop)));
end

function quality = assess_signal_quality(s, pulses, up)
% ASSESS_SIGNAL_QUALITY  Assesses the signal quality of each beat of the
% pulsatile signal.
% Inputs:
%       s           -  pulsatile signal, a structure containing s.v (a
%                       vector of values), and s.fs (sampling frequency in Hz).
%       pulse_inds  -  indices of the pulse peaks
%
% Outputs:
%       quality     -  the signal quality of each beat (1 indicates high
%                       quality, 0 low quality).
%
% Adapted from RRest
%
% Reference: This function uses an adaptation of the signal quality index
% for the photoplethysmogram described in:
%     Orphanidou, C. et al., 2015. Signal-quality indices for the electrocardiogram and photoplethysmogram: derivation and applications to wireless monitoring. IEEE Journal of Biomedical and Health Informatics, 19(3), pp.832–8. Available at: http://www.ncbi.nlm.nih.gov/pubmed/25069129.

% Skip if not required
if ~up.options.do_quality
    quality = true(length(pulses.peaks),1);    
    return
end

%% Skip if there is only a single pulse wave (or none detected)
if length(pulses.peaks) == 1
    quality(size(pulses.peaks)) = true;
    return
elseif isempty(pulses.peaks)
    quality = [];
    return
end

%% Setup
pulse_inds = pulses.onsets;
thresh_orig = 0.86; %threshold = 0.66 for ECG, 0.86 for PPG;    % ecg cross-correlation threshold value (for sqi)
thresh = 0.95; %1-(0.5*(1-thresh_orig));
s.t = [0:length(s.v)-1]/s.fs;
pulse_ms_inds = nan(length(pulse_inds)-1,1);
for pulse_no = 1 : length(pulse_inds)-1
    [~, temp] = max(diff(s.v(pulse_inds(pulse_no):pulse_inds(pulse_no+1))));
    pulse_ms_inds(pulse_no) = pulse_inds(pulse_no) -1 + temp; clear temp    
end
clear pulse_inds_temp pulse_no

%% Segment into windows of 10s duration
win_durn = 10;   % in secs
win_overlap = 2; % in secs
win.deb = s.t(1):(win_durn-win_overlap):s.t(end);
win.fin = win.deb + win_durn;
if win.fin(end)>s.t(end)
    win.fin(end) = s.t(end);
    win.deb(end) = win.fin(end)-win_durn;
end

high_quality_pulse_inds = [];
for win_no = 1 : length(win.deb)
    
    % identify data for this window
    
    rel_els = s.t >= win.deb(win_no) & s.t <= win.fin(win_no);
    first_el = find(rel_els,1);
    curr_sig.v = s.v(rel_els); 
    curr_sig.t = s.t(rel_els); clear rel_els
    curr_sig.t = curr_sig.t - curr_sig.t(1);
    curr_sig.pulse_ind_inds = find(s.t(pulse_ms_inds) >= win.deb(win_no) & s.t(pulse_ms_inds) <= win.fin(win_no));
    curr_pulse_inds = pulse_ms_inds(curr_sig.pulse_ind_inds) - first_el + 1;
    
    % find beat-to-beat interval
    
    beat_to_beat_interval = median(diff(curr_sig.t(curr_pulse_inds)));
    beat_to_beat_samples = round(beat_to_beat_interval*s.fs); clear beat_to_beat_interval
    
    % find a template beat
    ts = [];
    rel_els = curr_pulse_inds>beat_to_beat_samples/2 & ...
        curr_pulse_inds+floor(beat_to_beat_samples/2)<length(curr_sig.v);
    rel_pulse_inds = curr_pulse_inds(rel_els);
    curr_sig.used_pulse_ind_inds = curr_sig.pulse_ind_inds(rel_els);
    % find beat morphologies
    for rel_pulse_no = 1 : length(rel_pulse_inds)
        this_pulse = rel_pulse_inds(rel_pulse_no);
        t = curr_sig.v(this_pulse-floor(beat_to_beat_samples/2):this_pulse+floor(beat_to_beat_samples/2));
        tt = t/norm(t); tt = tt(:)';
        ts = [ts; tt]; clear tt t
    end
    clear k l j
    
    % find all templates in current window
    avtempl = mean(ts,1);
    
    % now calculate correlation for every beat in this window
    r2 = nan(size(ts,1),1);
    for k = 1:size(ts,1)
        r2(k) = corr2(avtempl,ts(k,:));
    end
    clear k
    
    high_quality_beats = r2 > thresh;
    
    % % calculate template using only high-quality beats
    do_repeat = 1;
    if do_repeat
        r2_orig = r2;
        avtempl_orig = avtempl;
        avtempl = mean(ts(r2>thresh,:),1);  
        r2 = nan(size(ts,1),1);
        for k = 1:size(ts,1)
            r2(k) = corr2(avtempl,ts(k,:));
        end
        clear k
        
        high_quality_beats = r2 > thresh;
        
    end
    
    high_quality_pulse_inds = [high_quality_pulse_inds, curr_sig.used_pulse_ind_inds(high_quality_beats)];
    
end

high_quality_pulse_inds = unique(high_quality_pulse_inds);
quality = false(length(pulse_inds),1);
for ind = 1 : length(quality)
    if intersect(high_quality_pulse_inds, ind)
        quality(ind) = true;
    end
end
end

function [sigs, pulses] = filter_signal(S, pulses, up)

%% Make output variable
sigs.fs = S.fs;
sigs.orig = S.v;

if up.options.do_filter
    
    %% Filter to remove high frequencies
    filt_characteristics = up.filtering.fiducial_points.elim_high_freqs;
    s_filt = elim_vhfs3(S, filt_characteristics);
    
    %% Filter to remove low frequencies
    if strcmp(up.analysis.no_of_pulses, 'multiple')
        filt_characteristics = up.filtering.fiducial_points.elim_low_freqs;
        s_filt = elim_vlfs(s_filt, filt_characteristics, up);
    else
        % remove baseline wander
        baseline = linspace(0, s_filt.v(end)-s_filt.v(1), length(s_filt.v)); baseline = baseline(:);
        s_filt.v = s_filt.v-baseline;
    end
    
    %% Output
    sigs.filt = s_filt.v;
    sigs.curr = sigs.filt;
    
    %% Adjust pulse onsets to be at minima on filtered signal
    sig.v = sigs.curr;
    sig.fs = sigs.fs;
    pulses = adjust_pulse_onsets(pulses, sig, up);
    
else
    % If not filtering
    sigs.curr = sigs.orig;
    pulses.curr_onsets = pulses.onsets;
end

end

function pulses = adjust_pulse_onsets(pulses, sig, up)

tol = round(up.paramSet.onset_tol*sig.fs);
for pw_no = 1 : length(pulses.onsets)
    start_el = pulses.onsets(pw_no)-tol;
    if start_el<1, start_el = 1; end
    end_el = pulses.onsets(pw_no)+tol;
    if end_el>length(sig.v), end_el = length(sig.v); end
    [~,min_el] = min(sig.v(start_el:end_el));
    pulses.onsets(pw_no) = start_el-1+min_el;
    clear min_el start_el end_el
end

end

function sigs = calib_signal(sigs, pulses, up)

% This function calibrates the signal to ensure that it occupies the
% required blood pressure range.

% Calibrate pulse wave
if up.options.calibrate_pw
    
    % setup
    sigs.calib = nan(size(sigs.filt));
    
    for pw_no = 1 : length(pulses.onsets)-1
        
        % identify current pulse wave, and signal portion to be calibrated
        temp.v_pw = sigs.curr(pulses.onsets(pw_no):pulses.onsets(pw_no+1));
        if pw_no == 1
            start_el = 1;
            end_el = pulses.onsets(pw_no+1);
        elseif pw_no == length(pulses.onsets)-1
            start_el = pulses.onsets(pw_no);
            end_el = length(sigs.curr);
        else
            start_el = pulses.onsets(pw_no);
            end_el = pulses.onsets(pw_no+1);
        end
        temp.v_to_calib = sigs.curr(start_el:end_el);
        
        % Ensure that pulse onset and end are at the same amplitude
        temp2 = interp1([1,length(temp.v_pw)], temp.v_pw([1,end]), 1:length(temp.v_pw));
        temp.v_pw = temp.v_pw - temp2(:);
        no_before = pulses.onsets(pw_no)-start_el;
        temp2a = interp1(temp2, 1-no_before:0, 'linear', 'extrap');
        no_after = end_el - pulses.onsets(pw_no+1);
        temp2b = interp1(temp2, length(temp2)+1:length(temp2)+no_after, 'linear', 'extrap');
        temp2 = [temp2a, temp2, temp2b];
        temp.v_to_calib = temp.v_to_calib-temp2(:);
        clear temp2 temp2a temp2b
        
        % Calibrate this pulse wave
        if sum(strcmp(fieldnames(up.options), 'calib_dbp')) && sum(strcmp(fieldnames(up.options), 'calib_sbp'))
        
            % Calibrate using systolic and diastolic pressures
            offset = -1*temp.v_pw(end);
            scale = max(temp.v_pw)-offset;
            temp.v_to_calib = (temp.v_to_calib+offset)/scale;
            temp.v_to_calib = temp.v_to_calib.*(up.options.calib_sbp-up.options.calib_dbp); % scale
            temp.v_to_calib = temp.v_to_calib + up.options.calib_dbp; % offset
            
        elseif sum(strcmp(fieldnames(up.options), 'calib_dbp')) && sum(strcmp(fieldnames(up.options), 'calib_map'))
            
            % Calibrate using diastolic and mean pressures
            offset = -1*temp.v_pw(end);
            scale = max(temp.v_pw)-offset;
                        
            mean_factor = mean( (temp.v_pw+offset)/max(temp.v_pw+offset) );
            temp.v_to_calib = ((temp.v_to_calib/scale)+offset)*((up.options.calib_map-up.options.calib_dbp)/mean_factor);
            temp.v_to_calib = temp.v_to_calib + up.options.calib_dbp;
            
        else
            error('The necessary calibration coefficients have not been specified.')
        end
        
        % Store the calibrated pulse wave
        sigs.calib(start_el:end_el) = temp.v_to_calib(1:end);
    end
    sigs.curr = sigs.calib; 
end

end

function [sigs, pulses] = generate_central_pressure(sigs, pulses, up)

if up.options.tran_func
    
    % setup variables
    old = sigs.curr;
    sig_type = up.options.sig_type;
    s.fs = sigs.fs;
    s.v = sigs.curr;
    
    % add on buffer at either end
    tol = floor(s.fs);
    s.v = [repmat(s.v(1), [tol,1]); s.v; repmat(s.v(end), [tol,1])];
    
    % transform signal
    temp = per_cen_TF(s, sig_type);
    
    % time align with original signal
    [r,lags] = xcorr(temp.centrABP.v, s.v);
    [~, lag_el] = max(r);
    lag = lags(lag_el); clear lags lag_el
    aligned_sig = temp.centrABP.v(tol+lag:tol+lag+length(old)-1);
    
%     % check
%     plot(old), hold on, plot(aligned_sig)
%     close all
    
    % store
    sigs.curr = aligned_sig;
    sigs.tf = sigs.curr;
    
    % Adjust pulse onsets to be at minima on filtered signal
    sig.v = sigs.curr;
    sig.fs = sigs.fs;
    pulses = adjust_pulse_onsets(pulses, sig, up);
    
end

end

function [sigs, pulses, up] = calculate_average_pw(sigs, pulses, up)

if up.options.calc_average_pw
    
    if ~strcmp(up.analysis.no_of_pulses, 'multiple')
        error('Couldn''t calculate average pulse wave as only a single pulse wave was detected')
    else
        
        % Extract relevant signal
        rel_sig.fs = sigs.fs;
        rel_sig.v = sigs.curr;
        if ~up.options.do_filter
            s_g_filter_len = up.filtering.derivatives.s_g_filt_len_no_filtering;
        else
            s_g_filter_len = up.filtering.derivatives.s_g_filt_len_no_filtered;
        end
        
        % calculate first derivative
        dt = 1/rel_sig.fs;
        first_d = savitzky_golay(rel_sig.v, 1, 9)./dt;
        
        % Identify ms points
        ms_inds = nan(length(pulses.quality)-1,1);
        for pw_no = 1 : length(pulses.quality)-1
            [~, temp_ms] = max(first_d(pulses.onsets(pw_no):pulses.onsets(pw_no+1)));
            ms_inds(pw_no) = temp_ms+pulses.onsets(pw_no)-1;
        end
        
        % Calculate median PW duration
        diffs = diff(pulses.onsets);
        diffs = diffs(pulses.quality(1:end-1));
        durn = round(mean(diffs)+1);
        if rem(durn,2)==0, durn=durn+1; end
                
        % extract data for each pulse wave, centred on the ms point
        pws = nan(length(pulses.quality)-1, durn);
        for pw_no = 1 : length(pulses.quality)-1
            start_ind = ms_inds(pw_no) - ((durn-1)/2);
            end_ind = ms_inds(pw_no) + ((durn-1)/2);
            if start_ind <1 || end_ind > length(rel_sig.v)
                continue
            end
            curr_pw = rel_sig.v(start_ind:end_ind);
            temp = interp1([1,length(curr_pw)], [curr_pw(1), curr_pw(end)], 1: length(curr_pw));
            curr_pw = curr_pw-temp(:)+curr_pw(1);
            pws(pw_no,:) = curr_pw;
            clear start_ind end_ind
        end
        
        % eliminate low quality pulse waves
        pws = pws(pulses.quality(1:end-1) & ~sum(isnan(pws'))',:);
        
        % align at ms point (this doesn't affect the average PW)
        pws = pws - pws(:,(durn-1)/2) + mean(pws(:,(durn-1)/2));
        
        % find average pw shape
        ave_pw.v = mean(pws);
        ave_pw.v = ave_pw.v(1:end-1);
        ave_pw.fs = sigs.fs;
        
        % start at pulse onset
        [ave_pw, align_el] = align_pulse(ave_pw, up);

        sigs.ave = ave_pw.v(:);
        sigs.curr = sigs.ave;
        
        up.analysis.no_of_pulses = 'single';
        pulses.ave.onsets = [1, length(sigs.ave)];
        if sum(isnan(sigs.ave)) == 0, pulses.ave.quality =1; else pulses.ave.quality = 0; end
        [~, pulses.ave.peaks] = max(sigs.ave);
        
        % store PWs
        sigs.pws = pws'; % - pws(:,(durn-1)/2)'; % + sigs.ave((durn-1)/2);
        sigs.pws = sigs.pws([align_el:end, 1:(align_el-1)],:);
        
    end    
end

end

function sigs = calculate_derivs(sigs, up, S)

% Skip if not calculating pulse wave indices
if up.options.calc_pw_inds == 0
    pts = [];
    return
end

%% Extract relevant signal (either filtered or not)
rel_sig.fs = sigs.fs;
rel_sig.v = sigs.curr;
% if ~up.options.do_filter
%     s_g_filter_len = up.filtering.derivatives.s_g_filt_len_no_filtering;
% else
%     s_g_filter_len = up.filtering.derivatives.s_g_filt_len_no_filtered;
% end

%% Repeat single pulse wave to avoid edge effects
if strcmp(up.analysis.no_of_pulses, 'single')
    rel_sig.v = repmat(rel_sig.v, [up.paramSet.no_pulse_repeats,1]);
end

%% Calculate derivatives
dt = 1/rel_sig.fs;
sigs.first_d = S.d1;%savitzky_golay(rel_sig.v, 1, 9)./dt;
sigs.second_d = S.d1;%savitzky_golay(sigs.first_d, 1, 9)./dt;
sigs.third_d = savitzky_golay(rel_sig.v, 3, 9)./dt;
% sigs.third_d = savitzky_golay(sigs.second_d, 1, 9)./dt;
sigs.fourth_d = savitzky_golay(sigs.third_d, 1, 9)./dt;
sigs.fifth_d = savitzky_golay(sigs.fourth_d, 1, 9)./dt;
sigs.sixth_d = savitzky_golay(sigs.fifth_d, 1, 9)./dt;
sigs.seventh_d = savitzky_golay(sigs.sixth_d, 1, 9)./dt;

%% Retain original pulse wave
if strcmp(up.analysis.no_of_pulses, 'single')
    % select derivative values corresponding to the original pulse
    orig_len = length(rel_sig.v)/up.paramSet.no_pulse_repeats;
    orig_els = (floor(up.paramSet.no_pulse_repeats/2)*orig_len)+1 : ((floor(up.paramSet.no_pulse_repeats/2)+1)*orig_len);
    sigs.first_d = S.d1;%sigs.first_d(orig_els);
    sigs.second_d = S.d2;%sigs.second_d(orig_els);
    sigs.third_d = sigs.third_d(orig_els);
end

end

function pts = identify_fiducial_point_indices(sigs, pulses, up)

% Skip if not calculating pulse wave indices
if up.options.calc_pw_inds == 0
    pts = [];
    return
end

%% Identify fiducial point indices

% use relevant pulses data
if up.options.calc_average_pw
    rel_pulses = pulses.ave;
else
    rel_pulses = pulses;
end

% setup variables to store fiducial point indices
fid_pt_names = {'a', 'b', 'c', 'd', 'e', 'f', 's', 'dia', 'dic', 'p1pk', 'p2pk', 'p1in', 'p2in', 'ms', 'f1', 'f2', 'ms2'};
for fid_pt_no = 1 : length(fid_pt_names)
    eval(['pts.ind.' fid_pt_names{fid_pt_no} ' = nan(length(rel_pulses.onsets)-1,1);'])
    eval(['pts.amp_norm.' fid_pt_names{fid_pt_no} ' = nan(length(rel_pulses.onsets)-1,1);'])
    eval(['pts.amp.' fid_pt_names{fid_pt_no} ' = nan(length(rel_pulses.onsets)-1,1);'])
    eval(['pts.t.' fid_pt_names{fid_pt_no} ' = nan(length(rel_pulses.onsets)-1,1);'])
end

% setup buffer zones
buffer_tol.deb = 0.005; % in secs
buffer_tol.fin = 0.2; % in proportion of beat
buffer_p1 = [0.1, 0.18]; % in secs

% cycle through each pulse wave
no_pulses = length(rel_pulses.onsets);

for pulse_no = 1 : no_pulses-1
    
    % extract data for this pulse wave
    curr_els = rel_pulses.onsets(pulse_no):rel_pulses.onsets(pulse_no+1);
    
    % (either filtered or not)
    rel_sig.fs = sigs.fs;
    curr.sig = sigs.curr(curr_els);
    curr.derivs.first = sigs.first_d(curr_els);
    curr.derivs.second = sigs.second_d(curr_els);
    curr.derivs.third = sigs.third_d(curr_els);
    
    %% Pre-process
    
    % Make sure the pulse wave ends at the same amplitude as its onset.
    old.v = curr.sig; old.fs = sigs.fs;
    temp = eliminate_low_freq_from_single_beat(old, up);
    
    % Normalise pulse wave
    if up.options.normalise_pw
        
        % Normalise to occupy a range of 0 to 1
        temp.v = (temp.v-min(temp.v))/range(temp.v);
        
        % Ensure that signal commences at start of systole
        %     temp = align_pulse(temp, up);
        
    end
    
    % Calibrate pulse wave
    if up.options.calibrate_pw
        if sum(strcmp(fieldnames(up.options), 'calib_dbp')) && sum(strcmp(fieldnames(up.options), 'calib_sbp'))
            temp.v = temp.v-min(temp.v);
            temp.v = (temp.v/max(temp.v));
            temp.v = temp.v.*(up.options.calib_sbp-up.options.calib_dbp); % scale
            temp.v = temp.v + up.options.calib_dbp; % offset
        elseif sum(strcmp(fieldnames(up.options), 'calib_dbp')) && sum(strcmp(fieldnames(up.options), 'calib_map'))
            temp.v = temp.v-min(temp.v);
            temp.v = (temp.v/max(temp.v));
            mean_factor = mean(temp.v);
            temp.v = temp.v*((up.options.calib_map-up.options.calib_dbp)/mean_factor);
            temp.v = temp.v + up.options.calib_dbp;
        else
            error('The necessary calibration coefficients have not been specified.')            
        end 
    end
    
    curr.sig = temp.v;
    
    %% Identify fiducial points
    
    % find buffers
    initial_buffer = floor(buffer_tol.deb * sigs.fs);  % in samples
    end_buffer = length(curr.sig) - ceil(buffer_tol.fin * length(curr.sig));  % in samples
    
    % find f1 and f2
    temp_f1 = 1;
    temp_f2 = length(curr.sig);
    
    % find s
    temp_s = identify_s(curr);
    
    % find ms
    temp_ms = identify_ms(curr);
    
    % find a
    temp_a = identify_a(curr, initial_buffer);
    
    if isempty(temp_a)
        continue
    end
    
    % find b
    temp_b = identify_b(curr, temp_a, up);
    
    if isempty(temp_b)
        continue
    end
    
    % find p1
    p1_buffer = floor(buffer_p1 .* sigs.fs);  % in samples
    temp_p1 = identify_p1(curr, temp_b, temp_ms, p1_buffer, up.analysis.no_of_pulses, sigs.fs, up);
    
    % find e
    redo_log = 0;
    temp_e = identify_e(curr, temp_s, temp_ms, temp_b, end_buffer, redo_log);
    
    if isempty(temp_e)
        [temp_f, temp_c, temp_d, temp_dic, temp_dia, temp_p2] = deal(nan);
    else
        
        % find c
        temp_c = identify_c(curr, temp_b, temp_e);
        
        if isempty(temp_c)
            redo_log = 1;
            
            % re-find e
            old_temp_e = temp_e;
            temp_e = identify_e(curr, temp_s, temp_ms, end_buffer, redo_log);
            
            % re-find c
            temp_c = identify_c(curr, temp_b, temp_e);
        end
        
        % find f
        temp_f = identify_f(curr, temp_e, end_buffer);
        
        % find dic
        temp_dic = identify_dic(curr,temp_e);
        
        % find dia
        temp_dia = identify_dia(curr, temp_dic, temp_e, end_buffer, up);
        
        if isempty(temp_c)
            [temp_d, temp_p2] = deal(nan);
        else
            % find d
            temp_d = identify_d(curr, temp_c, temp_e);
            
            % find p2
            temp_p2 = identify_p2(curr, temp_d, temp_p1, temp_dic);
            
        end
        
    end
    
    % retain timings of original p1 and p2 estimates
    temp_p1in = temp_p1; temp_p2in = temp_p2;
    
    % make p1 or p2 coincident with the systolic peak
    [~, rel_el] = min(abs(temp_s-[temp_p1,temp_p2]));
    if rel_el == 1
        temp_p1 = temp_s;
    else
        temp_p2 = temp_s;
    end
    
    if ~isnan(temp_p2) && ~isnan(temp_p1)
        % make sure p2 is at a peak if necessary
        pks = find_pks_trs(curr.sig, 'pk');
        cutoff = mean([temp_p1, temp_p2]);
        possible_pks = find(pks > cutoff & pks < temp_e & curr.sig(pks) > curr.sig(temp_p2));
        if ~isempty(possible_pks)
            [~, temp_el] = max(curr.sig(pks(possible_pks)));
            temp_p2 = pks(possible_pks(temp_el));
        end
        % make sure p1 is at a peak if necessary
        pks = find_pks_trs(curr.sig, 'pk');
        cutoff = mean([temp_p1, temp_p2]);
        possible_pks = find(pks < cutoff & pks > temp_ms & curr.sig(pks) > curr.sig(temp_p1));
        if ~isempty(possible_pks)
            [~, temp_el] = max(curr.sig(pks(possible_pks)));
            temp_p1 = pks(possible_pks(temp_el));
        end
    end
    
    % store p1pk and p2pk
    temp_p1pk = temp_p1;
    temp_p2pk = temp_p2;
    
    % identify ms2
    temp_ms2 = identify_ms2(curr, temp_p1in, temp_p2pk, up);
    
    % calculate scale factor
    curr_range = range(curr.sig);
    scale_factor = 1/curr_range;
    
    % store points
    for fid_pt_no = 1 : length(fid_pt_names)
        eval(['curr_temp_el = temp_' fid_pt_names{fid_pt_no} ';']);
        if ~isnan(curr_temp_el)
            
            % - index of fiducial point
            sig_ind = curr_temp_el + curr_els(1)-1;
            eval(['pts.ind.' fid_pt_names{fid_pt_no} '(pulse_no) = sig_ind;'])
            
            % - amplitude of fiducial point
            if sum(strcmp(fid_pt_names{fid_pt_no}, {'a','b','c','d','e','f'}))
                amp = curr.derivs.second(curr_temp_el);
            elseif sum(strcmp(fid_pt_names{fid_pt_no}, {'s','dia','dic','p1pk','p1in','p2pk','p2in','f1','f2'}))
                amp = curr.sig(curr_temp_el);
            elseif sum(strcmp(fid_pt_names{fid_pt_no}, {'ms','ms2'}))
                amp = curr.derivs.first(curr_temp_el);
            end
            eval(['pts.amp.' fid_pt_names{fid_pt_no} '(pulse_no) = amp;'])
            amp_norm = amp*scale_factor;
            eval(['pts.amp_norm.' fid_pt_names{fid_pt_no} '(pulse_no) = amp_norm;'])
            
            % - timing of fiducial point
            t = (curr_temp_el-1)/sigs.fs;
            eval(['pts.t.' fid_pt_names{fid_pt_no} '(pulse_no) = t;'])
            clear amp_norm sig_ind amp t
        end
    end
    
    
    clear curr temp* curr_els empty_log
    
end
clear pulse_no


%% Ensure only pts are provided for those pulses with all pts available
pt_names = fieldnames(pts.ind);
include_pulse = true(size(pts.ind.dia));
for pt_name_no = 1 : length(pt_names)
    eval(['curr_pt_measures = pts.ind.' pt_names{pt_name_no} ';']);
    include_pulse(isnan(curr_pt_measures)) = false;    
end
for pt_name_no = 1 : length(pt_names)
    if strcmp(pt_names{pt_name_no}, 'f1') || strcmp(pt_names{pt_name_no}, 'f2')
        continue
    end
    eval(['pts.ind.' pt_names{pt_name_no} '(~include_pulse) = nan;']);
    eval(['pts.amp_norm.' pt_names{pt_name_no} '(~include_pulse) = nan;']);
    eval(['pts.amp.' pt_names{pt_name_no} '(~include_pulse) = nan;']);
    eval(['pts.t.' pt_names{pt_name_no} '(~include_pulse) = nan;']);
end

end

function pw_inds = calculate_pw_inds(fid_pts, sigs, pulses, ht, up)

% Skip if not calculating pulse wave indices
if up.options.calc_pw_inds == 0
    pw_inds = [];
    return
end

%% Timings

% - Duration of pulse
pw_inds.T = fid_pts.t.f2-fid_pts.t.f1;

% - Instantaneous heart rate (bpm)
pw_inds.IHR = 60./pw_inds.T;

% - Time between systolic and diastolic peaks (secs)
pw_inds.delta_t = fid_pts.t.dia-fid_pts.t.s;

% - Crest time (secs)
pw_inds.CT = fid_pts.t.s-fid_pts.t.f1;

% - Crest time divided by height (s/m)
pw_inds.CT_div_ht = pw_inds.CT./ht;

% - Stiffness index (m/s)
pw_inds.SI = ht./pw_inds.delta_t;

% - Proportion of pulse wave duration which is spent on systolic upslope
pw_inds.prop_s = pw_inds.CT./pw_inds.T;

% - Duration of systole
pw_inds.t_sys = fid_pts.t.dic-fid_pts.t.f1;

% - Duration of diastole
pw_inds.t_dia = fid_pts.t.f2-fid_pts.t.dic;

% - Ratio of systolic to diastolic durations
pw_inds.t_ratio = pw_inds.t_sys./pw_inds.t_dia;

% - Proportion of pulse wave duration taken up by time between systolic and diastolic peaks
pw_inds.prop_delta_t = pw_inds.delta_t./pw_inds.T;

% - Time from p1 to diastolic peak
pw_inds.t_p1in_dia = fid_pts.t.dia-fid_pts.t.p1in;
pw_inds.t_p1pk_dia = fid_pts.t.dia-fid_pts.t.p1pk;

% - Time from p2 to diastolic peak
pw_inds.t_p2in_dia = fid_pts.t.dia-fid_pts.t.p2in;
pw_inds.t_p2pk_dia = fid_pts.t.dia-fid_pts.t.p2pk;

% - Time from 'b' to 'c'
pw_inds.t_b_c = fid_pts.t.c-fid_pts.t.b;

% - Time from 'b' to 'd'
pw_inds.t_b_d = fid_pts.t.d-fid_pts.t.b;

%% Amplitudes

% - Pulse amplitude
pw_inds.pulse_amp = fid_pts.amp.s-fid_pts.amp.f1;
pw_inds.pulse_amp_p1 = fid_pts.amp.p1in-fid_pts.amp.f1;
pw_inds.pulse_amp_p2 = fid_pts.amp.p2in-fid_pts.amp.f1;

% - Augmentation Pressure
pw_inds.AP = fid_pts.amp.p2pk-fid_pts.amp.p1in;

% - Agumentation Index
pw_inds.AI = 100*pw_inds.AP./pw_inds.pulse_amp;

% - Diastolic peak amplitude
pw_inds.dia_amp = fid_pts.amp.dia-fid_pts.amp.f1;

% - Reflection Index (calculated using systolic peak)
pw_inds.RI = pw_inds.dia_amp./pw_inds.pulse_amp;

% - Reflection Index (calculated using p1)
pw_inds.RI_p1 = pw_inds.dia_amp./pw_inds.pulse_amp_p1;

% - Reflection Index (calculated using p2)
pw_inds.RI_p2 = pw_inds.dia_amp./pw_inds.pulse_amp_p2;

% - Ratio of amplitudes of p2 and p1
pw_inds.ratio_p2_p1 = pw_inds.pulse_amp_p2./pw_inds.pulse_amp_p1;

%% Areas

% calculate area indices which require additional pulse wave analysis
for beat_no = 1:length(fid_pts.ind.f1)
    
    % skip if there isn't sufficient information for this pulse wave
    if isnan(fid_pts.ind.s(beat_no)), pw_inds.A1(beat_no) = nan; pw_inds.A2(beat_no) = nan; continue, end
    
    % find baseline of pulse wave (joining initial pulse onset to final pulse onset)
    baseline = linspace(sigs.curr(fid_pts.ind.f1(beat_no)), sigs.curr(fid_pts.ind.f2(beat_no)), fid_pts.ind.f2(beat_no) - fid_pts.ind.f1(beat_no)+1); baseline = baseline(:);
    
    % - systolic area
    rel_pts = fid_pts.ind.f1(beat_no) : fid_pts.ind.dic(beat_no);
    baseline_pts = rel_pts - fid_pts.ind.f1(beat_no) + 1;
    pw_inds.A1(beat_no) = sum(sigs.curr(rel_pts) - baseline(baseline_pts))/( sigs.fs*pw_inds.pulse_amp(beat_no));
    
    % - diastolic area
    rel_pts = fid_pts.ind.dic(beat_no) : fid_pts.ind.f2(beat_no);
    baseline_pts = rel_pts - fid_pts.ind.f1(beat_no) + 1;
    pw_inds.A2(beat_no) = sum(sigs.curr(rel_pts) - baseline(baseline_pts))/(sigs.fs*pw_inds.pulse_amp(beat_no));

end
clear beat_no rel_beats

% - Ratio of diastolic to systolic area (called Inflection point area)
pw_inds.IPA = pw_inds.A2 ./ pw_inds.A1;

%% First Derivative

% - Maximum slope
pw_inds.ms = fid_pts.amp.ms;

% - Maximum slope divided by the pulse amplitude
pw_inds.ms_div_amp = fid_pts.amp.ms./pw_inds.pulse_amp;

%% Second Derivative

% - Amplitude of 'b' relative to 'a'
pw_inds.b_div_a = fid_pts.amp.b./fid_pts.amp.a;

% - Amplitude of 'c' relative to 'a'
pw_inds.c_div_a = fid_pts.amp.c./fid_pts.amp.a;

% - Amplitude of 'd' relative to 'a'
pw_inds.d_div_a = fid_pts.amp.d./fid_pts.amp.a;

% - Amplitude of 'e' relative to 'a'
pw_inds.e_div_a = fid_pts.amp.e./fid_pts.amp.a;

% - Amplitude of 'a' relative to pulse amplitude
pw_inds.a_div_amp = fid_pts.amp.a./pw_inds.pulse_amp;

% - Amplitude of 'b' relative to pulse amplitude
pw_inds.b_div_amp = fid_pts.amp.b./pw_inds.pulse_amp;

% - Amplitude of 'c' relative to pulse amplitude
pw_inds.c_div_amp = fid_pts.amp.c./pw_inds.pulse_amp;

% - Amplitude of 'd' relative to pulse amplitude
pw_inds.d_div_amp = fid_pts.amp.d./pw_inds.pulse_amp;

% - Amplitude of 'e' relative to pulse amplitude
pw_inds.e_div_amp = fid_pts.amp.e./pw_inds.pulse_amp;

% - Ageing index: original
pw_inds.AGI = pw_inds.b_div_a - pw_inds.c_div_a - pw_inds.d_div_a - pw_inds.e_div_a;

% - Ageing index: informal
pw_inds.AGI_inf = pw_inds.b_div_a - pw_inds.e_div_a;

% - Ageing index: modified
pw_inds.AGI_mod = pw_inds.b_div_a - pw_inds.c_div_a - pw_inds.d_div_a;

% Calculate SIs from second derivative slopes
pt1.t = fid_pts.t.b;
pt1.v = fid_pts.amp.b;
pt2.t = fid_pts.t.c;
pt2.v = fid_pts.amp.c;
pw_inds.slope_b_c = ((pt2.v - pt1.v)./fid_pts.amp.a)./(pt2.t-pt1.t);
pt2.t = fid_pts.t.d;
pt2.v = fid_pts.amp.d;
pw_inds.slope_b_d = ((pt2.v - pt1.v)./fid_pts.amp.a)./(pt2.t-pt1.t);

%% Indices calculated from multiple derivatives

% - Ratio of diastolic to systolic area (called Inflection point area) plus d-peak
pw_inds.IPAD = pw_inds.IPA + pw_inds.d_div_a;

% - Stiffness constant
pw_inds.k = fid_pts.amp.s ./ ((fid_pts.amp.s - fid_pts.amp.ms ) ./ pw_inds.pulse_amp );

%% Calculate median values of each pulse wave index (using only high quality pulse waves)

pw_ind_names = fieldnames(pw_inds);
for pw_ind_no = 1 : length(pw_ind_names)
    
    % extract data for this pulse wave index
    curr_pw_ind = pw_ind_names{pw_ind_no};
    eval(['raw_data = pw_inds.' curr_pw_ind ';']);
    
    % If there are multiple pulse waves then find the median value of this index
    if length(raw_data)==1
        mod_data.v = raw_data;
    else
        mod_data.raw = raw_data;
        mod_data.v = nanmedian(raw_data(pulses.quality(1:end-1)));
    end
    
    % store data
    eval(['pw_inds.' curr_pw_ind ' = mod_data;']);
    
    clear mod_data
end


%% Remove current signal from "sigs"
sigs = rmfield(sigs, 'curr');

end

function deriv = savitzky_golay(sig, deriv_no, win_size)

%% assign coefficients
% From: https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter#Tables_of_selected_convolution_coefficients
% which are calculated from: A., Gorry (1990). "General least-squares smoothing and differentiation by the convolution (Savitzky?Golay) method". Analytical Chemistry. 62 (6): 570?3. doi:10.1021/ac00205a007.

switch deriv_no
    case 0
        % - smoothing
        switch win_size
            case 5
                coeffs = [-3, 12, 17, 12, -3];
                norm_factor = 35;
            case 7
                coeffs = [-2, 3, 6, 7, 6, 3, -2];
                norm_factor = 21;
            case 9
                coeffs = [-21, 14, 39, 54, 59, 54, 39, 14, -21];
                norm_factor = 231;
            otherwise
                error('Can''t do this window size')
        end
    case 1
        % - first derivative
        switch win_size
            case 5
                coeffs = -2:2;
                norm_factor = 10;
            case 7
                coeffs = -3:3;
                norm_factor = 28;
            case 9
                coeffs = -4:4;
                norm_factor = 60;
            otherwise
                error('Can''t do this window size')
        end
        
    case 2
        % - second derivative
        switch win_size
            case 5
                coeffs = [2,-1,-2,-1,2];
                norm_factor = 7;
            case 7
                coeffs = [5,0,-3,-4,-3,0,5];
                norm_factor = 42;
            case 9
                coeffs = [28,7,-8,-17,-20,-17,-8,7,28];
                norm_factor = 462;
            otherwise
                error('Can''t do this window size')
        end
        
    case 3
        % - third derivative
        switch win_size
            case 5
                coeffs = [-1,2,0,-2,1];
                norm_factor = 2;
            case 7
                coeffs = [-1,1,1,0,-1,-1,1];
                norm_factor = 6;
            case 9
                coeffs = [-14,7,13,9,0,-9,-13,-7,14];
                norm_factor = 198;
            otherwise
                error('Can''t do this window size')
        end
        
    case 4
        % - fourth derivative
        switch win_size
            case 7
                coeffs = [3,-7,1,6,1,-7,3];
                norm_factor = 11;
            case 9 
                coeffs = [14,-21,-11,9,18,9,-11,-21,14];
                norm_factor = 143;
            otherwise
                error('Can''t do this window size')
        end
        
    otherwise
        error('Can''t do this order of derivative')        
end

if rem(deriv_no, 2) == 1
    coeffs = -1*coeffs;
end

A = [1,0];
filtered_sig = filter(coeffs, A, sig);
s=length(sig);
half_win_size = floor(win_size*0.5);
deriv=[filtered_sig(win_size)*ones(half_win_size,1);filtered_sig(win_size:s);filtered_sig(s)*ones(half_win_size,1)];
deriv = deriv/norm_factor;

end

function temp_s = identify_s(curr)

[~, temp_s] = max(curr.sig);

end

function temp_a = identify_a(curr, initial_buffer)

[~,filtered_ms] = max(curr.derivs.first);

pks = find_pks_trs(curr.derivs.second, 'pk');
rel_pks = pks(pks > initial_buffer & pks < filtered_ms); 
if isempty(rel_pks) && sum(pks<=initial_buffer) > 0   % Added in PulseAnalyse5
    rel_pks = pks(find(pks <= initial_buffer, 1, 'last'));
end
[~, temp_el] = max(curr.derivs.second(rel_pks));
temp_a = rel_pks(temp_el);

end

function temp_e = identify_e(curr, temp_s, temp_ms, temp_b, end_buffer, redo_log)

% Find local maxima in the second derivative
pks = find_pks_trs(curr.derivs.second, 'pk');
% Set an upper bound of 60 % of the PW duration
upper_bound = 0.6*length(curr.sig);   % const from: https://en.wikipedia.org/wiki/QT_interval#/media/File:QT_interval_corrected_for_heart_rate.png
% Set a lower bound of 'ms'
lower_bound = temp_ms;
% Identify the highest local maximum in the second derivative between these two bounds
rel_pks = pks(pks >= lower_bound & pks <= upper_bound);
[~, max_el] = max(curr.derivs.second(rel_pks));
% If this is the first local maximum in this search region ...
if max_el == 1
    % work out whether this has detected the "c" wave
    % - find out if there's an inflection point between "b" and this
    temp_trs = find_pks_trs(curr.derivs.third, 'tr');
    no_infl = sum(temp_trs > temp_b & temp_trs < rel_pks(max_el));
    % - if not then take the next peak
    if no_infl == 0
        % if there is 1 peak in this search region ...
        if length(rel_pks) < max_el+1   % Added in PulseAnalyse5
            % then take the next peak (slightly above the upper bound
            orig_el = find(pks >= lower_bound & pks <= upper_bound);
            
            %%%%%%%% Added M. Goda, 10th of Jan, 2023
            if orig_el+1>length(pks)
                orig_el=length(pks)-1;
            end
            rel_pks = pks(orig_el:orig_el+1);
        end
        rel_pk = rel_pks(max_el+1);
    else
        rel_pk = rel_pks(max_el);
    end
else
    rel_pk = rel_pks(max_el);
end
temp_e = rel_pk;

end

function temp_f = identify_f(curr, temp_e, end_buffer)

lower_bound = temp_e;
upper_bound = end_buffer;
trs = find_pks_trs(curr.derivs.second, 'tr');
possible_els = trs(trs >= lower_bound & trs <= upper_bound);

if isempty(possible_els)    % Added in PulseAnalyse5
    possible_els = trs(find(trs >=lower_bound, 1));
end

if isempty(possible_els)
    temp_f = nan;
else
    temp_f = possible_els(1);
end

end

function temp_b = identify_b(curr, temp_a, up)

% find b (PulseAnalyse6 onwards)

% Find local minima in second derivative
trs = find_pks_trs(curr.derivs.second, 'tr');
% define an upper bound as 25% of the duration of the signal
upper_bound = 0.25*length(curr.sig);
% find local minima between 'a' and this upper bound
temp = find(trs > temp_a & curr.derivs.second(trs) < 0 & trs < upper_bound);
% Identify the lowest of these local minima
[~, rel_el] = min(curr.derivs.second(trs(temp)));
temp = temp(rel_el);
temp_b = trs(temp); clear temp

if isempty(temp_b)
    if up.options.verbose, fprintf('\n - ''b'' outside range'); end
    
    % find local minima after 'a'
    temp = find(trs > temp_a & curr.derivs.second(trs) < 0, 1, 'first');
    temp_b = trs(temp); clear temp

end

end

function temp_d = identify_d(curr, temp_c, temp_e)

% Identify "d" as the lowest minimum of the second deriv between "c" and "e"
trs = find_pks_trs(curr.derivs.second, 'tr');
possible_trs = find(trs > temp_c & trs < temp_e);
if ~isempty(possible_trs)
    temp = trs(possible_trs);
    [~, temp_el] = min(curr.derivs.second(temp));
    temp_d = temp(temp_el); clear temp
else
    % unless there isn't a minimum, in which case it's an inflection, and
    % "d" is the same as "c"
    temp_d = temp_c;
end

end

function temp_c = identify_c(curr, temp_b, temp_e)

% Identify C as the highest local maximum on the second derivative between "b" and "e"
pks = find_pks_trs(curr.derivs.second, 'pk');
temp = find(pks > temp_b & pks < temp_e);
[~, rel_tr_el] = max(curr.derivs.second(pks(temp)));
temp_c = pks(temp(rel_tr_el)); clear temp rel_tr_el pks

% If there aren't any peaks that satisfy this criterion ...
if isempty(temp_c)
    % then identify C as the lowest local minimum on the third derivative
    % after "b" and before "e"
    trs = find_pks_trs(curr.derivs.third, 'tr');
    temp = find(trs > temp_b & trs < temp_e);
    [~, rel_tr_el] = min(curr.derivs.third(trs(temp)));
    if ~isempty(rel_tr_el)
        temp_c = trs(temp(rel_tr_el)); clear temp rel_tr_el trs
    end
end

end

function temp_dic = identify_dic(curr,temp_e)

temp_dic = temp_e;

end

function temp_dia = identify_dia(curr, temp_dic, temp_e, end_buffer, up)

% if there is a diastolic peak, then use that:
%  -  first peak in signal after "dic"
pks = find_pks_trs(curr.sig, 'pks');
temp_dia = pks(find(pks > temp_dic & pks < end_buffer, 1));

% if not, then ...
% % I tried (i) take first peak on first derivative after "e"
if isempty(temp_dia)
    pks = find_pks_trs(curr.derivs.first, 'pks');
    temp_dia = pks(find(pks > temp_e & pks < end_buffer, 1));
end
% % But the problem is that the first derivative isn't necessarily a peak at
% % the diastolic peak - it can be an inflection point. So:
% (ii) the soonest out of (a) first peak on first derivative after "e"
%                         (b) first min on third derivative after "e"
% if isempty(temp_dia)
%     pks = find_pks_trs(curr.derivs.first, 'pks');
%     temp_dia1 = pks(find(pks > temp_e, 1));
%     trs = find_pks_trs(curr.derivs.third, 'trs');
%     temp_dia2 = trs(find(trs > temp_e, 1));
%     temp_dia = min([temp_dia1, temp_dia2]);
% end

if isempty(temp_dia)
    if up.options.verbose, fprintf('\n - ''dia'' outside range'); end
    pks = find_pks_trs(curr.derivs.first, 'pks');
    temp_dia = pks(find(pks > temp_e, 1, 'first'));
end   

end

function temp_ms = identify_ms(curr)

% find max slope in DPPG
[~, temp_ms] = max(curr.derivs.first);

end

function temp_ms2 = identify_ms2(curr, temp_p1in, temp_p2_pk, up)

% Identify peaks in DPPG
d_pks = find_pks_trs(curr.derivs.first, 'pk');

% See whether any are between p1in and p2pk
rel_pks = find(d_pks > temp_p1in & d_pks < temp_p2_pk);
if ~isempty(rel_pks)
    [~, rel_pk_el] = max(d_pks(rel_pks));
    temp_ms2 = d_pks(rel_pks(rel_pk_el));
else
    temp_ms2 = round(mean([temp_p1in, temp_p2_pk]));
    if up.options.verbose, fprintf("\n Fix this"); end
end

end

function temp_p1 = identify_p1(curr, temp_b, temp_ms, buffer_p1, no_of_pulses, fs, up)

% find p1

% find local minima in the first derivative
fd_trs = find_pks_trs(curr.derivs.first, 'tr');
% find local maxima in the second derivative
sd_pks = find_pks_trs(curr.derivs.second, 'pk');

% Find the first local minimum in the first derivative after 0.1 s
current_buffer = buffer_p1(1);
temp = find(fd_trs > current_buffer,1);
% Find the second local minimum (or failing that, the first) in the first derivative after 'b'
temp2 = find(fd_trs > temp_b, 2);
if length(temp2) > 1
    temp2 = temp2(2);
end
% Take whichever comes first:
if temp2 < temp
    temp = temp2;
end
temp_p1 = fd_trs(temp);

% If this value for p1 is after the buffer of 0.18 s ...
if temp_p1 > buffer_p1(2)
    curr.derivs.fourth = savitzky_golay(curr.derivs.third, 1, 9);
    % Then find the first local minimum in the fourth derivative after 0.1 s
    fd_trs = find_pks_trs(curr.derivs.fourth, 'tr');
    temp = find(fd_trs > current_buffer,1);
    temp_p1 = fd_trs(temp); clear temp
end

% If this value for p1 is after the buffer of 0.18 s ...
if temp_p1 > buffer_p1(2)
    if up.options.verbose, fprintf('\n - P1 outside range'); end
    % Then find the last local minimum in the first derivative before 0.18 s
    temp_p1 = fd_trs(find(fd_trs <= current_buffer,1,'last'));
   
    % If this doesn't find temp_p1, then extend the buffer
    if isempty(temp_p1)
        temp_p1 = fd_trs(find(fd_trs > current_buffer,1,'first'));
    end
end



end

function temp_p2 = identify_p2(curr, temp_d, temp_p1, temp_dic)

% Find "p2" from the minimum value of the third derivative immediately before "d"
td_trs = find_pks_trs(curr.derivs.third, 'tr');
temp = find(td_trs < temp_d,1,'last'); clear d_pks
temp_p2 = td_trs(temp);

% unless c=d, in which case p2 will now be before p1, so take the minimum
% value of the third derivative immediately after "d"
if temp_p2 < temp_p1
    temp_p2 = td_trs(find(td_trs<temp_dic,1,'last'));
end

% check to see if there is a maximum in the signal between the current
% estimate of p2, and temp_dic. If so, take this instead
pks = find_pks_trs(curr.sig, 'pk');
temp = find(pks> temp_p2 & pks < temp_dic);
if length(temp) == 1
    temp_p2 = pks(temp);
elseif length(temp) == 2
    temp_p2 = pks(temp(2));
elseif length(temp) > 1
    fprintf('\nCheck this')
end
clear pks

end

function trs = find_pks_trs(sig,type)

if strcmp(type, 'tr')
    sig = -sig;
end

temp1 = find(sig(2:end-1) > sig(1:(end-2)) & sig(2:end-1) > sig(3:end) );

temp2 = find(sig(2:end-2) > sig(1:(end-3)) & sig(2:end-2)== sig(3:(end-1)) & sig(3:end-1) > sig(4:end) );

temp = unique([temp1; temp2]);

trs = temp+1;

end

function S_elim_lf = eliminate_low_freq_from_single_beat(sig, up)

old = sig.v;

% Correct for low frequency baseline drift in a single beat
correction_line = linspace(sig.v(1), sig.v(end), length(sig.v));
S_elim_lf.v = sig.v - correction_line' + sig.v(1);
S_elim_lf.fs = sig.fs;

end

function [S_aligned, align_el] = align_pulse(sig, up)

% Ensure that signal commences at start of systole

[~, align_el] = min(sig.v);
S_aligned.v = sig.v([align_el:end, 1:(align_el-1)]);
S_aligned.fs = sig.fs;

% add on one additional point so that you can define a second onset
S_aligned.v(end+1) = S_aligned.v(1);

end

function make_plots(sigs, pulse_no, fid_pts, up, plot_inds, tol_samps)
% make plot of individual beat if needed

%% - setup
paper_size = 1.5*[300,1050];
if ~up.options.plot_third_deriv, paper_size(2) = paper_size(2)*0.75; end
if up.options.plot_pw_only, paper_size(2) = paper_size(2)*0.37; end
figure('Position', [50,50, paper_size])
ftsize = 1.5*12; lwidth = 2;
sigs.t = [0:length(sigs.v)-1]/sigs.fs;
sigs.t = sigs.t - (tol_samps/sigs.fs);

y_offset = 0.08;
if up.options.plot_third_deriv
    y_inc = 0.23;
    n_sub = 4;
elseif up.options.plot_pw_only
    y_inc = 0.8;
    n_sub = 1;
    y_offset = 0.17;
else
    y_inc = 0.31;
    n_sub = 3;
end

%% - plot sig
h_b(1) = subplot('Position', [0.21,y_offset+(n_sub-1)*y_inc,0.78,y_inc-0.01]);

% plot baseline curve
plot(sigs.t, sigs.v, 'b', 'LineWidth', lwidth); hold on,

if plot_inds && up.options.plot_areas
    
    % Areas: Systolic
    el = fid_pts.ind.dic(pulse_no) - fid_pts.ind.f1(pulse_no)+1;
    rel_els = 1:el-1;
    offset = 0.05*range(sigs.v);
    h = fill([sigs.t(rel_els), sigs.t(rel_els(end))], [sigs.v(rel_els); sigs.v(rel_els(1))], [1,0.8,0.8]);
    h.LineStyle = 'none';
    
    % Areas: Diastolic
    el2 = fid_pts.ind.f2(pulse_no) - fid_pts.ind.f1(pulse_no)+1;
    rel_els = el+1:el2;
    offset = 0.05*range(sigs.v);
    h = fill([sigs.t(rel_els), sigs.t(rel_els(1))], [sigs.v(rel_els); sigs.v(rel_els(end))], [0.8,0.8,1]);
    h.LineStyle = 'none';
    plot(sigs.t, sigs.v, 'b', 'LineWidth', lwidth)
    
end

% plot salient points
pt_names = {'s', 'dia', 'dic', 'f1', 'f2'}; %, 'p1pk', 'p2pk'};
pt_names = {'s', 'dia', 'dic', 'f1', 'f2', 'p1pk', 'p2pk'};
hspace0 = 0;
for pt_no = 1 : length(pt_names)
    
    curr_text = pt_names{pt_no};
    eval(['curr_pt.el = fid_pts.ind.' pt_names{pt_no} '(pulse_no) - fid_pts.ind.f1(pulse_no)+1+tol_samps;']);
    
    if isnan(curr_pt.el)
        eval(['sigs.pts.' curr_text ' = curr_pt.el;']);
        continue
    end
    
    % plot point
    curr_pt.v = sigs.v(curr_pt.el);
    curr_pt.t = sigs.t(curr_pt.el);
    plot(curr_pt.t, curr_pt.v, 'or')
    
    % annotate point
    vspace0 = 0.12*range(sigs.v);
    ftsize_to_use = ftsize;
    switch curr_text
        case {'dia','p2'}
            hspace0 = 0.06;
            if (strcmp(curr_text, 'p2pk') || strcmp(curr_text, 'p2in')) && curr_pt.el < (fid_pts.s(pulse_no) - fid_pts.f1(pulse_no)+1)
                hspace0 = -0.04;
            elseif (strcmp(curr_text, 'p1pk') || strcmp(curr_text, 'p1in'))
                hspace0 = 0;
            end
        case 'dic'
            hspace0 = -0.01;
            vspace0 = -1*vspace0;
        case {'p1','p1in'}
            hspace0 = -0.04;
        case 's'
            hspace0 = 0.07;
            vspace0 = 0.7*vspace0; 
        case {'f1', 'f2'}
            hspace0 = 0.03;  
            vspace0 = -0.6*vspace0;          
        case {'p1pk'}
            vspace0 = 0;          
            hspace0 = -0.06;
            %ftsize_to_use = ftsize - 4;
        case {'p2pk'}
            vspace0 = 1.0*vspace0;          
            hspace0 = 0.02;
            %ftsize_to_use = ftsize - 4;
    end
    
    lab_txt = strrep(curr_text, 'pk', '');
    if strcmp(lab_txt, 's'), lab_txt = 'sys'; end
    if strcmp(lab_txt, 'f1'), lab_txt = 'onset'; end
    if strcmp(lab_txt, 'f2'), lab_txt = 'end'; end
    
    text(sigs.t(curr_pt.el)+hspace0, sigs.v(curr_pt.el)+vspace0 , lab_txt,'FontSize', ftsize_to_use, 'Color', 'r', 'HorizontalAlignment', 'center');
    eval(['sigs.pts.' curr_text ' = curr_pt.el;']);
end

% set limits
curr_range = range(sigs.v);
if plot_inds
    factor1 = 0.4; factor2 = 0.2;
else
    factor1 = 0.15; factor2 = 0.15;
end
ylims = [min(sigs.v)-factor2*curr_range, max(sigs.v)+factor1*curr_range];

ylim(ylims)
curr_range = range(sigs.t);
xlims = [sigs.t(1)-0.05*curr_range, sigs.t(end)+0.1*curr_range];
if sigs.t(end) > 1.5
    xlim_labs = 0:0.5:0.25*ceil(sigs.t(end)*4);
else
    xlim_labs = 0:0.25:0.25*ceil(sigs.t(end)*4);
end
xlim(xlims)

% set labels
ylab = ylabel({'Pulse', 'Wave'}, 'FontSize', ftsize, 'Rotation', 0);
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(gca, 'FontSize', ftsize -2, 'XTick', xlim_labs, 'XGrid', 'on')
if ~up.options.normalise_pw
    ytick_vals = [round(min(sigs.v),3,'significant'), round(max(sigs.v),3,'significant')];
else
    ytick_vals = [0,1];
end
set(gca, 'YTick', ytick_vals)
box off
if up.options.plot_pw_only
    xlabel('Time (s)', 'FontSize', ftsize)
else
    set(gca, 'XTickLabel', {})
end

% Plot indices
if plot_inds
    
    color = 0.4;
    
    % - delta T
    if ~isnan(sigs.pts.s) && ~isnan(sigs.pts.dia)
        ht_annot = ylims(2)-0.11*range(ylims);
        plot(sigs.t(sigs.pts.s)*[1,1], [ht_annot, sigs.v(sigs.pts.s)], '--', 'color', color*ones(1,3))
        plot(sigs.t(sigs.pts.dia)*[1,1], [ht_annot, sigs.v(sigs.pts.dia)], '--', 'color', color*ones(1,3))
        normalised1  = coords_to_pos(sigs.t(sigs.pts.s), ht_annot);
        normalised2  = coords_to_pos(sigs.t(sigs.pts.dia), ht_annot);
        ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
        ah.Head1Length = 7; ah.Head2Length = 7;
        ah.Head1Width = 7; ah.Head2Width = 7;
        new_ht_annot= ht_annot +0.01*range(ylims);
        text(mean([sigs.t(sigs.pts.s),sigs.t(sigs.pts.dia)]), new_ht_annot, '\DeltaT','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
    
    % - Crest Time
    if ~isnan(sigs.pts.s)
        ht_annot = ylims(2)-0.11*range(ylims);
        plot(sigs.t(sigs.pts.f1)*[1,1], [ht_annot, sigs.v(1)], '--', 'color', color*ones(1,3))
        plot(sigs.t(sigs.pts.s)*[1,1], [ht_annot, sigs.v(sigs.pts.s)], '--', 'color', color*ones(1,3))
        normalised1  = coords_to_pos(sigs.t(sigs.pts.f1), ht_annot);
        normalised2  = coords_to_pos(sigs.t(sigs.pts.s), ht_annot);
        ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
        ah.Head1Length = 7; ah.Head2Length = 7;
        ah.Head1Width = 7; ah.Head2Width = 7;
        new_ht_annot= ht_annot +0.01*range(ylims);
        text(mean([sigs.t(sigs.pts.f1),sigs.t(sigs.pts.s)]), new_ht_annot, 'CT','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
    
    % - Reflection Index
    if ~isnan(sigs.pts.dia)
        
        % smaller amplitude
        t_annot = sigs.t(end)+0.03;
        plot([sigs.t(end), t_annot], [1,1]*sigs.v(end), '--', 'color', color*ones(1,3))
        plot([sigs.t(sigs.pts.dia), t_annot], [1,1]*sigs.v(sigs.pts.dia), '--', 'color', color*ones(1,3))
        normalised1  = coords_to_pos(t_annot, sigs.v(end));
        normalised2  = coords_to_pos(t_annot, sigs.v(sigs.pts.dia));
        ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
        ah.Head1Length = 7; ah.Head2Length = 7;
        ah.Head1Width = 7; ah.Head2Width = 7;
        new_ht_annot= mean([sigs.v(sigs.pts.dia), sigs.v(end)]);
        text(t_annot - 0.05, new_ht_annot, 'h1','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        
        % larger amplitude
        t_annot = sigs.t(end)+0.07;
        plot([sigs.t(end), t_annot], [1,1]*sigs.v(end), '--', 'color', color*ones(1,3))
        plot([sigs.t(sigs.pts.s), t_annot], [1,1]*sigs.v(sigs.pts.s), '--', 'color', color*ones(1,3))
        normalised1  = coords_to_pos(t_annot, sigs.v(end));
        normalised2  = coords_to_pos(t_annot, sigs.v(sigs.pts.s));
        ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
        ah.Head1Length = 7; ah.Head2Length = 7;
        ah.Head1Width = 7; ah.Head2Width = 7;
        temp1 = [sigs.v(sigs.pts.s), sigs.v(end)];
        new_ht_annot= min(temp1) + (2/3)*range(temp1);
        clear temp1
        text(t_annot - 0.05, new_ht_annot, 'h2','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        
        % RI calculation
        ht_annot = ylims(2)-0.11*range(ylims);
        new_ht_annot= ht_annot +0.01*range(ylims);
        text(t_annot - 0.15, new_ht_annot, 'RI = h1/h2','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
    
    do_systole_and_diastole = 0;
    if do_systole_and_diastole
        % - Systolic time
        if ~isnan(sigs.pts.dic)
            ht_annot = ylims(1)+0.02*range(ylims);
            plot(sigs.t(1)*[1,1], [ht_annot, sigs.v(1)], '--', 'color', color*ones(1,3))
            plot(sigs.t(sigs.pts.dic)*[1,1], [ht_annot, sigs.v(sigs.pts.dic)], '--', 'color', color*ones(1,3))
            normalised1  = coords_to_pos(sigs.t(1), ht_annot);
            normalised2  = coords_to_pos(sigs.t(sigs.pts.dic), ht_annot);
            ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
            ah.Head1Length = 7; ah.Head2Length = 7;
            ah.Head1Width = 7; ah.Head2Width = 7;
            new_ht_annot= ht_annot +0.01*range(ylims);
            text(mean([sigs.t(1),sigs.t(sigs.pts.dic)]), new_ht_annot, 'Systole','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
        
        % - Diastolic time
        if ~isnan(sigs.pts.dic)
            ht_annot = ylims(1)+0.02*range(ylims);
            plot(sigs.t(end)*[1,1], [ht_annot, sigs.v(1)], '--', 'color', color*ones(1,3))
            plot(sigs.t(sigs.pts.dic)*[1,1], [ht_annot, sigs.v(sigs.pts.dic)], '--', 'color', color*ones(1,3))
            normalised1  = coords_to_pos(sigs.t(sigs.pts.dic), ht_annot);
            normalised2  = coords_to_pos(sigs.t(end), ht_annot);
            ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
            ah.Head1Length = 7; ah.Head2Length = 7;
            ah.Head1Width = 7; ah.Head2Width = 7;
            new_ht_annot= ht_annot +0.01*range(ylims);
            text(mean([sigs.t(sigs.pts.dic),sigs.t(end)]), new_ht_annot, 'Diastole','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
    end
    
    if up.options.plot_areas
        % - Systolic area
        ht_annot = mean([sigs.v(1), sigs.v(sigs.pts.s)]);
        text(mean([sigs.t(sigs.pts.dic),sigs.t(1)]), ht_annot, 'A_s','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        
        % - Diastolic area
        ht_annot = mean([sigs.v(1), sigs.v(sigs.pts.s)]);
        text(0.71*sigs.t(sigs.pts.s)+0.29*sigs.t(end), ht_annot, 'A_d','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    
end

if ~up.options.plot_pw_only
    
    %% - plot first derivative
    h_b(2) = subplot('Position', [0.21,y_offset+(n_sub-2)*y_inc,0.78,y_inc - 0.01]);
    
    % Plot x-axis
    plot([-10, 10], [0,0], '--k'); hold on
    
    % plot baseline curve
    plot(sigs.t, sigs.first_d, 'b', 'LineWidth', lwidth); hold on,
    
    % plot salient points
    pt_names = {'ms'}; %, 'ms2', 'dia'};
    curr_range = range(sigs.first_d);
    vspace0 = 0.1*curr_range;
    hspace0 = 0.05;
    for pt_no = 1 : length(pt_names)
        
        curr_text = pt_names{pt_no};
        eval(['curr_fid_pt = fid_pts.ind.' pt_names{pt_no} '(pulse_no)+tol_samps;'])
        
        if isnan(curr_fid_pt)
            eval(['sigs.pts.' curr_text ' = nan;']);
            continue
        end
        curr_pt.el = curr_fid_pt - fid_pts.ind.f1(pulse_no)+1;
        if isnan(curr_pt.el)
            eval(['sigs.pts.' curr_text ' = nan;']);
            continue
        end
        
        % plot point
        curr_pt.v = sigs.first_d(curr_pt.el);
        curr_pt.t = sigs.t(curr_pt.el);
        plot(curr_pt.t, curr_pt.v, 'or')
        
        % annotate point
        text(sigs.t(curr_pt.el)+hspace0, sigs.first_d(curr_pt.el)+vspace0 , curr_text,'FontSize', ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
        eval(['sigs.pts.' curr_text ' = curr_pt.el;']);
        
    end
    
    % set limits
    curr_range = range(sigs.t);
    xlim(xlims)
    curr_range = range(sigs.first_d);
    ylim([min(sigs.first_d)-0.05*curr_range, max(sigs.first_d)+0.15*curr_range]); ylims = ylim;
    
    % set labels
    ylab = ylabel({'1st','deriv'}, 'FontSize', ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
    set(gca, 'FontSize', ftsize -2, 'YTick', [])
    set(gca, 'XTick', xlim_labs, 'XTickLabel', {}, 'XGrid', 'on')
    if ~up.options.normalise_pw
        ytick_vals = [0, round(max(sigs.first_d),3,'significant')];
        set(gca, 'YTick', ytick_vals)
    else
        set(gca, 'YTick', 0);
    end
    box off
    
    % Plot indices
    if plot_inds
        
        % - ms
        if ~isnan(sigs.pts.ms) && ~isnan(sigs.pts.dic)
            ht_annot = mean([sigs.first_d(sigs.pts.ms), min(sigs.first_d)]);
            [~, temp] = min(sigs.first_d(sigs.pts.ms:sigs.pts.dic));
            el = temp + sigs.pts.ms -1;
            plot([sigs.t(sigs.pts.ms), sigs.t(el)], sigs.first_d(sigs.pts.ms)*[1,1], '--', 'color', color*ones(1,3))
            %     plot([sigs.t(sigs.pts.ms), sigs.t(el)], sigs.first_d(el)*[1,1], '--', 'color', color*ones(1,3))
            normalised1  = coords_to_pos(sigs.t(el), 0);
            normalised2  = coords_to_pos(sigs.t(el), sigs.first_d(sigs.pts.ms));
            ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
            curr_range = range(sigs.t);
            text(sigs.t(el)+0.06*curr_range, ht_annot, 'ms','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
    end
    
    %% - plot Second derivative
    h_b(3) = subplot('Position', [0.21,y_offset+(n_sub-3)*y_inc,0.78,y_inc - 0.01]);
    
    % Plot x-axis
    plot([-10, 10], [0,0], '--k'); hold on
    
    % plot baseline curve
    curr_color = 0.4*0.5*[1,2,1];
    plot(sigs.t, sigs.second_d, 'b', 'LineWidth', lwidth); hold on,
    
    % plot salient points
    if up.options.plot_f_pt
        pt_names = {'a', 'b', 'c', 'd', 'e', 'f'};
    else
        pt_names = {'a', 'b', 'c', 'd', 'e'};
    end
    curr_range = range(sigs.second_d);
    vspace_const = 0.08;
    hspace0 = 0;
    for pt_no = 1 : length(pt_names)
        eval(['curr_pt.el = fid_pts.ind.' pt_names{pt_no} '(pulse_no) - fid_pts.ind.f1(pulse_no)+1+tol_samps;']);
        curr_text = pt_names{pt_no};
        if isnan(curr_pt.el)
            eval(['sigs.pts.' curr_text ' = nan;']);
            continue
        end
        
        % plot point
        curr_pt.v = sigs.second_d(curr_pt.el);
        curr_pt.t = sigs.t(curr_pt.el);
        plot(curr_pt.t, curr_pt.v, 'or')
        
        % annotate point
        switch curr_text
            case {'a','c', 'e'}
                vspace0 = (1.3*vspace_const)*curr_range;
            case {'b', 'd', 'f'}
                vspace0 = -1.3*vspace_const*curr_range;
        end
        text(sigs.t(curr_pt.el), sigs.second_d(curr_pt.el)+vspace0 , curr_text,'FontSize', ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
        eval(['sigs.pts.' curr_text ' = curr_pt.el;']);
    end
    
    % set limits
    curr_range = range(sigs.t);
    xlim(xlims)
    curr_range = range(sigs.second_d);
    ylim([min(sigs.second_d)-0.15*curr_range, max(sigs.second_d)+0.15*curr_range])
    
    % set labels
    ylab = ylabel({'2nd','deriv'}, 'FontSize', ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.65, 0]);
    set(gca, 'FontSize', ftsize -2, 'YTick', [], 'XTick', xlim_labs, 'XGrid', 'on')
    box off
    if ~up.options.plot_third_deriv
        xlabel('Time (s)', 'FontSize', ftsize)
    else
        set(gca, 'XTickLabel', {})
    end
    if ~up.options.normalise_pw
        ytick_vals = [round(min(sigs.second_d),3,'significant'), 0, round(max(sigs.second_d),3,'significant')];
        set(gca, 'YTick', ytick_vals)
    else
        set(gca, 'YTick', 0);
    end
    
    % Plot indices
    if plot_inds
        
        % - slope b-d
        if ~isnan(sigs.pts.b) && ~isnan(sigs.pts.d)
            plot(sigs.t([sigs.pts.b,sigs.pts.d]), sigs.second_d([sigs.pts.b, sigs.pts.d]), 'k', 'LineWidth', lwidth),
            curr_range = range(sigs.t);
            x = mean(sigs.t([sigs.pts.b,sigs.pts.d]));
            y = mean(sigs.second_d([sigs.pts.b,sigs.pts.b,sigs.pts.b,sigs.pts.d]));
            text(x+0.03*curr_range, y, 'slope_{b-d}','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        end
    end
    
    %% - plot Third derivative
    if up.options.plot_third_deriv
        h_b(4) = subplot('Position', [0.21,y_offset+(n_sub-4)*y_inc,0.78,y_inc-0.01]);
        
        % Plot x-axis
        plot([-10, 10], [0,0], '--k'); hold on
        
        % plot baseline curve
        curr_color = 0.4*0.5*[1,2,1];
        plot(sigs.t, sigs.third_d, 'b', 'LineWidth', lwidth); hold on,
        
        % plot salient points
        pt_names = {'p1pk', 'p2pk'};
        curr_range = range(sigs.third_d);
        vspace_const = 0.08;
        hspace0 = 0;
        [~, temp] = max(sigs.third_d);
        fid_pts.ind.p1pk(pulse_no) = temp + fid_pts.ind.f1(pulse_no) - 1 - tol_samps;
        for pt_no = 1 : length(pt_names)
            eval(['curr_pt.el = fid_pts.ind.' pt_names{pt_no} '(pulse_no) - fid_pts.ind.f1(pulse_no)+1+tol_samps;']);
            if isnan(curr_pt.el)
                continue
            end
            
            % plot point
            curr_pt.v = sigs.third_d(curr_pt.el);
            curr_pt.t = sigs.t(curr_pt.el);
            plot(curr_pt.t, curr_pt.v, 'or')
            
            % annotate point
            curr_text = pt_names{pt_no};
            switch curr_text
                case {'p1pk'}
                    vspace0 = vspace_const*curr_range;
                    curr_text = 'p1';
                case {'p2pk'}
                    vspace0 = -1*vspace_const*curr_range;
                    curr_text = 'p2';
            end
            text(sigs.t(curr_pt.el), sigs.third_d(curr_pt.el)+vspace0 , curr_text,'FontSize', ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
            
        end
        
        % set limits
        curr_range = range(sigs.t);
        xlim(xlims)
        ylims = ylim; curr_range = range(sigs.third_d);
        ylim([min(sigs.third_d)-0.05*curr_range, max(sigs.third_d)+0.15*curr_range])
        
        % set labels
        xlabel('Time (s)', 'FontSize', ftsize)
        ylab = ylabel({'3rd','deriv'}, 'FontSize', ftsize, 'Rotation', 0);
        set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
        set(gca, 'FontSize', ftsize -2, 'YTick', [], 'XTick', xlim_labs, 'XGrid', 'on')
        box off
    end
    
end

linkaxes(h_b, 'x')
shg

if ~isempty(up.options.save_folder)
    savepath = [up.options.save_folder, up.options.save_file];
    if ~isempty(up.options.save_file)
        savepath = [savepath, '_'];
    end
    if ~plot_inds
        savepath = [savepath, 'FidPts'];
    else
        savepath = [savepath, 'PWInds'];
    end
    set(gcf,'color','w');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize', [paper_size(1), paper_size(2)]./40);
    set(gcf,'PaperPosition',[0 0 paper_size(1) paper_size(2)]./40);
%     print(gcf,'-dpdf',savepath)
    print(gcf,'-depsc',savepath)
    print(gcf,'-dpng',savepath)
    print(gcf,'-dsvg',savepath)
end
shg

end

function normalised  = coords_to_pos(x_coord, y_coord)

pos = get(gca, 'Position');
normalised(1) = (x_coord - min(xlim))/diff(xlim) * pos(3) + pos(1);
normalised(2) = (y_coord - min(ylim))/diff(ylim) * pos(4) + pos(2);
 
end

function subtracted = subtract_baseline(sig)

baseline = linspace(sig(1),sig(end),length(sig));

subtracted = sig(:) - baseline(:); 

end

function subtracted = subtract_baseline_extended(sig, tol_samps)

baseline = linspace(sig(1+tol_samps),sig(end),length(sig)-(tol_samps));
baseline = interp1(1+tol_samps:length(sig), baseline, 1:length(sig),'linear','extrap');

subtracted = sig(:) - baseline(:); 

end

function [norm, scale_factor] = normalise(sig)

norm = sig - min(sig); 
scale_factor = max(norm);
norm = norm / scale_factor;

end

function plot_results(sigs, pulses, fid_pts, pw_inds, up)

% Skip if not calculating pulse wave indices
if up.options.calc_pw_inds == 0
    return
end

% skip if not plotting
if ~up.options.do_plot || sum(isnan(fid_pts.ind.a)) == length(fid_pts.ind.a)
    return
end

% extract data for the first high quality pulse wave
pulse_no = 1; tol_samps = 0;
if ~up.options.calc_average_pw
    if ~pulses.quality(pulse_no)
        temp = find(pulses.quality(1:end-1) & ~isnan(fid_pts.ind.a), 1);
        if ~isempty(temp)
            pulse_no = temp;
        end
    end
    % added to include a bit before and after the pulse wave
    if pulse_no > 1
        tol = 0.06;
        tol_samps = round(tol*sigs.fs);
        curr_els = pulses.onsets(pulse_no)-tol_samps:pulses.onsets(pulse_no+1);
        sigs.v = sigs.curr(curr_els);
        sigs.v = subtract_baseline_extended(sigs.v, tol_samps);
        fields = fieldnames(fid_pts.ind);
        for field_no = 1 : length(fields)
            curr_field = fields{field_no};
            %eval(['fid_pts.ind.' curr_field ' = fid_pts.ind.' curr_field '-tol_samps;']);
        end
    else
        curr_els = pulses.onsets(pulse_no):pulses.onsets(pulse_no+1);
        sigs.v = sigs.curr(curr_els);
        sigs.v = subtract_baseline(sigs.v);
    end
else
    curr_els = pulses.ave.onsets(1):pulses.ave.onsets(2);
    sigs.v = sigs.ave(curr_els);
    sigs.v = subtract_baseline(sigs.v);
end

sigs.first_d = sigs.first_d(curr_els);
sigs.second_d = sigs.second_d(curr_els);
sigs.third_d = sigs.third_d(curr_els);

% - subtract baseline from this pulse wave, and normalise (if specified)
if up.options.normalise_pw
    [sigs.v, scale_factor] = normalise(sigs.v);
    sigs.first_d = sigs.first_d./scale_factor;
    sigs.second_d = sigs.second_d./scale_factor;
    sigs.third_d = sigs.third_d./scale_factor;
end

% - Plot fiducial points
make_plots(sigs, pulse_no, fid_pts, up, 0, tol_samps)

% - Plot pulse wave indices
make_plots(sigs, pulse_no, fid_pts, up, 1, tol_samps)

end

function make_pre_processing_plot(sigs, pulses, up)

% Skip this step if not required
if ~up.options.do_plot
     return
end

% Setup figure
ftsize = 16;
lwidth = 2;
paper_size = [1600,900];
colors.all = 0.4*[1,1,1];
colors.used = [0,0,1];
colors.ave_uncalib = [1,0,0];
colors.smooth = [1,0,0];
colors.fill = 0.8*[1,1,1];
colors.beats = [1,0,0];
colors.tf = 0.6*[1,0,1];
offset.t = 0.6;
offset.v = 0.7;
figure('Position', [20,20,paper_size]);

%% raw signal plot
subplot(2,4,1:4)
sigs.t = [0:length(sigs.orig)-1]/sigs.fs;
h(1) = plot(sigs.t, sigs.orig, 'LineWidth', lwidth, 'Color', colors.all); hold on
counter = 1;
if sum(contains(fieldnames(pulses),'orig')), pulses_to_use = pulses.orig; else pulses_to_use = pulses; end
for pw_no = 1 : length(pulses_to_use.quality(1:end-1))
    if pulses_to_use.quality(pw_no)
        start_el = pulses_to_use.onsets(pw_no);
        end_el = pulses_to_use.onsets(pw_no+1);
        if counter == 1
            h(2) = plot(sigs.t(start_el:end_el), sigs.orig(start_el:end_el), 'Color', colors.used, 'LineWidth', lwidth);
        else
            plot(sigs.t(start_el:end_el), sigs.orig(start_el:end_el), 'Color', colors.used, 'LineWidth', lwidth);
        end
    end    
end
h(3) = plot(sigs.t(pulses_to_use.onsets), sigs.orig(pulses_to_use.onsets), 'or', 'Color', colors.beats, 'MarkerFaceColor', colors.beats, 'LineWidth', 2);
% tidy up
box off
xlabel('Time (s)', 'FontSize', ftsize)
ylabel('BP signal (au)', 'FontSize', ftsize)
ylim([min(sigs.orig)-0.03*range(sigs.orig), max(sigs.orig)+0.32*range(sigs.orig)])
set(gca, 'FontSize', ftsize, 'YTick', [])
if sum(pulses.quality)<2
    legend(h([1,3]), {'Raw signal', 'Pulse onsets'}, 'Location', 'northwest')
else
    legend(h, {'Raw signal', 'High quality PWs', 'Pulse onsets'}, 'Location', 'northwest')
end
title('(a) Identifying high quality pulse waves (PWs)', 'FontSize', ftsize)


%% Filtering plot

if up.options.do_filter
    
    % making smooth pulse wave plot
    subplot(2,4,5)
    rel_pulse = find(pulses.quality,1);
    pw.raw = sigs.orig(pulses.onsets(rel_pulse):pulses.onsets(rel_pulse+1));
    pw.filt = sigs.filt(pulses.onsets(rel_pulse):pulses.onsets(rel_pulse+1));
    pw.t = [0:length(pw.raw)-1]/sigs.fs;
    
    % plot portion to be blown up
    t_start =0.35; t_end = 0.43;
    rel_els = find(pw.t> t_start & pw.t < t_end);
    wav = pw.raw;
    x1 = pw.t(rel_els);
    y1 = wav(rel_els);
    wav = pw.raw;
    x2 = pw.t(rel_els);
    y2 = wav(rel_els);
    x = [min([x1(:); x2(:)]), min([x1(:); x2(:)]), max([x1(:); x2(:)]), max([x1(:); x2(:)])];
    y = [min([y1(:); y2(:)]), max([y1(:); y2(:)]), max([y1(:); y2(:)]), min([y1(:); y2(:)])];
    fill(x,y,colors.fill,'LineStyle', 'none'), hold on
    
    % plot whole PW
    wav = pw.raw;
    h_ave = plot(pw.t, wav, 'LineWidth', lwidth, 'Color', colors.used);
    hold on
    wav = pw.filt;
    h_sm = plot(pw.t, wav, 'LineWidth', lwidth, 'Color', colors.smooth);
    
    % plot blown up portion
    sc_factor = 2;
    wav = pw.raw;
    temp_offset.t = 0.3*(pw.t(end)-t_end);
    temp_offset.v = 0.1*(max(wav)-max(wav(rel_els)));
    x1 = pw.t(rel_els)+sc_factor*(pw.t(rel_els)-pw.t(rel_els(1)))+temp_offset.t;
    y1 = wav(rel_els)+sc_factor*(wav(rel_els)-min(wav(rel_els)))+temp_offset.v;
    wav = pw.filt;
    x2 = pw.t(rel_els)+sc_factor*(pw.t(rel_els)-pw.t(rel_els(1)))+temp_offset.t;
    y2 = wav(rel_els)+sc_factor*(wav(rel_els)-min(wav(rel_els)))+temp_offset.v;
    xa = [min([x1(:); x2(:)]), min([x1(:); x2(:)]), max([x1(:); x2(:)]), max([x1(:); x2(:)])];
    ya = [min([y1(:); y2(:)]), max([y1(:); y2(:)]), max([y1(:); y2(:)]), min([y1(:); y2(:)])];
    fill(xa,ya,colors.fill,'LineStyle', 'none')
    wav = pw.raw;
    plot(x1, y1, 'LineWidth', lwidth, 'Color', colors.used);
    wav = pw.filt;
    plot(x2, y2, 'LineWidth', lwidth, 'Color', colors.smooth);
    
    % plot lines between portions
    plot([x(2),xa(2)],[y(2),ya(2)],'Color', 0.7*colors.fill);
    plot([x(4),xa(4)],[y(4),ya(4)],'Color', 0.7*colors.fill);
    
    % tidy up
    box off
    % xlim([0, pw_t(end)])
    xlabel('Time (s)', 'FontSize', ftsize)
    % ylim([min_x-0.1*(max_x-min_x), max_x+0.1*(max_x-min_x)])
    ylabel('Pulse waves (au)', 'FontSize', ftsize)
    % ylim([min(s.v)-0.05*range(s.v), max(s.v)+0.3*range(s.v)])
    set(gca, 'FontSize', ftsize, 'YTick', [])
    legend([h_ave, h_sm], {'Original PW', 'Filtered PW'}, 'Location', 'northeast')
    ylim([min(wav)-0.05*range(wav), max(wav)+0.7*range(wav)])
    xlim([0 pw.t(end)])
    title('(b) Filtering', 'FontSize', ftsize)
    
end

%% Calibration plot
if up.options.calibrate_pw
    
    % making calibration plot
    subplot(2,4,6)
    
    % extract calibrated PW
    rel_pulse = find(pulses.quality,1);
    pw.calib = sigs.calib(pulses.onsets(rel_pulse):pulses.onsets(rel_pulse+1));
    
    % Extract calibration BPs
    if sum(strcmp(fieldnames(up.options),'calib_map'))
        mbp = up.options.calib_map;
        dbp = up.options.calib_dbp;
        cal_line = [mbp,mbp];
        sbp = max(pw.calib);
    else
        sbp = up.options.calib_sbp;
        dbp = up.options.calib_dbp;
        cal_line = [sbp,sbp];
        mbp = mean(pw.calib);
    end
    
    
    plot(pw.t, pw.calib, 'LineWidth', lwidth, 'Color', colors.smooth);
    hold on
    ylim([50, 140])
    
    % calibration lines
    plot([0, pw.t(end)],[dbp,dbp],'--','Color',colors.all,'LineWidth',lwidth)
    rel_el = floor(0.5*length(pw.calib));
    plot([0, pw.t(rel_el)], cal_line,'--','Color',colors.all,'LineWidth',lwidth)
    if sum(strcmp(fieldnames(up.options),'calib_map'))
        dim = [.36 .18 .1 .1];
        str = 'MBP';
        annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontSize', ftsize);
    else
        dim = [.41 .34 .1 .1];
        str = 'SBP';
        annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontSize', ftsize);
    end
    dim = [.41 .08 .1 .1];
    str = 'DBP';
    annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','FontSize', ftsize);
    
    % tidy up
    box off
    xlim([0, pw.t(end)])
    xlabel('Time (s)', 'FontSize', ftsize)
    ylim([dbp-0.1*(sbp-dbp), sbp+0.1*(sbp-dbp)])
    ylabel('Pulse wave (mmHg)', 'FontSize', ftsize)
    set(gca,'FontSize', ftsize)
    title('(c) Calibrating', 'FontSize', ftsize)
end

%% Transformation plot
if up.options.tran_func
    subplot(2,4,7)
    
    % extract pulse wave data
    rel_pulse = find(pulses.quality,1);
    pw.tf = sigs.tf(pulses.onsets(rel_pulse):pulses.onsets(rel_pulse+1));
    
    % Plot previous PW
    if up.options.calibrate_pw
        pw.tf = pw.tf - min(pw.tf);
        pw.tf = ( (mbp-dbp) * (pw.tf/mean(pw.tf)) ) + dbp;
        ref_pw = pw.calib;
        lab = 'Calibrated PW';
    elseif up.options.do_filter
        ref_pw = pw.filt;
        lab = 'Filtered PW';  
    else
        ref_pw = pw.raw;
        lab = 'Raw PW';  
    end
    plot(pw.t, ref_pw, 'LineWidth', lwidth, 'Color', colors.smooth);
    if ~up.options.calibrate_pw
        dbp = min(ref_pw); mbp = mean(ref_pw);
        pw.tf = pw.tf - min(pw.tf);
        pw.tf = ( (mbp-dbp) * (pw.tf/mean(pw.tf)) ) + dbp;
    end
    hold on
    plot(pw.t, pw.tf, 'LineWidth', lwidth, 'Color', colors.tf);
    ylim([50, 140])
    
    % tidy up
    box off
    xlim([0, pw.t(end)])
    xlabel('Time (s)', 'FontSize', ftsize)
    ylim([min(pw.tf)-0.1*(max(pw.tf)-min(pw.tf)), max(pw.tf)+0.4*(max(pw.tf)-min(pw.tf))])
    if up.options.calibrate_pw
        ylabel('Pulse wave (mmHg)', 'FontSize', ftsize)
    else
        ylabel('Pulse wave', 'FontSize', ftsize)
        set(gca, 'YTick', [])
    end
    set(gca,'FontSize', ftsize)
    title('(d) Transforming', 'FontSize', ftsize)
    legend({lab, 'Aortic PW'}, 'Location', 'northeast')
end

%% Average plot

if up.options.calc_average_pw

    subplot(2,4,8)
    clear h
    
    % plot PWs
    t = [0:size(sigs.pws,1)-1]/sigs.fs;
    h = plot(t, sigs.pws, 'Color', 1.2*colors.tf, 'LineWidth', 1);
    h = h(1);
    
%     if pw_no == length(rel_pulses)-1
%         h(1) = plot(t, pw.v, 'Color', 1.2*colors.tf, 'LineWidth', 1);
%     else
%         plot(t, pw.v, 'Color', 1.2*colors.tf, 'LineWidth', 1);
%     end
    hold on
    
    % plot ave PW (if there were any high quality PWs)
    h(2) = plot(t, sigs.ave, 'Color', 'k', 'LineWidth', 2);
        
    % Tidy up
    box off
    xlim([0, t(end)])
    xlabel('Time (s)', 'FontSize', ftsize)
    pws_lims = [min(min(sigs.pws)), max(max(sigs.pws))];
    ylim([pws_lims(1)-0.05*range(pws_lims), pws_lims(2)+0.4*range(pws_lims)])
    if up.options.calibrate_pw
        ylabel('Pulse wave (mmHg)', 'FontSize', ftsize)
    else
        ylabel('Pulse wave', 'FontSize', ftsize)
    end
    set(gca,'FontSize', ftsize)
    title('(e) Averaging', 'FontSize', ftsize)
    
    % Make legend
    
    if up.options.calibrate_pw
        lab = 'Calibrated';
    elseif up.options.do_filter
        lab = 'Filtered';  
    else
        lab = 'Raw';  
    end
    legend(h, {[lab, ' PWs'], 'Average PW'}, 'Location', 'northeast')
    
%     else
%         dim = [.8 .05 .2 .2];
%         str = {'No high quality', 'pulse waves'};
%         annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none', 'FontSize', ftsize);
%     end
    
end

% save figure
if ~isempty(up.options.save_folder)
    savepath = [up.options.save_folder, up.options.save_file];
    if ~isempty(up.options.save_file)
        savepath = [savepath, '_pre_processing'];
    end
    set(gcf,'color','w');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize', [paper_size(1), paper_size(2)]./40);
    set(gcf,'PaperPosition',[0 0 paper_size(1) paper_size(2)]./40);
    print(gcf,'-depsc',savepath)
    print(gcf,'-dpng',savepath)
    
    % save text file
    fid = fopen([savepath, '.txt'], 'w');
    p = mfilename('fullpath');
    p = strrep(p, '\', '\\');
    if exist('today') % requires Financial Toolbox
        fprintf(fid, ['Figures generated by:\n\n ' p '.m \n\n on ' datestr(today)]);
    else
        fprintf(fid, ['Figures generated by:\n\n ' p '.m \n\n']);
    end
    fclose all;
end

end

function [ peaks, sig_filt, peaks_filt, thres ] = PPG_pulses_detector( ppg, fsppg, pb, lpd_fp, lpd_fc, lpd_order, alpha, refract, taoRR, w_nA, plotflag )
% -- Source --
% This function was copied from the 'ecg-kit' repository by 'marianux' at: https://github.com/marianux/ecg-kit . See also: http://marianux.github.io/ecg-kit/ 
% This particular function was available at: https://github.com/marianux/ecg-kit/blob/master/common/ppg/PPG_pulses_detector.m
% It is available under the GNU General Public License v2.0: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
% ------------

% added by PC:
ignore_pb = 1; % skips all lines involving the 'pb' variable

%PPG_PULSES_DETECTOR     PPG signal pulses detector, based on low-pass
%                        differentiator filter.
%
% Created by Jesus Lazaro <jlazarop@unizar.es> in 2012
% adapted by Mariano Llamedo Soria to the ecg-kit project.
%--------
%   Sintax: [ peaks, sig_filt, peaks_filt, thres ] = PPG_pulses_detector( ppg, fsppg, lpd_fp, lpd_fc, lpd_order, alpha, refract, taoRR, w_nA, plotflag )
%   In:   ppg = PPG signal
%         fsppg = "ppg" sampling rate (Hz)
%         lpd_fp = last frequency of the pass-band for the
%                  low-pass-differentiator filter (Hz) [Default: 7.8]
%         lpd_fc = cut-off frequency for the low-pass-differentiator filter
%                  (Hz) [Default: 8]
%         lpd_order = filter order [Default: 3*fsppg]
%         alpha = multiplies previus amplitude of detected maximum in
%                 filtered signal for updating the threshold [Default: 0.2]
%         refract = refractary period for threshold (s) [Default: 150e-3]
%         taoRR = fraction of estimated RR where threshold reaches its
%                 minimum value (alpha*amplitude of previous SSF peak)
%                 [Default: 1]
%         w_nA = length of the window after a peak in filtered signal in
%                which is asociated peak in PPG will be searched (s)
%                [Default: 300e-3]
%         plotflag = if 1, plots a figure with PPG and SSF [Default: 0]
%
%   Out:  peaks = location of maximum of detected pulses (samples)
%         sig_filt = band-pass filtered signal
%         peaks_filt = location of peaks detected in filtered signal (samples)
%         thres = computed time varying theshold
    
    if nargin<2
        error('Not enough input arguments');
    end
    
    if nargin<3
        if ~ignore_pb
            pb = [];
        end
    end
    
    if nargin<4
        lpd_fp = 7.8;
    end
    
    if nargin<5
        lpd_fc = 8;
    end
    
    if nargin<6
        lpd_order = 3*fsppg;
    end
    
    if nargin<7
        alpha = 0.2;
    end
    
    if nargin<8
        refract = 150e-3;
    end
    
    if nargin<9
        taoRR = 1;
    end
    
    if nargin<10
         w_nA = 300e-3;
    end
    
    if nargin<11
        plotflag = 0;
    end
    
    
    ppg = ppg(:); %Ensure "ppg" is a column vector
    refract = round(refract*fsppg); %Seconds->samples
    
    %% Filtering:
    
    this_path = fileparts(mfilename('fullpath'));
    % default folder to look at
    cached_filter_filename = [this_path filesep sprintf( 'low_pass_differentiator_%dHz.mat', fsppg) ];
    if( exist(cached_filter_filename, 'file') )
        pepe = load(cached_filter_filename);
        bb = pepe.bb;
        clear pepe
    else
        % force an integer group delay
        if( rem(lpd_order,2) ~= 0 )
            lpd_order = lpd_order + 1;
        end
        d = fdesign.differentiator('n,fp,fst', lpd_order, lpd_fp*2/fsppg, lpd_fc*2/fsppg);
        % hd = design(d,'equiripple');
        hd = design(d, 'firls');
        bb = hd.Numerator*fsppg/(2*pi);
        save(cached_filter_filename, 'bb');
        clear d hd;
    end    
    
    %filter the signal
    delay = (numel(bb)-1)/2;
    sig_filt = filter(bb, 1, ppg);
    sig_filt = [sig_filt(1+delay:end); zeros(delay, 1)];
    
    
    %% Compute threshold for filtered signal:
    if ~ignore_pb
        pb.reset();
        pb.checkpoint('Compute threshold for filtered signal');
    end
    
    peaks_filt = [];
    thres_ini_w_ini = find(~isnan(sig_filt), 1, 'first');
    thres_ini_w_end = thres_ini_w_ini + round(10*fsppg);
    aux = sig_filt(thres_ini_w_ini:thres_ini_w_end);
    thres_ini = 3*mean(aux(aux>=0));
    thres = nan(size(sig_filt));
    t = 1:length(sig_filt);
    RR = round(60/80*fsppg);
    if (1+RR)<length(sig_filt)
        thres(1:1+RR) = thres_ini - (thres_ini*(1-alpha)/RR)*(t(1:RR+1)-1);
        thres(1+RR:end) = alpha*thres_ini;
    else
        thres(1:end) = thres_ini - (thres_ini*(1-alpha)/RR)*(t(1:end)-1);
    end
    
    if ~ignore_pb
        pb.Loops2Do = round(length(sig_filt) / RR);
    end
    
    kk=1;
    while true
        if ~ignore_pb
            pb.start_loop();
        end
        
        cross_u = kk-1 + find(sig_filt(kk:end)>thres(kk:end), 1, 'first'); %Next point to cross the actual threshold (down->up)
        if isempty(cross_u)
            % No more pulses -> end
            break;
        end
        
        cross_d = cross_u-1 + find(sig_filt(cross_u:end)<thres(cross_u:end), 1, 'first'); %Next point to cross the actual threshold (up->down)
        
        if isempty(cross_d)
            % No more pulses -> end
            break;
        end
        
        % Pulse detected:
        [vmax, imax] = max(sig_filt(cross_u:cross_d));
        p = cross_u-1+imax;
        peaks_filt = [peaks_filt, p];
        
        if ~ignore_pb
            pb.checkpoint([]);
        end
        
        % Update threshold
        N_RR_estimation = 3;
        N_ampli_est = 3;
        Npeaks = length(peaks_filt);
        if Npeaks>=N_RR_estimation+1;
            RR = round(median(diff(peaks_filt(end-N_RR_estimation:end))));
        elseif Npeaks>=2
            RR = round(mean(diff(peaks_filt)));
        end
        kk = min(p+refract, length(sig_filt));
        thres(p:kk) = vmax;
%         tao = 5/(taoRR*RR-refract);
%         thres(kk:end) = vmax*(1-alpha)*exp(-tao*(t(kk:end)-kk)) + vmax*alpha;
        
        if ~ignore_pb
            pb.checkpoint([]);
        end

        vfall = vmax*alpha;
        if Npeaks>=(N_ampli_est+1)
            ampli_est = median(sig_filt(peaks_filt(end-N_ampli_est:end-1)));
            if vmax>=(2*ampli_est)
                vfall = alpha*ampli_est;
                vmax = ampli_est;
            end
%             if vmax<=(0.5*ampli_est)
%                 vfall = alpha*ampli_est;
%                 vmax = ampli_est;
%             end
        end
        
        fall_end = round(taoRR*RR);
        if (kk+fall_end)<length(sig_filt)
            thres(kk:kk+fall_end) = vmax - (vmax-vfall)/fall_end*(t(kk:kk+fall_end)-kk);
            thres(kk+fall_end:end) = vfall;
        else
            thres(kk:end) = vmax - (vmax-vfall)/fall_end*(t(kk:end)-kk);
        end
        
        if ~ignore_pb
            pb.end_loop();
        end
        
    end
    
    
    %% Peak maximum search in PPG:
    if ~ignore_pb
        pb.reset();
        pb.checkpoint('Peak maximum search in PPG');
        pb.Loops2Do = length(peaks_filt);
    end
    
    peaks = nan(size(peaks_filt));
    w_nA_samples = round(fsppg*w_nA);
    for kk=1:length(peaks_filt)
        
        if ~ignore_pb
            pb.start_loop();
        end
        
%         w_int.begin = peaks_filt(kk) - w_nA_samples;
        w_int.begin = peaks_filt(kk);
        w_int.end = peaks_filt(kk) + w_nA_samples;

        if ~ignore_pb
            pb.checkpoint([]);
        end
        
        [aux, aux_t] = extract_interval(ppg, 1:length(ppg), w_int.begin, w_int.end);
        [~, pos] = max(aux);
        
        % Avoid possible detections at beggining or ending samples:
        if  ( aux_t(pos)==1 ) || ( ( aux_t(pos)==numel(ppg) ) )
            peaks(kk) = nan;
            continue;
        end
        
        if ( ppg(aux_t(pos))>=ppg(aux_t(pos)-1) ) && ( ppg(aux_t(pos))>=ppg(aux_t(pos)+1) )
            % Pulse is detected on a relative maximum -> right
            peaks(kk) = aux_t(pos);
        else
            % Pulse is not detected on a relative maximum -> wrong
            peaks(kk) = nan;
        end
        
        if ~ignore_pb
            pb.end_loop();
        end
        
    end
    peaks = peaks(~isnan(peaks));
    
    
    %% Figure:
    if plotflag==1
        figure;
        ax(1) = subplot(2,1,1); hold on;
        plot(ppg, 'b');
        plot(peaks_filt, ppg(peaks_filt), 'ro');
        plot(peaks, ppg(peaks), 'r*');
        title('PPG');
        ax(2) = subplot(2,1,2); hold on;
        plot(sig_filt, 'b');
        plot(thres, 'k');
        plot(peaks_filt, sig_filt(peaks_filt), 'ro');
        title('SSF');
        linkaxes(ax, 'x');
    end
    
end

function [ x_int, t_int, indexes ] = extract_interval( x, t, int_ini, int_end )
% -- Source --
% This function was copied from the 'ecg-kit' repository by 'marianux' at: https://github.com/marianux/ecg-kit . See also: http://marianux.github.io/ecg-kit/ 
% This particular function was available at: https://github.com/marianux/ecg-kit/blob/master/common/ppg/PPG_pulses_detector.m
% It is available under the GNU General Public License v2.0: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
% ------------

%EXTRACT_INTERVAL     Very simple function to extract an interval from a signal
%
% Created by Jesus Lazaro <jlazarop@unizar.es> in 2011
%--------
%   Sintax: [ x_int, t_int, indexes ] = extract_interval( x, t, int_ini, int_end )
%   In:   x = signal
%         t = time vector
%         int_ini = interval begin time (same units as 't')
%         int_end = interval end time (same units as 't')
%
%   Out:  x_int = interval [int_ini, int_end] of 'x'
%         t_int = interval [int_ini, int_end] of 't'
%         indexes = indexes corresponding to returned time interval

    if nargin<4
        error('Not enough input arguments');
    end
    
    indexes = find(t>=int_ini & t <=int_end);
    x_int = x(indexes);
    t_int = t(indexes);
end


%%%%%% qPPG code starting %%%%%%%%

function [idxPeaks] = qppg(data,fs,from,to)
% -- Source --
% This function is part of qppg.m, copied from the PhysioNet Cardiovascular Signal Toolbox 1.0.0 at: https://physionet.org/content/pcst/1.0.0/
% see also: https://doi.org/10.1088/1361-6579/aae021
% It is available under the GNU General Public License v3.0: https://www.gnu.org/licenses/gpl-3.0.en.html
% ------------

% This function is rewriten from wabp_pleth_new.c and wabp.c
% /* file wabp.c          Wei Zong       23 October 1998
%    			Last revised:   9 April 2010 (by G. Moody)
% -----------------------------------------------------------------------------
% wabp: beat detector for arterial blood presure (ABP) signal
% Copyright (C) 1998-2010 Wei Zong
% 
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
% PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 59 Temple
% Place - Suite 330, Boston, MA 02111-1307, USA.
% 
% You may contact the author by e-mail (wzong@mit.edu) or postal mail
% (MIT Room E25-505, Cambridge, MA 02139, USA).  For updates to this software,
% please visit PhysioNet (http://www.physionet.org/).
% ------------------------------------------------------------------------------
% 
% This program detects heart beats (pulse waveforms) in a continuous arterial
% blood pressure (ABP) signal.  This version of wabp works best with ABP signals
% sampled at 125 Hz, but it can analyze ABPs sampled at any frequency
% using on-the-fly resampling provided by the WFDB library.  'wabp' has been
% optimized for adult human  ABPs. For other ABPs, it may be necessary to
% experiment with the input sampling frequency and the time constants indicated
% below.
% 
% `wabp' can process records containing any number of signals, but it uses only
% one signal for ABP pulse detection (by default, the lowest-numbered signal
% labelled `ABP', `ART', or `BP';  this can be changed using the `-s' option, see
% below).
% 
% To compile this program under GNU/Linux, MacOS/X, MS-Windows, or Unix, use gcc:
%         gcc -o wabp wabp.c -lwfdb
% You must have installed the WFDB library, available at  
%         http://www.physionet.org/physiotools/wfdb.shtml
% gcc is standard with GNU/Linux and is available for other platforms from:
%         http://www.gnu.org/             (sources and Unix binaries)
%         http://fink.sourceforge.net     (Mac OS/X only)
%         http://www.cygwin.com/          (MS-Windows only)
%         
% For a usage summary, see the help text at the end of this file.  The input
% record may be in any of the formats readable by the WFDB library, and it may
% be anywhere in the WFDB path (in a local directory or on a remote web or ftp
% server).  The output of 'wabp' is an annotation file named RECORD.wabp (where
% RECORD is replaced by the name of the input record).  Within the output
% annotation file, the time of each NORMAL annotation marks an ABP pulse wave
% onset.
%
% 
%
% input :   data: PPG data
%           fs:   sampling frequency, default 125 Hz
%           from: begin point to analysis
%           to  : end point to analysis
% output:   beat: onset position of PPG beats in samples
%
%
%
% CHANGE LOG (should be moved to GitHub as commit messages...)
%
% 03 Aug 2010:  by Qiao Li
%   1) sample() function sometimes get WFDB_INVALID_SAMPLE value even when data is good.
%   2) Changed to mysample(), using isigsettime() and getvec() instead, but slow!!!
%   3) Add a mysetbuf() function, read the data to a databuf first, then mysample()
%      get data from the buf. It's much more fast than before.
%
%
% 04 Jul 2011 by Qiao Li
%   1) Set eye-closing period for PPG 0.34s, slope width 0.17s
%   2) Remove physical units to ADC units transform //Tm = physadu((unsigned)sig, Tm);
%   3) Add maxmin_2_3_threshold to modify the logic of eye-closing in order to minimize
%      the double beats
%	4) Change t += EyeClosing; to t = tpq+EyeClosing;
%
%
% 19 Jul 2011 by Qiao Li
% 
%   1) add: (before reading data,at line 480)    
%      isigsettime(from);
%      dataL=from;
%      to read data from 'from' setting
% 	2) Changed: (after learning period, at line 544)
%      (void)sample(sig, tpq);
%      if (sample_valid() == 0) break;
%       to:
%      if (dataend) break;
%
% 18 Oct 2016 by Qiao Li
%   1) add: input parameter fs for different sampling frequency data
%	2) add: re-scale data to ~ +/- 2000
%   3) add: find valley from the original data around 0.25s of tpq
%
% 03 Mar 2017 by Adriana Vest
%   Changed name of function to qppg to avoid confusion with wabp
%	Previous name: wabp_pleth_new.m


if nargin<3
    from=1;
    to=length(data);
end

if nargin<2
    fs=125;
end

global BUFLN ebuf lbuf tt_2  aet SLPwindow


idxPeaks=[];
beat_n=1;

sps=fs; % Sampling Frequency

BUFLN = 4096;           % /* must be a power of 2, see slpsamp() */
EYE_CLS = 0.34;         % /* eye-closing period is set to 0.34 sec (340 ms) for PPG */ 
LPERIOD  = sps*8;       % /* learning period is the first LPERIOD samples */
SLPW = 0.17;            % /* Slope width (170ms) for PPG */                        
NDP = 2.5;              % /* adjust threshold if no pulse found in NDP seconds */
TmDEF = 5;              % /* minimum threshold value (default) */
Tm = TmDEF;

BUFLN2 = (BUFLN*2);


INVALID_DATA=-32768;
if data(1)<=INVALID_DATA+10
    data(1)=mean(data);
end
inv=find(data<=INVALID_DATA+10);
for i=1:length(inv)
    data(inv(i))=data(inv(i)-1);
end

% re-scale data to ~ +/- 2000
if length(data)<5*60*sps
    data=(data-min(data))./(max(data)-min(data)).*4000-2000;
else
% find max/min every 5 minute for re-scaling data
    n=1;
    for i=1:5*60*sps:length(data)
        max_data(n)=max(data(i:min(i+5*60*sps-1,length(data))));
        min_data(n)=min(data(i:min(i+5*60*sps-1,length(data))));
        n=n+1;
    end
    data=(data-median(min_data))./(median(max_data)-median(min_data)).*4000-2000;
end

samplingInterval = 1000.0/sps;
spm = 60 * sps;
EyeClosing = round(sps * EYE_CLS);   % /* set eye-closing period */
ExpectPeriod = round(sps * NDP);	  % /* maximum expected RR interval */
SLPwindow = round(sps * SLPW);       % /* slope window size */
timer=0;

ebuf(1:BUFLN)=0;
lbuf=ebuf;
if from>BUFLN
    tt_2=from-BUFLN;
else
    tt_2=0;
end
aet=0;

t1=8*sps;
t1 = t1+from;
T0 = 0;
n=0;
for t = from:t1
    temp = slpsamp(t,data);
    if temp > INVALID_DATA+10
        T0 = T0+temp;
        n=n+1;
    end
end
T0 = T0/n; % T0=T0/(t1-from);
Ta = 3 * T0;

learning=1;
t=from;

%    /* Main loop */
t = from;
while t <= to
	
	if (learning) 

	    if (t > from + LPERIOD) 
    		learning = 0;
        	T1 = T0;
            t = from;	% /* start over */
        else
            T1 = 2*T0;
        end
    end
            
	temp = slpsamp(t,data);
    
	if (temp > T1)    % /* found a possible ABP pulse near t */ 
	    timer = 0; 
            % /* used for counting the time after previous ABP pulse */
	    maxd = temp;
        mind = maxd;
        tmax=t;
        for (tt = t + 1: t + EyeClosing-1)
            temp2=slpsamp(tt,data);
            if temp2 > maxd
                maxd=temp2;
                tmax=tt;
            end
        end
        if (maxd == temp)
            t=t+1;
            continue;
        end
        
        for (tt = tmax :-1: t - EyeClosing / 2 +1)
            temp2=slpsamp(tt,data);
            if temp2< mind
                mind=temp2;
            end
        end
        if maxd>mind+10
            onset=(maxd-mind)/100+2;
            tpq=t-round(0.04*fs);
            maxmin_2_3_threshold=(maxd-mind)*2.0/3;
            for tt=tmax:-1:t-EyeClosing/2+1
                temp2=slpsamp(tt,data);
                if temp2<maxmin_2_3_threshold
                    break;
                end
            end
            for tt=tt:-1:t - EyeClosing / 2 + round(0.024*fs)
                temp2=slpsamp(tt,data);
                temp3=slpsamp(tt-round(0.024*fs),data);
                if temp2-temp3<onset
                    tpq=tt-round(0.016*fs);
                    break;
                end
            end
            
            % find valley from the original data around 0.25s of tpq 
            valley_v = round(tpq);
            for valley_i=round(max(2,tpq-round(0.20*fs))):round(min(tpq+round(0.05*fs),length(data)-1))
                
                % If vally is too low, it cannot serve as an index, so move to the next time.
                if valley_v <= 0
                    t = t + 1;
                    continue;
                end
                
                if data(valley_v)>data(valley_i) && data(valley_i)<=data(valley_i-1) && data(valley_i)<=data(valley_i+1)
                    valley_v=valley_i;
                end
            end
            
            
            if (~learning) 
                
                % If we are looking for the first peak
                if beat_n == 1
                    
                    % If the proposed peak index > 0
                    if round(valley_v) > 0
                        idxPeaks(beat_n) = round(valley_v);
                        beat_n = beat_n + 1;
                    end
                else
                    % Check if rounded valley_v is greater than the prior beat index
                    if round(valley_v) > idxPeaks(beat_n-1)
                        idxPeaks(beat_n) = round(valley_v);
                        beat_n = beat_n + 1;
                    end
                end
            end
        

            % /* Adjust thresholds */
            Ta = Ta + (maxd - Ta)/10;
            T1 = Ta / 3;

            % /* Lock out further detections during the eye-closing period */
            t = tpq+EyeClosing;
        end
    else
        if (~learning) 
	    % /* Once past the learning period, decrease threshold if no pulse
	    %   was detected recently. */
            timer = timer+1;
	        if (timer > ExpectPeriod && Ta > Tm) 
                Ta=Ta-1;
                T1 = Ta / 3;
            end
        end
    end
    
    t=t+1;
    
end

% Discard first beat because algorithm always finds first minimum value, so trace-back logic
% will find a fir

end

function [beat1] = slpsamp(t,data) 
% -- Source --
% This function is part of qppg.m, copied from the PhysioNet Cardiovascular Signal Toolbox 1.0.0 at: https://physionet.org/content/pcst/1.0.0/
% see also: https://doi.org/10.1088/1361-6579/aae021
% It is available under the GNU General Public License v3.0: https://www.gnu.org/licenses/gpl-3.0.en.html
% ------------

global BUFLN ebuf lbuf tt_2  aet SLPwindow

    while (t > tt_2) 
        prevVal=0;
        
        if (tt_2>0) && (tt_2-1>0) && (tt_2<length(data)) && (tt_2-1<length(data))
            val2=data(tt_2 - 1);
            val1=data(tt_2);
        else
            val2=prevVal;
            val1=val2;
        end
    	prevVal=val2;
        dy =  val1-val2;
        if (dy < 0) 
            dy = 0;
        end
        tt_2=tt_2+1;
        M=round(mod(tt_2,(BUFLN-1))+1);
        et=dy;
        ebuf(M)=et;
%         M2=round(mod(tt_2-SLPwindow,(BUFLN-1))+1);
%         aet=aet+et-ebuf(M2);
        aet=0;
        for i=0:SLPwindow-1
            p=M-i;
            if p<=0
                p=p+BUFLN;
            end
            aet=aet+ebuf(p);
        end
        lbuf(M) = aet;

    end
    M3=round(mod(t,(BUFLN-1))+1);
    beat1=lbuf(M3);
end

%%%%%% qPPG code ended %%%%%%%%