%BlandAltman Draws a Bland-Altman and correlation graph for two datasets.
%
%   BlandAltman(Data1, Data2) - Data1 and Data2 have to be of the same size
%   and can be grouped for display purposes. 3rd dimension is encoded by
%   colors and 2nd dimension by symbols. The 1st dimension contains 
%   measurements within the groups.
%
%   BlandAltman(Data1,Data2,Label) - Names of data sets. Formats can be
%   - {'Name1'}
%   - {'Name1, 'Name2'}
%   - {'Name1, 'Name2', 'Units'}
%
%   BlandAltman(Data1,Data2,Label,Tit,GNames) - Specifies the names of the
%   groups for the legend.
%
%   BlandAltman(Fig, ...) - specify a figure handle in which to display 
%   figure in which the Bland-Altman and correlation will be displayed.
%
%   BlandAltman(ah, ...) - specify an axes which will be replaced by the 
%   Bland-Altman and correlation axes.
%
%   rpc = BlandAltman(...) - return the coefficient of reproducibility
%   (1.96 times the standard deviation of the differnces).
%
%   [rpc fig] = BlandAltman(...) - also return the figure handles.
%
%   [rpc fig sstruct] = BlandAltman(...) - also return the structure of
%   statistics for the analysis.
%
%   BlandAltman(..., GNames, Parameter, Value) - call with parameter and 
%   value pairs using the following parameters:
%
%    'corrInfo' - specifies what information to display on the correlation 
%     plot as a cell of string in order of top to bottom. The following codes 
%     are available:
%       - 'eq' - slope and intercept equation
%       - 'r' - Pearson r-value
%       - 'r2' - Pearson r-value squared
%       - 'p' - Pearson correlation p-value
%       - 'rho' - Spearman rho value
%       - 'rho (p)' - Spearman rho value and p-value
%       - 'SSE' - sum of squared error
%       - 'RMSE' - root mean squared error
%       - 'n' - number of data points used
%     {default = {'eq';'r2';'SSE';'n'} }
%
%    'baInfo' - specifies what information to display on the Bland-Altman plot
%     similar to corrInfo, but with the following codes:
%      - 'SD' - standard-deviation of the differences
%      - 'RPC' - reproducibility coefficient (1.96*SD)
%      - 'LOA' - limits of agreement (1.96*SD) - same as RPC but different labelling
%      - 'RPC(%)' - reproducibility coefficient and % of values
%      - 'LOA(%)' - limits of agreement and % of values
%      - 'CV' - coefficient of variation (SD of mean values in %)
%      - 'IQR' - interquartile range.
%      - 'RPCnp' - RPC estimate based on IQR (non-parametric statistics) where
%              RPCnp = 1.45*IQR ~ RPC (if distribution of differences is
%              normal). 
%              See: Peck, Olsen and Devore, Introduction to Statistics and 
%              Data Analysis. Nelson Education, 2011.
%      - 'ks' - Kolmogorov-Smirnov test that difference-data is Gaussian
%      - 'kurtosis' - Kurtosis test that difference-data is Gaussian
%      - 'skewness' - skewness test results
%    {default = {'RPC(%)';'CV'} }
%
%   'axesLimits' - specifies the axes limits:
%     - scalar - lower limit (eg. 0)
%     - [min max] - specifies minimum and maximum
%     - 'tight' - minimum and maximum of data.
%     - 'auto' - plot default. {default}
%     - 'auto0 - same as auto with 0,0 minma.
%
%    'colors' - specify the order of group colors. (eg. 'brg' for blue, 
%    red, green) or RGB columns. {default = 'rbgmcky'}
% 
%    'symbols' - specify the order of symbols. (eg. 'sod.' for squares, 
%    circles, diamonds, dots). Alternatively can be set to 'Num' to display 
%    the subject number. {default = 'sodp^v'};
%
%    'markerSize' - set the size of the symbols on the plot (or font size
%    if using 'Num' mode for symbols. {default is 4}
%
%    'data1Mode' - how to treat data set 1:
%      - 'Compare' - data sets 1 and 2 are being compared. Means of data1 
%         and data2 are used for x-coordinates on Bland-Altman. {default}
%      - 'Truth' - data set 1 is considered a true reference by which data2 
%         is being evaluated. Data 1 values are used for x-coordinates on 
%         Bland-Altman.
%
%    'forceZeroIntercept' - force the y-intercept of the linear fit on the 
%    correlation analysis to zero. {default is 'off'}
%
%    'showFitCI' - show fit line confidence intervals on correlation plot.
%    {default is 'off'};
%
%    'diffValueMode' - Units for differences:
%      - 'Absolute' - same units as the data {default}
%      - 'relative' - differences are normalized to the reference data 
%         (mean or data 1 depending on dataOneMode option).
%      - 'percent'  - same as relative, but in percent units.
%
%    'baYLimMode' - Mode for setting y-lim on BA axes.
%      - 'Auto' - Automatically fit to the data.
%      - 'Square' - Preserve 1:1 aspect ratio with x-axis and 0 is 
%        centered. {default}
%
%    'baStatsMode' - Statistical analysis mode for Bland-Altman (differnces).
%      - 'Normal' - normal (Gaussian) distributed statistics
%      - 'Gaussian' - same as 'Normal'.
%      - 'Non-parametric' - non-parametric statistics.
%         * NOTE: Gaussian distribution is tested using the Kolmogorov-
%         Smirnov test. If the data seems to violate the assumption of 
%         distribution type, a warning message is generated.
%
%    'legend' - show legend.
%      - 'On' - show legend {default}
%      - 'Off' - don't display the legend.
%
% See also correlationPlot
% By Ran Klein, The Ottawa Hospital, Department of Nuclear Medicine
% 2013-02-20  RK  Added unit labeling in SSE and RPC labels.
% 2013-02-20  RK  Switched to equal scaling on BA y axis.
% 2014-01-15  RK  Added colors and symbols input arguments.
% 2016-08-10  RK  Major overhaul with addition of parameter-value pairs to
%                 support the following new features: data1Mode, 
%                 forceZeroIntercept, showFitCI, diffValueMode, baYlimMode,
%                 baStatsMode.
% 2017-05-05  RK  Added support for data consisting of nan values.
% 2017-09-05  RK  Seperated out corelationPlot function so that can be
%                 called independently of BA (unsupprted feature for now).
%                 Correction of baYLimMode not working bug and improved
%                 labelling in absence of units.
% 2018-12-12  RK  Formally added correlationPlot function to only plot the
%                 correlation plot without the Bland-Altman.
%                 Text is moved to place if axes limits change.
% 2019-08-23  RK  Replaced ranksum with singnrank as data is paired
%                 (feedback from  Smilla on Mathworks Community credited).
function [rpc, fig, stats] = BlandAltman(varargin)
[fig, data, params] = ParseInputArguments(varargin{:});
[cAH, baAH, fig, drawAreaPos] = ConfigAxes(fig, params.processCorrelationOnly);
% Correlation plot
stats = CalcCorrelationStats(data, params);
PlotCorrelation(cAH, data, params);
params = FormatPlotAxes(cAH, data, params);
DisplayCorrelationStats(cAH, params, stats);
addTitle(cAH, baAH, drawAreaPos, params);
if strcmpi(params.Legend,'On')
	addLegend(cAH, drawAreaPos, params);
end
if strcmpi(params.processCorrelationOnly,'On')
	rpc = [];
	return
end
% Bland-Altman plot of differences plot
[stats, data, params] = CalcBAStats(stats, data, params);
params = PlotBA(baAH, data, stats, params);
DisplayBAStats(baAH, params, stats)
rpc = stats.rpc;
%% Helper functions
function [fig, data, params] = ParseInputArguments(varargin)
% optional 1st parameter is figure handle
if isscalar(varargin{1}) && isequal(size(varargin{1}),[1 1]) && ishandle(varargin{1})
	shift = 1;
	fig = varargin{1};
else
	shift = 0;
	fig = [];
end
% followed by two data sets of equal size
data.set1 = varargin{shift+1};
data.set2 = varargin{shift+2};
s = size(data.set1);
if ~isequal(s,size(data.set2))
	error('data1 and data2 must have the same size');
end
if nargin>=shift+3
	label = varargin{shift+3};
else
	label = '';params.baInfo
end
if nargin>=shift+4
	params.tit = varargin{shift+4};
else
	params.tit = '';
end
if nargin>=shift+5
	params.gnames = varargin{shift+5};
else
	params.gnames = '';
end
% default values
params.corrInfo = {'eq';'r2';'SSE';'n'};
params.baInfo = {'RPC(%)';'CV'};
params.defaultBaInfo = true;
params.axesLimits = 'auto';
params.colors = 'brgmcky';
params.symbols = 'sodp^v';
params.markerSize = 4;
params.data1TreatmentMode = 'Compare';
params.forceZeroIntercept = 'off';
params.showFitCI = 'off';
params.baYLimMode = 'Squared';
params.baStatsMode = 'Normal';
params.diffValueMode = 'Absolute';
params.processCorrelationOnly = 'Off';
params.Legend = 'On';
axesLimitsSpecified = false;
% parse parameter value pair options
i = shift+6;
while length(varargin)>i
	parameter = varargin{i};
	val = varargin{i+1};
	switch upper(parameter)
		case 'CORRINFO'
			if ischar(val)
				params.corrInfo = {val};
			else
				params.corrInfo = val;
			end
			
		case 'BAINFO'
			if ischar(val)
				params.baInfo = {val};
			else
				params.baInfo = val;
			end
			params.defaultBaInfo = false;
		case 'AXESLIMITS'
			params.axesLimits = val;
			axesLimitsSpecified = true;
		case 'COLORS', params.colors = val;
		case 'SYMBOLS', params.symbols = val;
		case 'MARKERSIZE', params.markerSize = val;
		case 'DATA1MODE', params.data1TreatmentMode = val; % use the 'Compare' mean of data1 and data2 or 'Truth' data1 
		case 'FORCEZEROINTERCEPT', params.forceZeroIntercept = val;
		case 'SHOWFITCI', params.showFitCI = val;
		case 'BASTATSMODE', params.baStatsMode = val;
		case 'DIFFVALUEMODE', params.diffValueMode = val;
		case 'BAYLIMMODE', params.baYLimMode = val;
		case 'PROCESSCORRELATIONONLY', params.processCorrelationOnly = val;
		case 'LEGEND', params.Legend = val;
        case 'RANGE', params.range = val;
		otherwise
			error(['Unknown parameter name ' parameter])
	end % of swich statement
	i = i+2;
end
% Default axes mode for a correlation plot only does not assume equal x and
% y axis limits
if ~axesLimitsSpecified && strcmpi(params.processCorrelationOnly, 'On')
	params.axesLimits = 'tightNonsquare';
end
switch length(s)
	case 1
		s = [s 1 1];
	case 2
		s = [s 1];
	case 3
	otherwise
		error('Data have too many dimension');
end
% reformat data as an array of elements and store grouping number
params.numElementsPerGroup = s(1); % number of elements in each group
params.numGroups = numel(data.set1)/params.numElementsPerGroup;
params.numGroupsBySymbol = s(2);
params.numGroupsByColor = s(3);
if ~ischar(params.colors)
	if size(params.colors,2)~=3
		if size(params.colors,1)==3
			params.colors = params.colors';
		else
			error('Colors must be specified in either character codes or RGB');
		end
	end
elseif size(params.colors,1)==1
	params.colors = params.colors';
end
if size(params.colors,1)<params.numGroupsByColor
	error('More groups than colors specified. Use the colors input variable to specify colors for each group.');
end
if ~strcmpi(params.symbols,'Num') && length(params.symbols)<params.numGroupsBySymbol
	error('More subgroups than symbols specified. Use the symbols input variable to specify symbols for each subgroup, or use the ''Num'' option.');
end
data.set1 = reshape(data.set1, [numel(data.set1),1]);
data.set2 = reshape(data.set2, [numel(data.set2),1]);
data.mask = isfinite(data.set1) & isnumeric(data.set1) & isfinite(data.set2) & isnumeric(data.set2);
data.maskedSet1 = data.set1(data.mask);
data.maskedSet2 = data.set2(data.mask);
AE=abs(data.maskedSet1-data.maskedSet2);
MAE=mean(AE);

% data.maskedSet1=data.maskedSet1(find(prctile(AE,99)>AE));
% data.maskedSet2=data.maskedSet2(find(prctile(AE,99)>AE));
AE2=abs(data.maskedSet1-data.maskedSet2);
MAE2=mean(AE2);


params = ResolveLabels(params,label);
%% Resolve labels and units
function params = ResolveLabels(params,label)
units = '';
if ~iscell(label)
	label = {label};
end
if length(label)==1
	params.d1Label = [label{1} '_1'];
	params.d2Label = [label{1} '_2'];
	if strcmpi(params.data1TreatmentMode,'Compare')
		params.meanLabel = ['Mean ' label{1}];
	else
		params.meanLabel = params.d1Label;
	end
	params.deltaLabel = ['\Delta ' label{1}];
elseif length(label)==2
	params.d1Label = label{1};
	params.d2Label = label{2};
	if strcmpi(params.data1TreatmentMode,'Compare')
		params.meanLabel = ['Mean ' label{1} ' & ' label{2}];
	else
		params.meanLabel = label{1};
	end
	params.deltaLabel = [label{2} ' - ' label{1}];
else % units also provided
	units = label{3};
	params.d1Label = [label{1} ' (' units ')'];
	params.d2Label = [label{2} ' (' units ')'];
	if strcmpi(params.data1TreatmentMode,'Compare')
		params.meanLabel = ['Mean of ' label{1} ' & ' label{2} ' (' units ')'];
	else
		params.meanLabel = [label{1} ' (' units ')'];
	end
	params.deltaLabel = [label{2} ' - ' label{1}];
end
switch upper(params.diffValueMode)
	case 'ABSOLUTE'
		diffUnits = units;
	case 'RELATIVE'
		diffUnits = '';
		params.baYLimMode = 'Auto';
	case 'PERCENT'
		diffUnits = '%';
		params.baYLimMode = 'Auto';
	otherwise
		error(['Unsupported diffValueMode ' params.diffValueMode])
end
if ~isempty(diffUnits)
	params.deltaLabel = [params.deltaLabel ' (' diffUnits ')'];
end
if isempty(units)
	params.unitsStr = '';
	params.diffUnitsStr = '';
else
	params.unitsStr = [' ' units];
	params.diffUnitsStr = [' ' diffUnits];
end
%% Initialize the axes (correlation and Bland-Altman) for display 
function [cAH, baAH, fig, drawAreaPos] = ConfigAxes(fig, processCorrelationOnly)
FontSize=15;
if ~strcmpi(processCorrelationOnly, 'On')
	if isempty(fig)
		fig = figure;
		set(fig,'Units','centimeters','Position',[1 3 36 18],'color','w');
		cAH = subplot(1,2,1);
        set(gca,'fontsize', FontSize)
		baAH = subplot(1,2,2);
        set(gca,'fontsize', FontSize)
		drawAreaPos = [0 0 1 1];
	elseif strcmpi(get(fig,'type'),'figure')
		cAH = subplot(1,2,1);
		baAH = subplot(1,2,2);
		drawAreaPos = [0 0 1 1];
	elseif strcmpi(get(fig,'type'),'axes')
		ah = fig;
		drawAreaPos = get(ah,'Position');
		fig = get(ah,'parent');
		delete(ah);
		cAH = axes('parent',fig,'Position',[drawAreaPos(1) drawAreaPos(2) drawAreaPos(3)/2 drawAreaPos(4)]);
		baAH = axes('parent',fig,'Position',[drawAreaPos(1)+drawAreaPos(3)/2 drawAreaPos(2) drawAreaPos(3)/2 drawAreaPos(4)]);
	else
		error('What in tarnations is the handle that was passed to Bland-Altman????')
	end
	set(cAH,'tag','Correlation Plot');
	set(baAH,'tag','Bland Altman Plot');
else
	baAH = [];
	if isempty(fig)
		fig = figure;
		set(fig,'Units','centimeters','Position',[3 3 10 10],'color','w');
		cAH = subplot(1,1,1);
		drawAreaPos = [0 0 1 1];
	elseif strcmpi(get(fig,'type'),'figure')
		cAH = subplot(1,1,1);
		drawAreaPos = [0 0 1 1];
	elseif strcmpi(get(fig,'type'),'axes')
		cAH = fig;
		drawAreaPos = get(cAH,'Position');
	else
		error('What in tarnation is the handle that was passed to Bland-Altman????')
	end
end
% Make room for title and legend (5% on top and bottom)
height = drawAreaPos(4) * 0.87;
bottom = drawAreaPos(2) + height * 0.07;
cPos = get(cAH,'Position');
cPos(2) = bottom;
cPos(4) = height;
set(cAH,'tag','Correlation Plot','OuterPosition',cPos);
if ~strcmpi(processCorrelationOnly, 'On')
	baPos = get(baAH,'Position');
	baPos(2) = bottom;
	baPos(4) = height;
	set(baAH,'tag','Bland Altman Plot','OuterPosition',baPos);
end
set(fig,'Units','normalized')
%% Plot the correlation graph
function PlotCorrelation(cAH, data, params)
hold(cAH,'on');
% FontSize = 15;
% set(gca,'fontsize', FontSize)
for groupi=1:params.numGroups
	if strcmpi(params.symbols,'Num')
		for i=1:params.numElementsPerGroup
			text(data.set1((groupi-1)*params.numElementsPerGroup+i),...
				 data.set2((groupi-1)*params.numElementsPerGroup+i),num2str(i),...
				 'parent',cAH,...
				 'fontsize',params.markerSize,...
				 'color',params.colors(floor((groupi-1)/params.numGroupsBySymbol)+1,:),...
				 'HorizontalAlignment','Center',...
				 'VerticalAlignment','Middle');
		end
	else
		if params.numGroupsByColor==1
			marker = params.symbols(1);
			color = params.colors(groupi,:);
		else
% 			marker = params.symbols(rem(groupi-1,params.numGroupsBySymbol)+1);
            marker = params.symbols(groupi);
			color = params.colors(floor((groupi-1)/params.numGroupsBySymbol)+1,:);
		end
		ph=plot(cAH, data.set1((groupi-1)*params.numElementsPerGroup+(1:params.numElementsPerGroup)), data.set2((groupi-1)*params.numElementsPerGroup+(1:params.numElementsPerGroup)),...
			marker,...
			'color',color);
		set(ph,'markersize',params.markerSize,'LineWidth',1);
	end
end
xlabel(cAH,params.d1Label); ylabel(cAH,params.d2Label);
%% Calculate the statistical results for correlation analysis.
function stats = CalcCorrelationStats(data, params)
% Linear regression
if strcmpi(params.forceZeroIntercept,'on')
	[stats.polyCoefs, stats.polyFitStruct] = polyfitZero(data.maskedSet1, data.maskedSet2, 1);
else
	[stats.polyCoefs, stats.polyFitStruct] = polyfit(data.maskedSet1, data.maskedSet2, 1);
end
[r, p] = corrcoef(data.maskedSet1,data.maskedSet2); 
stats.r=r(1,2);
stats.r2 = stats.r^2;
[stats.rho, stats.rhoP] = corr(data.maskedSet1,data.maskedSet2,'type','Spearman');

stats.corrP = p(1,2);
% stats.N = sum(data.mask);
stats.N=length(data.maskedSet1);
stats.SSE = sum((polyval(stats.polyCoefs,data.maskedSet1)-data.maskedSet2).^2);
% MG modified on 18/12/2023
stats.SSE = sum((data.maskedSet1-data.maskedSet2).^2);
stats.RMSE = sqrt(stats.SSE/(stats.N-2));
stats.slope = stats.polyCoefs(1);
stats.intercept = stats.polyCoefs(2);
%% Format the corelation plot
function params = FormatPlotAxes(cAH, data, params)
squareAxes = true;
if ischar(params.axesLimits)
	if strcmpi(params.axesLimits,'Auto') || strcmpi(params.axesLimits,'Auto0')
		% Workaround - Add invisible minimum and maximum point to fix Auto axes limits (text
		% does not count for axis('auto')
		if strcmpi(params.symbols,'Num')
			mindata = min( min(data.maskedSet1), min(data.maskedSet2) );
			maxdata = max( max(data.maskedSet1), max(data.maskedSet2) );
			ph = plot(cAH, [mindata maxdata], [mindata maxdata], '.', 'Visible','on');
		end
		if strcmpi(params.axesLimits,'Auto0')
			params.axesLimits = axis(cAH);
			params.axesLimits(1) = 0;
		else
			params.axesLimits = axis(cAH);
			params.axesLimits(1) = min(params.axesLimits(1),params.axesLimits(3));
		end
		params.axesLimits(2) = max(params.axesLimits(2),params.axesLimits(4));
		if strcmpi(params.symbols,'Num')
			delete(ph);
		end
	elseif strcmpi(params.axesLimits,'Tight')
		params.axesLimits = zeros(1,4);
		params.axesLimits(1) = min( min(data.maskedSet1), min(data.maskedSet2) );
		params.axesLimits(2) = max( max(data.maskedSet1), max(data.maskedSet2) );
	elseif strcmpi(params.axesLimits,'TightNonsquare')
		params.axesLimits = [min(data.maskedSet1),  max(data.maskedSet1), min(data.maskedSet2) , max(data.maskedSet2)];
		squareAxes = false;
	else
		error(['Unknown axes limit option (' params.axesLimits ') detected.']);
	end
else
	if length(params.axesLimits)==1
		a = axis(cAH);
		params.axesLimits(2) = max(a(2),a(4));
	elseif length(params.axesLimits)==2
		% Do nothing
	elseif length(params.axesLimits)==4
		squareAxes = false;
	else
		error('Axes limits must either be a string option or array of length 1,2 or 4 limits.');
	end
end
if squareAxes
	params.axesLimits(3) = params.axesLimits(1);
	params.axesLimits(4) = params.axesLimits(2);
	
	axis(cAH,params.axesLimits); axis(cAH,'square');
else
	axis(cAH,params.axesLimits); 
end

%% Display the correlation analyssis results
function DisplayCorrelationStats(cAH, params, stats)
x = linspace(params.axesLimits(1), params.axesLimits(2), 100);
[y, delta] = polyconf(stats.polyCoefs, x, stats.polyFitStruct,'simopt','on');
plot(cAH, x, y,	'-k');
if strcmpi(params.showFitCI,'on')
	plot(cAH, x, y+delta, '-', 'Color', 0.3*[1 1 1]);
	plot(cAH, x, y-delta, '-', 'Color', 0.3*[1 1 1]);
end
l = [min([params.axesLimits(1) params.axesLimits(3)]), max([params.axesLimits(2) params.axesLimits(4)])];
h = plot(cAH, l, l, ':'); set(h,'color',[0.6 0.6 0.6],'tag','IdentityLine');
if 0 % Add 95% CI lines
	xfit = params.axesLimits(1):(params.axesLimits(2)-params.axesLimits(1))/100:params.axesLimits(2);
	[yfit, delta] = polyconf(polyCoefs,xfit,S);
	h = [plot(cAH,xfit,yfit+delta);...
		plot(cAH,xfit,yfit-delta)];
	set(h,'color',[0.6 0.6 0.6],'linestyle','-');
end
corrtext = {};
for i=1:length(params.corrInfo)
	switch upper(params.corrInfo{i})
		case 'EQ'
			if ~strcmpi(params.forceZeroIntercept,'off')
				corrtext = [corrtext; ['y=' mynum2str(stats.slope,3,2) 'x']];
			elseif stats.intercept>=0
				corrtext = [corrtext; ['y=' mynum2str(stats.slope,3,2) 'x+' mynum2str(stats.intercept,3)]];
			else
				corrtext = [corrtext; ['y=' mynum2str(stats.slope,3,2) 'x' mynum2str(stats.intercept,3)]];
			end
		case 'R2', corrtext = [corrtext; ['r^2=' mynum2str(stats.r^2,4)]];
		case 'R', corrtext = [corrtext; ['r=' mynum2str(stats.r,4)]];
		case 'P', corrtext = [corrtext; ['p=' num2str(stats.corrP,2)]];
		case 'RHO', corrtext = [corrtext; ['rho=' mynum2str(stats.rho,4,4)]];
		case 'RHO (P)', corrtext = [corrtext; ['rho=' mynum2str(stats.rho,4,4) ' (p=' mynum2str(stats.rhoP,2) ')']];
		case 'SSE', corrtext = [corrtext; ['SSE=' mynum2str(stats.SSE,2) params.unitsStr]];
		case 'RMSE', corrtext = [corrtext; ['RMSE=' mynum2str(stats.RMSE,2) params.unitsStr]];
		case 'N', corrtext = [corrtext; ['n=' mynum2str(stats.N,4,0)]];
	end
end

FontSize = 15;
text(params.axesLimits(1),params.axesLimits(4),corrtext,'HorizontalAlignment','left', ...
     'VerticalAlignment','top','parent',cAH,'tag','CorrInfoText','FontSize',FontSize);
addlistener(cAH.XRuler,'MarkedClean',@(src,event)updateCorrAxesCallback(src,event) );

%% Callback function to keep results in place if the Correlation plot axis limits change
function updateCorrAxesCallback(ah, ~)
ah = get(ah, 'Parent');
a = axis(ah);
ch = get(ah,'children');
set(findobj(ch, 'flat', 'tag','CorrInfoText'),'Position', [a(1) a(4)]);
%% Calculate statistics for BA analysis
function [stats, data, params] = CalcBAStats(stats, data, params)
if strcmpi(params.data1TreatmentMode,'Compare')
	data.maskedBaRefData = mean([data.maskedSet1,data.maskedSet2],2);
	data.baRefData = mean([data.set1,data.set2],2);
else
	data.maskedBaRefData = data.maskedSet1; % previous version was calaculated as RPC/mean of data.
	data.baRefData = data.set1;
end
switch upper(params.diffValueMode)
	case 'ABSOLUTE'
		data.maskedDifferences = data.maskedSet2-data.maskedSet1;
		data.differences = data.set2-data.set1;
	case 'RELATIVE'
		data.maskedDifferences = (data.maskedSet2-data.maskedSet1) ./ data.maskedBaRefData;
		data.differences = (data.set2-data.set1) ./ data.baRefData;
	case 'PERCENT'
		data.maskedDifferences = (data.maskedSet2-data.maskedSet1) ./ data.maskedBaRefData*100;
		data.differences = (data.set2-data.set1) ./ data.baRefData*100;
end
stats.differenceSTD = std(data.maskedDifferences);
stats.differenceMean = mean(data.maskedDifferences);
stats.differenceMedian = median(data.maskedDifferences);
[~, stats.differenceMeanP] = ttest(data.maskedDifferences,0);
stats.differenceMedianP = signrank(data.set1,data.set2);
stats.rpc = 1.96*stats.differenceSTD;
stats.CV = 100*stats.differenceSTD/mean((data.maskedSet1+data.maskedSet2)/2);
stats.rpcPercent = 1.96*std(data.maskedDifferences ./ data.maskedBaRefData)*100; % previous version was calaculated as RPC/mean of data.
stats.IQR = iqr(data.maskedDifferences);
stats.rpcNP = stats.IQR * 1.45; % estimate of RPC if distribution was Gaussian: see: R. Peck, C. Olsen, and J. Devore, Introduction to Statistics and Data Analysis. Nelson Education, 2011.
[~, stats.ksp] = kstest((data.maskedDifferences-stats.differenceMean)/stats.differenceSTD);
stats.kurtosis = kurtosis(data.maskedDifferences);
stats.skewness = skewness(data.maskedDifferences,1);
%% Plot the BA plot
function params = PlotBA(baAH, data, stats, params)
set(baAH,'Units','normalized');
hold(baAH,'on');
for groupi=1:params.numGroups
	ref = data.baRefData((groupi-1)*params.numElementsPerGroup+(1:params.numElementsPerGroup));
	dif = data.differences((groupi-1)*params.numElementsPerGroup+(1:params.numElementsPerGroup));
	if strcmpi(params.symbols,'Num')
		for i=1:params.numElementsPerGroup
			text(ref(i), dif(i), num2str(i), 'parent',baAH,'fontsize',params.markerSize,'color',params.colors(floor((groupi-1)/params.numGroupsBySymbol)+1,:));
		end
	else
		if params.numGroupsByColor==1
			marker = params.symbols(1);
			color = params.colors(groupi,:);
		else
% 			marker = params.symbols(rem(groupi-1,params.numGroupsBySymbol)+1);
            marker = params.symbols(groupi);
			color = params.colors(floor((groupi-1)/params.numGroupsBySymbol)+1,:);
		end
		ph = plot(baAH,ref,dif,marker,'color',color);
		set(ph,'markersize',params.markerSize,'LineWidth',1);
	end
end
axis(baAH,'square')
xlabel(baAH,params.meanLabel); ylabel(baAH,params.deltaLabel);
% fix limits to +/- data limit
if strcmpi(params.baYLimMode,'Squared') 
	a = [params.axesLimits(1:2) [-1 1]*abs(params.axesLimits(2)-params.axesLimits(1))/2];
    %%LOA
    a = [params.axesLimits(1:2) [-1 1]*params.range];
	axis(baAH, a);
else
	a = axis(baAH);
end
fontsize = 15;
switch upper(params.baStatsMode)
	case {'NORMAL','GAUSSIAN'}
		plot(baAH,a(1:2),stats.differenceMean+[0 0],'k',"LineWidth",1.5)
		plot(baAH,a(1:2),stats.differenceMean+stats.rpc*[1 1],':k',"LineWidth",1)
		plot(baAH,a(1:2),stats.differenceMean-stats.rpc*[1 1],':k',"LineWidth",1)
		text(a(2),stats.differenceMean+stats.rpc, [mynum2str(stats.differenceMean+stats.rpc,2) ' (+1.96SD)'],'Parent',baAH,'HorizontalAlignment','left','VerticalAlignment','bottom','fontsize',fontsize,'Tag','ULimLabel');
		text(a(2),stats.differenceMean,[mynum2str(stats.differenceMean,2) ' [p=' mynum2str(stats.differenceMeanP,2) ']'],'Parent',baAH,'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fontsize,'Tag','MeanLabel');
		text(a(2),stats.differenceMean-stats.rpc, [mynum2str(stats.differenceMean-stats.rpc,2) ' (-1.96SD)'],'Parent',baAH,'HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize,'Tag','LLimLabel');
		if ~isGaussian(stats) 
			warning('Bland-Altman analysis is being performed using a Normal distribution assumptions, but the data does not appear to be normally distributed. Consider using a non-parametric analysis instead. See ''baStatsMode'' option for more details.')
		end
	case 'NON-PARAMETRIC'
		plot(baAH,a(1:2),stats.differenceMedian+[0 0],'k')
		plot(baAH,a(1:2),stats.differenceMedian+stats.rpcNP*[1 1],':k')
		plot(baAH,a(1:2),stats.differenceMedian-stats.rpcNP*[1 1],':k')
		text(a(2),stats.differenceMedian+stats.rpcNP, [mynum2str(stats.differenceMedian+stats.rpcNP,2) ' (+1.45IQR)'],'Parent',baAH,'HorizontalAlignment','left','VerticalAlignment','bottom','fontsize',fontsize,'Tag','ULimLabel');
		text(a(2),stats.differenceMedian,[mynum2str(stats.differenceMedian,2) ' [p=' mynum2str(stats.differenceMedianP,2) ']'],'Parent',baAH,'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fontsize,'Tag','MeanLabel');
		text(a(2),stats.differenceMedian-stats.rpcNP, [mynum2str(stats.differenceMedian-stats.rpcNP,2) ' (-1.45IQR)'],'Parent',baAH,'HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize,'Tag','LLimLabel');
		if isGaussian(stats) 
			warning('Bland-Altman analysis is being performed using a non-parametric distribution assumptions, but the data appears to be normally distributed. Consider using a Gaussian analysis instead. See ''baStatsMode'' option for more details.')
		end
		% Change BA summary data overlay to RPCnp if specific overlay was
		% not specified in he input arguments.
		if params.defaultBaInfo
			params.baInfo = {'RPCnp'};
		end
	otherwise
		error(['Unrecognized naStatsMode ' params.baStatsMode])
end
%% Display the stats of interest on the BA plot
function DisplayBAStats(baAH, params, stats)
BAtext = {};
for i=1:length(params.baInfo)
	switch upper(params.baInfo{i})
        case 'R2', BAtext = [BAtext; ['r^2=' mynum2str(stats.r^2,4)]];
		case 'R', BAtext = [BAtext; ['r=' mynum2str(stats.r,4)]];
		case 'P', BAtext = [BAtext; ['p=' num2str(stats.corrP,2)]];
		case 'RHO', BAtext = [BAtext; ['rho=' mynum2str(stats.rho,4,4)]];
		case 'RHO (P)', BAtext = [BAtext; ['rho=' mynum2str(stats.rho,4,4) ' (p=' mynum2str(stats.rhoP,2) ')']];
		case 'SSE', BAtext = [BAtext; ['SSE=' mynum2str(stats.SSE,2) params.unitsStr]];
		case 'RMSE', BAtext = [BAtext; ['RMSE=' mynum2str(stats.RMSE,2) params.unitsStr]];
		case 'N', BAtext = [BAtext; ['n=' mynum2str(stats.N,4,0)]];

		case 'SD', BAtext = [BAtext; ['{\bfSD: ' mynum2str(stats.differenceSTD,2) params.diffUnitsStr '}']];
		case 'RPC', BAtext = [BAtext; ['{\bfRPC: ' mynum2str(stats.rpc,2) params.diffUnitsStr '}']];
		case 'RPC(%)', BAtext = [BAtext; ['{\bfRPC: ' mynum2str(stats.rpc,2) params.diffUnitsStr '} (' mynum2str(stats.rpcPercent,2) '%)']];
		case 'LOA', BAtext = [BAtext; ['{\bfLOA: ' mynum2str(stats.rpc,2) params.diffUnitsStr '}']];
		case 'LOA(%)', BAtext = [BAtext; ['{\bfLOA: ' mynum2str(stats.rpc,2) params.diffUnitsStr '} (' mynum2str(stats.rpcPercent,2) '%)']];
		case 'CV', BAtext = [BAtext; ['CV: ' mynum2str(stats.CV,2) '%']];
		case 'RPCNP', BAtext = [BAtext; ['{\bfRPC_{np}: ' mynum2str(stats.rpcNP,2) params.diffUnitsStr '}']];
		case 'KS' % Kolmogorov-Smirnov test that difference-data is Gaussian
			BAtext = [BAtext; ['KS p-value: ' mynum2str(stats.ksp,3,3)]];
		case 'KURTOSIS' % Kurtosis test that difference-data is Gaussian
			BAtext = [BAtext; ['kurtosis: ' mynum2str(stats.kurtosis,2,2)]];
		case 'SKEWNESS'
			BAtext = [BAtext; ['skewness: ' mynum2str(stats.skewness,2,2)]];
	end
end
a = axis(baAH);
FontSize = 18;
text(a(1),a(4),BAtext,'interpreter','tex','HorizontalAlignment','left','VerticalAlignment','top','Parent',baAH,'tag','BAInfoText','FontSize',FontSize);
%STOP HERE
% addlistener(baAH.XRuler, 'MarkedClean', @(src,event)updateBAAxesCallback(src,event) );
function updateBAAxesCallback(ah, ~)
ah = get(ah,'Parent');
a = axis(ah);
ch = get(ah, 'Children');
set(findobj(ch, 'flat', 'tag','BAInfoText'),'Position',[a(2), a(4)]);
p = get(findobj(ch, 'flat', 'tag','ULimLabel'),'Position');
if ~isempty(p), set(findobj(ch, 'flat', 'tag','ULimLabel'),'Position', [a(2) p(2)]); end
p = get(findobj(ch, 'flat', 'tag','MeanLabel'),'Position');
if ~isempty(p), set(findobj(ch, 'flat', 'tag','MeanLabel'),'Position', [a(2) p(2)]); end
p = get(findobj(ch, 'flat', 'tag','LLimLabel'),'Position');
if ~isempty(p), set(findobj(ch, 'flat', 'tag','LLimLabel'),'Position', [a(2) p(2)]); end
%% Add legend to the plot
function addTitle(cAH, baAH, drawAreaPos, params)
if ~isempty(params.tit)
	fig = get(cAH,'Parent');
	units = get(fig,'Units');
	set(fig,'Units','normalized');
	cPos = get(cAH,'Position');
% 	if isempty(baAH)
% 		baPos = cPos;
% 	else
% 		baPos = get(baAH,'Position');
% 	end
	pos = [drawAreaPos(1), drawAreaPos(2)+drawAreaPos(4)*0.95, drawAreaPos(3), drawAreaPos(4)*0.05];
	titleAH = axes('Position',pos,'Xlim',[-1 1], 'YLim',[-1 1],'Visible','off');
    FontSize=20;
	text(0,0,params.tit,'HorizontalAlignment','Center','VerticalAlignment','Middle','FontSize',FontSize,'FontWeight','bold','Interpreter','Tex','Parent',titleAH)
	set(fig,'Units',units);
end
%% Add legend to the plot
function addLegend(cAH, drawAreaPos, params)
gnames = params.gnames;
if ~strcmpi(params.symbols,'Num') && ~isempty(gnames)
	lh = legend(cAH,'show');
	if iscell(gnames)
		if length(gnames)==2 
			if iscell(gnames{1}) 
				temp = cell(1,params.numGroups);
				for groupi=1:length(gnames{1})
					for j=1:length(gnames{2})
						temp{groupi+(j-1)*length(gnames{1})} = [gnames{1}{groupi} '-' gnames{2}{j}];
					end
				end	
				gnames = temp;
			elseif iscell(gnames{2})
				gnames = strcat(gnames{1}, '-', gnames{2});
            end
        else
            gnames = gnames{1};
		end
	end
	fig = get(cAH,'Parent');
	units = get(fig,'Units');
	set(fig,'Units','Normalized');
	set(lh,'string',gnames,'orientation','horizontal');
	drawnow;
	set(lh,'Units','Normalized');
	pos = get(lh,'Position');
	width = min(pos(3),drawAreaPos(3));
	set(lh,'Position',[drawAreaPos(1)+(drawAreaPos(3)-width)/2 drawAreaPos(2) width drawAreaPos(4)*0.05]);
	set(fig,'Units',units);
end
%% Is the BA stats data Gaussian?
function answer = isGaussian(stats)
answer = stats.ksp > 0.05;