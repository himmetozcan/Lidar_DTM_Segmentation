% VOTESEGDSM : probabilistic voting based DSM segmentation
%
%
% Syntax
% 
% 1) objs = votesegdsm(dsm)
% 2) objs = votesegdsm(dsm, ..., 'OptionName',OptionValue,...)
% 3) [objs, voteMap, Xe] = votesegdsm(...);
%
%
% 
% 1) objs = votesegdsm(dsm)
%        dsm: MxN digital surface model (DSM)
%        objs: MxN binary object class
%
% 2) objs = votesegdsm(dsm, ..., 'OptionName',OptionValue,...)
%        sets options OptionName to the specified OptionValue (see Options)
% 
% 3) [objs, voteMap, Xe] = votesegdsm(...);
%        returns object class, pdf and modes of pdf (local maximum locations)
% 
% 
% Description
%
% VOTESEGDSM is a probabilistic voting based DSM segmentation method. The
% output is the binary class of objects in DSM. Steps:
%  a) Extract edges in DSM
%  b) For each edge point (xe,ye), vote for possible object location
%  c) Gather votes and obtain a probability density function (pdf),
%     where modes (local maximums) of pdf give most probable object
%     locations (centers).
%  d) Extract modes (local maximums) of pdf 
%  e) Apply region growing (segmentation) for each mode and pixel wise
%     combine the results.
% 
%  Options 
% 
% The method is optimized for "0.5 meter" spatial resolution DSM.
% As described above, the method is applicable without entering any
% settings. All the parameters in the method is open to user if any setting
% want to be changed. The parameter definitions are given below: 
% 
% -(1)voting-
%
% cannysigma            : decrease for more edge points [0.1,5]
% removeedgesvariance   : decrease for more edge points [0.1,5]
% skipe                 : 1: all edge points, 2: 1:2:end edge points, and so on  
% votemaxwinsize        : maximum window size for voting on dsm
% voteheightdifference  : elevation height difference (m) for voting 
% kernelsigma           : bandwidth of Gaussian kernel (sigma of Gauss)
% 
% -(2)segmentation-
% 
% heightdifforinitialseg : elevation height difference (m) for initial segmentation
% segiteration           : max. iteration number at region growing
% thickstep              : size of thickenning at region growing
% segwinsize             : window size for segmentation (for faster run, try to choose larger than the maximum size of object) 
% muupper                : t_u, upper elev. threshold (m) for region growing
% mulower                : t_l, lower elev. threshold (m) for region growing
% mufinal                : d_mu, final elev. height dif. (m) for accepting the grown segment as object
% compareouters          : size of thickenning for comparing outer boundary
% 
% 
% Example(1):
% 
% Xkj = votesegdsm( Is ); % Is: MxN DSM, Xkj: binary object class
% 
% figure,imagesc(Xkj);axis image
% 
% 
% 
% Example(2):
% 
% [Xkj, kernelMap, Xe] = votesegdsm( Is, 'kernelsigma', 5, 'votemaxwinsize', 10 );
% 
% figure,imagesc(Xkj);axis image
% figure,imagesc(kernelMap);axis image
% figure,imagesc(Is);axis image
% hold on,plot(Xe(:,1),Xe(:,2),'k^','MarkerFaceColor','r') 
%
%
% Requirements:
% VOTESEGDSM requires 'imextrema.m' file for extracting local maximums. It
% is available via the Mathworks File Exchange Program at:
% http://www.mathworks.com/matlabcentral/fileexchange/41955-find-image-extrema
% 
% VOTESEGDSM uses parfor loop in parallel computing toolbox for faster run.
% If you dont have it, it should continue running as a regular for loop.
% 
%% 
% Author:
%     Abdullah H. Ozcan
%       PhD student at Department of Electrical and Electronics Engineering
%       Yeditepe University, Istanbul, Turkey.
%       himmetozcan@gmail.com
% 
% License:
% Copyright (c) 2016, Abdullah H. Ozcan
% All rights reserved.
% 
% 
%% 

function [objs, kernelMap, Xe] = votesegdsm(varargin)

[Is, cannySigma, removeEdgesVariance, skipE, voteMaxWinSize,...
    voteHeightDifference, kernelSigma, heightDifForInitialSeg, ...
    segIteration, thickStep, segWinSize, muUpper, muLower, muFinal,...
    compareOuters] = init(varargin{:});

%% a) Edge Detection: canny edge detector, then remove edges belong to non-objects

disp('..Edge Ext')

[edges, ~] = getclean_edges(Is, cannySigma, removeEdgesVariance);

%% b) Spatial Voting: pdf estimation and local maximums of pdf

disp('..Voting')

[edgeX, edgeY] = find(edges);

edgeX = edgeX(1 : skipE : end);
edgeY = edgeY(1 : skipE : end);

locMax = zeros(numel(edgeX), 2);
locMax(:, 1) = edgeY;
locMax(:, 2) = edgeX;

voteLocations = spatialvoting (Is, locMax, voteMaxWinSize, voteHeightDifference);

%% c) Gather votes and obtain a probability density function (pdf)

disp('..Kernel Est')

kernelMap = kernelestim (Is, voteLocations, kernelSigma);
kernelMap = ( kernelMap - min(kernelMap(:)) ) ./ ( max(kernelMap(:)) - min(kernelMap(:)) );

%% d) Extract modes (local maximums) of pdf

disp('..Loc Max')

[xm, ym, ~, cm] = imextrema(kernelMap);
xm = xm(cm == 1); ym = ym(cm == 1);
Xe(:,1) = xm(:);
Xe(:,2) = ym(:);

%% e) segmentation (morphology based region growing)

disp('..Segmentation')

objs = segmentation_process(Is,Is,xm,ym,segWinSize,heightDifForInitialSeg,...
    segIteration,thickStep, muUpper,muLower,muFinal, compareOuters);

end

function [Is, cannySigma, removeEdgesVariance, skipE, voteMaxWinSize,...
    voteHeightDifference, kernelSigma, heightDifForInitialSeg, segIteration,...
    thickStep, segWinSize, muUpper, muLower, muFinal, compareOuters] = init(varargin)

Is=varargin{1};
if nargin == 2
    if isstruct(varargin{2})
        inopts = varargin{2};
    else
        error('when using 2 arguments the first one is the analyzed DSM and the second one is a struct object describing the options')
    end
elseif nargin > 2
    try
        inopts = struct(varargin{2:end});
    catch
        error('bad argument syntax')
    end
end

opt_fields = {'cannysigma','removeedgesvariance','skipe','votemaxwinsize','voteheightdifference','kernelsigma','heightdifforinitialseg','segiteration','thickstep','segwinsize','muupper','mulower','mufinal','compareouters'};

%% Default parameter settings
% voting
defopts.cannysigma = 0.5; % decrease for more edge points [0.1,5]
defopts.removeedgesvariance = 0.5; % decrease for more edge points [0.1,5]
defopts.skipe = 1; % 1: all edge points, 2: 1:2:end edge points, and so on  
defopts.votemaxwinsize = 3; % maximum window size for voting on dsm
defopts.voteheightdifference = 1; % elevation height difference (m) for voting 
defopts.kernelsigma = 3; % bandwidth of Gaussian kernel (sigma of Gauss)

% segmentation
defopts.heightdifforinitialseg = 1; % elevation height difference (m) for initial segmentation
defopts.segiteration = 10; % max. iteration number at region growing
defopts.thickstep = 3; % size of thickenning at region growing
defopts.segwinsize = 150; % window size for segmentation (for faster run, try to choose larger than the maximum size of object) 
defopts.muupper = 1; % t_u: upper elev. threshold (m) for region growing
defopts.mulower = -5; % t_l: lower elev. threshold (m) for region growing
defopts.mufinal = 2; % deltaMu: final elev. height dif. (m) for accepting the grown segment as object
defopts.compareouters = 2; % size of thickenning for comparing outer boundary

opts = defopts;

if(nargin==1)
    inopts = defopts;
elseif nargin == 0
    error('not enough arguments')
end

names = fieldnames(inopts);
for nom = names'
    if ~any(strcmpi(char(nom), opt_fields))
        error(['bad option field name: ',char(nom)])
    end
    if ~isempty(eval(['inopts.',char(nom)])) % empty values are discarded
        eval(['opts.',lower(char(nom)),' = inopts.',char(nom),';'])
    end
end


cannySigma = opts.cannysigma;
removeEdgesVariance = opts.removeedgesvariance;
skipE = opts.skipe;
voteMaxWinSize = opts.votemaxwinsize;
voteHeightDifference = opts.voteheightdifference;
kernelSigma = opts.kernelsigma;

heightDifForInitialSeg = opts.heightdifforinitialseg;
segIteration = opts.segiteration;
thickStep = opts.thickstep;
segWinSize = opts.segwinsize;
muUpper = opts.muupper;
muLower = opts.mulower;
muFinal = opts.mufinal;
compareOuters = opts.compareouters;

if numel(size(Is)) > 2
    error('DSM must have only N rows and M columns')
end

end
