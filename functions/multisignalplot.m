function [fh,lh] = multisignalplot(d,varargin)
% -------- [fh,lh] = multisignalplot(d,varargin) ---------
% 
% Plots all of the signals in a data matrix in one figure, separated on the y-axis by a
% certain value. 
%
%               >>> INPUTS >>>
% Required:
%   d = column-oriented data matrix
% Optional:
%   Fs = sampling rate...if provided, will plot relative to time on x-axis
%   col = line color (default []...default matlab colors);
%   maxval = value to separate lines by...default is mean(range(d)); 
%
%               <<< OUTPUTS <<<
%   fh = handle to gca
%   lh = handle to each line separately
% scaler=gain to apply to signals
%
% By JMS, 11/04/2015
% ----------------------------------------------------------------------------
hold on;

if nargin>1; Fs = varargin{1}; 
else Fs = []; end
if nargin>2 && ~isempty(varargin{2}); col = varargin{2};
else col = []; end
if nargin>3 && ~isempty(varargin{3}); maxval = varargin{3};
else maxval = nanmean(range(d)); end
if maxval == 0
    maxval = 0.1;
end

% get time vec if supplied sampling rate
if ~isempty(Fs)
    time = (1:size(d,1))/Fs;
else
    time = (1:size(d,1));
end

n = size(d,2); 
if size(col,1) == 1
    col = repmat(col,n,1);
end
for i = 1:n
    yval(i) = maxval*(i-1);
    if ~isempty(col)
        lh(i) = plot(time,d(:,i)-yval(i),'color',col(i,:));
    else
        lh(i) = plot(time,d(:,i)-yval(i));
    end
end
fh = gca;
set(fh,'box','off','tickdir','out','ytick',-fliplr(yval),'yticklabel',fliplr([1:n]));
axis tight
ylim([-maxval*(n) maxval]);
    
end