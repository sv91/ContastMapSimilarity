function ContrastMapComparator(varargin)
%% ContrastMapComparator
% Compare Contrast maps to show in which condition is the strongest at each
% voxel and save the result into "compare.mat" file.
%
% Arguments:
%   'folder'   : Folder containing the contrast maps. 
%                       Default : pwd.
%   'regex'    : Regex of the files to use. 
%                       Default : '*mat'.
%   'thr_coef' : Coefficient to use to multiply the threshold. 
%                       Default : 1.    
%   'thr_p'    : Threshold p value, in case of -log(p) threshold.
%                       Default : 0.01.    
%   'save'     : Path where to save.
%                       Default : pwd.
%   'type'     : Type of comparison to do.
%                       -1 : How is the most negative.
%                        0 : How has the highest absolute value. (Default)
%                        1 : How is the most positive.
%
% Example:
%   ContrastMapSimilarities('folder','example/folder','thr_coef', 10, 'min_sym', 5);

%% Managing the inputs
if nargin > 0
   for i=1:2:nargin
       switch varargin{i}
           case 'folder'
               folder        = fullfile(pwd,varargin{i+1});
           case 'regex'
               rgx           = varargin{i+1};
           case 'thr_coef'
               thr_coef_mult = varargin{i+1};
           case 'thr_p'
               thr_p         = varargin{i+1};
           case 'type'
               type         = varargin{i+1};
           case 'save'
               save_file     = strcat(varargin{i+1},'/compare.mat');
           otherwise
               fprintf('Unknown argument "%s"\n',varargin{i});
       end
   end
end

%% Setting to default all the unset variables.
if notDefined('folder')   
    folder       = pwd;
end

if notDefined('rgx')   
    rgx          = '*.mat';
end

if notDefined('thr_coef_mult')   
    thr_coef_mult = 1;
end

if notDefined('thr_p')   
    thr_p       = 0.01;
end

if notDefined('type')   
    type       = 0;
end

if notDefined('save')   
    save_file = strcat(pwd,'/compare.mat');
    counter = 1;
    while exist(save_file, 'file') == 2
        save_file = strcat(pwd,'/compare',num2str(counter) ,'.mat');
        counter = counter + 1;
    end
end


%% Calculations
% Getting all the files coresponding to the regex.
tmp = dir(fullfile(folder,rgx));
nb_files = numel(tmp);
for ii = 1:nb_files
    mat_file{ii} = fullfile(folder, tmp(ii).name); 
    assert(exist(mat_file{ii}, 'file')>0)
end

% Calculating the threshold.
threshold = -log(thr_p)*thr_coef_mult;

% Preloading to get size
L = load(mat_file{1});
t = L.map{1};
[a,b,c] = size(t);
temp     = zeros(a,b,c);

fprintf('=====\nComputing similarities\n=====\n');

% Keeping only the needed values.
for ii = 1:nb_files
   L = load(mat_file{ii});
   t = L.map{1};
   [a,b,c] = size(t);
   for i=1:a
       for j=1:b
           for k=1:c
               if (abs(t(i,j,k))>threshold) 
                   change = false;
                   switch type
                       case -1
                           change = (t(i,j,k)<temp(i,j,k));
                           mapName = 'Contrast Maps with most Negative Value';
                       case 0
                           change = (abs(t(i,j,k))>temp(i,j,k));
                           mapName = 'Contrast Maps with Highest Absolute Value';
                       case 1
                           change = (t(i,j,k)>temp(i,j,k));
                           mapName = 'Contrast Maps with most Positive Value';
                   end
                   if change
                        temp(i,j,k) = ii;
                   end
               end
           end
       end
   end
end

%% Saving the resulting matrix
% Renaming
map = cell(1);
map{1} = temp;
co = L.co;
mapUnits = 'Contrast Map ID';

save(save_file,'map','mapName','mapUnits','co');
