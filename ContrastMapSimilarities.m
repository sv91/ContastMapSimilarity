function ContrastMapSimilarities(varargin)
%% ContrastMapSimilarities
% Compute and show the similarities between the contrast maps provided.
%
% Arguments:
%   'folder'    : Folder containing the contrast maps. 
%                       Default : pwd.
%   'regex'     : Regex of the files to use. 
%                       Default : '*mat'.
%   'thr_coef'  : Coefficient to use to multiply the threshold. 
%                       Default : 1.    
%   'thr_type'  : Type of threshold.
%                       0 : mean * coef. (Default)
%                       1 : -log(p) * coef.
%   'thr_p'     : Threshold p value, in case of -log(p) threshold.
%                       Default : 0.01.
%   'min_sym'   : Minimal number of appearance a similarity have to have to
%                 appear in the graph.
%                       Default : 1.
%   'clr_type'  : Coloring type to apply.
%                       0 : Color depends of the number of times a
%                       similarity appears. (Default)
%                       1 : Different colors are assigned to each contrast
%                       map. 
%   'clr_scheme': Color scheme to apply.
%                       Default : [1 0 0;0 1 0;0 0 1;1 1 0;0 1 1;1 0 1].
%   'show_brain': Show the middle inplane slice.
%                       0 : Not showing.
%                       1 : Show. (Default)
%   'inp_file'  : Path to the inplane picture.
%                       Default : "inplane.png".
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
           case 'thr_type'
               thr_type      = varargin{i+1};
           case 'thr_p'
               thr_p         = varargin{i+1};
           case 'min_sym'
               limit         = varargin{i+1};
           case 'clr_type'
               clr_type      = varargin{i+1};
           case 'clr_scheme'
               clr_scheme    = varargin{i+1};
           case 'show_brain'
               show_brain    = varargin{i+1};
           case 'inp_file'
               inp_file      = fullfile(pwd,varargin{i+1});
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

if notDefined('thr_type')   
    thr_type    = 0;
end

if notDefined('thr_p')   
    thr_p       = 0.01;
end

if notDefined('limit')   
    limit       = 1;
end

if notDefined('clr_type')   
    clr_type    = 0;
end

if notDefined('clr_scheme')   
    clr_scheme  = [1 0 0;0 0 1;0 1 0;0 1 1;1 1 0;1 0 1];
end

if notDefined('show_brain')   
    show_brain  = 1;
end

if notDefined('inp_file')   
    inp_file    = 'inplane.png';
end

%% Calculations
% Getting all the files coresponding to the regex.
tmp = dir(fullfile(folder,rgx));
nb_files= numel(tmp);
for ii = 1:nb_files
    mat_file{ii} = fullfile(folder, tmp(ii).name); 
    assert(exist(mat_file{ii}, 'file')>0)
end

% Preallocting all the needed matrices.
temp     = zeros(128,128,54);
clr_temp = zeros(128,128,54,3);
maxV     = zeros(0,nb_files);
minV     = zeros(0,nb_files);
avgV     = zeros(0,nb_files);
absV     = zeros(0,nb_files);

% Checking the values of each file.
for ii = 1:nb_files
   L          = load(mat_file{ii});
   maxV(1,ii) = max(max(max(L.map{1})));
   minV(1,ii) = min(min(min(L.map{1})));
   avgV(1,ii) = mean(mean(mean(L.map{1})));
   absV(1,ii) = mean(mean(mean(abs(L.map{1}))));
end

fprintf('Maximal values:\n');
disp(maxV);
fprintf('Minimal values:\n');
disp(minV);
fprintf('Mean values:\n');
disp(avgV);
fprintf('Absolute mean values:\n');
disp(absV);

% Calculating the threshold depending of the type.
switch thr_type
    case 0
        threshold = mean(absV)*thr_coef_mult;
        thr_val   = thr_coef_mult;
    case 1
        threshold = -log(thr_p)*thr_coef_mult;
        thr_val   = thr_p;
end

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
                   temp(i,j,k) = temp(i,j,k) + 1;
                   for color_i=1:3
                        clr_temp(i,j,k,color_i) = (clr_temp(i,j,k,color_i)/temp(i,j,k))*(temp(i,j,k)-1)+clr_scheme(ii,color_i)/temp(i,j,k);
                   end
               end
           end
       end
   end
end


figure
hold on;
%% Color the brain.
fprintf('=====\nGeneration started\n=====\n');
[a,b,c] = size(temp);
for i=1:a
    if mod(i,10)==0
        fprintf('Generated %d slices out of %d.\n',i,a);
    end
    for j=1:b
        for k=1:c
            % Coloring type 0
            if temp(i,j,k)>=limit && clr_type==0
                col='';
                alp = 0.7;
                switch temp(i,j,k)
                    case nb_files-5
                        col = clr_scheme(6,:);
                    case nb_files-4
                        col = clr_scheme(5,:);
                    case nb_files-3
                        col = clr_scheme(4,:);
                    case nb_files-2
                        col = clr_scheme(3,:);
                    case nb_files-1
                        col = clr_scheme(2,:);
                    case nb_files
                        col = clr_scheme(1,:);
                end
                voxel([i j k],[1 1 1],col,alp);
            end
            % Coloring type 1
            if temp(i,j,k)>=1 && clr_type==1
                voxel([i j k],[1 1 1],clr_temp(i,j,k,:),0.8);
            end
        end
    end
end
fprintf('=====\nGeneration finished\n=====\n');
rotate3d on
axis([0 a 0 b 0 c]);
fldrs = strsplit(folder,'/');
title(strcat('Subject ',fldrs(6),'; Selection: ',rgx,'; ThrCoef: ',num2str(thr_val),'; Showing over: ',num2str(limit)));
%% Creating middle pictures
if show_brain
    I = imread(inp_file);
    [pic_a,pic_b] = size(I);
    pic_s = 0;
    good = false;
    while ~good & pic_s<pic_a-1
        pic_s = pic_s + 1;
        if(I(pic_s,:) == (ones(1,pic_b)*255) & mod(pic_a,pic_s)==0)
            good = true;
        end
    end
    pic_x = pic_a/pic_s;
    pic_y = pic_b/pic_s;
    
    pic_pos_y = floor((c/2)/pic_y);
    pic_pos_x = mod(c/2,pic_x)-1;
    
    mid_picture_I = I(pic_pos_x*pic_s+1:(pic_pos_x+1)*pic_s-1,pic_pos_y*pic_s+1:(pic_pos_y+1)*pic_s-1);
    mid_picture_I = fliplr(mid_picture_I);
    
    [X, Y, Z] = meshgrid(1:a,1:b,c/2);
    w = warp(Y,X,Z,mid_picture_I);
    set(w,'FaceAlpha',0.2);
end
hold off;
