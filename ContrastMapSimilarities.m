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

%%
% Specify the MAT files to use.
folder = pwd;
% Regex of the files to analyse.
rgx = '*V0.mat';
% Coefficient to use in case of a threshold defined compared to the
% average.
thr_coef_mult = 2;
% Type of threshold.
%   0: mean * coef.
%   1: -log(p) * coef.
thr_type = 1;
% Thresholp p value:
thr_p = 0.01;
% Minimal number of similarities for which are shown on the graph.
limit = 3;

% Coloring type to apply.
% 0: Color depends of the number of times a similarity appears.
% 1: Color depends of what contast map are here.
clr_type = 1;

% Color scheme to apply.
clr_sceme = [1 0 0;0 0 1;0 1 0;0 1 1;1 1 0;1 0 1];

% Getting all the files coresponding to the regex.
tmp = dir(fullfile(folder,rgx));
nb_files= numel(tmp);
for ii = 1:nb_files
    mat_file{ii} = fullfile( tmp(ii).name); 
    assert(exist(mat_file{ii}, 'file')>0)
end

% Preallocting all the needed matrices.
temp = zeros(128,128,54);
clr_temp = zeros(128,128,54,3);
maxV = zeros(0,nb_files);
minV = zeros(0,nb_files);
avgV = zeros(0,nb_files);
absV = zeros(0,nb_files);

% Checking the values of each file.
for ii = 1:nb_files
   L = load(mat_file{ii});
   maxV(1,ii) = max(max(max(L.map{1})));
   minV(1,ii) = min(min(min(L.map{1})));
   avgV(1,ii) = mean(mean(mean(L.map{1})));
   absV(1,ii) = mean(mean(mean(abs(L.map{1}))));
end

disp('max')
disp(maxV)
disp('minV')
disp(minV)
disp('avgV')
disp(avgV)
disp('absV')
disp(absV)

% Calculating the threshold depending of the type.
switch thr_type
    case 0
        threshold = mean(absV)*thr_coef_mult;
        thr_val = thr_coef_mult;
    case 1
        threshold = -log(thr_p)*thr_coef_mult;
        thr_val = thr_p;
end

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
                        clr_temp(i,j,k,color_i) = (clr_temp(i,j,k,color_i)/temp(i,j,k))*(temp(i,j,k)-1)+clr_sceme(ii,color_i)/temp(i,j,k);
                   end
               end
           end
       end
   end
end

%Find head borders
% temp_c = L.co{1};
% [x,y,z] = size(temp_c);
% for i =1:z
%     for j=1:y
%         found_start = false;
%         found_end = false;
%         counter = 1;
%         while ~found_start && counter < x
%             if temp_c(counter,j,i) > 0.015
%                 found_start = true;
%                 temp(counter,j,i) = 100;
%             end
%             counter = counter +1;
%         end
%         counter = x;
%         while ~found_end && counter > 0
%             if temp_c(counter,j,i) > 0.015
%                 found_end = true;
%                 temp(counter,j,i) = 100;
%             end
%             counter = counter -1;
%         end
%     end
% end

%color the brain
[a,b,c] = size(temp);
figure
for i=1:a
    if mod(i,10)==0
        disp(i)
    end
    for j=1:b
        for k=1:c
            if temp(i,j,k)>=limit && clr_type==0
                col='';
                alp = 0.7;
                switch temp(i,j,k)
                    case nb_files-5
                        col = clr_sceme(6,:);
                    case nb_files-4
                        col = clr_sceme(5,:);
                    case nb_files-3
                        col = clr_sceme(4,:);
                    case nb_files-2
                        col = clr_sceme(3,:);
                    case nb_files-1
                        col = clr_sceme(2,:);
                    case nb_files
                        col = clr_sceme(1,:);
                end
                voxel([i j k],[1 1 1],col,alp);
            end
            if temp(i,j,k)>=1 && clr_type==1
                voxel([i j k],[1 1 1],clr_temp(i,j,k,:),0.8);
            end
        end
    end
end
rotate3d on
axis([0 a 0 b 0 c]);
fldrs = strsplit(folder,'/');
title(strcat('Subject ',fldrs(6),'; Selection: ',rgx,'; ThrCoef: ',num2str(thr_val),'; Showing over: ',num2str(limit)));
%% Creating middle pictures
[I,pic_map] = imread('Images/inplane.png');
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
%mid_picture_map = pic_map(pic_pos_x*pic_s+1:(pic_pos_x+1)*pic_s,pic_pos_y*pic_s+1:(pic_pos_y+1)*pic_s,:);
%%
[X,Y] = meshgrid(0:a,0:b);
Z = c/2;
warp(X,Y,Z,mid_picture_I);
