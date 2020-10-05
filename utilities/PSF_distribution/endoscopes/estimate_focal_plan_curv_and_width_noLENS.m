function [curv_radius, r_z, z_depth] = estimate_focal_plan_curv_and_width_noLENS

if ispc
    file_list = dir('..\..\data\fluorescent films per simulazione\500NoLens\');
elseif isunix
    file_list = dir('../../data/fluorescent films per simulazione/500NoLens/');
end
curv_radius = zeros(1,length(file_list)-2);
r_z = zeros(1,length(file_list)-2);
um_per_px_Z = 1;
um_per_px_X = 0.91;
um_per_px_Y = 0.91;
for id_list = 3:length(file_list)
filename = [file_list(id_list).folder '/' file_list(id_list).name];%'E:\data\endoscopes\fluorescent films per simulazione\500Lens\ZSeries-05162016-1202-2553.tif';

info = imfinfo(filename);
intensity_profile = uint16(zeros(info(1).Width, info(1).Height));
for i = 1:length(info)
    intensity_profile(:,:,i) = imread(filename,i);
end


% xz_proj = flipud(squeeze(mean(intensity_profile,2))');
xz_proj = flipud(squeeze(intensity_profile(:,round(size(intensity_profile,2)/2),:))');
% xz_proj_pad = padarray(xz_proj,[500 100]); 
figure; subplot(4,1,1);
imagesc(xz_proj); colormap(gray); axis image;
thr1 = prctile(xz_proj(:),91)';
xz_bin = double(xz_proj>=thr1);
subplot(4,1,2);
imagesc(xz_bin); colormap(gray); axis image;
if mod(size(xz_bin,2),2)==0
    z_min = mean([mean(find(xz_bin(:,size(xz_bin,2)/2))) mean(find(xz_bin(:,size(xz_bin,2)/2 + 1)))]);
else
    z_min = mean(find(xz_bin(:,(size(xz_bin,2)+1)/2)));
end
% min_single_col = zeros(1,size(xz_bin,2));
r_z_temp = zeros(1,size(xz_bin,2));
for i = 1:size(xz_bin,2)
    if ~isempty(find(xz_bin(:,i)))
%         min_single_col(i) = mean(find(xz_bin(:,i)));
        r_z_temp(i) = length(find(xz_bin(:,i)));
    else
%         min_single_col(i) = NaN;
        r_z_temp(i) = NaN;
    end
end
% [~,x_min] = max(min_single_col);
% z_min = mean(find(xz_bin(:,x_min)));
z_depth(id_list-2) = z_min*um_per_px_Z;
% r_z_temp = zeros(1,size(xz_bin,2));
% for i = 1:size(xz_bin,2)
%     if ~isempty(find(xz_bin(:,i)))
%         r_z_temp(i) = length(find(xz_bin(:,i)));
%     else
%         r_z_temp(i) = NaN;
%     end
% end
r_z(id_list-2) = nanmean(r_z_temp)*um_per_px_Z/2;
xz_proj2=xz_proj;
% xz_proj2(find(xz_bin))=0;
% xz_proj2=xz_proj2(1:z_min,:);
xz_proj2=xz_proj(1:50,:);
subplot(4,1,3);
imagesc(xz_proj2); colormap(gray); axis image;
thr2 = prctile(xz_proj2(:),80)';
xz_bin = double(xz_proj2>=thr2);
subplot(4,1,4);
imagesc(xz_bin); colormap(gray); axis image;
% dy = 45;
% row = round(z_min)-dy;
row = 3;
dy = round(z_min)-row;
if mod(size(xz_bin,2),2)==0
    dx = find(xz_bin(row,:));
    dx1 = dx(1);
    dx1 = size(xz_bin,2)/2 - dx1;
    dx2 = dx(end);
    dx2 = dx2-size(xz_bin,2)/2;
else
    dx = find(xz_bin(row,:));
    dx1 = dx(1);
    dx1 = (size(xz_bin,2)-1)/2 - dx1;
    dx2 = dx(end);
    dx2 = dx2-(size(xz_bin,2)+1)/2;
end
dx = round(nanmax([dx1 dx2]))*um_per_px_X;
% dy = (z_min-row)*um_per_px_Z;
dy = dy*um_per_px_Z;
close;

figure;
xx= 1:1:size(xz_proj,2); xx = xx*um_per_px_X;
zz= 1:1:size(xz_proj,1);
imagesc(xx,zz,xz_proj); colormap(gray); axis image;
hold on;
if mod(size(xz_bin,2),2)==0
    scatter(size(xz_bin,2)*um_per_px_X/2,z_min,30,'r');
%     scatter(x_min*um_per_px_X,z_min,30,'r');
    hold on;
    scatter(size(xz_bin,2)*um_per_px_X/2 - dx ,row,30,'r');
%     scatter(x_min*um_per_px_X - dx ,row,30,'r');
    hold on;
    scatter(size(xz_bin,2)*um_per_px_X/2 + dx ,row,30,'r');
%     scatter(x_min*um_per_px_X + dx ,row,30,'r');
else
    scatter((size(xz_bin,2)+1)*um_per_px_X/2,z_min,30,'r');
    hold on;
    scatter( (size(xz_bin,2)+1)*um_per_px_X/2 - dx ,row, 30,'r');
    hold on;
    scatter((size(xz_bin,2)+1)*um_per_px_X/2 + dx ,row,30,'r');
end
% close;

alpha = acos(dy/dx);
beta = pi/2-alpha;
gamma = alpha-beta;
teta = pi/2-gamma;
curv_radius(id_list-2) = dx/cos(teta);
end

close all;

curv_radius = nanmean(curv_radius);
r_z = nanmean(r_z);
z_depth = nanmean(z_depth);