function place_field_matrix = define_place_field(firing_rate_map,imageX,imageY)
%Written by Seth Koenig 9/19/2016 taken from various other codes.
%code defines place field by clustering firing rate map pixels into groups
%using kmeans clustering. The highest cluster is defined as the place
%field.
%
% Inputs:
%   1) Firing_rate_map: binned firing rate map for place field
%   2) imageX: horizontal size of display e.g. 800
%   3) imageY: vertical size of display e.g. 600
%
% Outputs:
%   1) place_field_matrix: pixel by pixel maxtrix of the place field
%   with 1's denoting the place field and 0's denoting outside place field


%---Normalize to remove outliers at the very high end---%
fr = firing_rate_map(:);
fr(isnan(fr)) = [];
maxfr = prctile(fr,95); %95 percentile of firing rate
fr(fr > maxfr) = maxfr; %so doesn't screw up fitting of low firing rate data

%---Define Place Field using Kmeans Clustering---%
sil = zeros(1,3); %determines the number of clusters by comparing the ratio
%of intercluster and intracluster distances, faster mod of silhouette
%try 2 or 3 Clusters depending on the distribution may have a moderate firing
%rate cluster that's in between the place field and outside the place field
for numclusts = 2:3
    T = kmeans(fr,numclusts,'replicate',5);
    [silh] = InterVSIntraDist(fr,T);
    sil(numclusts) = mean(silh);
end
numclusters = find(sil == max(sil));
T = kmeans(fr,numclusters(end),'replicate',25);

region_mean_fr = zeros(1,numclusters);
for t = 1:numclusters
    region_mean_fr(t) = mean(fr(T == t));
end
[~,place_field_cluster] = max(region_mean_fr);
min_firing_rate = min(fr(T == place_field_cluster));

%---Make place_field_matrix---%
[r,c] = find(firing_rate_map > min_firing_rate);
firing_ind = sub2ind(size(firing_rate_map),r,c);
place_field_matrix = zeros(size(firing_rate_map));
place_field_matrix(firing_ind) = 1;

area = sum(sum(place_field_matrix));
%clean up place field matrix with missing pixels should not affect analysis
%since missing pixels occur when no eye data but will improve estimate of
%area/size of field and may help with sequence task analysis if item falls
%on missing spot
if area < 0.66*numel(place_field_matrix) %will not work well for very large field
    %---Fill in Holes due to low coverage---%
    %do upper left first
    pfm = padarray(place_field_matrix,[1 1],1,'pre');
    pfm = imfill(pfm,'holes');
    pfm = pfm(2:end,2:end);
    if sum(sum(pfm)) < 1.25*sum(sum(place_field_matrix)) %so not a whole ton was added
        place_field_matrix = pfm;
    end
    %then do lower right
    pfm = padarray(place_field_matrix,[1 1],1,'post');
    pfm = imfill(pfm,'holes');
    pfm = pfm(1:end-1,1:end-1);
    if sum(sum(pfm)) < 1.25*sum(sum(place_field_matrix)) %so not a whole ton was added
        place_field_matrix = pfm;
    end
else
    disp('big field')
end

%---Resize to Make Same as Image Size---%
place_field_matrix = imresize(place_field_matrix,[imageY,imageX],'method','nearest');%upsample to images size in pixel coordinates

end