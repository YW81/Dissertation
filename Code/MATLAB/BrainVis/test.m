clear all
addpath('/Users/Runze/Documents/GitHub/LLG/Code/MATLAB/BrainVis/NII')

vertexVec = [34,38];

% Load atlas
mask = load_nii('desikan.nii');
mask.img = double(mask.img);
maskCopy = mask.img;

% Get brain background
brainBackground = 1*((~ismember(mask.img, vertexVec))&(mask.img~=0));
brainBackground(brainBackground==0)=NaN;

mask.img((~ismember(mask.img, vertexVec))&(mask.img~=0))=0;

figure(1);

FV = {};
FV{1} = isosurface(double(maskCopy~=0), .5);
colors = distinguishable_colors(21,'w');
p = patch(FV{1});
isonormals(double(maskCopy~=0),p)
p.FaceColor = colors(1,:);
p.EdgeColor = 'none';
p.FaceAlpha = .15;
daspect([1,1,1])
view([1,0,0]); axis tight
lighting gouraud
material dull
light('Position', [100 100 100], 'Style', 'local')
set(gca,'Visible','off')
hold on

for i = 1:length(vertexVec)
    V = double(mask.img == vertexVec(i));
    FV{i+1} = isosurface(V, .5);
    p = patch(FV{i+1});
    isonormals(double(maskCopy~=0),p)
    p.FaceColor = colors(i+1,:);
    p.EdgeColor = 'none';
    p.FaceAlpha = 1;
    lighting gouraud
    material dull
end

%%

fileName = '../../../Data/desikan/centroid.csv';
centroids = csvread(fileName);
centroids = reshape(centroids, 3, length(centroids)/3)';
centroids = centroids(1:2, :);

centroids = zeros(70,3);
for i = 1:70
    for x = 1:size(maskCopy,1)
        for y = 1:size(maskCopy,2)
            for z = 1:size(maskCopy,3)
                if (maskCopy(x,y,z) == i)
                    centroids(i,:) = centroids(i,:)+[x,y,z];
                end
            end
        end
    end
    centroids(i,:) = centroids(i,:)/sum(sum(sum(maskCopy==i)));
end


edgeVec = [34,38];

FV = {};
FV{1} = isosurface(double(maskCopy~=0), .5);
colors = distinguishable_colors(21,'w');
p = patch(FV{1});
isonormals(double(maskCopy~=0),p)
p.FaceColor = colors(1,:);
p.EdgeColor = 'none';
p.FaceAlpha = .15;
daspect([1,1,1])
view([1,0,0]); axis tight
lighting gouraud
material dull
light('Position', [100 100 100], 'Style', 'local')
set(gca,'Visible','off')

hold on;

% for i = 1:length(edgeVec)
%     cen = [centroids((edgeVertex==edgeVec(i,1)),:);
%         centroids((edgeVertex==edgeVec(i,2)),:)];
%     plot3(cen(:,1),cen(:,2), cen(:,3),'-o','LineWidth',1.5)
% end
for i = 1:size(edgeVec,1)
    cen = [centroids((edgeVec(i,1)),:);
        centroids((edgeVec(i,2)),:)];
    plot3(cen(:,2),cen(:,1), cen(:,3),'-o','LineWidth',1.5)
end
hold off
