clear all
% addpath('E:/GitHub/LLG/Code/MATLAB/BrainVis/NII')
addpath('./NII')

fileName = '../../../Result/eigenvector.csv';
xHat = reshape(csvread(fileName), 70, 2);
xHat = xHat(:, 2);
max(xHat)
min(xHat)

% Load atlas
mask = load_nii('desikan.nii');
mask.img = double(mask.img);
maskCopy = mask.img;

% Get brain background
brainBackground = 1*((mask.img~=0));
brainBackground(brainBackground==0)=NaN;

mask.img((mask.img~=0))=0;

%%
n = length(unique(maskCopy))-1;
centroids = zeros(n,3);
for x = 1:size(maskCopy,1)
    for y = 1:size(maskCopy,2)
        for z = 1:size(maskCopy,3)
            if (maskCopy(x,y,z)~=0)
                centroids(maskCopy(x,y,z),:) = ...
                    centroids(maskCopy(x,y,z),:)+[x,y,z];
            end
        end
    end
end
for i = 1:n
    centroids(i,:) = centroids(i,:)/sum(sum(sum(maskCopy==i)));
end


%%

figure(1);

p1 = subplot(1,3,1);

FV = {};
FV{1} = isosurface(double(maskCopy~=0), .5);
% colors = distinguishable_colors(71,'w');
p = patch(FV{1});
isonormals(double(maskCopy~=0),p)
p.FaceColor = [0, 0, 0];
p.EdgeColor = 'none';
p.FaceAlpha = .15;
daspect([1,1,1])
view([1,0,0]); axis tight
lighting gouraud
material dull
light('Position', [100 100 100], 'Style', 'local')
set(gca,'Visible','off')
hold on

for i = 1:70
    V = double(maskCopy == i);
    FV{i+1} = isosurface(V, .5);
    p = patch(FV{i+1});
    isonormals(double(maskCopy~=0),p)
    p.FaceColor = [1, 1, (xHat(i)+1)/2];
    p.EdgeColor = 'none';
    p.FaceAlpha = 0.2;
    lighting gouraud
    material dull
end
% colorbar

hold off





p2 = subplot(1,3,2);

FV = {};
FV{1} = isosurface(double(maskCopy~=0), .5);
% colors = distinguishable_colors(71,'w');
p = patch(FV{1});
isonormals(double(maskCopy~=0),p)
p.FaceColor = [0, 0, 0];
p.EdgeColor = 'none';
p.FaceAlpha = .15;
daspect([1,1,1])
view([0,1,0]); axis tight
lighting gouraud
material dull
light('Position', [100 100 100], 'Style', 'local')
set(gca,'Visible','off')
hold on

for i = 1:70
    V = double(maskCopy == i);
    FV{i+1} = isosurface(V, .5);
    p = patch(FV{i+1});
    isonormals(double(maskCopy~=0),p)
    p.FaceColor = [1, 1, (xHat(i)+1)/2];
    p.EdgeColor = 'none';
    p.FaceAlpha = 0.2;
    lighting gouraud
    material dull
end
% colorbar

hold off




p3 = subplot(1,3,3);

FV = {};
FV{1} = isosurface(double(maskCopy~=0), .5);
% colors = distinguishable_colors(71,'w');
p = patch(FV{1});
isonormals(double(maskCopy~=0),p)
p.FaceColor = [0, 0, 0];
p.EdgeColor = 'none';
p.FaceAlpha = .15;
daspect([1,1,1])
view([0,0,1]); axis tight
lighting gouraud
material dull
light('Position', [100 100 100], 'Style', 'local')
set(gca,'Visible','off')
hold on

for i = 1:70
    V = double(maskCopy == i);
    FV{i+1} = isosurface(V, .5);
    p = patch(FV{i+1});
    isonormals(double(maskCopy~=0),p)
    p.FaceColor = [1, 1, (xHat(i)+1)/2];
    p.EdgeColor = 'none';
    p.FaceAlpha = 0.2;
    lighting gouraud
    material dull
end
% colorbar

hold off




% p1Pos = get(p1,'position');
% p2Pos = get(p2,'position');
% p2Pos(2:4) = p1Pos(2:4);
% set(p2, 'position', p2Pos);