clear all
addpath('E:/GitHub/LLG/Code/MATLAB/BrainVis/NII')

% model = 1; Between
% model = 2; Abar
% model = 3; Phat
model = 1;

% Vertices to visualize
if (model == 1)
    fileName = '../../../Result/Vertex_Diff_Between_desikan.csv';
elseif (model == 2)
    fileName = '../../../Result/Vertex_Diff_Abar_desikan.csv';
else
    fileName = '../../../Result/Vertex_Diff_Phat_desikan.csv';
end
vertexVec = csvread(fileName);

% Load atlas
mask = load_nii('desikan.nii');
mask.img = double(mask.img);
maskCopy = mask.img;

% Get brain background
brainBackground = 1*((~ismember(mask.img, vertexVec))&(mask.img~=0));
brainBackground(brainBackground==0)=NaN;

mask.img((~ismember(mask.img, vertexVec))&(mask.img~=0))=0;

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

if (model == 1)
    fileName = '../../../Result/Edge_Diff_Between_desikan.csv';
elseif (model == 2)
    fileName = '../../../Result/Edge_Diff_Abar_desikan.csv';
else
    fileName = '../../../Result/Edge_Diff_Phat_desikan.csv';
end
edgeVec = csvread(fileName);
edgeVec = reshape(edgeVec, 3, length(edgeVec)/3)';
edgeMagnitude = edgeVec(:,3);
edgeVec = edgeVec(:,1:2);

minMag = min(abs(edgeMagnitude));
maxMag = max(abs(edgeMagnitude));
minWidth = 1;
maxWidth = 5;

edgeMagnitude = sign(edgeMagnitude).* ...
    ((maxWidth - minWidth)/(maxMag - minMag)* ...
    (abs(edgeMagnitude) - minMag) + minWidth);

%%

figure(1);

p1 = subplot(1,3,1);

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

hold on;

for i = 1:size(edgeVec,1)
    cen = [centroids((edgeVec(i,1)),:);
        centroids((edgeVec(i,2)),:)];
    if (edgeMagnitude(i) > 0)
        plot3(cen(:,2), cen(:,1), cen(:,3), '-r', ...
            'LineWidth', edgeMagnitude(i));
    else
        plot3(cen(:,2), cen(:,1), cen(:,3), '-b', ...
            'LineWidth', -edgeMagnitude(i));
    end
end
hold off




p2 = subplot(1,3,2);

FV = {};
FV{1} = isosurface(double(maskCopy~=0), .5);
colors = distinguishable_colors(21,'w');
p = patch(FV{1});
isonormals(double(maskCopy~=0),p)
p.FaceColor = colors(1,:);
p.EdgeColor = 'none';
p.FaceAlpha = .15;
daspect([1,1,1])
view([0,1,0]); axis tight
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

hold on;

for i = 1:size(edgeVec,1)
    cen = [centroids((edgeVec(i,1)),:);
        centroids((edgeVec(i,2)),:)];
    if (edgeMagnitude(i) > 0)
        pp1 = plot3(cen(:,2), cen(:,1), cen(:,3), '-r', ...
            'LineWidth', edgeMagnitude(i));
    else
        pp2 = plot3(cen(:,2), cen(:,1), cen(:,3), '-b', ...
            'LineWidth', -edgeMagnitude(i));
    end
end
hold off

p3 = subplot(1,3,3);

FV = {};
FV{1} = isosurface(double(maskCopy~=0), .5);
colors = distinguishable_colors(21,'w');
p = patch(FV{1});
isonormals(double(maskCopy~=0),p)
p.FaceColor = colors(1,:);
p.EdgeColor = 'none';
p.FaceAlpha = .15;
daspect([1,1,1])
view([0,0,1]); axis tight
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

hold on;

for i = 1:size(edgeVec,1)
    cen = [centroids((edgeVec(i,1)),:);
        centroids((edgeVec(i,2)),:)];
    if (edgeMagnitude(i) > 0)
        pp3 = plot3(cen(:,2), cen(:,1), cen(:,3), '-r', ...
            'LineWidth', edgeMagnitude(i));
    else
        pp4 = plot3(cen(:,2), cen(:,1), cen(:,3), '-b', ...
            'LineWidth', -edgeMagnitude(i));
    end
end
hold off

hL = legend([pp1,pp2], {'Overestimate','Underestimate'}, ...
    'Location','South');
set(hL,'FontSize',20)

p1Pos = get(p1,'position');
p2Pos = get(p2,'position');
p2Pos(2:4) = p1Pos(2:4);
set(p2, 'position', p2Pos);

% if (model == 1)
%     figureName = '../../../Result/Diff_Between_desikan.png';
% elseif (model == 2)
%     figureName = '../../../Result/Diff_Abar_desikan.png';
% else
%     figureName = '../../../Result/Diff_Phat_desikan.png';
% end
% 
% r = 150; % pixels per inch
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 1080 480]/r);
% print(gcf,'-dpng',sprintf('-r%d',r), figureName);

