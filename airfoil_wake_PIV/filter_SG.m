SourceFolder = '../../Experiment_wing3/OUT_TRPIV/';
OutputFile = 'V_SG.mat';
ROI = [41 110 41 150];         % subdomain of PIV field
                               % [first-last row, first-last column]
Mm_per_px_ratio = 0.1197;      % spatial resolution of photos
Sample_rate = 30;              % sample rate of PIV
Vector_Spacing = 10;           % vector per pixel in PIV

AFrame = 11:530;               % frames to be processed
[~,M_SG,~] = MatSG(5,1,1);     % Local Polynomial (Savitzky-Golay) filter

load([SourceFolder, 'Grid_Wing.mat'], 'X','Y');
% entire domain coordinates
X_etr = X;                     Y_etr = Y;
tic

disp('loading data...');
U_3D_array = zeros([size(X_etr), length(AFrame)]);
V_3D_array = U_3D_array;
iCount = 0;
for iFrame = AFrame
    iCount = iCount + 1;
    load([SourceFolder,'Wing_',num2str(iFrame,'%06u'),'.mat'], 'U','V');
    U_3D_array(:, :, iCount) = V;    V_3D_array(:, :, iCount) = U;
    % swapping U and V
end
clear iFrame iCount U V
toc

disp('smoothing fields using Savitzky-Golay filter...');
weight = convn(ones(size(U_3D_array)), M_SG, 'same');
U_3D_array = convn(U_3D_array, M_SG, 'same')./weight;
V_3D_array = convn(V_3D_array, M_SG, 'same')./weight;
U = reshape(U_3D_array(:,:,:), [numel(X_etr), length(AFrame)]);
V = reshape(V_3D_array(:,:,:), [numel(X_etr), length(AFrame)]);
toc

% save
U = U.*Sample_rate.*Mm_per_px_ratio.*1e-3;
V = V.*Sample_rate.*Mm_per_px_ratio.*1e-3;
save(OutputFile, 'U','V','AFrame');

%% functions
function [MM,M,Mdt] = MatSG(w,dx,dt)
% w: window size; dx: spatial increment in grid; dt: temporal ...

polynomial =@(i,j,k) [1,...
                      i*dx, j*dx, k*dt,...
                      i*j*dx^2, i*k*dx*dt, j*k*dx*dt,...
                      i*i*dx^2, j*j*dx^2, k*k*dt^2];

w2=floor(w/2);
range = -w2:w2;
M = zeros(numel(range)^3,numel(polynomial(1,1,1)));

cont=0;
for k=-range
    for j=range
        for i=range
            cont=cont+1;
            M(cont,:) = polynomial(i,j,k);
        end
    end
end

MM  = (M'*M)\(M');%inv(M'*M)*M';

convcoef =@(col,w) reshape(col,[w w w]);

M    = convcoef(MM(1,:),w);
Mdt  = convcoef(MM(4,:),w);
Mddx = convcoef(MM(9,:),w);
Mddy = convcoef(MM(8,:),w);
end