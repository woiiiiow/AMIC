SourceFolder = '../../Experiment_wing3/OUT_TRPIV/';
ROI = [41 110 41 150];         % subdomain of PIV field
                               % [first-last row, first-last column]
Mm_per_px_ratio = 0.1197;      % spatial resolution of photos
Sample_rate = 30;              % sample rate of PIV
Vector_Spacing = 10;           % vector per pixel in PIV

AFrame = 21:520;               % frames for statistics
ShowFrame = 24;                % single frame to show

load([SourceFolder, 'Grid_Wing.mat'], 'X','Y');
% entire domain coordinates
X_etr = X;                     Y_etr = Y;
% map of subdomain
map = false(size(X_etr));
map(ROI(1):ROI(2), ROI(3):ROI(4)) = true;
% subdomain coordinates normalized by the chord length
X = X_etr(ROI(1):ROI(2),ROI(3):ROI(4));
Y = Y_etr(ROI(1):ROI(2),ROI(3):ROI(4));
X = (X - X(1)).*Mm_per_px_ratio./80 + 0.70;
Y = (Y - Y(1)).*Mm_per_px_ratio./80 - 0.53;

V_PIV.U = zeros(numel(X_etr), length(AFrame));   V_PIV.V = V_PIV.U;
iCount = 0;
for iFrame = AFrame
    iCount = iCount + 1;
    load([SourceFolder,'Wing_',num2str(iFrame,'%06u'),'.mat'], 'U','V');
    V_PIV.U(:, iCount) = V(:);    V_PIV.V(:, iCount) = U(:);
    % swapping U and V
end
V_PIV.U = V_PIV.U.*Sample_rate.*Mm_per_px_ratio.*1e-3;
V_PIV.V = V_PIV.V.*Sample_rate.*Mm_per_px_ratio.*1e-3;
clear iFrame iCount U V

%%
u = reshape(V_PIV.U(:,ShowFrame), size(X_etr));
u = reshape(u(map(:)), size(X));
figure; hold on; pcolor(X, Y, u); shading flat; axis equal;
colormap jet; caxis([0 0.08]); colorbar('northoutside');
title('u piv'); AddMark;

u = reshape(V_PIV.V(:,ShowFrame), size(X_etr));
u = reshape(u(map(:)), size(X));
figure; hold on; pcolor(X, Y, u); shading flat; axis equal;
colormap jet; caxis([-0.02 0.02]); colorbar('northoutside');
title('v piv'); AddMark;

%%
V_SG = load('V_SG.mat', 'U','V','AFrame');
u = reshape(V_SG.U(:,ShowFrame+AFrame(1)-V_SG.AFrame(1)), size(X_etr));
u = reshape(u(map(:)), size(X));
figure; hold on; pcolor(X, Y, u); shading flat; axis equal;
colormap jet; caxis([0 0.08]); colorbar('northoutside');
title('u sg'); AddMark;

u = reshape(V_SG.V(:,ShowFrame+AFrame(1)-V_SG.AFrame(1)), size(X_etr));
u = reshape(u(map(:)), size(X));
figure; hold on; pcolor(X, Y, u); shading flat; axis equal;
colormap jet; caxis([-0.02 0.02]); colorbar('northoutside');
title('v sg'); AddMark;

%%
V_POD = load('V_POD.mat', 'U','V','AFrame');
u = reshape(V_POD.U(:,ShowFrame+AFrame(1)-V_POD.AFrame(1)), size(X_etr));
u = reshape(u(map(:)), size(X));
figure; hold on; pcolor(X, Y, u); shading flat; axis equal;
colormap jet; caxis([0 0.08]); colorbar('northoutside');
title('u pod'); AddMark;

u = reshape(V_POD.V(:,ShowFrame+AFrame(1)-V_POD.AFrame(1)), size(X_etr));
u = reshape(u(map(:)), size(X));
figure; hold on; pcolor(X, Y, u); shading flat; axis equal;
colormap jet; caxis([-0.02 0.02]); colorbar('northoutside');
title('v pod'); AddMark;

%%
V_AMIC = load('V_AMIC.mat', 'U','V','AFrame');
u = reshape(V_AMIC.U(:,ShowFrame+AFrame(1)-V_AMIC.AFrame(1)), size(X_etr));
u = reshape(u(map(:)), size(X));
figure; hold on; pcolor(X, Y, u); shading flat; axis equal;
colormap jet; caxis([0 0.08]); colorbar('northoutside');
title('u amic'); AddMark;

u = reshape(V_AMIC.V(:,ShowFrame+AFrame(1)-V_AMIC.AFrame(1)), size(X_etr));
u = reshape(u(map(:)), size(X));
figure; hold on; pcolor(X, Y, u); shading flat; axis equal;
colormap jet; caxis([-0.02 0.02]); colorbar('northoutside');
title('v amic'); AddMark;

%%
V_AE = load('V_AE.mat', 'U','V','AFrame');
u = reshape(V_AE.U(:,ShowFrame+AFrame(1)-V_AE.AFrame(1)), size(X));
figure; hold on; pcolor(X, Y, u); shading flat; axis equal;
colormap jet; caxis([0 0.08]); colorbar('northoutside');
title('u ae'); AddMark;

u = reshape(V_AE.V(:,ShowFrame+AFrame(1)-V_AE.AFrame(1)), size(X));
figure; hold on; pcolor(X, Y, u); shading flat; axis equal;
colormap jet; caxis([-0.02 0.02]); colorbar('northoutside');
title('v ae'); AddMark;

%%
function AddMark
NACA0018 = [
1.0000     0.00189
0.9500     0.01210
0.9000     0.02172
0.8000     0.03935
0.7000     0.05496
0.6000     0.06845
0.5000     0.07941
0.4000     0.08705
0.3000     0.09003
0.2500     0.08912
0.2000     0.08606
0.1500     0.08018
0.1000     0.07024
0.0750     0.06300
0.0500     0.05332
0.0250     0.03922
0.0125     0.02841
0.0000     0.00000
0.0125     -0.02841
0.0250     -0.03922
0.0500     -0.05332
0.0750     -0.06300
0.1000     -0.07024
0.1500     -0.08018
0.2000     -0.08606
0.2500     -0.08912
0.3000     -0.09003
0.4000     -0.08705
0.5000     -0.07941
0.6000     -0.06845
0.7000     -0.05496
0.8000     -0.03935
0.9000     -0.02172
0.9500     -0.01210
1.0000     -0.00189
];
x_wing = NACA0018(:,1)' - 1;
y_wing = NACA0018(:,2)';
phi = -10/180*pi;
tmp = [cos(phi),-sin(phi);sin(phi),cos(phi)]*[x_wing; y_wing];
x_wing = tmp(1,:);
y_wing = tmp(2,:);

fill(x_wing, y_wing, 'k');

xlim([-1.03 2.33]);
ylim([-0.5 0.5]);
end
