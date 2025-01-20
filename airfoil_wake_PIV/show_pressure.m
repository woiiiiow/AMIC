SourceFolder = '../../Experiment_wing3/OUT_TRPIV/';
load('Pressure.mat');
ROI = [41 110 41 150];         % subdomain of PIV field
                               % [first-last row, first-last column]
Mm_per_px_ratio = 0.1197;      % spatial resolution of photos
Sample_rate = 30;              % sample rate of PIV
Vector_Spacing = 10;           % vector per pixel in PIV

ShowFrame = 44;                % single frame to show (44)
Caxis = [-0.4 0.4];

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

%%
p = reshape(P_PIV(:,ShowFrame-AFrame_P(1)+1), size(X));
figure; hold on; pcolor(X, Y, p); shading flat; axis equal;
colormap jet; caxis(Caxis); colorbar('northoutside');
title(['p piv at frame ', num2str(ShowFrame)]); AddMark;

%%
p = reshape(P_SG(:,ShowFrame-AFrame_P(1)+1), size(X));
figure; hold on; pcolor(X, Y, p); shading flat; axis equal;
colormap jet; caxis(Caxis); colorbar('northoutside');
title(['p sg at frame ', num2str(ShowFrame)]); AddMark;

%%
p = reshape(P_POD(:,ShowFrame-AFrame_P(1)+1), size(X));
figure; hold on; pcolor(X, Y, p); shading flat; axis equal;
colormap jet; caxis(Caxis); colorbar('northoutside');
title(['p pod at frame ', num2str(ShowFrame)]); AddMark;

%%
p = reshape(P_AE(:,ShowFrame-AFrame_P(1)+1), size(X));
figure; hold on; pcolor(X, Y, p); shading flat; axis equal;
colormap jet; caxis(Caxis); colorbar('northoutside');
title(['p ae at frame ', num2str(ShowFrame)]); AddMark;

%%
p = reshape(P_AMIC(:,ShowFrame-AFrame_P(1)+1), size(X));
figure; hold on; pcolor(X, Y, p); shading flat; axis equal;
colormap jet; caxis(Caxis); colorbar('northoutside');
title(['p amic at frame ', num2str(ShowFrame)]); AddMark;

%% functions
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