% Woii 240616
% for additional wing foil LES

SourceFolder = '../../../airfoilLES/airfoilLES_midspan/';
AFrame = 1:320;                % Frame to adding noise
OutputFile = 'V_NOISE.mat';    % Output file name
NoiseLevel = 0.25;             % std of noise
FieldStd = 0.07;               % std of velocity field (estimated)
tic

% loading LES grid
FileGrid = '../../../airfoilLES/airfoilLES_grid.h5';
x0 = h5read(FileGrid, '/x');
y0 = h5read(FileGrid, '/y');
xa = h5read(FileGrid, '/xa'); % coordinates of the airfoil
ya = h5read(FileGrid, '/ya');

% interpolation grid
dx = 0.01;
xb = 1.5:dx:3;
yb = -0.5:dx:0.5;
[X, Y] = meshgrid(xb, yb);
toc

rng(19260817);                 % random number seed

disp('processing flow field data...');
U_clean = zeros(numel(X),length(AFrame));
V_clean = U_clean; W_clean = U_clean;
U_noisy = U_clean; V_noisy = U_clean; W_noisy = U_clean;
iCount = 0;
for iFrame = AFrame
    iCount = iCount + 1;
    % loading data
    sub_file = [SourceFolder,...
        'airfoilLES_t', num2str(iFrame,'%05u'), '.h5'];
    u0 = h5read(sub_file, '/ux');
    v0 = h5read(sub_file, '/uy');
    % w0 = h5read(sub_file, '/uz');
    F = scatteredInterpolant(x0, y0, double(u0)); u = F(X, Y);
    F = scatteredInterpolant(x0, y0, double(v0)); v = F(X, Y);
    % F = scatteredInterpolant(x0, y0, double(w0)); w = F(X, Y);

    U_clean(:,iCount) = u(:);
    V_clean(:,iCount) = v(:);
    % W_clean(:,iCount) = w(:);

    % adding noise
    u = u + FieldStd*NoiseLevel*randn(size(u));
    v = v + FieldStd*NoiseLevel*randn(size(v));
    % w = w + FieldStd*NoiseLevel*randn(size(w));

    U_noisy(:,iCount) = u(:);
    V_noisy(:,iCount) = v(:);
    % W_noisy(:,iCount) = w(:);
end
toc

% saving
disp('saving data...');
save(OutputFile, 'xb','yb','AFrame',...
    'U_clean','V_clean','U_noisy','V_noisy', '-v7.3');
    % 'U_clean','V_clean','W_clean','U_noisy','V_noisy','W_noisy', '-v7.3');
toc

%%
v = reshape(U_clean(:,1), size(X));
figure; pcolor(X, Y, v); shading flat; colormap jet; axis equal; colorbar;
hold on; fill(xa, ya, 'k');
v = reshape(U_noisy(:,1), size(X));
figure; pcolor(X, Y, v); shading flat; colormap jet; axis equal; colorbar;
hold on; fill(xa, ya, 'k');
FieldStd = std([U_clean-mean(U_clean,2);V_clean-mean(V_clean,2)],1,'all');
disp(['Standard deviation of the flow: ', num2str(FieldStd)]);