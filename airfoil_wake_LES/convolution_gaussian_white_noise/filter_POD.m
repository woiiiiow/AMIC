% Woii 240627
% for additional wing foil LES

SourceFolder = '../../../airfoilLES/airfoilLES_midspan/';
PathField = 'V_NOISE.mat';
FileGrid = '../../../airfoilLES/airfoilLES_grid.h5';
AFrame_ext = 1:3200;           % Frame for POD spectrum
OutputFile = 'V_POD.mat';      % Output file name
NoiseLevel = 0.25*3.5;             % std of noise
FieldStd = 0.07;               % std of velocity field (estimated)
Th_sigma = sqrt(0.99);         % filter parameter
tic

% loading
load(PathField);
[X, Y] = meshgrid(xb, yb);
x0 = h5read(FileGrid, '/x');
y0 = h5read(FileGrid, '/y');

rng(19260817);                 % random number 
disp('loading and processing flow field data...');
U_ext = zeros(numel(X),length(AFrame_ext));
V_ext = U_ext; W_ext = U_ext;
iCount = 0;
for iFrame = AFrame_ext
    iCount = iCount + 1;
    % loading data
    sub_file = [SourceFolder,...
        'airfoilLES_t', num2str(iFrame,'%05u'), '.h5'];
    u0 = h5read(sub_file, '/ux');
    v0 = h5read(sub_file, '/uy');
    F = scatteredInterpolant(x0, y0, double(u0)); u = F(X, Y);
    F = scatteredInterpolant(x0, y0, double(v0)); v = F(X, Y);

    % adding noise
    u = u + FieldStd*NoiseLevel*imgaussfilt(randn(size(u)),1);
    v = v + FieldStd*NoiseLevel*imgaussfilt(randn(size(v)),1);

    U_ext(:,iCount) = u(:);
    V_ext(:,iCount) = v(:);
end
toc

disp('POD and truncation...');
Um = mean(U_ext, 2);
Vm = mean(V_ext, 2);
[PsiU, SigmaU, PhiU] = svd([U_ext-Um; V_ext-Vm]', 'econ');
sigma = diag(SigmaU);
b_sigma = sigma(2:end)./sigma(1:end-1) < Th_sigma;
b_sigma = circshift(b_sigma,-1) | b_sigma | circshift(b_sigma, 1);
mode_end = find(~b_sigma,1);
disp(['truncating mode: ', num2str(mode_end)]);
psi = [U_noisy-Um; V_noisy-Vm]'*PhiU/SigmaU;
u_filt = psi(:,1:mode_end)*SigmaU(1:mode_end,1:mode_end)*...
    transpose(PhiU(:,1:mode_end));
U = transpose(u_filt(:,           1:1*numel(X))) + Um;
V = transpose(u_filt(:,  numel(X)+1:2*numel(X))) + Vm;
toc

% saving
save(OutputFile, 'U','V','AFrame');

%% display
figure; iFrame = 24;
subplot(3,2,1); pcolor(X, Y, reshape(U_clean(:,iFrame),size(X)));
shading flat; axis equal; colormap jet; caxis([0.6 1.1]);
subplot(3,2,3); pcolor(X, Y, reshape(U_noisy(:,iFrame),size(X)));
shading flat; axis equal; colormap jet; caxis([0.6 1.1]);
subplot(3,2,5); pcolor(X, Y, reshape(U      (:,iFrame),size(X)));
shading flat; axis equal; colormap jet; caxis([0.6 1.1]);
subplot(3,2,2); pcolor(X, Y, reshape(V_clean(:,iFrame),size(X)));
shading flat; axis equal; colormap jet; caxis([-0.25 0.25]);
subplot(3,2,4); pcolor(X, Y, reshape(V_noisy(:,iFrame),size(X)));
shading flat; axis equal; colormap jet; caxis([-0.25 0.25]);
subplot(3,2,6); pcolor(X, Y, reshape(V      (:,iFrame),size(X)));
shading flat; axis equal; colormap jet; caxis([-0.25 0.25]);