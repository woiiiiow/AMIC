SourceFolder = '../../Experiment_wing3/OUT_TRPIV/';
OutputFile = 'V_POD.mat';
ROI = [41 110 41 150];         % subdomain of PIV field
                               % [first-last row, first-last column]
Mm_per_px_ratio = 0.1197;      % spatial resolution of photos
Sample_rate = 30;              % sample rate of PIV
Vector_Spacing = 10;           % vector per pixel in PIV

AFrame = 11:530;               % frames to be saved
AFrame_train = 11:3610;        % larger dataset for POD
Th_sigma = sqrt(0.99);         % filter parameter

load([SourceFolder, 'Grid_Wing.mat'], 'X','Y');
% entire domain coordinates
X_etr = X;                     Y_etr = Y;
% AFrame_save in AFrame
AFrame_sub = AFrame - AFrame_train(1) + 1;
tic

disp('loading data...');
U_array = zeros(numel(X_etr), length(AFrame_train));   V_array = U_array;
iCount = 0;
for iFrame = AFrame_train
    iCount = iCount + 1;
    load([SourceFolder,'Wing_',num2str(iFrame,'%06u'),'.mat'], 'U','V');
    U_array(:, iCount) = V(:);    V_array(:, iCount) = U(:);
    % swapping U and V
end
U_array = U_array.*Sample_rate.*Mm_per_px_ratio.*1e-3;
V_array = V_array.*Sample_rate.*Mm_per_px_ratio.*1e-3;
clear iFrame iCount U V
toc

disp('POD and truncation...');
Um = mean(U_array, 2);
Vm = mean(V_array, 2);
[PsiU, SigmaU, PhiU] = svd([U_array-Um; V_array-Vm]', 'econ');
sigma = diag(SigmaU);
b_sigma = sigma(2:end)./sigma(1:end-1) < Th_sigma;
b_sigma = circshift(b_sigma,-1) | b_sigma | circshift(b_sigma, 1);
mode_end = find(~b_sigma,1);
disp(['truncating mode: ', num2str(mode_end)]);
psi = [U_array(:,AFrame_sub)-Um; V_array(:,AFrame_sub)-Vm]'*PhiU/SigmaU;
u_filt = psi(:,1:mode_end)*SigmaU(1:mode_end,1:mode_end)*...
    transpose(PhiU(:,1:mode_end));
U = transpose(u_filt(:,           1:1*numel(X))) + Um;
V = transpose(u_filt(:,  numel(X)+1:2*numel(X))) + Vm;
toc

% save
save(OutputFile, 'U','V','AFrame');
