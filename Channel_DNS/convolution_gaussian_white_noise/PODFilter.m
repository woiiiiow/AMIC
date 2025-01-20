SourceFile = 'V_NOISE.mat';
OutputFile = 'V_POD.mat';
AFrame = 1:120;
% grid
dt = 0.0065;   dx = 0.0114;
xb = 0:dx:1;   yb = 0:dx:1;   zb = 0:dx:0.5;   zb = zb(1:12);
[X,Y,Z] = meshgrid(xb,yb,zb);
% filter parameter
Th_sigma = sqrt(0.99);
tic

load(SourceFile, 'U_array','V_array','W_array',...
    'Um','Vm','Wm','SigmaU','PhiU');
disp('smoothing fields on POD spectrum...');
sigma = diag(SigmaU);
b_sigma = sigma(2:end)./sigma(1:end-1) < Th_sigma;
b_sigma = circshift(b_sigma,-1) | b_sigma | circshift(b_sigma, 1);
mode_end = find(~b_sigma,1);
disp(['ending mode: ', num2str(mode_end)]);
psi = [U_array(:,AFrame)-Um; V_array(:,AFrame)-Vm;...
    W_array(:,AFrame)-Wm]'*PhiU/SigmaU;
u_filt = psi(:,1:mode_end)*SigmaU(1:mode_end,1:mode_end)*...
    transpose(PhiU(:,1:mode_end));
U = transpose(u_filt(:,           1:1*numel(X))) + Um;
V = transpose(u_filt(:,  numel(X)+1:2*numel(X))) + Vm;
W = transpose(u_filt(:,2*numel(X)+1:3*numel(X))) + Wm;
toc

% save
save(OutputFile, 'xb','yb','zb','dt', 'U','V','W');