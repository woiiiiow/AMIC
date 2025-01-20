SourceFolder = '../../Channel/volumn_edge/Fields_testing/';
iFrame = 34;
iSlice = 6;
dx = 0.0114;  xb = 0:dx:1;   yb = 0:dx:1;   zb = 0:dx:0.5;   zb = zb(1:12);
[X,Y,Z] = meshgrid(xb,yb,zb);

%%
sub_name = ['Field_', num2str(iFrame,'%06u'), '.mat'];
load([SourceFolder, sub_name], 'u', 'v', 'w');
figure; imagesc(u(:,:,iSlice)'); axis equal; title('u');
colormap jet; caxis([0 1.2]); colorbar;

%%
V_NOISE = load('V_NOISE.mat', 'U_noise','V_noise','W_noise');
u = reshape(V_NOISE.U_noise(:,iFrame), size(X));
figure; imagesc(u(:,:,iSlice)'); axis equal; title('u noise');
colormap jet; caxis([0 1.2]); colorbar;

%%
V_SG = load('V_SG.mat', 'U','V','W');
u = reshape(V_SG.U(:,iFrame), size(X));
figure; imagesc(u(:,:,iSlice)'); axis equal; title('u sg');
colormap jet; caxis([0 1.2]); colorbar;

%%
V_POD = load('V_POD.mat', 'U','V','W');
u = reshape(V_POD.U(:,iFrame), size(X));
figure; imagesc(u(:,:,iSlice)'); axis equal; title('u pod');
colormap jet; caxis([0 1.2]); colorbar;

%%
V_AMIC = load('V_AMIC.mat', 'U','V','W');
u = reshape(V_AMIC.U(:,iFrame), size(X));
figure; imagesc(u(:,:,iSlice)'); axis equal; title('u amic');
colormap jet; caxis([0 1.2]); colorbar;

%% error
ROI = [1 88 1 88 1 12];
AFrame = 11:110;
iCount = 0;
U_ref = zeros(numel(X),length(AFrame));
V_ref = U_ref;                 W_ref = U_ref;
load('V_NOISE.mat', 'Um','Vm','Wm')
for iFrame = AFrame
    sub_name = ['Field_', num2str(iFrame,'%06u'), '.mat'];
    % loading field
    load([SourceFolder, sub_name], 'u', 'v', 'w');
    % croping
    u = u(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6));
    v = v(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6));
    w = w(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6));
    iCount = iCount + 1;
    U_ref(:,iCount) = u(:);
    V_ref(:,iCount) = v(:);
    W_ref(:,iCount) = w(:);
end
% V_error = @(U,V,W) std([U(:,AFrame)-U_ref;...
%     V(:,AFrame)-V_ref;W(:,AFrame)-W_ref], 1, 'all')/...
%     std([U_ref-Um;V_ref-Vm;W_ref-Wm], 1, 'all');
V_error = @(U,V,W) std([U(:,AFrame)-U_ref;...
    V(:,AFrame)-V_ref;W(:,AFrame)-W_ref], 1, 'all')/1;
fprintf(['velocity error from field with noise: ',...
    num2str(V_error(V_NOISE.U_noise,V_NOISE.V_noise,...
    V_NOISE.W_noise)),'\n']);
fprintf(['velocity error from SG filtered field: ',...
    num2str(V_error(V_SG.U,V_SG.V,V_SG.W)),'\n']);
fprintf(['velocity error from POD filtered field: ',...
    num2str(V_error(V_POD.U,V_POD.V,V_POD.W)),'\n']);
fprintf(['velocity error from AMIC filtered field: ',...
    num2str(V_error(V_AMIC.U,V_AMIC.V,V_AMIC.W)),'\n']);
% fprintf(['velocity error from AE filtered field: ',...
%     num2str(V_error(V_AE.U,V_AE.V,V_AE.W)),'\n']);
% fprintf(['velocity error from U-Net filtered field: ',...
%     num2str(V_error(V_UNET.U,V_UNET.V,V_UNET.W)),'\n']);

%% cosine similarity
V_simi = @(U,V,W) mean(sum([U(:,AFrame);V(:,AFrame);W(:,AFrame)].*...
    [U_ref;V_ref;W_ref])./...
    sqrt(sum([U(:,AFrame);V(:,AFrame);W(:,AFrame)].^2).*...
    sum([U_ref;V_ref;W_ref].^2)));
fprintf(['velocity similarity from field with noise to reference: ',...
    num2str(V_simi(V_NOISE.U_noise,V_NOISE.V_noise,...
    V_NOISE.W_noise)),'\n']);
fprintf(['velocity similarity from SG filtered to reference: ',...
    num2str(V_simi(V_SG.U,V_SG.V,V_SG.W)),'\n']);
fprintf(['velocity similarity from POD filtered to reference: ',...
    num2str(V_simi(V_POD.U,V_POD.V,V_POD.W)),'\n']);
fprintf(['velocity similarity from AMIC filtered to reference: ',...
    num2str(V_simi(V_AMIC.U,V_AMIC.V,V_AMIC.W)),'\n']);
% fprintf(['velocity similarity from AE filtered to reference: ',...
%     num2str(V_simi(V_AE.U,V_AE.V,V_AE.W)),'\n']);
% fprintf(['velocity similarity from U-Net filtered to reference: ',...
%     num2str(V_simi(V_UNET.U,V_UNET.V,V_UNET.W)),'\n']);