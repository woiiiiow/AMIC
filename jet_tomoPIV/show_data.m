load('JetAR1_02_conAMIC_5.mat');

%%
% iSlice = 17;
% iFrame = 999;
% 
% v = reshape(V_array(:, iFrame), size(X));
% figure; pcolor(X(:,:,iSlice), Y(:,:,iSlice), v(:,:,iSlice));
% shading flat; colorbar('northoutside'); colormap(jet); axis equal;
% clim([-0.06 0.06]);
% 
% U_mean = reshape(mean(U_array, 2), size(X));
% figure; pcolor(X(:,:,iSlice), Y(:,:,iSlice), U_mean(:,:,iSlice));
% shading flat; colorbar('northoutside'); colormap(jet); axis equal;
% clim([-0.02 0.16]);
% 
% % U_p = U_array - mean(U_array, 2);
% % U_pp = reshape(mean(U_p.*U_p, 2), size(X));
% % figure; pcolor(X(:,:,iSlice), Y(:,:,iSlice), U_pp(:,:,iSlice));
% % shading flat; colorbar('northoutside'); colormap(jet); axis equal;
% % % clim([-0.02 0.16]);
% 
% V_p  = V_array - mean(V_array, 2);
% V_pp = reshape(mean(V_p.*V_p, 2), size(X));
% figure; pcolor(X(:,:,iSlice), Y(:,:,iSlice), V_pp(:,:,iSlice));
% shading flat; colorbar('northoutside'); colormap(jet); axis equal;
% clim([0 1e-3]);
% 
% P_p  = P_array - mean(P_array, 2);
% P_pp = reshape(mean(P_p.*P_p, 2), size(X));
% figure; pcolor(X(:,:,iSlice), Y(:,:,iSlice), P_pp(:,:,iSlice));
% shading flat; colorbar('northoutside'); colormap(jet); axis equal;
% clim([0 2.5]);

%%
iSlice = 17;
iFrame = 999;
figure; set(gcf,'outerposition',get(0,'screensize'));

v = reshape(V_array(:, iFrame), size(X));
subplot(2, 2, 1);
pcolor(X(:,:,iSlice), Y(:,:,iSlice), v(:,:,iSlice));
shading flat; colorbar('northoutside'); colormap(jet); axis equal;
clim([-0.06 0.06]);

U_mean = reshape(mean(U_array, 2), size(X));
subplot(2, 2, 3);
pcolor(X(:,:,iSlice), Y(:,:,iSlice), U_mean(:,:,iSlice));
shading flat; colorbar('northoutside'); colormap(jet); axis equal;
clim([-0.02 0.16]);

% U_p = U_array - mean(U_array, 2);
% U_pp = reshape(mean(U_p.*U_p, 2), size(X));
% figure; pcolor(X(:,:,iSlice), Y(:,:,iSlice), U_pp(:,:,iSlice));
% shading flat; colorbar('northoutside'); colormap(jet); axis equal;
% % clim([-0.02 0.16]);

V_p  = V_array - mean(V_array, 2);
V_pp = reshape(mean(V_p.*V_p, 2), size(X));
subplot(2, 2, 2);
pcolor(X(:,:,iSlice), Y(:,:,iSlice), V_pp(:,:,iSlice));
shading flat; colorbar('northoutside'); colormap(jet); axis equal;
clim([0 1e-3]);

P_p  = P_array - mean(P_array, 2);
P_pp = reshape(mean(P_p.*P_p, 2), size(X));
subplot(2, 2, 4);
pcolor(X(:,:,iSlice), Y(:,:,iSlice), P_pp(:,:,iSlice));
shading flat; colorbar('northoutside'); colormap(jet); axis equal;
clim([0 2.5]);