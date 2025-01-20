load('Pressure.mat');
iFrame = 24;
iSlice = 6;
dx = 0.0114;  xb = 0:dx:1;   yb = 0:dx:1;   zb = 0:dx:0.5;   zb = zb(1:12);
[X,Y,Z] = meshgrid(xb,yb,zb);

%%
p = reshape(P_ref(:,iFrame), size(X));
figure; imagesc(p(:,:,iSlice)'); axis equal; title('p');
colormap jet; caxis([-0.02 0.02]); colorbar;

%%
% p = reshape(P_vari(:,iFrame), size(X));
% figure; imagesc(p(:,:,iSlice)'); axis equal; title('p vari');
% colormap jet; caxis([-0.02 0.02]); colorbar;

%%
p = reshape(P_noise(:,iFrame), size(X));
figure; imagesc(p(:,:,iSlice)'); axis equal; title('p on noise');
colormap jet; caxis([-0.02 0.02]); colorbar;

%%
p = reshape(P_sg(:,iFrame), size(X));
figure; imagesc(p(:,:,iSlice)'); axis equal; title('p on sg');
colormap jet; caxis([-0.02 0.02]); colorbar;

%%
p = reshape(P_pod(:,iFrame), size(X));
figure; imagesc(p(:,:,iSlice)'); axis equal; title('p on pod');
colormap jet; caxis([-0.02 0.02]); colorbar;

%%
p = reshape(P_amic(:,iFrame), size(X));
figure; imagesc(p(:,:,iSlice)'); axis equal; title('p on amic');
colormap jet; caxis([-0.02 0.02]); colorbar;

%% error
% P_error = @(P) std(P-P_ref, 1, 'all')/std(P_ref, 1, 'all');
P_error = @(P) std(P-P_ref, 1, 'all')/(0.5*1*1);
fprintf(['pressure error from field with noise: ',...
    num2str(P_error(P_noise)),'\n']);
fprintf(['pressure error from SG filtered field: ',...
    num2str(P_error(P_sg)),'\n']);
fprintf(['pressure error from POD filtered field: ',...
    num2str(P_error(P_pod)),'\n']);
fprintf(['pressure error from AMIC filtered field: ',...
    num2str(P_error(P_amic)),'\n']);
% fprintf(['pressure error from AE filtered field: ',...
%     num2str(P_error(P_ae)),'\n']);

%% cosine similarity
P_simi = @(P) mean(sum(P.*P_ref)./sqrt(sum(P.^2).*sum(P_ref.^2)));
fprintf(['pressure similarity from field with noise to reference: ',...
    num2str(P_simi(P_noise)),'\n']);
fprintf(['pressure similarity from SG filtered field to reference: ',...
    num2str(P_simi(P_sg)),'\n']);
fprintf(['pressure similarity from POD filtered field to reference: ',...
    num2str(P_simi(P_pod)),'\n']);
fprintf(['pressure similarity from AMIC filtered field to reference: ',...
    num2str(P_simi(P_amic)),'\n']);
% fprintf(['pressure similarity from AE filtered field to reference: ',...
%     num2str(P_simi(P_ae)),'\n']);