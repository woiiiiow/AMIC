% Woii 240619
% for additional wing foil LES
% show result

iFrame = 42;

% loading grid
FileGrid = '../../../airfoilLES/airfoilLES_grid.h5';
% x0 = h5read(FileGrid, '/x');
% y0 = h5read(FileGrid, '/y');
global xa ya
xa = h5read(FileGrid, '/xa'); % coordinates of the airfoil
ya = h5read(FileGrid, '/ya');
load('V_NOISE.mat', 'xb','yb');
[X, Y] = meshgrid(xb, yb);

%% velocity
load('V_NOISE.mat');
field = reshape(U_clean(:,iFrame), size(X));
showimage(X, Y, field);
caxis([0.4 1.4]); title('U REF');
field = reshape(V_clean(:,iFrame), size(X));
showimage(X, Y, field);
caxis([-0.5 0.5]); title('V REF');

field = reshape(U_noisy(:,iFrame), size(X));
showimage(X, Y, field);
caxis([0.4 1.4]); title('U NOISY');
field = reshape(V_noisy(:,iFrame), size(X));
showimage(X, Y, field);
caxis([-0.5 0.5]); title('V NOISY');

U_AMIC = load('V_AMIC.mat');

field = reshape(U_AMIC.U(:,iFrame), size(X));
showimage(X, Y, field);
caxis([0.4 1.4]); title('U AMIC');
field = reshape(U_AMIC.V(:,iFrame), size(X));
showimage(X, Y, field);
caxis([-0.5 0.5]); title('V AMIC');




%% functions
function showimage(X, Y, field)
global xa ya
figure; hold on;
pcolor(X, Y, field); shading flat; colormap jet;
fill(xa, ya, 'k');
axis equal; box on;
xlim([-0.1 3]); ylim([-0.5 0.5]);
set(gca, 'fontsize', 24);
set(gcf, 'position', [10 10 1080 480]);
end