function I = hole_conc(a,L,H,S,out,cmap,px,fname)
% HOLE_CONC - Contour plot of stress concentration of circular hole
%
% COPYRIGHT (C) 2017 Russell Maguire

if nargin < 8
    fname = 'report/img/_.png';
end
if nargin < 7
    px = 1000;
end
if nargin < 6
    cmap = @jet;
end
if nargin < 5
    out = 's2';
end
if nargin < 4
    S = 1;
end
if nargin < 3
    H = 200e-3;
end
if nargin < 2
    L = 300e-3;
end
if nargin < 1
    a = 25e-3;
end

% Coordinates
xv = linspace(-L/2, L/2, px);
yv = linspace(-H/2, H/2, round(length(xv)*H/L));
[x, y] = meshgrid(xv, yv);
[th, r] = cart2pol(x, y);

% Radial stress
srr =   0.5*S*(1 - a^2./r.^2)...
       +0.5*S*(1 - 4*a^2./r.^2 + 3*a^4./r.^4).*cos(2*th);
srr(r < a) = nan;
% Transverse stress
sthth = 0.5*S*(1 + a^2./r.^2)...
       -0.5*S*(1 + 3*a^4./r.^4).*cos(2*th);
sthth(r < a) = nan;
% Radial--Transverse shear stress
srth = -0.5*S*(1 + 2*a^2./r.^2 - 3*a^4./r.^4).*sin(2*th);
srth(r < a) = nan;

% X stress
sxx = (srr + sthth)/2 + (srr - sthth)/2.*sin(2*th) + srth.*cos(2*th);
% Y stress
syy = (srr + sthth)/2 - (srr - sthth)/2.*sin(2*th) - srth.*cos(2*th);
% XY shear stress
sxy = -(srr - sthth)/2.*sin(2*th) + srth.*cos(2*th);

% Principle stresses
sa = (srr + sthth)/2 - sqrt((srr - sthth).^2/4 + srth.^2);
sb = (srr + sthth)/2 + sqrt((srr - sthth).^2/4 + srth.^2);
s2 = arrayfun(@max, sa, sb);
s1 = arrayfun(@min, sa, sb);

% Principle shear stress
t = s2 - s1;

I = eval(out);
imshow(I, [], 'Colormap', cmap());
colorbar();

if ~isempty(fname)
    % Normalise image
    img = I - min(I(:));
    img = 1 + ceil(255 * img / max(img(:)));
    img = ind2rgb(img, cmap(256));

    % Write to file
    imwrite(img, fname, 'png');
end
end