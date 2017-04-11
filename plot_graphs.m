function plot_graphs(hole_mat, fillet_mat)
% PLOT_GRAPHS - plot graphs and process data for report

% Value/error format.
err2latex = @(x, se) sprintf('%.3f\\pm%.3f', x, se);

% Load data from file.
data = load(hole_mat);
e = data.strain;

data = load(fillet_mat);
s = data.stress';
r_d = data.normalised_radius;
H_d = data.normalised_height;
d_r = 1./r_d;
d_H = 1./H_d;

%% 1a. Find Poisson's Ratio
e_1 = e(:,1,:);
e_1 = e_1(:);
e_2 = e(:,2,:);
e_2 = e_2(:);
[nu, nu_int, ~, ] = regress(e_1, [ones(size(e_2)) e_2]);
nu = -nu(end);
nu_se = (nu_int(end) - nu_int(1))/(2*1.96);
t = (nu - 0.33)/nu_se;
p = tcdf(t, length(e_1)-1);

%% 1b. Calculate strain concentrations.
% Handle repeats.
e_se = std(e, [], 3)/sqrt(size(e,3));
e = mean(e, 3);
e_rel = e_se./e;

% Calculate strain concentrations.
Ke = bsxfun(@rdivide, e, e(:,2));
Ke_se = abs(Ke).*hypot(e_rel, e_rel(:,2));
Ke_se(:,2) = 0;

% Calculate a weighted mean.
mKe = sum(Ke.*Ke_se.^-2,1)./sum(Ke_se.^-2,1);
mKe(:,2) = mean(Ke(:,2));
mKe_se = sqrt(1./sum(Ke_se.^-2,1));

%% 2a. Calculate Stress Concentration
Ks = s/100;
for i = 1:length(H_d)
    y = log(Ks(i,~isnan(Ks(i,:))))';
    X = log(r_d(~isnan(Ks(i,:))))';
    X = [ones(size(X)) X];
    [b,~,~,~,stats] = regress(y, X);
    R2 = stats(1);
end

%% 3a. Display Poisson's Ratio
tab = {err2latex(nu, nu_se) p};
tab = cell2table(tab, 'Variablenames', {'nu', 'p'});
disp(tab);

%% 3b. Tabulate Strain Concentrations
tab = cell(size(Ke)+1);
tab(1:end-1,1) = arrayfun(err2latex, e(:,7), e_se(:,7), 'UniformOutput', false);
tab{end,1} = 'mean';
tab(1:end-1,2:end) = arrayfun(err2latex, Ke, Ke_se, 'UniformOutput', false);
tab(end,2:end) = arrayfun(err2latex, mKe, mKe_se, 'UniformOutput', false);

names = cell(1,size(tab,2));
names{1} = 'e7';
names(2:end) = arrayfun(@(i){sprintf('Ke%d', i)}, 1:size(names,2)-1);
tab = cell2table(tab, 'VariableNames', names);
disp(tab);

%% 3c. Plot Strains and Trendlines
sz = [250 500];
sc = 0.7;
f = figure();
f.Position(3:4) = sz/sc;
ax = axes(f);
ax.Position = [flip(sz)/max(sz)*sc*(1-sc) sc sc];
hold(ax, 'on');
line_colour = lines(size(e,2));

% Make sure legend is in the right order
[~,order] = sort(e(end,:), 'descend');

% Plot trendlines.
for i = order
    b = regress(e(:,i), [ones(size(e(:,2))) e(:,2)]);
    L = plot(ax, e(:,2), [ones(size(e(:,2))) e(:,2)]*b);
    L.Color = line_colour(i,:);
    L.Marker = 'none';
    L.LineStyle = '-';
    L.LineWidth = 1;
end

% Plot error bars.
for i = order
    eb = errorbar(ax, e(:,2), e(:,i), 3*e_se(:,i), 3*e_se(:,i), 3*e_se(:,2), 3*e_se(:,2));
    eb.Color = 0.8 * line_colour(i, :);
    eb.Marker = 'x';
    eb.MarkerSize = 4;
    eb.MarkerEdgeColor = 0.8 * line_colour(i, :);
    eb.LineStyle = 'none';
    eb.LineWidth = 0.8;
end

% Label axes.
names = cell(size(order));
for i = 1:length(order)
    names{i} = sprintf('e%d', order(i));
end
xlabel(ax, 'e2');
ylabel(ax, 'ei');
g = legend(ax, [names names]); % Produce legends for points and trendlines
g.Location = 'northwest';
end