function plot_graphs(mat_file)
%%

% Load data from file.
data = load(mat_file);
e = data.strain;

%% Find Poisson's Ratio
e_1 = e(:,1,:);
e_1 = e_1(:);
e_2 = e(:,2,:);
e_2 = e_2(:);
[nu, nu_int] = regress(e_1, [ones(size(e_2)) e_2]);
nu = nu(end);
nu_err = (nu_int(end) - nu_int(1))/2;

%% Remove repeats.
e_se = std(e, [], 3)/sqrt(size(e,3));
e = mean(e, 3);
e_rel = e_se./e;

%% Calculate strain concentrations.
Ke = bsxfun(@rdivide, e, e(:,2));
Ke_se = abs(Ke).*hypot(e_rel, e_rel(:,2));
Ke_se(:,2) = 0;

% Calculate a weighted mean.
mKe = sum(Ke.*Ke_se.^-2,1)./sum(Ke_se.^-2,1);
mKe(:,2) = mean(Ke(:,2));
mKe_se = sqrt(1./sum(Ke_se.^-2,1));

%% Tabulate Strain Concentrations
fun = @(x, se) sprintf('%.3f&%.3f', x, 1.96*se);
tab = cell(size(Ke)+1);
tab(1:end-1,1) = arrayfun(fun, e(:,7), e_se(:,7), 'UniformOutput', false);
tab{end,1} = 'mean';
tab(1:end-1,2:end) = arrayfun(fun, Ke, Ke_se, 'UniformOutput', false);
tab(end,2:end) = arrayfun(fun, mKe, mKe_se, 'UniformOutput', false);

names = cell(1,size(tab,2));
names{1} = 'e7';
names(2:end) = arrayfun(@(i){sprintf('Ke%d', i)}, 1:size(names,2)-1);
tab = cell2table(tab, 'VariableNames', names);
disp(tab);
return

%% Plot strains and trendlines
f = figure();
ax = axes(f);
hold(ax, 'on');
line_colour = lines(size(e_mean,2));

% Plot trendlines.
for i = 1:size(e_mean, 2)
    b = regress(e(:,i), [ones(size(e,2)), e]);
    y_hat = X*b;
    L = plot(ax, e(:,2), y_hat);
    L.Color = line_colour(i,:);
    L.Marker = 'none';
    L.LineStyle = '-';
    L.LineWidth = 1;
end

% Plot error bars.
for i = 1:size(e, 2)
    eb = errorbar(ax, e(:, 2), e(:, i), ...
                  1.96*e_se(:, i), 1.96*e_se(:, i),...
                  1.96*e_se(:, 2), 1.96*e_se(:, 2));
    eb.Color = 0.8 * [1 1 1];
    eb.Marker = 'x';
    eb.MarkerEdgeColor = 0.8 * line_colour(i, :);
    eb.LineStyle = 'none';
    eb.LineWidth = 0.8;
end

% Label axes.
names = cell(size(e_mean, 2));
for i = 1:size(e_mean, 2)
    names{i} = sprintf('i=%d', i);
end
xlabel(ax, 'epsilon 2');
ylabel(ax, 'epsilon i');
g = legend(ax, names{1:size(e_mean, 2)});
g.Location = 'northwest';
end