function plot_graphs(hole_mat, fillet_mat)
% PLOT_GRAPHS - plot graphs and process data for report

% Value/error format.
err2latex = @(x, e) sprintf('%.4f\\pm%.4f', x, e);

% Load data from file.
data = load(hole_mat);
e = data.strain;
e_all = zeros(size(e,2), size(e,1)*size(e,3));
for i = 1:size(data.strain,2)
    e_all(i,:) = reshape(e(:,i,:), 1, []);
end

data = load(fillet_mat);
s = data.stress';
rd = data.normalised_radius;
Hd = data.normalised_height;

%% 1a. Calculate strain concentrations by fitting gradient.
Ke = zeros(1, size(e_all,1));
Ke_e = zeros(1, size(e_all,1));
b = zeros(size(e_all,1), 2);
X = [e_all(2,:)'];
for i = 1:size(e_all,1)
    [b(i,:), b_ci] = regress(e_all(i,:)', X);
    Ke(i) = b(i,end);
    Ke_e(i) = (b_ci(end,end) - b_ci(end,1))/2;
end

%% 1b. Handle repeats
e_e = tpdf(0.975,size(e,3))*std(e, [], 3)/sqrt(size(e,3));
e_e = shiftdim(e_e, 1);
e = mean(e, 3);
e = shiftdim(e, 1);

%% 2a. Tabulate Strain Concentrations
tab = arrayfun(err2latex, Ke, Ke_e, 'UniformOutput', false);
names = arrayfun(@(i){sprintf('Ke%d', i)}, 1:length(Ke));
tab = cell2table(tab, 'VariableNames', names);
disp(tab);

%% 2b. Plot Strains and Trendlines
sz = [250 500];
sc = 0.7;
f = figure();
f.Position(3:4) = sz/sc;
ax = axes(f);
ax.Position = [sqrt(sc)*(1-sc) sqrt(sc)*(1-sc) sc sc];
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
hold(ax, 'on');
line_colour = lines(size(e,1));

% Make sure legend is in the right order
[~,order] = sort(e(:,end)', 'descend');

% Plot data points.
for i = order
    P = scatter(ax, e(2,:), e(i,:));
    P.Marker = 'x';
    P.MarkerEdgeColor = 0.8*line_colour(i,:);
    P.LineWidth = 1;
end

% Plot trendlines.
for i = order
    L = plot(ax, e(2,[1 end]), [ones(1,2)' e(2,[1 end])']*b(i,:)');
    L.Color = line_colour(i,:);
    L.Marker = 'none';
    L.LineStyle = '-';
    L.LineWidth = 1;
end

% Label axes.
names = cell(size(order));
for i = 1:length(order)
    names{i} = sprintf('e%d', order(i));
end
xlabel(ax, 'e2');
ylabel(ax, 'ei');
g = legend(ax, names); % Produce legends for points and trendlines
g.Location = 'northwest';

%% 3a. Calculate Stress Concentration
Ks = s/100;

%% 3b. Fit Ks = b_0(r/d)^b_1
b = zeros(2,length(Hd));
b_e = zeros(2,length(Hd));
R2 = zeros(1,length(Hd));
p = zeros(1,length(Hd));
for i = 1:length(Hd)
    y = log10(Ks(i,~isnan(Ks(i,:))))';
    X = log10(rd(~isnan(Ks(i,:))))';
    X = [ones(size(X)) X];
    [b(:,i),b_ci,~,~,stats] = regress(y, X);
    b_e(:,i) = (b_ci(:,2) - b_ci(:,1))/2;
    R2(i) = stats(1);
    p(i) = stats(end);
end
b(1,:) = 10.^b(1,:);
b_e(1,:) = log(10)*b(1,:).*b_e(1,:);

%% 4a. Tabulate stress concentration plots
tab = cell(length(Hd),5);
for i = 1:length(Hd)
    tab{i,1} = sprintf('%.1f', Hd(i));
    tab{i,2} = err2latex(b(1,i), b_e(1,i));
    tab{i,3} = err2latex(b(2,i), b_e(2,i));
    tab{i,4} = sprintf('%.5f', R2(i));
    tab{i,5} = sprintf('%.3g', p(i));
end

names = {'H_d', 'a', 'b', 'R2', 'p'};
tab = cell2table(tab, 'Variablenames', names);
disp(tab);


%% 4b. Plot stress concentration against r/d on a log-log scale.
sz = [250 250];
sc = 0.7;
f = figure();
f.Position(3:4) = sz/sc;
ax = axes(f);
ax.Position = [sqrt(sc)*(1-sc) sqrt(sc)*(1-sc) sc sc];
ax.XScale = 'log';
ax.YScale = 'log';
hold(ax, 'on');
line_colour = lines(length(Hd));

% Make sure legend is in the right order
[~,order] = sort(Ks(:,1), 'descend');
order = order';

% Plot datapoints.
for i = order
    P = scatter(ax, rd, Ks(i,:));
    P.Marker = 'x';
    P.MarkerEdgeColor = 0.8 * line_colour(i, :);
    P.LineWidth = 1;
end

% Plot trendlines.
for i = order
    x = rd(~isnan(Ks(i,:)));
    x = [x(1) x(end)];
    L = plot(ax, x, b(1,i)*x.^b(2,i));
    L.Color = line_colour(i,:);
    L.Marker = 'none';
    L.LineStyle = '-';
    L.LineWidth = 1;
end

% Label axes.
names = cell(size(order));
for i = 1:length(order)
    names{i} = sprintf('H/d=%.1f', Hd(order(i)));
end
xlabel(ax, 'r/d');
ylabel(ax, 'Ks');
g = legend(ax, names);
g.Location = 'northeast';


%% 4b. Plot stress concentration against r/d on a linear scale.
sz = [250 250];
sc = 0.7;
f = figure();
f.Position(3:4) = sz/sc;
ax = axes(f);
ax.Position = [sqrt(sc)*(1-sc) sqrt(sc)*(1-sc) sc sc];
hold(ax, 'on');
line_colour = lines(length(Hd));

% Make sure legend is in the right order
[~,order] = sort(Ks(:,1), 'descend');
order = order';

% Plot datapoints.
for i = order
    P = scatter(ax, rd, Ks(i,:));
    P.Marker = 'x';
    P.MarkerEdgeColor = 0.8 * line_colour(i, :);
    P.LineWidth = 1;
end

% Plot trendlines.
lim = 0.5*(Hd-1);
for i = order
    x = linspace(0, lim(i), round(100*lim(i)/max(lim)));
    L = plot(ax, x, b(1,i)*x.^b(2,i));
    L.Color = line_colour(i,:);
    L.Marker = 'none';
    L.LineStyle = '-';
    L.LineWidth = 1;
end

% Label axes.
names = cell(size(order));
for i = 1:length(order)
    names{i} = sprintf('H/d=%.1f', Hd(order(i)));
end
xlabel(ax, 'r/d');
ylabel(ax, 'Ks');
g = legend(ax, names);
g.Location = 'northeast';
end

