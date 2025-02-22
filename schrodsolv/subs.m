
figure;

% First subplot
subplot(2, 2, 1); % 2x2 grid, first subplot
plot(fort((1*a:2*a),1),fort((1*a:2*a),2)')

pbaspect([1 1 1])

xlim([-100 100])
ylim([0 0.15])
title('T=30k');

% Second subplot
subplot(2, 2, 2); % 2x2 grid, second subplot
plot(fort((5*a:6*a),1),fort((6*a:7*a),2)')

pbaspect([1 1 1])

xlim([-100 100])
ylim([0 0.15])
title('T=24k');

% Third subplot
subplot(2, 2, 3); % 2x2 grid, third subplot
plot(fort((11*a:12*a),1),fort((11*a:12*a),2)')
pbaspect([1 1 1])

xlim([-100 100])
ylim([0 0.15])
title('T=14k');

% Fourth subplot
subplot(2, 2, 4); % 2x2 grid, fourth subplot
plot(fort((14*a:15*a),1),fort((14*a:15*a),2)')
pbaspect([1 1 1])

xlim([-100 100])
ylim([0 0.15])
title('T=4k');


