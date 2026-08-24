figure

global xi

normalized_xi = xi/xi(end-1);

plot(normalized_xi, u(end, :), '-')

xlabel('Normalized x')
ylabel('Velocity (m a^{-1})')

xlim([0 1])
ylim([0 max(max(u)) + 20])

set(gcf, 'Position', [465 220 511 260])