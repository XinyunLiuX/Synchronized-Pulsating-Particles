folder = 'RES';
filename = dir(, ".mat");
[up, lo] = envelope(x(:,1), 10, 'peak');
secondPart = up-lo;
secondPart = secondPart(floor(length(secondPart)/2):end);
plot(secondPart)
mean(secondPart)
std(secondPart)
hold off