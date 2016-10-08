### Interval arithmetic class for matlab
A simple method to perform numerical calculation on regions (bounds; [intervals](https://en.wikipedia.org/wiki/Interval_arithmetic)) instead on values.
This concept is suitable for a variaty of problems, the most common one - handle uncertainties.

####Example #1 - general uncertainty calculation:
```matlab
% under assumption that the number of paint cans to paint
% the living room is calculated using this equatlon:
% paint cans = 2 * (room width + room length) * room height / (paint per can * paint efficiency)
% lets measure the room and calculate how many cans of paint we need:
room_width = interval(5.4, 6.6); % [m] (bounded interval of room width)
room_length = interval(3.6, 4.4);% [m] (bounded interval of room length)
room_height = interval(2.98, 3.1); % [m] bounded interval of room height)
paint_efficiency = interval(0.7, 1.3); % efficiency changes between painters (amature vs. professional)
paint_per_can = interval(39.8, 40.2); % [liter] (bounded region specifying amount of paint in a can)
paint_cans = 2 * (room_width + room_length) * room_height / (paint_efficiency * paint_per_can);
% > we need between 1.02 to 2.448 cans of paint, which means we should buy 3 cans!
```

####Example #2 - set calculations with intervals:
```matlab
% intervals
a = interval(-3, 5);
b = interval(-5, 3);
c = interval(0, 1);
d = interval(-8, -6);

% set operations
a_and_b = a & b; % [-3, 3]
a_or_b = a | b; % [-5, 5]
a_xor_b = xor(a, b); % [-5, 5]
is_c_in_a = in(c, a); % true
b_diff_a = diff(b, a); % [-5, -3]

% various distances between sets
hd = hausdorffDist(a, b); % hausedorff distance = 2
id = innerDist(d, a); % inner distance = 3
sd = signedDist(d, b); % signed distance = -1
```

####Example #3 - array uncertainty calculation (operations on vectors in which an element is a region rather then a value):
```matlab
% a measurement, bounded by uncertainty.
% each measurement is given in a resoultion or with uncertainty of 0.1
a = interval(0 : 0.3 : pi, 0.1 : 0.3 : pi+0.1);

% calculation based upon uncertainty
b = sin(a.^2) - 2 * sinh(0.5*sqrt(a));

% visualize
figure;
subplot(2,1,1);
plot(a.infimum, 'b'); hold on; plot(a.supremum, 'r');
grid on;
legend('measurement lower bound', 'measurement upper bound');
ylabel('[rad]');
title('measurement (bounded by uncertainty)');
subplot(2,1,2);
plot(b.infimum, 'b'); hold on; plot(b.supremum, 'r');
grid on;
legend('function lower bound', 'function upper bound');
title('function operated on measurement (f = sin(x^2) - 2 * sinh(sqrt(x) / 2))');
```

####example usage #4 - interval matrix (operations on matrix in which an element is a region rather then a value):
```matlab
% define an interval matrix
%     [[8, 9], [1, 2],  [6, 7]
% A =  [3, 4], [5, 6],  [7, 8]
%      [4, 5], [9, 10], [2, 3]]
A = interval(magic(3), magic(3)+1);

% what is it's inverse?
Ainv = inv(A);
%           [[0.0743,   0.2161], [-0.2556, -0.0373],  [-0.0215,  0.1453]
% inv(A) =   [-0.1006, -0.0256], [-0.0375,  0.0779],  [0.0595,   0.1476]
%            [-0.0843,  0.0414], [0.0902,   0.2835],  [-0.1787, -0.0309]]

% what is the determinant region?
Adet = det(A); % [-852.125, -142.5098]

% it's upper triangular decomposition?
[l, u, p] = lu(A);
%     [[8, 9], [1, 2],      [6, 7]
% u =  [0, 0], [4, 5.6667], [3.5, 6]
%      [0, 0], [0, 0],      [-16.7083, -4.4534]]

% if A is bounded region of linear system coefficients,
% what is the region of the solution A*x = [1; 1; 1] ?
b = [1; 1; 1];
x = linsolve(A, b);
%    [[-0.1438, 0.1271]
% x = [-0.0632, 0.1183]
%     [0.0077,  0.2734]]
```


####Example #5 - root of a function defined upon an interval:
```matlab
% lets find the zero of function: 'f' in interval 'a'
f = @(x) x.^3 - 2*x - 5;
a = interval(-20, 20);

% define solver parameters
opt.maxIter = 20;
opt.tolX = 0.1;
opt.tolFun = 0.1;

% the function zero inside the given inteval:
x = fzero(f, a, [], opt);
```

####Example #6 - root of an interval function with a supplied gradient:
```matlab
% lets find the zero of function: 'f' with gradient 'g', in interval 'a'
f = @(x) cos(x);
g = @(x) -sin(x);
a = interval(-9, -5);

% define solver parameters
opt.maxIter = 20;
opt.tolX = 0.1;
opt.tolFun = 0.1;

% the function zero inside the given inteval:
x = fzero(f, a, g, opt);
```
