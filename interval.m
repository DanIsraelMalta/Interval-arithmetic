%{
% Interval Arithmetic Class
%
% properties:
% --------------------
% - infimum - interval lower bound
% - supremum - interval upper bound
%
% construction methods:
% --------------------------------------
% - full interval construction:
%   >> a = interval(); % a.infimum= -Inf, a.supremum= Inf
% - scalar intervals:
%   >> a = interval(1); % a.infimum = 1, a.supremum = 1
%   >> a = interval(1, 2); % a.infimum = 1, a.supremum = 2
% - vector interval:
%   >> a = interval([1, 2]); % a.infimum = [1, 2], a.supremum = [1, 2]
%   >> a = interval({1, 2}); % a.infimum = [1, 2], a.supremum = [1, 2]
%   >> a = interval([1, 2], [2, 3]); % a.infimum = [1, 2], a.supremum = [2, 3]
%   >> a = interval(1:5, 2:5); % a.infimum = [1, 2, 3, 4, 5], a.supremum = [2, 3, 4, 5, 6]
%   >> a=linspace(interval(1,2), interval(3,6), 3); % a.infimum = [1, 2, 3], a.supremum = [2, 4, 6]
%   >> a=logspace(interval(1,2), interval(3,6), 3); % a.infimum = [10, 100, 100], a.supremum = [100, 10000, 1e6]
% - matrix interval:
%   >>  a = interval(magic(3), magic(3) + 1);
%          % a.infimum = [8, 1, 6; 3, 5, 7; 4, 9, 2]
%          % a.supremum = [9, 2, 7; 4, 6, 8; 5, 10, 3]
%   >> a = meshgrid(interval(1:3), interval(2:4))
%         % a.infimum = [1, 2, 3; 1, 2, 3;1, 2, 3]
%         % a.supremum = [1, 2, 3; 1, 2, 3;1, 2, 3]
% - notes:
%   1) a NaN shall be transformed into an 'empty' interva - [Inf, -Inf].
%   2) a reversed input argument (one whos infimum is larger then its
%       supremum shall be transformed into an 'empty' interval - [-Inf, Inf].
%
% overlaoded operators:
% ---------------------------------------
% +, -, *, .*, /, ./, .^, .', >, >=, <, <=, ==, ~=, and, or, xor
%
% algebric methods:
% -------------------------------
% min, max, abs, sign, round, ceil, floor, fix, sqrt, log, log2,
% log10, exp, pow2, sum, dot, prod, hypot, polyval
%
% trigonometric ahd hyperbolic methods:
% --------------------------------------------------------------------
% sin, cos, tan, atan, atan2, sinh, cosh, tanh
%
% matrix relates methods:
% -----------------------------------------
% transpose, diag, inv, lu, det, linsolve
%
% visualization methods:
% ----------------------------------------
% plot
%
% optimization methods:
% --------------------------------------
% fzero
%
% special methods:
% ----------------------------
% - middle - return the middle point of an interval
% - width - return the width of an interval
% - radius - return the radius of an interval
% - mignitude - return the mignitude (minimum of absolute of infimum/supremum pairs of interval)
% - magnitude - return the maximum of absolute value of infimum/supremum pairs of interval
% - isempty - return true if inteval is empty (Inf, -Inf)
% - issingle - return true if an interval is actualy 1D, e.g - infimum is also the supremum
% - ismember - return true if a given value is within a given interval
% - next - increase an interval boundary (lower and upper) to the next value
% - in - return true if one interval is contained within another interval
% - intersect - return an interval holding the intersection of two intervals
% - union - return an interval holding the union of two intervals
% - diff - return the relative complement of first interval in second interval
% - hausdorffDist - return the hausdorff distance between two intervals
% - innerDist - return the inner distance between two intervals
% - signedDist - return the signed distance between two intervals
% 
%
% remarks:
% --------------
% 1) mldivide (\) and mrdivide (/) numerical accuracy is correct only
%    for 'not to irregular' matrix (matrix with minimal condition number).
% 2) power operator can be performed for the following cases:
%     a. base interval is positive, and required power is positive.
%          in  this case the power can be an interval as well.
%     b. base interval include a zero or is negative, and the required
%          power is a positive integer (not interval).
% 3) matrix power is not supported.
% 4) complex numbers are not supported.
% 5) if an euclidean norm is needed, use the 'hypot' method.
% 6) fzero method uses a recrusive newton method, so:
%     a. its run time might be high for complex problems.
%     b. it might not converge if the zeros are to close to each other.
%
%
% example usage #1 - general uncertainty calculation:
% -------------------------------------------------------------------------------------------

% under assumption that the number of paint cans to paint
% the living room is calculated using this equatlon:
% paint cans = 2 * (room width + room length) * room height / (paint per can * paint efficiency)
% lets measure the room and calculate how many cans of paint we need:
room_width = interval(5.4, 6.6); % [m]
room_length = interval(3.6, 4.4);% [m]
room_height = interval(2.98, 3.1); % [m]
paint_efficiency = interval(0.7, 1.3); % efficiency changes between painters
paint_per_can = interval(39.8, 40.2); % [liter]
paint_cans = 2 * (room_width + room_length) * room_height / (paint_efficiency * paint_per_can);
% > we need between 1.02 to 2.448 cans of paint

% example usage #2 - set calculations with intervals:
% -----------------------------------------------------------------------------------------

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

% example usage #3 - numeric uncertainty calculation:
% -------------------------------------------------------------------------------------------

% a measurement, bounded by uncertainty
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

% example usage #4 - root of an interval function:
% ----------------------------------------------------------------------------------

% lets find the zero of function: 'f' in interval 'a'
f = @(x) x.^3 - 2*x - 5;
a = interval(-20, 20);

% define solver parameters
opt.maxIter = 20;
opt.tolX = 0.1;
opt.tolFun = 0.1;

% the function zero inside the given inteval:
x = fzero(f, a, [], opt)

%
% example usage #5 - root of an interval function with a supplied gradient:
% ------------------------------------------------------------------------------------------------------------------------------

% lets find the zero of function: 'f' with gradient 'g', in interval 'a'
f = @(x) cos(x);
g = @(x) -sin(x);
a = interval(-9, -5);

% define solver parameters
opt.maxIter = 20;
opt.tolX = 0.1;
opt.tolFun = 0.1;

% the function zero inside the given inteval:
x = fzero(f, a, g, opt)

%
% Dan I. Malta (malta.dan@gmail.com)
%}
classdef interval < handle
    % properties
    properties (SetAccess = protected)
        infimum = [];       % infimum boundery
        supremum = [];       % supremum boundery
    end % properties

    % constructors and core methods
    methods

        % constructor
        function xo_interval = interval(xi_inf, xi_sup)
            switch nargin
                case 0 % "empty" interval
                    xo_interval.infimum = -Inf;
                    xo_interval.supremum = Inf;
                case 1
                    if isa(xi_inf, 'interval')       % input is interval - then clone it
                        xo_interval.infimum = xi_inf.infimum;
                        xo_interval.supremum = xi_inf.supremum;
                    elseif iscell(xi_inf)   % input is cell of structure {{val},... , {vals}} or {{inf1, sup1}, ..., {infn, supn}}
                        for i = 1 : max(size(xi_inf))
                            if isa(xi_inf{i}, 'double')
                                temp = xi_inf{i};
                                if isa(temp, 'double')
                                    if isnan(temp)
                                        xo_interval.infimum(i) = -Inf;
                                        xo_interval.supremum(i) = Inf;
                                    else
                                        xo_interval.infimum(i) = temp;
                                        xo_interval.supremum(i) = temp;
                                    end
                                else
                                    error(INTERVAL_CTOR, 'interval constructor whos argument is a cell can not hold non numerical values.');
                                end
                            elseif iscell(xi_inf{i})
                                temp = cell2mat(xi_inf{i});
                                if numel(temp) == 1
                                    if isa(temp, 'double')
                                        if isnan(temp)
                                            xo_interval.infimum(i) = -Inf;
                                            xo_interval.supremum(i) = Inf;
                                        else
                                            xo_interval.infimum(i) = temp;
                                            xo_interval.supremum(i) = temp;
                                        end
                                    else
                                        error(INTERVAL_CTOR, 'interval constructor whos argument is a cell can not hold non numerical values.');
                                    end
                                elseif numel(temp) == 2
                                    if temp(1) > temp(2)
                                        % transpose to "empty set"
                                        xo_interval.infimum(i) = Inf;
                                        xo_interval.supremum(i) = -Inf;
                                    end
                                    if isa(temp(1), 'double')
                                        if isnan(temp(1))
                                            xo_interval.infimum(i) = -Inf;
                                        else
                                            xo_interval.infimum(i) = temp(1);
                                        end
                                    end
                                    if isa(temp(2), 'double')
                                        if isnan(temp(2))
                                            xo_interval.supremum(i) = Inf;
                                        else
                                            xo_interval.supremum(i) = temp(2);
                                        end
                                    end
                                else
                                    error(INTERVAL_CTOR, 'interval constructor input argument of type cell shoule be of the following form: {{inf, sup}, ..., {inf, sup}} or {{value}, ..., {value}}.');
                                end
                            end
                        end
                    elseif isa(xi_inf, 'double')   % input is double
                        [row, col] = size(xi_inf);
                        for j = 1 : row
                            for i = 1 : col
                                if isnan(xi_inf(j, i)) % NaN are transformed to empty intervals
                                    xo_interval.infimum(j, i) = -Inf;
                                    xo_interval.supremum(j, i) = Inf;
                                else
                                    xo_interval.infimum(j, i) = xi_inf(j, i);
                                    xo_interval.supremum(j, i) = xi_inf(j, i);
                                end
                            end % column iteration
                        end % row iteration
                    else    % do not handle other type of inputs
                        error(INTERVAL_CTOR, 'interval constructor input argument can only be of type ''double'', ''cell'' or ''interval''.');
                    end
                case 2 % interval given as a pair of infimum & supremum
                    if isa(xi_inf, 'double') && isa(xi_sup, 'double')
                        if numel(xi_inf) ~= numel(xi_sup)
                            error(INTERVAL_CTOR, 'if interval constructor accepts two input arguments, they must be of equal length.');
                        else
                            [row, col] = size(xi_inf);
                            for j = 1 : row
                                for i = 1 : col
                                    if xi_inf(j, i) > xi_sup(j, i)
                                        % transpose to "empty set"
                                        xi_inf(j, i) = Inf;
                                        xi_sup(j, i) = -Inf;
                                    end
                                    if isnan(xi_inf(j, i)) % NaN are transformed to empty intervals
                                        xo_interval.infimum(j, i) = -Inf;
                                    else
                                        xo_interval.infimum(j, i) = xi_inf(j, i);
                                    end
                                    if isnan(xi_sup(j, i)) % NaN are transformed to empty intervals
                                        xo_interval.supremum(j, i) = Inf;
                                    else
                                        xo_interval.supremum(j, i) = xi_sup(j, i);
                                    end
                                end % column iteration
                            end % row iteration
                        end
                    else
                        error(INTERVAL_CTOR, 'if interval constructor accepts two input arguments, they must be of type double.');
                    end
                otherwise
                    error(INTERVAL_CTOR, 'interval constructor can maximum number of input arguments is two.');
            end % switch
        end

        % linearly spaced intervals constructor
        function xo_c = linspace(xi_min, xi_max, xi_len)
            % number of argument check
            if nargin ~= 3
                error(INTERVAL_LINSPACE, 'interval::linspace - method accepts three input argument.');
            end
            
            % linspace
            if isa(xi_min, 'interval') && ~isa(xi_max, 'interval')
                xi_max = interval(xi_max);
                xo_c = interval(linspace(xi_min.infimum, xi_max.infimum, xi_len),...
                    linspace(xi_min.supremum, xi_max.supremum, xi_len));
            elseif isa(xi_max, 'interval') && ~isa(xi_min, 'interval')
                xi_min = interval(xi_min);
                xo_c = interval(linspace(xi_min.infimum, xi_max.infimum, xi_len),...
                    linspace(xi_min.supremum, xi_max.supremum, xi_len));
            elseif isa(xi_max, 'interval') && isa(xi_min, 'interval')
                xo_c = interval(linspace(xi_min.infimum, xi_max.infimum, xi_len),...
                    linspace(xi_min.supremum, xi_max.supremum, xi_len));
            else
                xo_c = interval();
            end
        end

        % logarithmically spaced intervals constructor
        function xo_c = logspace(xi_min, xi_max, xi_len)
            % number of argument check
            if nargin ~= 3
                error(INTERVAL_LOGSPACE, 'interval::logspace - method accepts three input argument.');
            end
            
            % logspace
            if isa(xi_min, 'interval') && ~isa(xi_max, 'interval')
                xi_max = interval(xi_max);
                xo_c = interval(logspace(xi_min.infimum, xi_max.infimum, xi_len),...
                    logspace(xi_min.supremum, xi_max.supremum, xi_len));
            elseif isa(xi_max, 'interval') && ~isa(xi_min, 'interval')
                xi_min = interval(xi_min);
                xo_c = interval(logspace(xi_min.infimum, xi_max.infimum, xi_len),...
                    logspace(xi_min.supremum, xi_max.supremum, xi_len));
            elseif isa(xi_max, 'interval') && isa(xi_min, 'interval')
                xo_c = interval(logspace(xi_min.infimum, xi_max.infimum, xi_len),...
                    logspace(xi_min.supremum, xi_max.supremum, xi_len));
            else
                xo_c = interval();
            end
        end

        % generate X, Y interval arrays for 3D plots
        function [xo_x, xo_y] = meshgrid(xi_x, xi_y)
            % number of argument check
            if nargin ~= 2
                error(INTERVAL_MESHGRID, 'interval::meshgrid - method accepts only two input argument.');
            end
            
            % meshgrid
            if isa(xi_x, 'interval') && ~isa(xi_y, 'interval')
                xi_y = interval(xi_y);
                [xlower, ylower] = meshgrid (xi_x.infimum, xi_y.infimum);
                [xupper, yupper] = meshgrid (xi_x.supremum, xi_y.supremum);
                xo_x = interval(xlower, xupper);
                xo_y = interval(ylower, yupper);
            elseif isa(xi_y, 'interval') && ~isa(xi_x, 'interval')
                xi_x = interval(xi_x);
                [xlower, ylower] = meshgrid (xi_x.infimum, xi_y.infimum);
                [xupper, yupper] = meshgrid (xi_x.supremum, xi_y.supremum);
                xo_x = interval(xlower, xupper);
                xo_y = interval(ylower, yupper);
            elseif isa(xi_x, 'interval') && isa(xi_y, 'interval')
                [xlower, ylower] = meshgrid (xi_x.infimum, xi_y.infimum);
                [xupper, yupper] = meshgrid (xi_x.supremum, xi_y.supremum);
                xo_x = interval(xlower, xupper);
                xo_y = interval(ylower, yupper);
            end
        end

        % size of interval
        function xo_size = size(xi_interval)
            xo_size = size(xi_interval.infimum);
        end
        
        % reshape interval array
        function xo_b = reshape(xi_a, xi_m, xi_n)
            % input argument manager
            if nargin == 2
                if numel(xi_m) ~= 2
                    error (INTERVAL_RESHAPE, 'interval::reshape - second input argument must include two elements.');
                else
                    xi_n = xi_m(2);
                    xi_m = xi_m(1);
                end
            elseif nargin == 3
                % best possible input type
            else
                error (INTERVAL_RESHAPE, 'interval::reshape - method requires two or three input arguments.');
            end

            % input type
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end

            % reshape
            lower = reshape(xi_a.infimum, xi_m, xi_n);
            upper = reshape(xi_a.supremum, xi_m, xi_n);

            % output
            xo_b = interval(lower, upper);
        end
        
        % subscript assignment overload
        function xo_interval = subsasgn(xi_a, xi_s, xi_b)
            % verify string as '()'
            if strcmp(xi_s, '()')
                error(INTERVAL_SUBSASGN, 'interval::subsasgn object requires subscripts with ''()''.');
            end
            
            % type transformation
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end
            
            % assign
            xo_interval = interval(subsasgn(xi_a.infimum, xi_s, xi_b.infimum),...
                                                                    subsasgn(xi_a.supremum, xi_s, xi_b.supremum));
        end

        % subscript reference overload
        function varagout = subsref(xi_a, xi_s)
            % string type
            switch xi_s(1).type
                case '()'
                    infi = subsref(xi_a.infimum, xi_s(1));
                    supre = subsref(xi_a.supremum, xi_s(1));
                    varagout = interval(infi, supre);

                case '.'
                    if any(strcmp(xi_s(1).subs, methods(xi_a)))
                        varagout = feval(xi_s(1).subs, xi_a);
                    elseif strcmp(xi_s(1).subs, 'infimum')
                        varagout = xi_a.infimum;
                    elseif strcmp(xi_s(1).subs, 'supremum')
                        varagout = xi_a.supremum;
                    else
                        error(INTERVAL_SUBSREF, 'interval::subsref - invalid property/method.');
                    end
                   
                otherwise
                    error(INTERVAL_SUBSREF, 'interval::subsref is used with an invalid subscript type..');
            end
            
            % longer string
            if numel(xi_s) > 1
                varagout = subsref(varagout, xi_s(2 : end));
            end
        end
        
    end % constructors and core

    % interval specific methods
    methods

        % return the middle point of an interval
        function xo_mid = middle(xi_interval)
            % number of argument check
            if nargin ~= 1
                error(INTERVAL_MIDDLE, 'interval::middle - method accepts only one input argument.');
            end
            
            % output
            xo_mid = (xi_interval.infimum / 2) + (xi_interval.supremum / 2);
            xo_mid(xi_interval.infimum == -Inf & xi_interval.supremum == Inf) = 0;
            xo_mid(isempty(xi_interval)) = NaN;
        end

        % return the width of an interval
        function xo_width = width(xi_interval)
            % number of argument check
            if nargin ~= 1
                error(INTERVAL_WIDTH, 'interval::width - method accepts only one input argument.');
            end
            
            % output
            xo_width = abs(xi_interval.supremum - xi_interval.infimum);
            xo_width(isempty(xi_interval)) = NaN;
        end

        % return the radius of an interval
        function xo_radius = radius(xi_interval)
            % number of argument check
            if nargin ~= 1
                error(INTERVAL_RADIUS, 'interval::radius - method accepts only one input argument.');
            end
            
            % output
            mid = middle(xi_interval);
            xo_radius = max(bsxfun(@minus, mid, xi_interval.infimum),...
                                                 bsxfun(@minus, xi_interval.supremum, mid));
            xo_radius(isempty(xi_interval)) = NaN;
        end

        % return the mignitude (minimum of absolute of infimum/supremum pairs of interval
        function xo_mig = mignitude(xi_interval)
            % number of argument check
            if nargin ~= 1
                error(INTERVAL_MIGNITUDE, 'interval::mignitude - method accepts only one input argument.');
            end
            
            % output
            xo_mig = min(abs(xi_interval.infimum), abs(xi_interval.supremum));
            xo_mig(abs(xo_mig) == Inf) = 0;
            xo_mig(isempty(xi_interval)) = NaN;
        end

        % return the maximum of absolute value of infimum/supremum pairs of interval
        function xo_mag = magnitude(xi_interval)
            % number of argument check
            if nargin ~= 1
                error(INTERVAL_MAGNITUDE, 'interval::magnitude - method accepts only one input argument.');
            end
            
            % output
            xo_mag = max(abs(xi_interval.infimum), abs(xi_interval.supremum));
            xo_mag(isempty(xi_interval)) = NaN;
        end

        % return true if inteval is empty (Inf, -Inf)
        function xo_empty = isempty(xi_interval)
            xo_empty = xi_interval.infimum > xi_interval.supremum;
        end

        % return true if an interval is actualy 1D, e.g - infimum is also the supremum
        function xo_single = issingle(xi_interval)
            % number of argument check
            if nargin ~= 1
                error(INTERVAL_ISSINGLE, 'interval::issingle - method accepts only one input argument.');
            end
            
            % output
            xo_single = xi_interval.infimum == xi_interval.supremum;
        end

        % check if a given value is within interval
        function xo_member = ismember(xi_value, xi_interval)
             % number of argument check
            if nargin ~= 2
                error(INTERVAL_ISMEMBER, 'interval::ismembere - method requires two input arguments.');
            end
            
            % input check
            if ~isscalar(xi_value) || isnan(xi_value)
                error(INTERVAL_ISMEMBER, 'interval::ismembere - first input argument must be a scalar.');
            end
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % test membership
            if isnumeric(xi_value)
                xo_member = ~isinf(xi_value) & xi_value >= xi_interval.infimum & xi_value <= xi_interval.supremum;
            else
                xo_member = in(interval(xi_value), xi_interval);
            end
        end
        
        % return true if first interval is contained within second interval
        function xo_in = in(xi_a, xi_b)
            % number of argument check
             if nargin ~= 2
                error(INTERVAL_IN, 'interval::in - method accepts two input argument.');
             end
             
            % housekeeping
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end

            % calculation
            xo_in = ((xi_b.infimum <= xi_a.infimum | (xi_b.infimum == -inf & xi_a.infimum == -inf)) & ...
                              (xi_a.supremum <= xi_b.supremum | (xi_a.supremum == inf & xi_b.supremum == inf)));
            xo_in(isempty(xi_a) & isempty(xi_b)) = 1;
        end

        % return an interval holding the intersection of two intervals
        % and a logical index showing where intersection did not occure.
        % if no intersection - then the interval is reversed - {Inf, -Inf}
        function [xo_c, xo_no_intersection] = intersect(xi_a, xi_b)
            % number of argument check
             if nargin ~= 2
                error(INTERVAL_INTERSECT, 'interval::intersect - method accepts two input argument.');
             end
             
            % housekeeping
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end

            % intersection
            lower = max(xi_a.infimum, xi_b.infimum);
            upper = min(xi_a.supremum, xi_b.supremum);
            xo_no_intersection = lower > upper;
            lower(xo_no_intersection) = Inf;
            upper(xo_no_intersection) = -Inf;
            xo_c = interval(lower, upper);
        end

        % return an interval holding the union of two intervals
        function xo_c = union(xi_a, xi_b)
            % number of argument check
             if nargin ~= 2
                error(INTERVAL_UNION, 'interval::union - method accepts two input argument.');
             end
             
            % housekeeping
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end

            % union
            xo_c = interval(min(xi_a.infimum, xi_b.infimum), max(xi_a.supremum, xi_b.supremum));
        end

        % return the relative complement of first interval in second interval
        function xo_c = diff(xi_a, xi_b)
            % number of argument check
             if nargin ~= 2
                error(INTERVAL_DIFF, 'interval::diff - method accepts two input argument.');
             end
            
            % housekeeping
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end
            
            % dimension adjustment
            if isscalar(xi_a.infimum) ~= isscalar(xi_b.infimum)
                xi_a.infimum = ones(size(xi_b.infimum)) .* xi_a.infimum;
                xi_a.supremum = ones(size (xi_b.infimum)) .* xi_a.supremum;
                xi_b.infimum = ones(size (xi_a.infimum)) .* xi_b.infimum;
                xi_b.supremum = ones(size (xi_a.infimum)) .* xi_b.supremum;
            end
            
            % initialization (lower <- xi_a)
            lower = xi_a.infimum;
            upper = xi_a.supremum;
            
            % [xi_a, xi_b)
            rightSide = (xi_b.supremum >= xi_a.supremum) && (xi_b.infimum > xi_a.infimum);
            upper(rightSide) = min(upper(rightSide), xi_b.infimum(rightSide));

            % (xi_a, xi_b]
            leftSide = (xi_b.infimum <= xi_a.infimum) && (xi_b.supremum < xi_a.supremum);
            lower(leftSide) = max(lower(leftSide), xi_b.supremum (leftSide));
            
            % in(xi_a, xi_b) <- empty
            emptyInterval = (xi_b.infimum <= xi_a.infimum) && (xi_b.supremum >= xi_a.supremum);
            lower(emptyInterval) = Inf;
            upper(emptyInterval) = -Inf;
            
            % output
            xo_c = interval(lower, upper);
        end
        
        % increase an interval boundary (lower and upper) to the next value
        function xo_interval = next(xi_interval)
            % number of argument check
            if nargin ~= 1
                error(INTERVAL_NEXT, 'interval::next - method accepts only one input argument.');
            end
            
            % output
            xo_interval = interval(xi_interval.infimum - eps, xi_interval.supremum + eps);
        end
        
        % return the hausdorff distance between two intervals
        function xo_c = hausdorffDist(xi_a, xi_b)
             % number of argument check
             if nargin ~= 2
                error(INTERVAL_HAUSDORFDIST, 'interval::hausdorffDist - method accepts two input argument.');
             end
            
            % housekeeping
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end
            
            % dimension adjustment
            if isscalar(xi_a.infimum) ~= isscalar(xi_b.infimum)
                xi_a.infimum = ones(size(xi_b.infimum)) .* xi_a.infimum;
                xi_a.supremum = ones(size(xi_b.infimum)) .* xi_a.supremum;
                xi_b.infimum = ones(size(xi_a.infimum)) .* xi_b.infimum;
                xi_b.supremum = ones(size(xi_a.infimum)) .* xi_b.supremum;
            end
            
            % housekeeping
            xo_c = zeros(size(xi_a.infimum));
            
            % xi_a.infimum < xi_b.infimum
            select = xi_a.infimum < xi_b.infimum;
            if any(any (select))
                xo_c(select) = xi_b.infimum(select) - xi_a.infimum(select);
            end
            
            % xi_a.infimum > xi_b.infimum
            select = xi_a.infimum > xi_b.infimum;
            if any(any(select))
                xo_c(select) = max(xo_c (select), xi_a.infimum(select) - xi_b.infimum(select));
            end
            
            % xi_a.supremum < xi_b.supremum
            select = xi_a.supremum < xi_b.supremum;
            if any(any(select))
                xo_c(select) = max(xo_c(select), xi_b.supremum(select) - xi_a.supremum(select));
            end
            
            % xi_a.supremum > xi_b.supremum
            select = xi_a.supremum > xi_b.supremum;
            if any (any (select))
                xo_c(select) = max(xo_c(select), xi_a.supremum(select) - xi_b.supremum(select));
            end

            % empty set
            xo_c(isempty(xi_a) || isempty(xi_b)) = NaN;
        end
        
        % inner distance
        function xo_c = innerDist(xi_a, xi_b)
             % number of argument check
             if nargin ~= 2
                error(INTERVAL_INNERDIST, 'interval::innerDist - method accepts two input argument.');
             end
            
            % housekeeping
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end
            
            % output dimension determination
            if isscalar(xi_a.infimum) ~= isscalar(xi_b.infimum)
                xi_a.infimum = ones(size(xi_b.infimum)) .* xi_a.infimum;
                xi_a.supremum = ones(size(xi_b.infimum)) .* xi_a.supremum;
                xi_b.infimum = ones(size(xi_a.infimum)) .* xi_b.infimum;
                xi_b.supremum = ones(size(xi_a.infimum)) .* xi_b.supremum;
            end

            % output initialization
            xo_c = zeros(size(xi_a.infimum));
            
            % "right boundary" distance
            condition = xi_a.supremum < xi_b.infimum;
            if any(any(condition))
                xo_c(condition) = xi_b.infimum(condition) - xi_a.supremum(condition);
            end
            
            % "left boundery" distance
            condition = xi_a.infimum > xi_b.supremum;
            if any(any(condition))
                xo_c(condition) = max(xo_c(condition), xi_a.infimum(condition) - xi_b.supremum (condition));
            end

            % empty intervals have no distance
            xo_c(isempty(xi_a) | isempty(xi_b)) = NaN;
        end
        
        % signed distance
        function xo_c = signedDist(xi_a, xi_b)
            % number of argument check
             if nargin ~= 2
                error(INTERVAL_SIGNEDDIST, 'interval::signedDist - method accepts two input argument.');
             end
            
            % housekeeping
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end
            
            % output dimension determination
            if isscalar(xi_a.infimum) ~= isscalar(xi_b.infimum)
                xi_a.infimum = ones(size(xi_b.infimum)) .* xi_a.infimum;
                xi_a.supremum = ones(size(xi_b.infimum)) .* xi_a.supremum;
                xi_b.infimum = ones(size(xi_a.infimum)) .* xi_b.infimum;
                xi_b.supremum = ones(size(xi_a.infimum)) .* xi_b.supremum;
            end

            % output initialization
            xo_c = zeros(size(xi_a.infimum));
            
            % "right boundary" distance
            condition = xi_a.supremum < xi_b.infimum;
            if any(any(condition))
                xo_c(condition) = xi_a.supremum(condition) - xi_b.infimum(condition);
            end
            
            % "left boundery" distance
            condition = xi_a.infimum > xi_b.supremum;
            if any(any(condition))
                xo_c(condition) = max(xo_c(condition), xi_a.infimum(condition) - xi_b.supremum (condition));
            end

            % empty intervals have no distance
            xo_c(isempty(xi_a) | isempty(xi_b)) = NaN;
        end
        
    end % interval specific methods

    % operator overloading
    methods

        % '+'
        function xo_c = plus(xi_a, xi_b)
            % check
            if nargin ~= 2
                error(INTERVAL_PLUS, 'interval::plus - two arguments are required.');
            end

            % type adjustment
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end

            % housekeeping
            aSize = size(xi_a);
            bSize = size(xi_b);
            if (aSize(1) ~= bSize(1)) || (aSize(2) ~= bSize(2))
                error(INTERVAL_PLUS, 'interval::plus - input arguments should be of equal dimensions.');
            end
            xo_c = interval(xi_a);

            % add column-wise
            for i = 1 : aSize(2)
                % addition
                tinfimum = bsxfun(@plus, xo_c.infimum(:, i), xi_b.infimum(:, i));
                tsupremum = bsxfun(@plus, xo_c.supremum(:, i), xi_b.supremum(:, i));

                % equality test
                equality = bsxfun(@eq, tinfimum, tsupremum);
                if any(equality)
                    tinfimum(equality) = tinfimum(equality) - eps;
                    tsupremum(equality) = tsupremum(equality) + eps;
                end

                % empty set
                empty = isempty(xi_a) | isempty(xi_b);
                tinfimum(empty) = Inf;
                tsupremum(empty) = -Inf;

                % output
                xo_c.infimum(:, i) = tinfimum;
                xo_c.supremum(:, i) = tsupremum;
            end
        end

        % '-'
        function xo_c = minus(xi_a, xi_b)
            % check
            if nargin ~= 2
                error(INTERVAL_MINUS, 'interval::minus - two arguments are required.');
            end

            % type adjustment
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end

            % housekeeping
            aSize = size(xi_a);
            bSize = size(xi_b);
            if (aSize(1) ~= bSize(1)) || (aSize(2) ~= bSize(2))
                error(INTERVAL_MINUS, 'interval::minus - input arguments should be of equal dimensions.');
            end
            xo_c = interval(xi_a);

            % add column-wise
            for i = 1 : aSize(2)
                % addition
                tinfimum = bsxfun(@minus, xo_c.infimum(:, i), xi_b.supremum(:, i));
                tsupremum = bsxfun(@minus, xo_c.supremum(:, i), xi_b.infimum(:, i));

                % equality test
                equality = bsxfun(@eq, tinfimum, tsupremum);
                if any(equality)
                    tinfimum(equality) = tinfimum(equality) - eps;
                    tsupremum(equality) = tsupremum(equality) + eps;
                end

                % empty set
                empty = isempty(xi_a) | isempty(xi_b);
                tinfimum(empty) = Inf;
                tsupremum(empty) = -Inf;

                % output
                xo_c.infimum(:, i) = tinfimum;
                xo_c.supremum(:, i) = tsupremum;
            end
        end

        % '*'
        function xo_c = times(xi_a, xi_b)
            % check
            if nargin ~= 2
                error(INTERVAL_TIMES, 'interval::times - two arguments are required.');
            end

            % type adjustment
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end

            %interval
            minMin = xi_a.infimum * xi_b.infimum;
            minMax = xi_a.infimum * xi_b.supremum;
            maxMin = xi_a.supremum * xi_b.infimum;
            maxMax = xi_a.supremum * xi_b.supremum;
            lower = min(min(min(minMin, minMax), maxMin), maxMax);
            upper = max(max(max(minMin, minMax), maxMin), maxMax);

            % full set
            temp = xi_a.infimum == -inf | xi_a.supremum == inf;
            lower(temp) = -Inf;
            upper(temp) = Inf;

            % if zero
            temp = (xi_a.infimum == 0 & xi_a.supremum == 0) |...
                (xi_b.infimum == 0 & xi_b.supremum == 0);
            lower(temp) = 0;
            upper(temp) = 0;

            % empty set
            temp = isempty(xi_a) | isempty(xi_b);
            lower(temp) = Inf;
            upper(temp) = -Inf;

            % output
            xo_c = interval(lower, upper);
        end

        % '.*'
        function xo_c = mtimes(xi_a, xi_b)
            % check
            if nargin ~= 2
                error(INTERVAL_TIMES, 'interval::mtimes - two arguments are required.');
            end

            % type adjustment
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end

            % if an input argument is scalar?
            if isscalar(xi_a.infimum) || isscalar(xi_b.infimum)
                xo_c = times(xi_a, xi_b);
                return;
            end

            % matrix multiplication using the "seven point multiplication" method
            aMiddle = middle(xi_a);
            bMiddle = middle(xi_b);
            aRadius = radius(xi_a);
            bRadius = radius(xi_b);
            aRho = sign(aMiddle) .* min(abs(aMiddle), aRadius);
            bRho = sign(bMiddle) .* min(abs(bMiddle), bRadius);
            cRadius = abs(aMiddle) * bRadius +...
                aRadius * (abs(bMiddle) + bRadius) -...
                abs(aRho) * abs(bRho);
            xo_c = interval(aMiddle * bMiddle + aRho * bRho - cRadius,...
                aMiddle * bMiddle + aRho * bRho + cRadius);
        end

        % './' (divide each entry of xi_a by the corresponding entry of xi_b)
        function xo_c = rdivide(xi_a, xi_b)
            % check
            if nargin ~= 2
                error(INTERVAL_RDIVIDE, 'interval::rdivide - two arguments are required.');
            end

            % type adjustment
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end

            % if xi_a == 1
            if all(all(xi_a.infimum == 1 & xi_a.supremum == 1))
                infInv = 1 ./ xi_b.infimum;
                supInv = 1 ./ xi_b.supremum;
                xo_c = interval(min(infInv, supInv), max(infInv, supInv));

                %rescale to matrix form
                if isscalar(xi_a.infimum) ~= isscalar(xi_b.infimum)
                    xo_c.infimum = xo_c.infimum .* ones(size(xi_a.infimum));
                    xo_c.supremum = xo_c.supremum .* ones(size(xi_a.infimum));
                end

                return;
            end

            % resize element to equal dimensions
            if isscalar(xi_a.infimum) ~= isscalar(xi_b.infimum)
                xi_a.infimum = xi_a.infimum .* ones(size(xi_b.infimum));
                xi_a.supremum = xi_a.supremum .* ones(size(xi_b.infimum));
                xi_b.infimum = xi_b.infimum .* ones(size(xi_a.infimum));
                xi_b.supremum = xi_b.supremum .* ones(size(xi_a.infimum));
            end

            % all possible interval slicing
            slice1 = xi_a.supremum <= 0 & xi_b.supremum < 0;
            slice2 = xi_a.supremum <= 0 & xi_b.infimum > 0;
            slice3 = xi_a.supremum <= 0 & xi_b.supremum == 0;
            slice4 = xi_a.supremum <= 0 & xi_b.supremum > 0 & xi_b.infimum == 0;
            slice5 = xi_a.infimum >= 0 & xi_a.supremum > 0 & xi_b.supremum < 0;
            slice6 = xi_a.infimum >= 0 & xi_a.supremum > 0 & xi_b.infimum > 0;
            slice7 = xi_a.infimum >= 0 & xi_a.supremum > 0 & xi_b.supremum == 0;
            slice8 = xi_a.infimum >= 0 & xi_a.supremum > 0 & xi_b.supremum > 0 & xi_b.infimum == 0;
            slice9 = xi_a.infimum < 0 & 0 < xi_a.supremum & xi_b.supremum < 0;
            slice10 = xi_a.infimum < 0 & 0 < xi_a.supremum & xi_b.infimum > 0;

            % infimum & supremum according to slicing
            lower = zeros(size(xi_a.infimum));
            upper = zeros(size(xi_a.infimum));

            lower(slice1) = xi_a.supremum (slice1) ./ xi_b.infimum (slice1);
            lower(slice2) = xi_a.infimum (slice2) ./ xi_b.infimum (slice2);
            lower(slice3) = xi_a.supremum (slice3) ./ xi_b.infimum (slice3);
            lower(slice4) = -inf;
            lower(slice5) = xi_a.supremum (slice5) ./ xi_b.supremum (slice5);
            lower(slice6) = xi_a.infimum (slice6) ./ xi_b.supremum (slice6);
            lower(slice7) = -inf;
            lower(slice8) = xi_a.infimum (slice8) ./ xi_b.supremum (slice8);
            lower(slice9) = xi_a.supremum (slice9) ./ xi_b.supremum (slice9);
            lower(slice10) = xi_a.infimum (slice10) ./ xi_b.infimum (slice10);
            upper(slice1) = xi_a.infimum (slice1) ./ xi_b.supremum (slice1);
            upper(slice2) = xi_a.supremum (slice2) ./ xi_b.supremum (slice2);
            upper(slice3) = inf;
            upper(slice4) = xi_a.supremum (slice4) ./ xi_b.supremum (slice4);
            upper(slice5) = xi_a.infimum (slice5) ./ xi_b.infimum (slice5);
            upper(slice6) = xi_a.supremum (slice6) ./ xi_b.infimum (slice6);
            upper(slice7) = xi_a.infimum (slice7) ./ xi_b.infimum (slice7);
            upper(slice8) = inf;
            upper(slice9) = xi_a.infimum (slice9) ./ xi_b.supremum (slice9);
            upper(slice10) = xi_a.supremum (slice10) ./ xi_b.infimum (slice10);

            % handle intervals which include a zero
            zeroEncircle = (xi_b.infimum < 0 & xi_b.supremum > 0) |...
                (xi_a.infimum < 0 & xi_a.supremum > 0 &...
                (xi_b.infimum == 0 | xi_b.supremum == 0));
            lower (zeroEncircle) = -inf;
            upper (zeroEncircle) = inf;

            % handle results which are equal to zero
            isZero = xi_a.infimum == 0 & xi_a.supremum == 0;
            lower(isZero) = 0;
            upper(isZero) = 0;

            % handle empty intervals
            emptyInterval = isempty (xi_a) | isempty (xi_b) | (xi_b.infimum == 0 & xi_b.supremum == 0);
            lower (emptyInterval) = inf;
            upper (emptyInterval) = -inf;

            % output
            xo_c = interval(lower, upper);
        end

        % '.\' (divide each entry of xi_b by the corresponding entry of xi_a)
        function xo_c = ldivide(xi_a, xi_b)
             % check
            if nargin ~= 2
                error(INTERVAL_LDIVIDE, 'interval::ldivide - two arguments are required.');
            end
            
            % output
            xo_c = rdivide(xi_b, xi_a);
        end

        % matrix '\' (matrix left divison)
        % if xi_a is square xo_c = inv(xi_a) * xi_b
        % if xi_a is square and xi_b is colmn xo_c = norm(xi_a  * X - xi_b)
        function xo_c = mldivide(xi_a, xi_b)
            % check
            if nargin ~= 2
                error(INTERVAL_MLDIVIDE, 'interval::mldivide - two arguments are required.');
            end

            % type adjustment
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end
            
            % housekeeping
            [arow, acol] = size(xi_a.infimum);
            [brow, bcol] = size(xi_b.infimum);
            
            %input argument compatibility
            if isscalar(xi_a.infimum) % xi_a is scalar
                xo_c = ldivide(xi_a, xi_b);
                return;
            elseif acol ~= arow % xi_a is not square
                error(INTERVAL_MLDIVIDE, 'interval::mldivide - first input argument is not a square matrix.');
            elseif arow ~= brow % xi_a & xi_b rows are not equal
                error(INTERVAL_MLDIVIDE, 'interval::mldivide - input arguments are not dimension compatible.');
            elseif isempty(xi_a) % xi_a empty interval
                xo_c = interval(zeros(0, bcol));
                return;
            end
            
            % housekeeping #2
            conditionTol = 1e-3;
            aMiddle = middle(xi_a);
            bMiddle = middle(xi_b);
            
            % invert xi_a
            aInv = inv(aMiddle);
            condNum = cond(aInv, 1) - 1;
            if abs(condNum < conditionTol)
                xo_c = linsolve(xi_a, xi_b);
                return;
            end
            aInvSize = size(aInv);
            aInvRow = aInvSize(1);
            aInvCol = aInvSize(2);
            
            % solution for standard numeric type
            xAprox = aInv * bMiddle;
            xAprox = xAprox + aInv * (bMiddle - aMiddle * xAprox);
            
            % residual in interval notation
            res = mtimes(aInv, xi_b - mtimes(xi_a, xAprox));
            c = eye(aInvRow) - mtimes(aInv, xi_a);
            
            % is this solution good enough?
            [res, goodSol] = intervalRefine(res, c);
            if goodSol
                xo_c = xAprox + res;
                return;
            end
            
            % precise solution
            res2 = zeros(size(aInvSize));
            for i = 1 : aInvRow
                for j = 1 : aInvCol
                    res2(i, j) = dot(aInv(i, :), aMiddle(:, j));
                end
            end
            
            % is this solution good enough?
            res2 = inv(res2);
            invFail = cond(res2, 1) - 1;
            if abs(invFail) < conditionTol
                xo_c = linsolve(xi_a, xi_b);
                return;
            end

            % solution is too sensitive to error in suplied data
            error(INTERVAL_MLDIVIDE, 'interval::mldivide - numerical operation too sensitive to error in first input argument.');
        end

        % matrix '/' (matrix right division)
        % if xi_a is square xo_c = xi_a * inv(xi_b)
        function xo_c = mrdivide(xi_a, xi_b)
             % check
            if nargin ~= 2
                error(INTERVAL_MRDIVIDE, 'interval::mrdivide - two arguments are required.');
            end

            % type adjustment
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end
            
            % scalar division
            if isscalar(xi_a) || isscalar(xi_b)
                xo_c = rdivide(xi_a, xi_b);
                return;
            end
            
            % division
            xo_c = transpose(mldivide(transpose(xi_b), transpose(xi_a)));
            %xo_c = xi_a * inv((mignitude(xi_b)));
        end

         % '.^'
        function xo_c = power(xi_a, xi_b)
            % number of arguments check
            if nargin ~= 2
                error(INTERVAL_POWER, 'interval::power - two arguments are required.');
            end
            
            % base can not be a matrix
            if min(size(xi_a)) > 1
                error(INTERVAL_POWER, 'interval::power - matrix base is not supported.');
            end
            
             % type adjustment
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            
            % housekeeping
            xo_c = interval(zeros(size(xi_a)));
            
            % iterate along base
            for i = 1 : numel(xi_a)
                % positive intervals
                if xi_a.infimum(i) >= 0
                    if isscalar(xi_b) &  xi_b >= 0 |...
                          isa(xi_b, interval) & xi_b.supremum >= 0 & xi_b.infimum >= 0 %#ok<OR2,AND2>
                      xo_c(i) = positivePower(xi_a(i), xi_b);
                    else
                        error(INTERVAL_POWER, 'interval::power - for positive interval base, only positive powers are allowed.');
                    end
                else % non positive intervals
                    if isnumeric(xi_b) && fix(xi_b) == xi_b
                        xo_c(i) = powerN(xi_a(i), xi_b);
                    else
                        error(INTERVAL_POWER, 'interval::power - for negative interval base, only positive integer powers are allowed.');
                    end
                end
            end
        end
        
        % unary minus overload
        function xo_interval = uminus(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_UMINUS, 'interval::uminus - one input arguments is required.');
            end
            
            % housekeeping
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % output
            xo_interval = interval(-xi_interval.supremum, -xi_interval.infimum);
        end

        % ==
        function xo_c = eq(xi_a, xi_b)
             % check input arguments
            if nargin ~= 2
                error(INTERVAL_EQ, 'interval::eq - two arguments are required.');
            end
            
            % housekeeping
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end

            % calculation
            xo_c = (xi_a.infimum == xi_b.infimum) & (xi_a.supremum == xi_b.supremum);
        end

        % ~=
        function xo_c = neq(xi_a, xi_b)
            % check input arguments
            if nargin ~= 2
                error(INTERVAL_NEQ, 'interval::neq - two arguments are required.');
            end
            
            % neq
            xo_c = ~eq(xi_a, xi_b);
        end

        % <
        function xo_c = lt(xi_a, xi_b)
            % check input arguments
            if nargin ~= 2
                error(INTERVAL_LT, 'interval::lt - two arguments are required.');
            end
            
            % housekeeping
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end

            % calculation
            xo_c = (((xi_a.infimum < xi_b.infimum) | (xi_a.infimum == -Inf & xi_b.infimum == -Inf)) &...
                ((xi_a.supremum < xi_b.supremum) | (xi_a.supremum == Inf & xi_b.supremum == Inf)));
            xo_c(isempty(xi_a) & isempty(xi_b)) = 1;
        end

        % >
        function xo_c = gt(xi_a, xi_b)
            % check input arguments
            if nargin ~= 2
                error(INTERVAL_GT, 'interval::gt - two arguments are required.');
            end
            
            % gt
            xo_c = lt(xi_b, xi_a);
        end

        % <=
        function xo_c = le(xi_a, xi_b)
            % check input arguments
            if nargin ~= 2
                error(INTERVAL_LE, 'interval::le - two arguments are required.');
            end
            
            % housekeeping
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end

            % calculation
            xo_c = xi_a.infimum <= xi_b.infimum & xi_a.supremum <= xi_b.supremum;
        end

        % >=
        function xo_c = ge(xi_a, xi_b)
            % check input arguments
            if nargin ~= 2
                error(INTERVAL_GE, 'interval::ge - two arguments are required.');
            end
            
            % ge
            xo_c = le(xi_b, xi_a);
        end

        % and (set intersection))
        function xo_c = and(xi_a, xi_b)
            % check input arguments
            if nargin ~= 2
                error(INTERVAL_AND, 'interval::and - two arguments are required.');
            end
            
            % housekeeping
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end

            % and
            [xo_c, dummy] = intersect(xi_a, xi_b); %#ok<NASGU>
        end

        % or (set union)
        function xo_c = or(xi_a, xi_b)
            % check input arguments
            if nargin ~= 2
                error(INTERVAL_OR, 'interval::or - two arguments are required.');
            end
            
            % housekeeping
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end

            % and
            xo_c = union(xi_a, xi_b);
        end

        % xor
        % non xor get the reversed interval notation {Inf, -Inf}
        function xo_c = xor(xi_a, xi_b)
            % check input arguments
            if nargin ~= 2
                error(INTERVAL_XOR, 'interval::xor - two arguments are required.');
            end
            
            % housekeeping
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end

            % size update (in case one is a lonely interval and the other a matrix)
            if isscalar(xi_a.infimum) ~= isscalar(xi_b.infimum)
                len = size(xi_b.infimum);
                xi_a.infimum = ones(len) .* xi_a.infimum;
                xi_a.supremum = ones(len) .* xi_a.supremum;

                len = size(xi_a.infimum);
                xi_b.infimum = ones(len) .* xi_b.infimum;
                xi_b.supremum = ones(len) .* xi_b.supremum;
            end

            % all possible bounderies
            lower = min(xi_a.infimum, xi_b.infimum);
            maxMin = max(xi_a.infimum, xi_b.infimum);
            minMax = min(xi_a.supremum, xi_b.supremum);
            upper = max(xi_a.supremum, xi_b.supremum);

            % bounderies borders
            [maxMin, minMax] = deal(min(maxMin, minMax), max(maxMin, minMax));

            % equal infinite infimum
            equal = xi_a.infimum == xi_b.infimum;
            lower(equal) = minMax(equal);

            % equal supremum
            equal = xi_a.supremum == xi_b.supremum;
            upper(equal) = maxMin (equal);

            % equal interval's (xor is an entire interval
            equal = (xi_a.infimum == xi_b.infimum) & (xi_a.supremum == xi_b.supremum);
            lower(equal) = inf;
            upper(equal) = -inf;

            % xor
            xo_c = interval(lower, upper);
        end

    end
    
    % general algebric & numerical methods
    methods
        
        % return the minimal interval enclosure
        function xo_min = min(xi_interval1, xi_interval2)
            if nargin == 1
                lower = min(xi_interval1.infimum);
                lower(any(isempty(xi_interval1))) = Inf;
                xo_min = interval(lower, min(xi_interval1.supremum));
            else
                lower = min(xi_interval1.infimum, xi_interval2.infimum);
                lower(isempty(xi_interval1) | isempty(xi_interval2)) = inf;
                xo_min = interval(lower, min(xi_interval1.supremum, xi_interval2.supremum));
            end
        end

        % return the maximal interval enclosure
        function xo_max = max(xi_interval1, xi_interval2)
            if nargin == 1
                upper = max(xi_interval1.supremum);
                upper(any(isempty(xi_interval1))) = -Inf;
                xo_max = interval(max(xi_interval1.infimum), upper);
            else
                upper = max(xi_interval1.supremum, xi_interval2.supremum);
                upper(any(isempty(xi_interval1) | isempty(xi_interval2))) = -Inf;
                xo_max = interval(max(xi_interval1.infimum, xi_interval2.infimum), upper);
            end
        end

        % hypot
        function xo_c = hypot(xi_a, xi_b)
             % input argument check
            if nargin ~= 2
                error(INTERVAL_HYPOT, 'interval::hypot - two input arguments are required.');
            end

            % housekeeping
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end
            
            % housekeeping
            aMig = mignitude(xi_a);
            aMag = magnitude(xi_a);
            bMig = mignitude(xi_b);
            bMag = magnitude(xi_b);
            
            % interval bounderies
            lower = sqrt(aMig * aMig + bMig * bMig);
            upper = sqrt(aMag * aMag + bMag * bMag);
            
            % empty interval
            empty = isempty(xi_a) | isempty(xi_b);
            lower(empty) = Inf;
            upper(empty) = -Inf;
            
            % output
            xo_c = interval(lower, upper);
        end

        % absolute
        function xo_interval = abs(xi_interval)
            % input argument check
            if nargin ~= 1
                error(INTERVAL_ABS, 'interval::abs - can not handle more then one input argument.');
            end

            % housekeeping
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % initialization is good for positive intervals
            xo_interval = interval(xi_interval);

            % negative supremum
            negSup = bsxfun(@le, xi_interval.supremum, 0);
            xo_interval.infimum(negSup) = -xi_interval.supremum(negSup);
            xo_interval.supremum(negSup) = -xi_interval.infimum(negSup);

            % intervals surrounding a zero
            includeZero = bsxfun(@le, xi_interval.infimum, 0) & ~negSup;
            xo_interval.infimum(includeZero) = 0;
            xo_interval.supremum(includeZero) = max(-xi_interval.infimum(includeZero), xi_interval.supremum(includeZero));
        end

        % sign
        function xo_interval = sign(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_SIGN, 'interval::sign - one input arguments is required.');
            end
            
            % housekeeping
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % sign
            lower = sign(xi_interval.infimum);
            upper = sign(xi_interval.supremum);

            % empty interval
            empty = isempty(xi_interval);
            lower(empty) = Inf;
            upper(empty) = -Inf;

            % output
            xo_interval = interval(lower, upper);
        end

        % round
        function xo_interval = round(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_ROUND, 'interval::round - one input arguments is required.');
            end
            
            % housekeeping
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % output
            xo_interval = interval(round(xi_interval.infimum), round(xi_interval.supremum));
        end
        
        % ceil
        function xo_interval = ceil(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_CEIL, 'interval::ceil - one input arguments is required.');
            end
            
            % housekeeping
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % output
            xo_interval = interval(ceil(xi_interval.infimum), ceil(xi_interval.supremum));
        end
        
        % floor
        function xo_interval = floor(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_FLOOR, 'interval::floor - one input arguments is required.');
            end
            
            % housekeeping
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % output
            xo_interval = interval(floor(xi_interval.infimum), floor(xi_interval.supremum));
        end
        
        % fix
        function xo_interval = fix(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_FIX, 'interval::fix - one input arguments is required.');
            end
            
            % housekeeping
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % output
            xo_interval = interval(fix(xi_interval.infimum), fix(xi_interval.supremum));
        end
                                
        % sqrt
        function xo_interval = sqrt(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_SQRT, 'interval::sqrt - one input arguments is required.');
            end
            
            % housekeeping
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % '0' manage
            lower = sqrt(max(0, xi_interval.infimum));
            upper = sqrt(max(0, xi_interval.supremum));
            
            % empty/negative interval
            negInterval = isempty(xi_interval) | (xi_interval.supremum < 0);
            lower(negInterval) = Inf;
            upper(negInterval) = -Inf;

            % output
            xo_interval = interval(lower, upper);
        end
        
        % log
        function xo_interval = log(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_LOG, 'interval::log - one input arguments is required.');
            end
            
            % housekeeping
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % transform to [0, Inf]
            [xi_interval, dummy] = intersect(xi_interval, interval(0, Inf)); %#ok<NASGU>

            % log
            lower = log(xi_interval.infimum);
            upper = log(xi_interval.supremum);

            % zero check
            supZero = xi_interval.supremum == 0;
            lower(supZero) = Inf;
            upper(supZero) = -Inf;

            % empty interval
            upper(isempty(xi_interval)) = -Inf;

            % assignment
            xo_interval = interval(lower, upper);
        end

        % log2
        function xo_interval = log2(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_LOG2, 'interval::log2 - one input arguments is required.');
            end
            
            % housekeeping
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % transform to [0, Inf]
            [xi_interval, dummy] = intersect(xi_interval, interval(0, Inf)); %#ok<NASGU>

            % log2
            lower = log2(xi_interval.infimum);
            upper = log2(xi_interval.supremum);

            % zero check
            supZero = xi_interval.supremum == 0;
            lower(supZero) = Inf;
            upper(supZero) = -Inf;

            % empty interval
            upper(isempty(xi_interval)) = -Inf;

            % assignment
            xo_interval = interval(lower, upper);
        end

        % log10
        function xo_interval = log10(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_LOG10, 'interval::log10 - one input arguments is required.');
            end
            
            % housekeeping
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % transform to [0, Inf]
            [xi_interval, dummy] = intersect(xi_interval, interval(0, Inf)); %#ok<NASGU>

            % log10
            lower = log10(xi_interval.infimum);
            upper = log10(xi_interval.supremum);

            % zero check
            supZero = xi_interval.supremum == 0;
            lower(supZero) = Inf;
            upper(supZero) = -Inf;

            % empty interval
            upper(isempty(xi_interval)) = -Inf;

            % assignment
            xo_interval = interval(lower, upper);
        end

        % exp
        function xo_interval = exp(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_EXP, 'interval::exp - one input arguments is required.');
            end
            
            % housekeeping
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % exp
            lower = exp(xi_interval.infimum);
            upper = exp(xi_interval.supremum);
            upper(isempty(xi_interval)) = -inf;

            % assignment
            xo_interval = interval(lower, upper);
        end

        %pow2 (c's ldexp function)
        function xo_interval = pow2(xi_interval)
            % input argument check
            if nargin ~= 1
                error(INTERVAL_POW2, 'interval::pow2 - only one input argument is expected.');
            end
            
            % interval type
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % pow2
            lower = pow2(xi_interval.infimum);
            upper = pow2(xi_interval.supremum);
            
            % empty set
            upper(isempty(xi_interval)) = Inf;
            
            % output
            xo_interval = interval(lower, upper);
        end
        
        % sum
        function xo_sum = sum(xi_interval, xi_dim)
            % check input arguments
            if nargin > 2
                error(INTERVAL_SUM, 'interval::sum - one or two input arguments are required.');
            end
            
            % interval type
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % dimension not supplied
            if nargin < 2
                xi_dim = find(size(xi_interval.infimum) > 1, 1);
                if isempty (xi_dim)
                    xi_dim = 1;
                end
            end
            
            % one element interval
            if size (xi_interval.infimum, xi_dim) == 1
                xo_sum = xi_interval;
                return;
            end
            
            % output size determination
            if xi_dim == 1 % array
                outDim = [1, max(1, size(xi_interval.infimum, 2))];
            elseif xi_dim == 2 % matrix
                outDim = [max(1, size(xi_interval.infimum, 1)), 1];
            else
                error(INTERVAL_SUM, 'interval::sum - interval class does not support dimesnions higher then 2.');
            end
            
            % housekeeping
            lower = zeros (outDim);
            upper = zeros (outDim);

            % iterate over vectors
            for i = 1 : size (xi_interval.infimum, 3 - xi_dim)
                % prepare type and dimension
                ident.type = '()';
                ident.subs = cell(1, 2);
                ident.subs{xi_dim} = '.';
                ident.subs{3 - xi_dim} = i;

                %pick vector
                if size(xi_interval.infimum, 3 - xi_dim) == 1
                     vector.x = xi_interval;
                else
                    vector.x = subsref(xi_interval, ident);
                end
                
                % sum
                if any(isempty (vector.x))
                    lower(i) = Inf;
                    upper(i) = -Inf;
                else
                    lower(i) = sum(vector.x.infimum);
                    upper(i) = sum(vector.x.supremum);
                end
            end
            
            % output
            xo_sum = interval(lower, upper);
        end
        
        % dot product
        function xo_interval = dot(xi_a, xi_b, xi_dim)
            % input argument check
            if nargin == 1 || nargin > 3
                error(INTERVAL_DOT, 'interval::dot - two or three input argument are expected.');
            end
            
            % interval variable reauired
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end
            
            % dimension adaptation
            if nargin < 3
                if isvector(xi_a.infimum) && isvector(xi_b.infimum)
                    xi_dim = 1;
                    xi_a.infimum = xi_a.infimum(:);
                    xi_a.supremum = xi_a.supremum(:);
                    xi_b.infimum = xi_b.infimum(:);
                    xi_b.supremum = xi_b.supremum(:);
                else
                    xi_dim = find(any(min(size(xi_a.infimum)), min(size(xi_b.infimum))) > 1, 1);
                    if (isempty (xi_dim))
                        xi_dim = 1;
                    end
                end
            end
            
            % empty interval matrix
            if isempty(xi_a.infimum) || isempty(xi_b.infimum)
                xo_interval = interval(zeros(min(size(xi_a.infimum), size(xi_b.infimum))));
                return;
            end
            
            % dimension check
            if ((min(size(xi_a.infimum, 1), size(xi_b.infimum, 1)) > 1 && ...
                    size(xi_a.infimum, 1) ~= size(xi_b.infimum, 1)) || ...
                   (min(size(xi_a.infimum, 2), size(xi_b.infimum, 2)) > 1 && ...
                    size(xi_a.infimum, 2) ~= size(xi_b.infimum, 2)))
                error (INTERVAL_DOT, 'interval::dot - input arguments size missmatch.');
            end
            
            % housekeeping
            outputDim = max(size(xi_a.infimum), size(xi_b.infimum));
            outputDim(xi_dim) = 1;
            lower = zeros(outputDim);
            upper = zeros(outputDim);

            % dot product by dimension
            for n = 1 : numel (lower)
                ident.type = '()';
                ident.subs = cell(1, 2);
                ident.subs{xi_dim} = ':';
                ident.subs{3 - xi_dim} = n;
                
                % column/row selection
                if size(xi_a.infimum, 3 - xi_dim) == 1
                    left = xi_a;
                else
                    left = subsref(xi_a, ident);
                end
                if size(xi_b.infimum, 3 - xi_dim) == 1
                    right = xi_b;
                else
                    right = subsref(xi_b, ident);
                end
                
                % dot
                if isscalar(left) && isscalar(right)
                    temp = left .* right;
                    lower(n) = temp.infimum;
                    upper(n) = temp.supremum;
                else
                    lower(n) = sum(left.infimum .* right.infimum);
                    upper(n) = sum(left.supremum .* right.supremum);
                end
            end
            
            % output
            xo_interval = interval(lower, upper);
        end
        
        % product
        function xo_interval = prod(xi_interval, xi_dim)
            % input argument check
            if nargin ~= 1 && nargin ~= 2
                error(INTERVAL_PROD, 'interval::prod - one or two input argument are expected.');
            end
            
            % dimension determination
            if nargin < 2
                xi_dim = find(size(xi_interval.infimum) > 1, 1);
                if isempty (xi_dim)
                    xi_dim = 1;
                end
            end
            
            % output dimension
            if xi_dim == 1 % array
                xo_interval = interval(ones(1, max(1, size (xi_interval.infimum, 2))));
            elseif xi_dim == 2 % matrix
                xo_interval = interval(ones(max(1, size(xi_interval.infimum, 1)), 1));
            else
                 error(INTERVAL_PROD, 'interval::prod - interval class does not support dimesnions higher then 2.');
            end
            
            % special sets (empty/full)
            emptyInterval = any(isempty(xi_interval), xi_dim);
            xo_interval.infimum(emptyInterval) = Inf;
            xo_interval.supremum(emptyInterval) = -Inf;
            
            zeroInterval = ~(emptyInterval) & any(xi_interval.infimum == 0 & xi_interval.supremum == 0, xi_dim);
            xo_interval.infimum(zeroInterval) = 0;
            xo_interval.supremum(zeroInterval) = 0;
            
            fullInterval = ~(emptyInterval | zeroInterval) & any(xi_interval.infimum == -Inf & xi_interval.supremum == Inf, xi_dim);
            xo_interval.infimum(fullInterval) = -Inf;
            xo_interval.supremum(fullInterval) = Inf;

            % product
            ident.type = '()';
            ident.subs = {':', ':'};
            ident.subs{3 - xi_dim} = ~(emptyInterval || zeroInterval || fullInterval);
            if any(ident.subs{3 - xi_dim})
                ident.subs{xi_dim} = 1;
                outcome = subsref(xo_interval, ident);
                for i = 1 : size(xi_interval.infimum, xi_dim)
                    ident.subs{xi_dim} = i;
                    outcome = times(outcome, subsref(xi_interval, ident));
                end
                ident.subs{xi_dim} = 1;
                xo_interval = subsasgn(xo_interval, ident, outcome);
            end
        end

    end % general algebric methods
    
    % trigonometric
    methods
        
        % sin (input is in radians)
        function xo_interval = sin(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_SIN, 'interval::sin - one input arguments is required.');
            end
            
            % input argument type
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % housekeeping
            lower = zeros(size(xi_interval.infimum));
            upper = lower;
            cosSignLower = lower;
            cosSignUpper = lower;

            % |xi_interval| >= 2*pi
            wid = width(xi_interval);
            twoPi = wid >= 2*pi;
            lower(twoPi) = -1;
            upper(twoPi) = 1;

            % |xi_interval| < 2*pi
            inRegion = ~twoPi;
            lower(inRegion) = min(sin(xi_interval.infimum (inRegion)), sin(xi_interval.supremum(inRegion)));
            upper(inRegion) = max(sin(xi_interval.infimum (inRegion)), sin(xi_interval.supremum (inRegion)));

            % derivative @ < 2*pi.
            cosSignLower(inRegion) = sign(cos(xi_interval.infimum(inRegion)));
            cosSignUpper(inRegion) = sign(cos(xi_interval.supremum (inRegion)));

            % crossover points.
            cosSignLower(cosSignLower == 0) = sign(lower(cosSignLower == 0));
            cosSignUpper(cosSignUpper == 0) = -1 * sign(upper(cosSignUpper == 0));

            % {-Inf[    ]Inf}
            divergence = inRegion & ((cosSignLower == -1 & cosSignUpper == 1) | ...
                                                                      (cosSignLower == cosSignUpper & wid >= pi));
            lower (divergence) = -1;

            divergence = inRegion & ((cosSignLower == 1 & cosSignUpper == -1) | ...
                                                                      (cosSignLower == cosSignUpper & wid >= pi));
            upper (divergence) = 1;

            % empty interval
            emptyInterval = isempty(xi_interval);
            lower(emptyInterval) = Inf;
            upper(emptyInterval) = -Inf;

            % output
            xo_interval = interval(lower, upper);
        end

        % cos (input in radians)
        function xo_interval = cos(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_COS, 'interval::cos - one input arguments is required.');
            end
            
            % input argument type
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % housekeeping
            lower = zeros(size(xi_interval.infimum));
            upper = lower;
            sinSignLower = lower;
            sinSignUpper = lower;

            % |xi_interval| >= 2*pi
            wid = width(xi_interval);
            twoPi = wid >= 2*pi;
            lower (twoPi) = -1;
            upper (twoPi) = 1;

             % |xi_interval| < 2*pi
             inRegion = ~twoPi;
             lower(inRegion) = min(cos(xi_interval.infimum(inRegion)), cos(xi_interval.supremum(inRegion)));
             upper(inRegion) = max(cos(xi_interval.infimum(inRegion)), cos(xi_interval.supremum(inRegion)));

                % derivative @ < 2*pi
                sinSignLower(inRegion) = sign(sin(xi_interval.infimum(inRegion)));
                sinSignUpper(inRegion) = sign(sin(xi_interval.supremum(inRegion)));

                 % crossover points
                sinSignLower(sinSignLower == 0) = -1 * sign(lower(sinSignLower == 0));
                sinSignUpper(sinSignUpper == 0) = sign(upper(sinSignUpper == 0));

                % {-Inf[    ]Inf}
                divergence = inRegion & ((sinSignLower == -1 & sinSignUpper == 1) | ...
                                                                          (sinSignLower == sinSignUpper & wid >= pi)) & ne(0, xi_interval);
                lower(divergence) = -1;

                divergence = inRegion & ((sinSignLower == 1 & sinSignUpper == -1) | ...
                                                                          (sinSignLower == sinSignUpper & wid >= pi));
                upper(divergence) = 1;

                % empty interval
                emptyInterval = isempty(xi_interval);
                lower(emptyInterval) = inf;
                upper(emptyInterval) = -inf;

                % output
                xo_interval = interval(lower, upper);
        end

        % tan (input in radians)
        function xo_interval = tan(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_TAN, 'interval::tan - one input arguments is required.');
            end
            
            % input argument type
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % housekeeping
            lower = zeros (size (xi_interval.infimum));
            upper = lower;

            % |xi_interval| < pi
            wid = width(xi_interval);
            overPi = wid >= pi;
            inRegion = ~overPi;
            lower(inRegion) = tan(xi_interval.infimum(inRegion));
            upper(inRegion) = tan(xi_interval.supremum(inRegion));

            % divergence
            divergence = overPi | (lower > upper) | (wid > 2 & (sign(lower) == sign(upper) | max(abs(lower), abs(upper)) < 1));
            lower(divergence) = -Inf;
            upper(divergence) = Inf;

            % empty interval
            emptyInterval = isempty(xi_interval);
            lower(emptyInterval) = Inf;
            upper(emptyInterval) = -Inf;

            % output
            xo_interval = interval(lower, upper);
        end

        % atan
        function xo_interval = atan(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_ATAN, 'interval::atan - one input arguments is required.');
            end
            
            % input argument type
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % atan
            lower = atan(xi_interval.infimum);
            upper = atan(xi_interval.supremum);
            
            % empty interval
            empty = isempty(xi_interval);
            lower(empty) = Inf;
            upper(empty) = -Inf;
            
            % output
            xo_interval = interval(lower, upper);
        end
        
        % atan2
        function xo_c = atan2(xi_a, xi_b)
            % check input arguments
            if nargin ~= 2
                error(INTERVAL_ATAN2, 'interval::atan2 - two input arguments are required.');
            end
            
            % input argument type
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end
            
            % output dimension determination
            if isscalar(xi_b.infimum) ~= isscalar(xi_a.infimum)
                xi_b.infimum = ones(size(xi_a.infimum)) .* xi_b.infimum;
                xi_b.supremum = ones(size(xi_a.infimum)) .* xi_b.supremum;
                xi_a.infimum = ones(size(xi_b.infimum)) .* xi_a.infimum;
                xi_a.supremum = ones(size(xi_b.infimum)) .* xi_a.supremum;
            end

            % what portions of input arguments is in what unit circle quadrant
            xLeftTop = intersect(xi_b, interval (-Inf, 0));            
            yLeftTop = intersect(xi_a, interval (0, Inf));
            xRightTop = intersect(xi_b, interval (0, Inf));
            yRightBottom = intersect(xi_a, interval (-Inf, 0));
            xLeftBottom = xLeftTop;
            yRightTop = yLeftTop;
            xRightBottom = xRightTop;
            yLeftBottom = yRightBottom;

            % unit circle quadrants
            leftTop = ~(isempty(xLeftTop) | isempty(yLeftTop)) & (xLeftTop.infimum < 0 | yLeftTop.supremum > 0);
            rightTop = ~(isempty(xRightTop) | isempty(yRightTop)) & (xRightTop.supremum > 0 | yRightTop.supremum > 0);
            rightBottom = ~(isempty(xRightBottom) | isempty(yRightBottom)) & (xRightBottom.supremum > 0 | yRightBottom.infimum < 0);
            leftBottom = ~(isempty(xLeftBottom) | isempty(yLeftBottom)) & (xLeftBottom.infimum < 0 | yLeftBottom.infimum < 0);

            % atan2(0, ~)
            condition = leftTop & yLeftTop.supremum == 0;
            xLeftTop.infimum(condition) =  -1;
            xLeftTop.supremum (condition) = -1;
            
            condition = rightTop & yRightTop.supremum == 0;
            xRightTop.infimum(condition) = 1;
            xRightTop.supremum(condition) = 1;
            
            condition = rightBottom & yRightBottom.infimum == 0;
            xRightBottom.infimum(condition) = 1;
            xRightBottom.supremum(condition) = 1;
            
            leftBottom(yLeftBottom.infimum == 0) = 0;

            % atan2(~, 0)
            condition = leftTop & xLeftTop.infimum == 0;
            yLeftTop.infimum(condition) =  1;
            yLeftTop.supremum(condition) = 1;
            
            condition = rightTop & xRightTop.supremum == 0;
            yRightTop.infimum(condition) =  1;
            yRightTop.supremum(condition) = 1;
            
            condition = rightBottom & xRightBottom.supremum == 0;
            yRightBottom.infimum(condition) =  -1;
            yRightBottom.supremum(condition) = -1;
            
            condition = leftBottom & xLeftBottom.infimum == 0;
            yLeftBottom.infimum(condition) =  -1;
            yLeftBottom.supremum(condition) = -1;

            % atan2(0, <0)
            yLeftTop.infimum(leftTop & yLeftTop.infimum == 0) = +0;
            yLeftBottom.supremum(leftBottom & yLeftBottom.supremum == 0) = -0;

            % housekeeping
            sze = ones(size(leftTop));
            lower = sze * Inf;
            upper = -Inf * sze;
            
            % lower boundery
            select = leftBottom;
            lower(select) = atan2(yLeftBottom.supremum(select), xLeftBottom.infimum(select));
            
            select = rightBottom & ~leftBottom;
            lower(select) = atan2(yRightBottom.infimum(select), xRightBottom.infimum(select));
            
            select = rightTop & ~(rightBottom | leftBottom);
            lower(select) = atan2(yRightTop.infimum(select), xRightTop.supremum(select));
            
            select = leftTop & ~(rightTop | rightBottom | leftBottom);
            lower(select) = atan2(yLeftTop.supremum(select), xLeftTop.supremum(select));

            % upper boundary
            select = leftTop;
            upper(select) = atan2(yLeftTop.infimum(select), xLeftTop.infimum(select));
            
            select = rightTop & ~leftTop;
            upper(select) = atan2(yRightTop.supremum(select), xRightTop.infimum(select));
            
            select = rightBottom & ~(leftTop | rightTop);
            upper(select) = atan2(yRightBottom.supremum(select), xRightBottom.supremum(select));
            
            select = leftBottom & ~(leftTop | rightTop | rightBottom);
            upper(select) = atan2(yLeftBottom.infimum(select), xLeftBottom.supremum(select));

            % output
            xo_c = interval(lower, upper);
        end
        
    end % trigonometric
    
    % hyperbolic
    methods
        
        % sinh
        function xo_interval = sinh(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_SINH, 'interval::sinh - one input arguments is required.');
            end
            
            % input argument type
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % atan
            lower = sinh(xi_interval.infimum);
            upper = sinh(xi_interval.supremum);
            
            % empty interval
            empty = isempty(xi_interval);
            lower(empty) = Inf;
            upper(empty) = -Inf;
            
            % output
            xo_interval = interval(lower, upper);
        end
        
        % cosh
        function xo_interval = cosh(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_COSH, 'interval::cosh - one input arguments is required.');
            end
            
            % input argument type
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % atan
            lower = cosh(xi_interval.infimum);
            upper = cosh(xi_interval.supremum);
            
            % empty interval
            empty = isempty(xi_interval);
            lower(empty) = Inf;
            upper(empty) = -Inf;
            
            % output
            xo_interval = interval(lower, upper);
        end
        
        % tanh
        function xo_interval = tanh(xi_interval)
            % check input arguments
            if nargin ~= 1
                error(INTERVAL_TANH, 'interval::tanh - one input arguments is required.');
            end
            
            % input argument type
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % atan
            lower = tanh(xi_interval.infimum);
            upper = tanh(xi_interval.supremum);
            
            % empty interval
            empty = isempty(xi_interval);
            lower(empty) = Inf;
            upper(empty) = -Inf;
            
            % output
            xo_interval = interval(lower, upper);
        end
        
    end  % hyperbolic
    
    % matrix methods
    methods
        % transpose
        function xo_interval = transpose(xi_interval)
            % input argument check
            if nargin ~= 1
                error(INTERVAL_TRANSPOSE, 'interval::transpose - only one input argument is expected.');
            end
            
            % housekeeping
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % calculation
            xo_interval = interval(transpose(xi_interval.infimum),...
                                                                    transpose(xi_interval.supremum));
        end
        
        % digonal element
        function xo_diag = diag(xi_interval)
            % input argument check
            if nargin ~= 1
                error(INTERVAL_DIAG, 'interval::diag - only one input argument is expected.');
            end
            
            % housekeeping
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % output
            xo_diag = interval(diag(xi_interval.infimum), diag(xi_interval.supremum));
        end
        
        % LU decomposition of interval matrix
        % third input argument (xo_p), the permutation matrix, can be excluded
        function [xo_l, xo_u, xo_p] = lu(xi_a)
            % input argument check
            if nargin ~= 1
                error(INTERVAL_LU, 'interval::lu - only one input argument is expected.');
            end
            
            % interval variable reauired
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            
            % housekeeping
            aSize = size(xi_a);
            arow = aSize(1);
            acol = aSize(2);
            
            % check that input argument is square
            if acol ~= arow
                error(INTERVAL_LU, 'interval::lu - input argument must be a square matrix.');
            end
            
            % input is scalar            
            if isscalar(xi_a.infimum)
                xo_l = eye(arow);            
                xo_p = xo_l;
                xo_u = matrixPermute(xo_p, xi_a);
                return;
            end
            
            % permutation matrix preallocation
            xo_p = zeros(arow);
            
            % pivot manufactures
            aMig = mignitude(xi_a);
            aMag = magnitude(xi_a);
            aMig(isnan(aMig)) = 0;
            aMag(isnan(aMag)) = 0;
            
            % permutation matrix (pivot-style)
            for i = 1 : arow
                % columns with minimal number of positive mignitude
                colScore = sum(aMig > 0, 1);
                colScore(max(aMig) == Inf) = Inf;
                
                % column with minimal score
                colPossible = colScore == min(colScore);
                col = find(colPossible, 1);
                
                % zero containing intervals should not be used
                if colScore(col) >= 1
                    row = aMig(:, col) > 0;
                else
                    % choose interval which other than zero includes a number
                    row = (aMig(:, col) > 0) & (aMag(:, col) > 0);
                    
                    % all intervals contain zero
                    if ~max(row)
                        row = aMig(:, col) >= 0;
                    end
                end
                
                % find row with minimal number of positive mignitude
                if sum(row) == 1
                    row = find(row, 1);
                else
                    rowScore = sum(aMig > 0, 2);
                    rowScore = rowScore + sum(aMig(:, col) > 0, 2) / 2;
                    rowScore(~row) = Inf;
                    row = find(rowScore == min(rowScore), 1);
                end
                
                % stamp the chosen row and column
                xo_p(row, col) = 1;
                
                % mark used rows/columns from further usage
                aMig(row, :) = -Inf;
                aMig(isnan(aMig)) = Inf;
                aMig(:, col) = Inf;
            end
            
            % L & U preallocation
            xo_l = interval(eye(arow));
            xo_u = matrixPermute(xo_p, xi_a);
            
            % L & U calculation
            v.type = '()';
            u.type = '()';
            urow.type = '()';
            r.type = '()';            
            for i = 1 : arow - 1
                % subsref
                v.subs = {i, i};
                u.subs = {i, i : arow};
                
                % iterate row-wise
                for k = i + 1 : arow
                    r.subs = {k, i};
                    
                    % lower matrix
                    l = mulRevToPair(subsref(xo_u, v), subsref(xo_u, r));
                    xo_l = subsasgn(xo_l, r, l);
                    
                    % remaining columns
                    urow.subs = {k, i : arow};
                    
                    % upper matrix
                    uEnd = subsref(xo_u, urow);
                    sEnd = l .* subsref(xo_u, u);
                    xo_u = subsasgn(xo_u, urow, uEnd - sEnd);
                end
            end
            
            % force U to be upper triangular
            xo_u.infimum = triu(xo_u.infimum);
            xo_u.supremum = triu(xo_u.supremum);
        end
        
        % invert an interval matrix
        function xo_interval = inv(xi_interval)
            % input argument check
            if nargin ~= 1
                error(INTERVAL_INV, 'interval::inv - only one input argument is expected.');
            end
            
            % housekeeping
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % housekeeping
            len = length(xi_interval.infimum);
            
            % scalar
            if len <= 1
                xo_interval = rdivide(1, xi_interval);
            else
                xo_interval = mldivide(xi_interval, eye(len));
            end
        end

        % return the determinant of a matrix
        function xo_det = det(xi_interval)
             % input argument check
            if nargin ~= 1
                error(INTERVAL_DET, 'interval::det - only one input argument is expected.');
            end
            
            % housekeeping
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            
            % input argument must be square
            inSize = size(xi_interval.infimum);
            if inSize(1) ~= inSize(2)
                error(INTERVAL_DET', 'interval::det - input argument must be square.');
            end
            
            % empty interval
            if any(any(isempty(xi_interval)))
                xo_det = interval();
                return;
            end
            
            % scalar input
            if isscalar(xi_interval)
                xo_det = xi_interval;
                return;
            end
            
            % LU decomposition
            [L, U, P] = lu(xi_interval);
            xo_det = times(det(P), prod(diag (U)));
        end
        
        % solve linear system equation (xi_a*xo_x = xi_b)
        function xo_x = linsolve(xi_a, xi_b)
            % check
            if nargin ~= 2
                error(INTERVAL_LINSOLVE, 'interval::linsolve - two arguments are required.');
            end

            % type adjustment
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end
            
            % one of input arguments is scalar
            if isscalar(xi_a.infimum) || isscalar(xi_b.infimum)
                xo_x = rdivide(xi_b, xi_a);
                return;
            end
            
            % is xi_a square?
            aSize = size(xi_a);
            arow = aSize(1);
            acol = aSize(2);
            bSize = size(xi_b);
            brow = bSize(1);
            bcol = bSize(2);
            if acol ~= arow
                error(INTERVAL_LINSOLVE, 'interval::linsolve - first input argument must be a square matrix.');
            end
            
            % are xi_a & xi_b compatible?
            if arow ~= brow
                error(INTERVAL_LINSOLVE, 'interval::linsolve - input arguments dimenstions are incompatible.');
            end
            
            % housekeeping #2
            n = length(xi_a.infimum);
            m = bcol;
            
            % LU decomposition of xi_a
            [L, U, P] = lu(xi_a);
            
            % Forward substitute
            s = matrixPermute(inv(P), xi_b);
            curr.type = '()';
            prevv.type = '()';
            LRow.type = '()';
            
            % iterate over columns
            for i = 1 : m
                % iterate over rows
                for k = 2 : n
                    curr.subs = {k, i};
                    prevv.subs = {1 : k, i};
                    LRow.subs = {k, 1 : k};
                    vCol = subsref(s, prevv);
                    Lrow = subsref(L, LRow);
                    
                    Lrow.infimum(k) = -1;
                    Lrow.supremum(k) = -1;
                    
                    s = subsasgn(s, curr, -dot(Lrow, vCol));
                end
            end
            
            % Backward substitution
            z = s;
            uRow1.type = '()';
            uRow2.type = '()';
            
            % iterate over columns
            for i = 1 : m
                curr.subs = {n, i};
                uRow1.subs = {n, n};
                
                % z = s / U
                z = subsasgn (z, curr,  mulRevToPair(subsref (U, uRow1), subsref (z, curr)));
                for k = n - 1 : -1 : 1
                    curr.subs = {k, i};
                    uRow1.subs = {k, k};
                    prevv.subs = {k : n, i};
                    uRow2.subs = {k, k : n};
                    vCol = subsref(z, prevv);
                    Urow = subsref (U, uRow2);
                    Urow.infimum(1) = -1;
                    Urow.supremum(1) = -1;
                    z = subsasgn(z, curr, mulRevToPair(subsref (U, uRow1), -dot (Urow, vCol)));
                end
            end
            
            % gauss-seidel (tricky)
            aIdx.type = '()';
            bIdx.type = '()';
            zIdx.type = '()';
            aMig = mignitude(xi_a);
            aMig(isnan(aMig)) = 0;
            
            % iterate over columns
            for k = 1 : m
                zIdx.subs = {1 : n, k};
                zCol = subsref (z, zIdx);
                
                % iterate over rows
                for j = n : -1 : 1
                    zEelement = interval(zCol.infimum(j), zCol.supremum(j));
                    if isempty(zEelement) || issingle(zEelement)
                        continue;
                    end
                    
                    i = find(aMig(:, j) == max(aMig(:, j)), 1);
                    aIdx.subs = {i, 1 : n};
                    a = subsref(xi_a, aIdx);
                    if a.infimum(j) < 0 && a.supremum(j) > 0
                        continue;
                    end
                    
                    aElement = interval(a.infimum(j), a.supremum(j));
                    bIdx.subs = {i, k};
                    b = subsref(xi_b, bIdx);
                    
                    a.infimum(j) = b.infimum;
                    a.supremum(j) = b.supremum;
                    zCol.infimum(j) = -1;
                    zCol.supremum(j) = -1;
                    zEelement = mulRevToPair(aElement, -dot (a, zCol), zEelement);
                    
                    zCol.infimum(j) = zEelement.infimum;
                    z.infimum(j, k) = zEelement.infimum;
                    zCol.supremum(j) = zEelement.supremum;
                    z.supremum(j, k) = zEelement.supremum;
                end
            end
            
            xo_x = z;
        end
        
    end % matrix methods
   
    % optimization methods
    methods
        
        % fzero
        % xi_option.maxIter - maximum iterations (a good number is 200)
        % xi_option.tolX - convergence tolerance on x
        % xi_option.tolFun - convergence tolerance on function evaluation
        function xo_interval = fzero(xi_func, xi_x0, xi_gradient, xi_option)
            % input arguments
            if nargin ~= 4
                error(INTERVAL_FZERO, 'interval::fzero - four input arguments are required.');
            end

            % parmeter type assurance
            if ~isa(xi_x0, 'interval')
                xi_x0 = interval(xi_x0);
            end
            if ~isscalar(xi_x0)
                error(INTERVAL_FZERO, 'interval::fzero - initial starting point (x0) must be a scalar.');
            end
            if isempty(xi_x0)
                error(INTERVAL_FZERO, 'interval::fzero - initial starting point (x0) is empty.');
            end
            
            if ~isa(xi_func, 'function_handle') && ~ischar(xi_func)
                error(INTERVAL_FZERO, 'interval::fzero - supplied function is not a function handle.');
            end
            if ~isempty(xi_gradient) && ~isa(xi_gradient, 'function_handle') && ~ischar(xi_gradient)
                    error(INTERVAL_FZERO, 'interval::fzero - supplied gradient is not a function handle.');
            end

            % newton solution
            [lower, upper] = newton(xi_func, xi_gradient, xi_x0, 0, xi_option);

            % output
            xo_interval = interval(lower, upper);
        end
    end % optimization methods
    
    % visualization methods
    methods

        % plot one/two intervals
        % - if only one interval is given - it is drawn against it self
        %    where a line connects the infimum and a shaded rectangele filles
        %    the rectange [infimum, infimum, supremum, supremum]
        % - if two intervals are given -they are drawn such that the first
        %    interval is along the X axis and the second interval is along
        %    Y axis. as before - a line connects the infimum's and a shaded
        %    rectangle filles the intervals themself.
        function plot(xi_x, xi_y)
            % check input arguments
            if nargin > 2
                error(INTERVAL_PLOT, 'interval::plot - one or two input arguments are required.');
            end

            % maintain drawn data
            hold on;

            % draw only one interval
            if nargin == 1
                % define counterpart
                xi_y = xi_x;

                % reshape
                xi_x = interval(reshape(1 : numel(xi_y.infimum), size(xi_y.infimum)));
            end

            % plot type using GEM trick
            empty = isempty(xi_x) | isempty(xi_y);
            notInterval = issingle(xi_x) & issingle(xi_y);
            curve = xor(issingle(xi_x), issingle(xi_y)) & ~empty;
            rectangle = ~(empty | notInterval | curve);

            % non-interval's plot
            if any(any(notInterval))
                scatter(xi_x.infimum(notInterval), xi_y.infimum(notInterval), 3, 'k');
            end

            % curve plot
            if any(any(curve))
                x = [xi_x.infimum(curve), xi_x.supremum(curve)]';
                y = [xi_y.infimum(curve), xi_y.supremum(curve)]';
                plot(x, y, 'k');
            end

            % interval area plot
            if any(any(rectangle))
                % 'middle' line
                x = [xi_x.infimum(rectangle), xi_x.supremum(rectangle)];
                y = [xi_y.infimum(rectangle), xi_y.supremum(rectangle)];
                plot(x, y, 'k');

                % interval 'area'
                i = 1 : numel(x) - 1;
                xx = [x(i); x(i+1); x(i+1); x(i)];
                yy = [y(i); y(i); y(i+1); y(i+1)];
                patch(xx,yy, [0.95, 0.95, 0.95], 'FaceAlpha', 0.75);
            end
        end

        function plo3t()
        end

    end % % visualization methods

    % internal methods
    methods (Access = protected)
        
        % interval xi_x infinitly inflated by xi_eps
        function xo_y = inflate(xi_x, xi_eps)
            xo_y = next((1 + xi_eps) .* xi_x - xi_eps .* xi_x);
        end

        % iteratively iterate along an interval until it includes a given value
        % xo_x = final interval
        % xo_flag = 1 if interval refinment was a succesful
        % xi_x = initial interval
        % xi_c = value to include
        function [xo_x, xo_flag] = intervalRefine(xi_x, xi_c)
            % housekeeping
            xo_flag = 0;
            xo_x = xi_x;

            % iterative refinment
            for i = 1 : 5
                % inflate
                y = inflate(xo_x, 0.5);

                % advance
                xo_x = xi_x + mtimes(xi_c, y);

                % verify
                xo_flag = all(all(in(xo_x, y)));
                if xo_flag
                    break;
                end
            end

            % verification
            if xo_flag
                % iterative refinment
                for i = 1 : 5
                    y = xo_x;
                    [xo_x, dummy] = intersect(xi_x + mtimes(xi_c, xo_x), xo_x); %#ok<NASGU>

                    % refinment stoppage criteria
                    wid = max(abs(xo_x.infimum - y.infimum), abs(xo_x.supremum - y.supremum));
                    if max(max(wid)) <= 1e-5
                        break;
                    end
                end
            end
        end
        
        % perform fast multiplication of a permutation matrix by interval matrix
        % B = P * A (P - non numeric matrix)
        function xo_b = matrixPermute(xi_p, xi_a)
            % housekeeping
            aSize = size(xi_a);
            aRow = aSize(1);
            xo_b = interval(xi_a.infimum, xi_a.supremum);

             % row wise permutation
            for i = 1 : aRow
                row = find(xi_p(i, :) == 1, 1);
                xo_b.infimum(row, :) = xi_a.infimum(i, :);
                xo_b.supremum(row, :) = xi_a.supremum(i, :);
            end
        end
        
        % mulRevToPair (IEEE Std) function implementation in matlab
        % Notice: this method is only used when REALLY needed and is not optimized.
        %                    it should be tested and treated better. tricky function to implement.
        function [xo_negative, xo_positive] = mulRevToPair(xi_a, xi_b, xi_c)
            % handle only two input arguments
            if nargin == 2
                xi_c = interval(-Inf, Inf);
            end
            
            % type assurence
            if ~isa(xi_a, 'interval')
                xi_a = interval(xi_a);
            end
            if ~isa(xi_b, 'interval')
                xi_b = interval(xi_b);
            end
            if ~isa(xi_c, 'interval')
                xi_c = interval(xi_c);
            end
            
            % housekeepings
            aSize = size(xi_a);
            arow = aSize(1);
            acol = aSize(2);
            bSize = size(xi_b);
            brow = bSize(1);
            bcol = bSize(2);
            
            % if xi_a & xi_b are of different dimensions
            if arow ~= brow
                xi_a.infimum = ones(brow, acol) .* xi_a.infimum;
                xi_a.supremum = ones(brow, acol) .* xi_a.supremum;
                xi_b.infimum = ones(arow, bcol) .* xi_b.infimum;
                xi_b.supremum = ones(arow, bcol) .* xi_b.supremum;
            end
            if acol ~= bcol
                xi_a.infimum = ones(arow, bcol) .* xi_a.infimum;
                xi_a.supremum = ones(arow, bcol) .* xi_a.supremum;
                xi_b.infimum = ones(brow, acol) .* xi_b.infimum;
                xi_b.supremum = ones(brow, acol) .* xi_b.supremum;
            end
            
            % preallocation (empty intervals)
            xo_negative.infimum = Inf(size(xi_a.infimum));
            xo_positive.infimum = xo_negative.infimum;
            xo_negative.supremum = -Inf(size(xi_a.infimum));
            xo_positive.supremum = xo_negative.supremum;
            
            % empty interval
            empty = (xi_a.infimum == Inf) || (xi_b.infimum == Inf) ||...
                ((xi_a.infimum == 0) && (xi_a.supremum == 0) && ((xi_b.supremum < 0) || (xi_b.infimum > 0)));
            
            % seperated intervals
            seperated = ~empty  && (xi_a.infimum < 0) && (xi_a.supremum > 0) && ((xi_b.supremum < 0) || (xi_b.infimum > 0));
            
            % integrated intervals (as one)
            joined = ~seperated && ~empty;
            
            % mark seperated intervals
            xo_negative.infimum(seperated) = -Inf;
            xo_positive.supremum(seperated) = Inf;
            
            % seperated with 0 inclusion
            temp = seperated && xi_b.infimum <= 0 && xi_b.supremum >= 0;
            xo_negative.supremum(temp) = 0;
            xo_positive.infimum(temp) = 0;
            
            % seperated and to the righ of 0
            temp = seperated && xi_b.infimum > 0;
            xo_negative.supremum(temp) = rdivide(xi_b.infimum(temp), xi_a.infimum(temp));
            xo_positive.infimum(temp) = rdivide(xi_b.infimum(temp), xi_a.supremum(temp));

            % seperated and to the left of 0
            temp = seperated && xi_b.supremum < 0;
            xo_negative.supremum(temp) = rdivide(xi_b.supremum(temp), xi_a.supremum(temp));
            xo_positive.infimum(temp) = rdivide(xi_b.supremum(temp), xi_a.infimum(temp));
            
            % joined and to the right (or including) 0
            temp = joined && xi_a.infimum >= 0 && xi_b.infimum >= 0;
            xi_a.infimum(temp && xi_a.infimum == 0) = 0; % should eps becpme 0
            xi_b.infimum(temp && xi_b.infimum == 0) = 0;
            xo_negative.infimum(temp) = rdivide(xi_b.infimum(temp), xi_a.supremum(temp));
            xo_negative.supremum(temp) = rdivide(xi_b.supremum(temp), xi_a.infimum(temp));
    
            % joined and inversely include 0
            temp = joined && xi_a.supremum <= 0 && xi_b.infimum >= 0;
            xi_a.supremum(temp && xi_a.supremum == 0) = 0;
            xi_b.infimum(temp && xi_b.infimum == 0) = 0;
            xo_negative.infimum(temp) = rdivide(xi_b.supremum(temp), xi_a.supremum(temp));
            xo_negative.supremum(temp) = rdivide(xi_b.infimum(temp), xi_a.infimum(temp));
    
            temp = joined && xi_a.infimum >= 0 && xi_b.supremum <= 0;
            xi_a.infimum(temp && xi_a.infimum == 0) = 0;
            xi_b.supremum(temp && xi_b.supremum == 0) = 0;
            xo_negative.infimum(temp) = rdivide(xi_b.infimum(temp), xi_a.infimum(temp));
            xo_negative.supremum(temp) = rdivide(xi_b.supremum(temp), xi_a.supremum(temp));

            % joined and to the left (or including) 0
            temp = joined && xi_a.supremum <= 0 && xi_b.supremum <= 0;
            xi_a.supremum(temp && xi_a.supremum == 0) = 0;
            xi_b.supremum(temp && xi_b.supremum == 0) = 0;
            xo_negative.infimum(temp) = rdivide(xi_b.supremum(temp), xi_a.infimum(temp));
            xo_negative.supremum(temp) = rdivide(xi_b.infimum(temp), xi_a.supremum(temp));
            
            % joined and xi_b include 0 and xi_a is to the right of 0
            temp = joined && xi_b.infimum < 0 && xi_b.supremum > 0 && xi_a.infimum > 0;
            xo_negative.infimum(temp) = rdivide(xi_b.infimum(temp), xi_a.infimum (temp));
            xo_negative.supremum(temp) = rdivide(xi_b.supremum (temp), xi_a.infimum (temp));

            % joined and xi_b include 0 and xi_a is to the left of 0
            temp = joined && xi_b.infimum < 0 && xi_b.supremum > 0 && xi_a.supremum < 0;
            xo_negative.infimum(temp) = rdivide(xi_b.supremum (temp), xi_a.supremum (temp));
            xo_negative.supremum(temp) = rdivide(xi_b.infimum (temp), xi_a.supremum (temp));
            
            % joined and both xi_a and xi_b include 0
            temp = joined && xi_a.infimum <= 0 && xi_a.supremum >= 0 && xi_b.infimum <= 0 && xi_b.supremum >= 0;
            xo_negative.infimum(temp) = -Inf;
            xo_negative.supremum(temp) = Inf;
            
            % intersect (xo_negative || xo_positive, xi_c)
            xo_negative.infimum = max(xo_negative.infimum, xi_c.infimum);
            xo_negative.supremum = min(xo_negative.supremum, xi_c.supremum);
            xo_positive.infimum = max(xo_positive.infimum, xi_c.infimum);
            xo_positive.supremum = min(xo_positive.supremum, xi_c.supremum);

            % some applications require only te unique negative part
            if nargout == 1
                % union(negative, positive)
                xo_negative.infimum(seperated) = min(xo_negative.infimum (seperated), xo_positive.infimum (seperated));
                xo_negative.supremum(seperated) = max(xo_negative.supremum (seperated), xo_positive.supremum (seperated));
                
                % empty interval
                empty = xo_negative.infimum > xo_negative.supremum;
                xo_negative.infimum(empty) = Inf;
                xo_negative.supremum(empty) = -Inf;
                
                % output
                xo_negative = interval(xo_negative.infimum, xo_negative.supremum);
            else
                % empty negative interval
                empty = xo_negative.infimum > xo_negative.supremum;
                xo_negative.infimum(empty) = Inf;
                xo_negative.supremum(empty) = -Inf;
                
                % empty positive interval
                empty = xo_positive.infimum > xo_positive.supremum;
                xo_positive.infimum(empty) = Inf;
                xo_positive.supremum(empty) = -Inf;
                
                % output as intervals
                xo_negative = interval(xo_negative.infimum, xo_negative.supremum);
                xo_positive = interval(xo_positive.infimum, xo_positive.supremum);
                
                % seperated intervals with a joined group (it could happens in obscure situations)
                swap = seperated && isempty(xo_negative) && ~(isempty(xo_positive));
                [xo_negative.infimum(swap),...
                 xo_negative.supremum(swap),...
                 xo_positive.infimum(swap),...
                 xo_positive.supremum(swap)] = deal(xo_positive.infimum (swap),...
                                                                                                  xo_positive.supremum (swap),...
                                                                                                  xo_negative.infimum (swap),...
                                                                                                  xo_negative.supremum (swap));                                                                  
            end
        end
        
        % calculate the square value of each interval element
        function xo_interval = square(xi_interval)
            % square
            lower = mignitude(xi_interval) .* mignitude(xi_interval);
            upper = magnitude(xi_interval) .* magnitude(xi_interval);
            
            % empty interval
            empty = isempty(xi_interval);
            lower(empty) = Inf;
            upper(empty) = -Inf;
            
            % output
            xo_interval = interval(lower, upper);
        end
        
        % calculate the power of positive intervals
        function xo_interval = positivePower(xi_interval, xi_power)
             % type assurence
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end
            if ~isa(xi_power, 'interval')
                xi_power = interval(xi_power);
            end
            
            % output dimension determination
                if isscalar (xi_interval.infimum) ~= isscalar (xi_power.infimum)
                    xi_interval.infimum = ones(size(xi_power.infimum)) .* xi_interval.infimum;
                    xi_interval.supremum = ones(size(xi_power.infimum)) .* xi_interval.supremum;
                    xi_power.infimum = ones(size(xi_interval.infimum)) .* xi_power.infimum;
                    xi_power.supremum = ones(size(xi_interval.infimum)) .* xi_power.supremum;
                end

                % positive set
                xi_interval = intersect(xi_interval, interval(0, Inf));
                
                % remove '0' and empty interval power
                xi_power.infimum(xi_interval.supremum == 0) = max(0, xi_power.infimum(xi_interval.supremum == 0));
                empty = isempty(xi_power);
                xi_power.supremum(empty) = -Inf;
                xi_power.infimum(empty) = Inf;

                % power calculation
                lower = min(min(min(power(abs(xi_interval.infimum), xi_power.infimum), ...
                                                             power(abs(xi_interval.infimum), xi_power.supremum)), ...
                                                             power(xi_interval.supremum, xi_power.infimum)), ...
                                                             power(xi_interval.supremum, xi_power.supremum));
                upper = max(max(max(power(abs(xi_interval.infimum), xi_power.infimum), ...
                                                             power(abs(xi_interval.infimum), xi_power.supremum)), ...
                                                             power(xi_interval.supremum, xi_power.infimum)), ...
                                                             power(xi_interval.supremum, xi_power.supremum));

                 %  empty interval
                empty = isempty(xi_interval) | isempty(xi_power) | (xi_interval.supremum == 0 & xi_power.supremum == 0);
                lower(empty) = Inf;
                upper(empty) = -Inf;

                % 0^(+) = (+)
                upper(xi_interval.supremum == 0) = 0;

                % output
                xo_interval = interval(lower, upper);
        end
        
        % power with a positive integer exponent
        function xo_interval = powerN(xi_interval, xi_power)
             % power argument check
            if ~isnumeric(xi_power) || any(any(fix(xi_power) ~= xi_power))
                error (INTERVAL_POWERN, 'interval::powerN - exponent is not an integer.');
            end
            
             % type assurence
            if ~isa(xi_interval, 'interval')
                xi_interval = interval(xi_interval);
            end

            % output dimension determination
            if isscalar(xi_interval.infimum) ~= isscalar(xi_power)
                xi_interval.infimum = ones(size(xi_power)) .* xi_interval.infimum;
                xi_interval.supremum = ones(size(xi_power)) .* xi_interval.supremum;
                xi_power = ones(size(xi_interval.infimum)) .* xi_power;
            end

            % housekeeping
            xo_interval = xi_interval;

            % powerN(~, 0)
            pow0 = (xi_power == 0 & ~isempty(xi_interval));
            xo_interval.infimum(pow0) = 1;
            xo_interval.supremum(pow0) = 1;

            % powerN(~, 2)
            tempSub.type = '()';
            tempSub.subs = {(xi_power == 2)};
            if any(any(tempSub.subs{1}))
                xo_interval = subsasgn(xo_interval, tempSub, square(subsref(xi_interval, tempSub)));
            end

            % powerN(~, -1)
            tempSub.subs = {(xi_power == -1)};
            if any(any(tempSub.subs{1}))
                xo_interval = subsasgn (xo_interval, tempSub, 1 ./ subsref(xi_interval, tempSub));
            end

            % even powers (not zero or two)
            tempSub.subs = {(rem (xi_power, 2) == 0 & xi_power ~= 2 & xi_power ~= 0)};
            if any(any(tempSub.subs{1}))
                % mignitude & magnitude
                x_mig = mignitude(subsref(xi_interval, tempSub));
                x_mag = magnitude(subsref(xi_interval, tempSub));
                x_mig(isnan(x_mig)) = Inf;
                x_mag(isnan(x_mag)) = -Inf;

                % power
                xi_interval.infimum = subsasgn(xi_interval.infimum, tempSub, x_mig);
                xi_interval.supremum = subsasgn(xi_interval.supremum, tempSub, x_mag);

                xo_interval = subsasgn(xo_interval, tempSub, positivePower(subsref(xi_interval, tempSub), subsref(xi_power, tempSub)));
            end

            % odd powers
            tempSub.subs = {(rem (xi_power, 2) ~= 0 & xi_power ~= -1)};
            if any (any (tempSub.subs{1}))
                intervalIndex = subsref(xi_interval, tempSub);
                powerIndex = interval(subsref (xi_power, tempSub));
                xo_interval = subsasgn(xo_interval, tempSub, ...
                                                                        union(positivePower(intervalIndex, powerIndex), ...
                                                                                       -positivePower(-intervalIndex, powerIndex)));
            end

            % zero interval with positive power
            zeroInterval = xi_power > 0 & xi_interval.infimum == 0 & xi_interval.supremum == 0;
            xo_interval.infimum(zeroInterval) = -0;
            xo_interval.supremum(zeroInterval) = +0;
        end
        
        % bisect an interval such that the pivot
        % element is at pow2(middle(log2(abs (interval))))
        function [xo_lower, xo_upper] = bisect(xi_interval)
            % housekeeping
            m = zeros(size(xi_interval.infimum));
            eps64 = pow2(-1074);
            subIndex.type = '()';

            % positive interval bisection
            posInterval = xi_interval.infimum >= 0 & xi_interval.supremum > 0;
            if any(any(posInterval))
                subIndex.subs = {posInterval};
                bounded = intersect(interval(eps64, realmax), subsref(xi_interval, subIndex));
                m(posInterval) = min(realmax, pow2(middle(log2(bounded))));
            end

            % negative interval bisection
            negative = xi_interval.infimum < 0 & xi_interval.supremum <= 0;
            if (any (any (negative)))
                subIndex.subs = {negative};
                bounded = intersect (interval(-realmax, -eps64), subsref(xi_interval, subIndex));
                m(negative) = uminus(min(realmax, pow2(middle(log2(uminus(bounded))))));
            end

            % zero encompesing interval bisection, positive part
            both_signs = xi_interval.infimum < 0 & xi_interval.supremum > 0;
            both_signs_p = both_signs & xi_interval.supremum > -xi_interval.infimum;
            if any(any(both_signs_p))
                subIndex.subs = {both_signs_p};
                [bounded_n, dummy] = intersect(interval(-realmax, -eps64), subsref(xi_interval, subIndex));
                [bounded_p, dummy] = intersect(interval(eps64, realmax), subsref(xi_interval, subIndex));
                m(both_signs_p) = min(realmax, pow2(middle(log2(bounded_p)) - middle(log2(uminus(bounded_n))) - 1074));
            end

            % zero encompesing interval bisection, negative part
            both_signs_n = both_signs & xi_interval.supremum < -xi_interval.infimum;
            if any(any(both_signs_n))
                subIndex.subs = {both_signs_n};
                [bounded_n, dummy] = intersect(interval(-realmax, -eps64), subsref(xi_interval, subIndex));
                [bounded_p, dummy] = intersect(interval(eps64, realmax), subsref(xi_interval, subIndex));
                m(both_signs_n) = uminus(min(realmax, pow2(middle(log2(uminus(bounded_n))) - middle(log2(bounded_p)) - 1074)));
            end

            % empty interval
            m(isempty(xi_interval)) = -Inf;
            xo_lower = interval(xi_interval.infimum, min(m, xi_interval.supremum));
            m(isempty(xi_interval)) = Inf;
            xo_upper = interval(max (m, xi_interval.infimum), xi_interval.supremum);
        end
        
        % recrusive Newton (gradient) method
        % xi_option.maxIter - maximum iterations
        % xi_option.tolX - convergence tolerance on x
        % xi_option.tolFun - convergence tolerance on function evaluation
        function [xo_lower, xo_upper] = newton(xi_func, xi_gradient, xi_x0, xi_iteration, xi_option)
            % housekeeping
            xo_lower = [];
            xo_upper = [];

            % Try the newton step, if derivative is known
            if ~isempty(xi_gradient)
                mid = interval(middle(xi_x0));
                [a, b] = mulRevToPair(feval(xi_gradient, xi_x0), feval(xi_func, mid));
                if isempty(a) && isempty(b)
                    a = xi_x0;
                else
                    a = intersect(xi_x0, mid - a);
                    b = intersect(xi_x0, mid - b);
                    if isempty(a)
                        [a, b] = deal(b, a);
                    end
                end
            else
                a = xi_x0;
                b = interval();
            end

            % newton iteration fail, bisect the interval
            if (eq(xi_x0, a) || isempty(b)) && (a.infimum ~= a.supremum) && ~isempty(a)                
                %mid = middle(a); 
                %a = interval(a.infimum, mid);
                %b = interval(mid, a.supremum);
                %mid = pow2(middle(log2(abs(a))));
                [a, b] = bisect(a);
            elseif (b < a)
                [a, b] = deal(b, a);
            end

            % iterate along bisectioned region
            for x1 = {a, b}
                % sample
                x1 = x1{1}; %#ok<FXSET>
                f_x1 = feval(xi_func, x1);
                
                % no roots
                if  ~ismember(0, f_x1)
                    continue;
                end
                
                % very slow convrgence, minimize iteration number
                if (f_x1.infimum == -Inf && f_x1.supremum == Inf) || (width(f_x1) / max(realmin, width(x1)) < pow2(-20))
                    xi_option.maxIter = xi_option.maxIter / 1.5;
                end

                % stopage criteria (tolerance)
                if eq(x1, xi_x0) || (xi_iteration >= xi_option.maxIter) ||...
                        (width(x1) <= xi_option.tolX) || (width(f_x1) <= xi_option.tolFun)
                    [newLower, newUpper] = deal(x1.infimum, x1.supremum);
                else
                    [newLower, newUpper] = newton(xi_func, xi_gradient, x1, xi_iteration + 1, xi_option);
                end
                
                % good interval
                if ~isempty(newLower)
                    if isempty(xo_lower)
                        xo_lower = newLower;
                        xo_upper = newUpper;
                    elseif xo_upper(end) == newLower(1)
                        xo_upper(end) = newUpper(1); %#ok<AGROW>
                        xo_lower = [xo_lower; newLower(2 : end, 1)]; %#ok<AGROW>
                        xo_upper = [xo_upper; newUpper(2 : end, 1)]; %#ok<AGROW>
                    else
                        xo_lower = [xo_lower; newLower]; %#ok<AGROW>
                        xo_upper = [xo_upper; newUpper]; %#ok<AGROW>
                    end
                end
            end
        end
        
    end % internal

end  % interval class
