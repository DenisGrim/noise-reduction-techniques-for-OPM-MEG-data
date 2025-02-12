function set_axis(h)

if nargin < 1
    h = gca;
end

set(h, 'box', 'off', 'tickdir', 'out');
