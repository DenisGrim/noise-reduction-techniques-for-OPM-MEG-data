function set_font(fontsize)
% SET FONTSIZE OF FIGURES
% 
% > set_font(fontsize)

if nargin < 1 | isempty(fontsize)  % Default MATLAB Setting
    set(0, 'defaultTextFontSize', 10);
    set(0, 'defaultTextFontName', 'helovetica');
    set(0, 'defaultAxesFontSize', 10);
    set(0, 'defaultAxesFontName', 'helovetica');
else

    set(0, 'defaultTextFontSize', fontsize);
    set(0, 'defaultTextFontName', 'Arial');
    set(0, 'defaultAxesFontSize', fontsize);
    set(0, 'defaultAxesFontName', 'Arial');
end



