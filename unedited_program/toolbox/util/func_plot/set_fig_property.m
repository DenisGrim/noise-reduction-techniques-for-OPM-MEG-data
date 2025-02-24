function set_fig_property(width_inch, height_inch, tfont_size, afont_size, alinew, llinew)
% Set figure property
%
% - Input
%  width_inch : Width of figure [inch]
%  height_inch : Height of figure [inch]
%
% Y. Takeda 2018-10-02
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

if nargin<=5, llinew=2; end 
if nargin<=4, alinew=1.5; end 
if nargin<=3, afont_size=30; end 
if nargin<=2, tfont_size=25; end 

dpi = 300;
set(0, 'defaultFigureUnit', 'Pixel')
scsz = get(0, 'ScreenSize');
set(0, 'defaultFigurePosition', [1 scsz(4)-height_inch*dpi width_inch*dpi height_inch*dpi]);
set(0, 'defaultFigureColor', [1 1 1]);
set(0, 'defaultTextFontSize', tfont_size);
set(0, 'defaultAxesFontSize', afont_size);
set(0, 'defaultAxesLineWidth', alinew);
set(0, 'defaultLineLineWidth', llinew);