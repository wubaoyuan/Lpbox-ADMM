function rpm_demo(cmd_str)
% Robust point matching demo. 
% 
% Purpose: 
%     1. A small GUI.
%     2. load a few non-rigid point matching examples.
%     3. run RPM and ICP.
%
% Usage: [] = rpm_demo;
%     1. Click the buttons in figure(2) to run the examples.
%     2. The matching process and results are displayed in figure(1).
%     3. You can resize figure(1), if it is too small.
%
% Notes: There are a total of 5 examples included to demonstrate 
%        the non-rigid point matching algorithm.
%
% Input
%   cmd_str  -  'init' | 'load_ex1' | ...
%
% History
%   create  -  Anand Rangarajan	(anand@noodle.med.yale.edu), 04-27-2000
%   modify  -  Feng Zhou (zhfe99@gmail.com), 01-05-2012

% global variables
global x y frac T_init T_finalfac disp_flag m_method lam1 lam2 perTmaxit
global c d m

% init the command
if nargin < 1
    figure(1); delete(1); % Clean up previously opened figures.
    figure(2); delete(2); 
    cmd_str = 'init';
end;

if strcmp (cmd_str, 'init')
    % init the figure windows.
    h1 = figure(1); 
    set(gcf,'position', [10 110 600 500],'menubar','none');
    set(gcf, 'color', [0 0 0]);
    
    h2 = figure(2); 
    set(gcf,'position', [10 10 600 40], 'menubar','none'); 
    
    % Init the command buttons.
    h_fig = h2;
    col_ex0 = 10; 
    row_ex0 = 5;
    col_ex1 = 10+100;
    row_ex1 = row_ex0;
    col_ex = 30;
    row_ex = 30;
    h_jnk = uicontrol(h_fig, 'style', 'text', ...
                      'position', [col_ex0 row_ex0 100 row_ex], ...
                      'string', 'Load Data:','fontsize',15);
    h_ex1 = uicontrol(h_fig, 'style', 'pushbutton', ...
                      'callback', 'rpm_demo(''load_ex1'');',...
                      'position', [col_ex1 row_ex0 col_ex row_ex], ...
                      'string', '1');
    h_ex2 = uicontrol(h_fig, 'style', 'pushbutton', ...
                      'callback', 'rpm_demo(''load_ex2'');',...
                      'position', [col_ex1+(col_ex*1) row_ex0 col_ex row_ex], ...
                      'string', '2');
    h_ex3 = uicontrol(h_fig, 'style', 'pushbutton', ...
                      'callback', 'rpm_demo(''load_ex3'');',...
                      'position', [col_ex1+(col_ex*2) row_ex0 col_ex row_ex], ...
                      'string', '3');
    h_ex4 = uicontrol(h_fig, 'style', 'pushbutton', ...
                      'callback', 'rpm_demo(''load_ex4'');',...
                      'position', [col_ex1+(col_ex*3) row_ex0 col_ex row_ex], ...
                      'string', '4');
    h_ex5 = uicontrol(h_fig, 'style', 'pushbutton', ...
                      'callback', 'rpm_demo(''load_ex5'');',...
                      'position', [col_ex1+(col_ex*4) row_ex0 col_ex row_ex], ...
                      'string', '5');

    col_run0 = 10+100+30*5+20; 
    row_run0 = 5;
    col_run1 = col_run0+100;
    row_run1 = row_run0;
    col_run  = 80;
    row_run  = 30;
    h_jnk = uicontrol(h_fig, 'style', 'text', ...
                      'position', [col_run0 row_run0 100 row_run], ...
                      'string', 'Run:','fontsize',15);
    h_run_rpm = uicontrol(h_fig, 'style', 'pushbutton', ...
                          'callback', 'rpm_demo(''run_rpm'');',...
                          'position', [col_run1 row_run0 col_run row_run], ...
                           'string', 'RPM');
    h_run_icp = uicontrol(h_fig, 'style', 'pushbutton', ...
                          'callback', 'rpm_demo(''run_icp'');',...
                          'position', [col_run1+(col_run*1) row_run0 col_run row_run], ...
                          'string', 'ICP');

    % Init all the parameters:
    load demodata_ex1; 
    x = x1; 
    y = y1; 
    frac       = 1;
    T_init     = 0.5;
    T_finalfac = 500;
    disp_flag  = 1;
    m_method   = 'mix-rpm';
    lam1       = 1;
    lam2       = 0.01;
    perTmaxit  = 3;
    
% Load data
elseif strcmp(cmd_str, 'load_ex1')
    rpm_demo('reset_all');
    load demodata_ex1;
    x = x1; y = y1; 
    demo_disp;
    
elseif strcmp(cmd_str, 'load_ex2') 
    load demodata_ex2;
    x = x2; y = y2; 
    demo_disp;
    
elseif strcmp (cmd_str, 'load_ex3') 
    load demodata_ex3;
    x = x3; y = y3; 
    demo_disp;
    
elseif strcmp (cmd_str, 'load_ex4') 
    load demodata_ex4;
    x = x1 * 2; y = y2a * 2; 
    demo_disp;
    
elseif strcmp (cmd_str, 'load_ex5') 
    load demodata_ex5;
    x = x1; y = y2a; 
    demo_disp;
  
% run
elseif strcmp(cmd_str, 'run_rpm')
    [c, d, m] = cMIX(x, y, frac, T_init, T_finalfac);
    disp('RPM point matching done ...');
  
elseif strcmp(cmd_str, 'run_icp')
    [c, d, m] = cMIX(x, y, frac, T_init, T_finalfac, 1, 'icp3');
    disp('ICP point matching done ...');
end

%%%%%%%%%%%%%%%%%%
function demo_disp
% Dispaly two point sets.

global x y 
global axis_save 

fig(1); clf; 
cplot(x,y); hold on; 

jnk  = axis; fac = 10;
xmin = jnk(1); xmax = jnk(2); ymin = jnk(3); ymax = jnk(4);
xmin = xmin - (xmax-xmin)/fac;
xmax = xmax + (xmax-xmin)/fac;
ymin = ymin - (ymax-ymin)/fac;
ymax = ymax + (ymax-ymin)/fac; axis_save = [xmin xmax ymin ymax];

[jnk,itmp] = max(x(:,2)); 
pttmp = x(itmp,:); 
htmp = text (pttmp(1,1), pttmp(1,2)+(ymax-ymin)/fac/2, 'Template Point Set'); 
set(htmp, 'color', 'g', 'fontsize', 15);
[jnk,itmp] = max(y(:,2)); 
pttmp = y(itmp,:); 
htmp = text (pttmp(1,1), pttmp(1,2)+(ymax-ymin)/fac/2, 'Data Point Set'); 
set(htmp, 'color', 'r', 'fontsize', 15);

axis(axis_save); 
axis('off');
