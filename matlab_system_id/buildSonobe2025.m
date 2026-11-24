function file = buildSonobe2025
% Build an inspectable implementation of Sonobe & Miura 2025, Eqs. 6-8.
name = 'sonobe2025_explorer';
file = fullfile(fileparts(mfilename('fullpath')),[name '.slx']);
assert(~isfile(file),'Model already exists. Open the existing SLX to preserve edits.');
load_system('simulink');
new_system(name);
w = get_param(name,'ModelWorkspace');
% Illustrative parameters, NOT identified subject parameters.
values = struct('mb',60,'lb',0.9,'Jb',8,'g',9.80665, ...
    'Kp',10.53,'Kd',3.97,'td',0.155,'Tv',1,'gamma',0.02);
fields = fieldnames(values);
for i=1:numel(fields), assignin(w,fields{i},values.(fields{i})); end
assignin(w,'A',values.Jb/values.lb+values.mb*values.lb);
% State realization exposes position and velocity without numerical derivative.
add('Sources/Sine Wave','Support acceleration',[30 200 100 235], ...
    'Amplitude','0.02','Frequency','2*pi*0.2');
add('Continuous/Integrator','Support velocity',[140 70 180 105]);
add('Continuous/Transfer Fcn','Visual target Eq8',[230 65 365 110], ...
    'Numerator','-gamma','Denominator','[Tv 1]');
add('Math Operations/Sum','Position error',[420 135 445 170],'Inputs','+-');
add('Math Operations/Gain','Position gain',[485 135 555 170],'Gain','-Kp*A');
add('Math Operations/Gain','Velocity gain',[485 280 555 315],'Gain','-Kd*A');
add('Math Operations/Sum','Controller torque',[595 175 620 215],'Inputs','++');
add('Continuous/Transport Delay','Reaction delay',[655 175 745 215],'DelayTime','td');
add('Math Operations/Gain','Inertial forcing',[155 200 245 235],'Gain','-mb*lb');
add('Math Operations/Gain','Gravity',[650 350 730 385],'Gain','mb*g');
add('Math Operations/Sum','Net torque',[795 200 820 245],'Inputs','+++');
add('Math Operations/Gain','Inverse inertia',[855 200 925 240],'Gain','1/A');
add('Continuous/Integrator','COM velocity',[965 200 1000 240]);
add('Continuous/Integrator','COM displacement',[1050 200 1085 240]);
add('Sinks/Scope','COM scope',[1170 195 1220 235]);
add('Sinks/To Workspace','Save displacement',[1150 280 1270 315], ...
    'VariableName','COM_displacement','SaveFormat','Timeseries');
wire('Support acceleration/1','Support velocity/1');
wire('Support velocity/1','Visual target Eq8/1');
wire('Visual target Eq8/1','Position error/2');
wire('COM displacement/1','Position error/1');
wire('Position error/1','Position gain/1');
wire('Position gain/1','Controller torque/1');
wire('COM velocity/1','Velocity gain/1');
wire('Velocity gain/1','Controller torque/2');
wire('Controller torque/1','Reaction delay/1');
wire('Reaction delay/1','Net torque/1');
wire('Support acceleration/1','Inertial forcing/1');
wire('Inertial forcing/1','Net torque/2');
wire('COM displacement/1','Gravity/1');
wire('Gravity/1','Net torque/3');
wire('Net torque/1','Inverse inertia/1');
wire('Inverse inertia/1','COM velocity/1');
wire('COM velocity/1','COM displacement/1');
wire('COM displacement/1','COM scope/1');
wire('COM displacement/1','Save displacement/1');
Simulink.Annotation(name,sprintf(['Sonobe & Miura 2025: Eqs.6-8 / Fig.3\n' ...
    'Support-surface acceleration -> relative COM displacement.\n' ...
    'Visual target uses -gamma/(Tv*s+1), following Eq.8 and Eq.9.\n' ...
    'Fig.3 appears to show the opposite visual-path sign.\n' ...
    'Illustrative parameters only. Edit Model Explorer > Model Workspace.\n' ...
    'A=Jb/lb+mb*lb; recompute A when changing body parameters.\n' ...
    'VR-only input requires a separate adaptation; do not feed VR acceleration here.']));
set_param(name,'StopTime','40','Solver','ode45','MaxStep','0.005');
save_system(name,file);
fprintf('Created: %s\n',file);
    function add(lib,label,pos,varargin)
        add_block(['simulink/' lib],[name '/' label],'Position',pos,varargin{:});
    end
    function wire(src,dst)
        add_line(name,src,dst,'autorouting','on');
    end
end
