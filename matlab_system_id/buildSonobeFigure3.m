function file = buildSonobeFigure3
% Paper-style layout, with explicit routing and editable symbolic parameters.
name = 'sonobe2025_figure3';
folder = fileparts(mfilename('fullpath'));
file = fullfile(folder,[name '.slx']);
assert(~isfile(file),'Existing model preserved. Open the SLX instead.');
load_system('simulink'); new_system(name);
w = get_param(name,'ModelWorkspace');
p = struct('mb',60,'lb',0.9,'Jb',8,'g',9.80665, ...
    'kp',10.53*(8/0.9+60*0.9),'kd',3.97*(8/0.9+60*0.9), ...
    'td',0.155,'Tv',1,'gamma_v',0.02,'visual_sign',-1);
for f=string(fieldnames(p))', assignin(w,f,p.(f)); end
add('Sources/Sine Wave','Support acceleration',[45 460 95 500], ...
    'Amplitude','0.02','Frequency','2*pi*0.2','ShowName','off');
add('Continuous/Integrator','Velocity',[175 120 220 175],'ShowName','off');
add('Continuous/Transfer Fcn','Visual dynamics',[400 115 565 180], ...
    'Numerator','visual_sign*gamma_v','Denominator','[Tv 1]','ShowName','off');
icon('Visual dynamics','text(0.5,0.68,''visual_sign * gamma_v''); text(0.5,0.50,''--------------------''); text(0.5,0.30,''Tv s + 1'');');
add('Math Operations/Sum','Position error',[810 240 840 270], ...
    'Inputs','-+','IconShape','round','Orientation','left','ShowName','off');
add('Math Operations/Gain','Position gain',[635 230 700 280], ...
    'Gain','kp','Orientation','left','ShowName','off');
add('Continuous/Derivative','Velocity feedback',[735 335 790 385], ...
    'Orientation','left','ShowName','off');
icon('Velocity feedback','text(0.5,0.5,''s'');');
add('Math Operations/Gain','Velocity gain',[635 335 700 385], ...
    'Gain','kd','Orientation','left','ShowName','off');
add('Math Operations/Sum','Negative feedback',[555 290 585 320], ...
    'Inputs','--','IconShape','round','Orientation','left','ShowName','off');
add('Continuous/Transport Delay','Reaction delay',[365 280 480 330], ...
    'DelayTime','td','Orientation','left','ShowName','off');
icon('Reaction delay','text(0.5,0.5,''exp(-td s)'');');
add('Math Operations/Gain','Inertial force',[175 455 250 505], ...
    'Gain','mb*lb','ShowName','off');
add('Math Operations/Sum','Torque balance',[300 465 330 495], ...
    'Inputs','-+','IconShape','round','ShowName','off');
add('Continuous/Transfer Fcn','Body dynamics',[430 445 715 515], ...
    'Numerator','1','Denominator','[Jb/lb+mb*lb 0 -mb*g]','ShowName','off');
icon('Body dynamics','text(0.5,0.73,''1''); text(0.5,0.52,''--------------------------------''); text(0.5,0.28,''(Jb/lb + mb lb) s^2 - mb g'');');
add('Sinks/Scope','COM scope',[965 460 1010 500]);
add('Sinks/To Workspace','Save COM',[960 570 1055 605], ...
    'VariableName','COM_displacement','SaveFormat','Timeseries');
wire('Support acceleration',1,'Inertial force',1,[],'support acceleration [m/s^2]');
wire('Support acceleration',1,'Velocity',1,[130 480;130 148],'');
wire('Velocity',1,'Visual dynamics',1,[],'support velocity [m/s]');
wire('Visual dynamics',1,'Position error',1,[825 148],'target COM [m]');
wire('Position error',1,'Position gain',1,[],'');
wire('Position gain',1,'Negative feedback',1,[570 255],'');
wire('Velocity feedback',1,'Velocity gain',1,[],'COM velocity');
wire('Velocity gain',1,'Negative feedback',2,[570 360],'');
wire('Negative feedback',1,'Reaction delay',1,[],'');
wire('Reaction delay',1,'Torque balance',2,[315 305],'ankle torque Na');
wire('Inertial force',1,'Torque balance',1,[],'');
wire('Torque balance',1,'Body dynamics',1,[],'');
wire('Body dynamics',1,'COM scope',1,[],'COM displacement xb [m]');
wire('Body dynamics',1,'Position error',2,[890 480;825 480],'');
wire('Body dynamics',1,'Velocity feedback',1,[890 480;890 360],'');
wire('Body dynamics',1,'Save COM',1,[920 480;920 588],'');
note('VISUAL: equilibrium position shift',[270 65 680 95],'blue');
note('DELAYED STATE FEEDBACK',[380 210 720 230],'red');
note('BODY: single inverted pendulum',[440 530 760 550],'black');
note(sprintf(['Sonobe & Miura (2025), Fig.3 / Eqs.6-8 | EXAMPLE parameters\n' ...
    'Edit: Model Explorer > Model Workspace. Run, then open COM scope.\n' ...
    'visual_sign=-1 follows Eqs.8-9; +1 follows the printed Fig.3 sign.\n' ...
    'Input is MOVING SUPPORT acceleration. VR-only data are not connected.']), ...
    [70 640 1040 735],'black');
set_param(name,'StopTime','40','Solver','ode45','MaxStep','0.005', ...
    'Location',[60 60 1250 900]);
save_system(name,file);
open_system(name); set_param(name,'ZoomFactor','FitSystem');
print(['-s' name],'-dpng','-r130',fullfile(folder,'sonobe2025_figure3_preview.png'));
fprintf('Created %s\n',file);
    function add(lib,label,pos,varargin)
        add_block(['simulink/' lib],[name '/' label],'Position',pos, ...
            'FontSize',12,varargin{:});
    end
    function icon(label,display)
        m=Simulink.Mask.create([name '/' label]); m.IconUnits='normalized';
        m.Display=strrep(display,''');',''',''horizontalAlignment'',''center'');');
    end
    function wire(src,sp,dst,dp,via,label)
        a=get_param([name '/' src],'PortHandles'); b=get_param([name '/' dst],'PortHandles');
        start=get_param(a.Outport(sp),'Position'); finish=get_param(b.Inport(dp),'Position');
        h=add_line(name,[start;via;finish]);
        if ~isempty(label), set_param(h,'Name',label); end
    end
    function note(txt,pos,color)
        a=Simulink.Annotation(name,txt); a.Position=pos;
        a.ForegroundColor=color; a.FontSize=12;
    end
end
