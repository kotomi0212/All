function auditSonobeFigure3
% Read the saved SLX, check real block settings and wiring, then check Eq.9.
folder=fileparts(mfilename('fullpath'));
n='sonobe2025_figure3'; load_system(fullfile(folder,[n '.slx']));
check('Visual dynamics','BlockType','TransferFcn');
check('Visual dynamics','Numerator','visual_sign*gamma_v');
check('Visual dynamics','Denominator','[Tv 1]');
check('Body dynamics','BlockType','TransferFcn');
check('Body dynamics','Numerator','1');
check('Body dynamics','Denominator','[Jb/lb+mb*lb 0 -mb*g]');
check('Velocity feedback','BlockType','Derivative');
check('Velocity','BlockType','Integrator');
check('Position error','Inputs','-+');
check('Negative feedback','Inputs','--');
check('Torque balance','Inputs','-+');
check('Position gain','Gain','kp'); check('Velocity gain','Gain','kd');
check('Inertial force','Gain','mb*lb');
check('Reaction delay','BlockType','TransportDelay');
check('Reaction delay','DelayTime','td');
edges={ 'Velocity',1,'Support acceleration'; 'Visual dynamics',1,'Velocity'; ...
    'Position error',1,'Visual dynamics'; 'Position error',2,'Body dynamics'; ...
    'Position gain',1,'Position error'; 'Velocity feedback',1,'Body dynamics'; ...
    'Velocity gain',1,'Velocity feedback'; 'Negative feedback',1,'Position gain'; ...
    'Negative feedback',2,'Velocity gain'; 'Reaction delay',1,'Negative feedback'; ...
    'Inertial force',1,'Support acceleration'; 'Torque balance',1,'Inertial force'; ...
    'Torque balance',2,'Reaction delay'; 'Body dynamics',1,'Torque balance'};
for k=1:size(edges,1)
    ports=get_param([n '/' edges{k,1}],'PortHandles');
    line=get_param(ports.Inport(edges{k,2}),'Line');
    source=get_param(line,'SrcBlockHandle');
    assert(strcmp(getfullname(source),[n '/' edges{k,3}]),'Unexpected wiring.');
end
w=get_param(n,'ModelWorkspace');
v=@(key)getVariable(w,key);
assert(v('visual_sign')==-1,'Equation 9 check requires visual_sign=-1.');
s=1i*2*pi*0.2; A=v('Jb')/v('lb')+v('mb')*v('lb');
delay=exp(-v('td')*s);
expected=(-v('mb')*v('lb')*v('Tv')*s^2-v('mb')*v('lb')*s ...
    -v('gamma_v')*v('kp')*delay)/(s*(v('Tv')*s+1)* ...
    (A*s^2+(v('kd')*s+v('kp'))*delay-v('mb')*v('g')));
in=Simulink.SimulationInput(n); in=in.setModelParameter('StopTime','100');
out=sim(in); y=out.COM_displacement;
keep=y.Time>=60; t=y.Time(keep); yy=y.Data(keep);
% Fit sin/cos plus a constant to remove the input-integrator DC component.
B=[sin(imag(s)*t) cos(imag(s)*t) ones(size(t))]; c=B\yy;
observed=(c(1)+1i*c(2))/0.02;
relativeError=abs(observed-expected)/abs(expected);
assert(relativeError<0.02,'Response differs from Eq.9 by >=2%%.');
fprintf('PASS: saved block types, settings and 14 connections.\n');
fprintf('PASS: Eq.9 at 0.20 Hz, relative complex response error %.5g%%.\n',100*relativeError);
fprintf('Parameters are illustrative, not a reproduction of a fitted subject.\n');
    function check(block,property,value)
        assert(strcmp(get_param([n '/' block],property),value), ...
            'Unexpected setting: %s / %s',block,property);
    end
end
