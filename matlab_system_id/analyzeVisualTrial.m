function result = analyzeVisualTrial(folder, opts)
% One-command visual acceleration -> COM acceleration identification.
arguments
    folder (1,1) string = ""
    opts.Window (1,2) double = [10 90]
    opts.SampleTime (1,1) double {mustBePositive} = 0.01
    opts.ForceColumn (1,1) string = "sensor_Fy_N"
    opts.ForceSign (1,1) double {mustBeMember(opts.ForceSign,[-1 1])} = 1
    opts.MassKg (1,1) double = NaN
end
if folder == ""
    selected = uigetdir(pwd,'Select trial folder');
    if isequal(selected,0), result = []; return; end
    folder = string(selected);
end
assert(exist('procest','file') == 2,'System Identification Toolbox is required.');
[~,stem] = fileparts(folder);
q = readtable(fullfile(folder,stem+"_quest.csv"));
f = readtable(fullfile(folder,stem+"_combalance.csv"));
mass = opts.MassKg;
if isnan(mass), mass = q.body_mass_kg(1); end
assert(isfinite(mass) && mass > 0,'Invalid body mass.');
assert(opts.Window(2) > opts.Window(1),'Invalid analysis window.');
qt = q.pc_mapped_time_s; qe = q.stimulus_elapsed_s;
ft = f.pc_estimated_sample_time_s;
assert(all(isfinite([qt;qe;ft])),'Nonfinite timestamps.');
assert(all(diff(qt)>0) && all(diff(qe)>0) && all(diff(ft)>0), ...
    'Timestamps must be strictly increasing; inspect synchronization.');
if opts.Window(2)>qe(end) && opts.Window(2)-qe(end)<0.1
    opts.Window(2) = qe(end);
end
assert(opts.Window(1)>=qe(1) && opts.Window(2)<=qe(end), ...
    'Requested window is outside the stimulus recording. Specify Window explicitly.');
bounds = interp1(qe,qt,opts.Window);
assert(bounds(1)>=ft(1) && bounds(2)<=ft(end),'Force recording does not cover window.');
assert(max(diff(qt))<0.1 && max(diff(ft))<0.1,'Recording gap >= 100 ms.');
t = (bounds(1):opts.SampleTime:bounds(2)-opts.SampleTime)';
u = interp1(qt,q.ap_acceleration_m_s2,t);
y = opts.ForceSign*interp1(ft,f.(opts.ForceColumn),t)/mass;
assert(all(isfinite([u;y])),'Missing or nonfinite signal values.');
assert(std(u)>1e-6,'No usable visual acceleration input.');
% Contiguous holdout. Estimate offsets from training samples only.
n = floor(numel(t)*0.6);
assert(n>=100 && numel(t)-n>=100,'Insufficient samples.');
offset = [mean(u(1:n)) mean(y(1:n))];
u = u-offset(1); y = y-offset(2);
data = iddata(y,u,opts.SampleTime,'TimeUnit','seconds');
data.InputName = {'visual acceleration'};
data.OutputName = {'COM acceleration'};
train = data(1:n); validation = data(n+1:end);
initial = idproc('P1');
initial.Structure.Kp.Value = -0.05;
initial.Structure.Kp.Minimum = -100;
initial.Structure.Kp.Maximum = -1e-8;
initial.Structure.Tp1.Value = 0.25;
initial.Structure.Tp1.Minimum = 0.001;
initial.Structure.Tp1.Maximum = 1000;
model = procest(train,initial,procestOptions('Display','off'));
[predicted,fit] = compare(validation,model);
yp = predicted.OutputData;
zeroRMSE = sqrt(mean(validation.OutputData.^2));
modelRMSE = sqrt(mean((validation.OutputData-yp).^2));
status = "REVIEW_REQUIRED";
if fit>50 && modelRMSE<zeroRMSE && model.Tp1>0.0011 && model.Tp1<990
    status = "PROVISIONAL_FIT";
end
result = struct('Model',model,'Data',data,'Validation',validation, ...
    'Tv',model.Tp1,'Gamma',-model.Kp,'FitPercent',fit, ...
    'ModelRMSE',modelRMSE,'ZeroRMSE',zeroRMSE,'Status',status, ...
    'Options',opts,'Offsets',offset,'Trial',stem);
out = string(tempname(folder))+"_system_id";
mkdir(out);
save(fullfile(out,'result.mat'),'result');
writetable(table(model.Tp1,-model.Kp,fit,modelRMSE,zeroRMSE,status, ...
    'VariableNames',{'Tv_s','Gamma','ValidationFitPercent','ModelRMSE', ...
    'ZeroBaselineRMSE','Status'}),fullfile(out,'summary.csv'));
writetable(table(t-t(1),u,y,'VariableNames', ...
    {'Time_s','VisualAcceleration','COMAcceleration'}),fullfile(out,'aligned.csv'));
fig = figure('Name','Held-out validation');
compare(validation,model); grid on;
exportgraphics(fig,fullfile(out,'validation.png'));
fig = figure('Name','Validation residual diagnostics');
resid(validation,model);
exportgraphics(fig,fullfile(out,'residuals.png'));
result.OutputFolder = out;
fprintf('Tv=%.4g s, gamma=%.4g, validation fit=%.2f%%, %s\n%s\n', ...
    result.Tv,result.Gamma,fit,status,out);
end
