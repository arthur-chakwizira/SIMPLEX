%load spin transition history

close all

s_fn = "..\output\stat_calib.bin";

if isfile(s_fn)
fileID = fopen(s_fn);
Npart = fread(fileID, [1,1], 'int32');
T =  fread(fileID, [1,1], 'single');
Nt = fread(fileID,[1,1], 'int32');
S = fread(fileID, [Npart*Nt, 1], 'int32');
fclose(fileID);
else
   error("File " + s_fn + " not found.") 
end

%in St, 1 = intra, 0 = extra
St = zeros(Npart, Nt);
for c = 1:Npart
   idx = ((c-1)*Nt+1):c*Nt;
   St(c,:) = S(idx);
end

%separate intra from extra
intra = zeros(Nt, 1);
extra = zeros(Nt, 1);
for c = 1:Nt
   tmp_s = St(:,c);
   intra(c) = sum(tmp_s==1);
   extra(c) = sum(tmp_s == 0);
end

%plot
time = linspace(0, T, Nt);
figure('Color', 'w');
plot(time*1e3, (intra), 'LineWidth', 2)
hold on
plot(time*1e3, extra, 'LineWidth', 2)
xlabel("time [ms]")
ylabel("N")
legend(["Intra", "Extra"])
ylim([0 Npart])



%do some fits
where = 1:100; %crop time

%fit N(t) = N(0)*exp(-k*t)
y = log(intra(where));
x = time(where);
p = polyfit(x,y, 1);
disp(p)

figure;
plot(x,y, 'LineWidth', 2);
hold on;
plot(x, polyval(p, x), '--', 'LineWidth', 2);
legend(["Data", "Fit"])
movegui center


k = abs(p(1));
V = 2.3896e+06;
S = 1.6810e+06;
S_V =  1e6*S/V;
% disp(" Surface to volume = " + num2str(S_V))
S_V_voxelised = 1.2273e+06;
kappa = 0.0135e-3;
k_theory = S_V*kappa;
k_theory_voxelised = S_V_voxelised*kappa;
disp("k sim = " + num2str(k))
disp("k theory = " + num2str(k_theory))
disp("k theory voxelised = " + num2str(k_theory_voxelised))
Vtot =  45320796;
Vin =  2.3490e+06;
Vout = Vtot-Vin;
figure;
plot(intra/Vin)
hold on
plot(extra/Vout)
disp("Difference theory = " + num2str(100*(k - k_theory)/k_theory) + "%")
disp("Difference voxelised = " + num2str(100*(k - k_theory_voxelised)/k_theory_voxelised) + "%")


% do Lizzie's analysis
mfs = fit_model(intra(where)/Vin, time(where)')
hold on
plot(time(where)*1e3, mfs.P, 'k--')

function mfs = fit_model(P, time)

t_ub      = [1e6   1e6  500]; %f(t=0), f(t=inf), k
t_lb      = [0   0   0];

t_0 = [0.6        0.5    5];

fun = @(x)get_target(x); %objective

options = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 1e6,  'MaxIterations', 1e6,...
    'Display', 'off');
t = lsqnonlin(fun,t_0,t_lb, t_ub, options);

% Objective
    function target = get_target(t)
        f_0 = t(1); f_inf = t(2); 
        k = t(3);
        P_pred = f_inf + (f_0 - f_inf)*exp(-time*k);
        target = (P)-(P_pred);
    end

% Evaluate final results
mfs.P = get_target(t)+ P;
mfs.params = t;
mfs.par_names = {'f(t=0)', 'f(t = inf)', 'k [/s]'};

end
