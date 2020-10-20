%% prelims
clear all 
clc 
close all 


%% read material database
opts = spreadsheetImportOptions("NumVariables", 13);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A3:M11";

% Specify column names and types
opts.VariableNames = ["names", "price", "density", "youngsModulus", "poissonsRatio", "yieldStrength", "tensileStrength", "productionEnergy", "productionCO2", "disposalEnergy", "disposalCO2", "eolEnergy", "eolCO2"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "names", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "names", "EmptyFieldRule", "auto");

% Import the data
materialsData = readtable("/Users/khashi03/Desktop/Doktorand_/materialData.xlsx", opts, "UseExcel", false);

%clear opts; 

%% initiate model of the panel
%model= initModel;
  
%% select a material


%% determine height to fit constraint and the resulting mass
%INSERT CODE HERE!
thickness = []; 
vikt = [];
cost = []; 
energy_consumption= [];
CO2_footprint= [];
[i,j] = size(materialsData); 

for m = 1:i 
    material=table2struct(materialsData(m,:));
    model= initModel(material);
    model = computeOptimalVariable(model,material);
    thickness(m) = model.h; 
    vikt(m)=computeMass(model);
    NAMN(m) = material.names;
    cost(m) = material.price * vikt(m);
    [energy_consumption(m), CO2_footprint(m)] = Product_Environmental_Footprint(material,vikt(m));
end 
%figure(3)
%scatter3(Ip,cost,prod_energy,'filled');
[fitresult, gof] = createFit(cost, CO2_footprint, energy_consumption); 


figure(2), clf
subplot(2,1,1)
bargraph(struct('title','Thickness','xlabel',{{''}},'ylabel','[m]','legend',NAMN,'values',[thickness]))
subplot(2,1,2)
bargraph(struct('title','Mass','xlabel',{{''}},'ylabel','[kg]','legend',NAMN,'values',[vikt])) 



for m = 1:i
    
    material=table2struct(materialsData(m,:));
    material.names=materialsData(m,1);
    names=table2struct(material.names); % this is only for the legend in plot
    model(m).materialNamn = (materialsData(m,1));
 
    model(m).F=1e5;
    model(m).delta=1e-3;
    model(m).L=1;
    model(m).b=1;
    model(m).h = sqrt((3/2) .* (model(m).F .* model(m).L)/(material.tensileStrength .*model(m).b));
    model(m).E= material.youngsModulus;
    model(m).rho=material.density;
    model(m).mass = model(m).rho .* (model(m).b .* model(m).L .* model(m).h);
    
    model(m).Ip_energy=material.productionEnergy*model(m).mass; %[J] In Production
    model(m).Ip_CO_2 = material.productionCO2 * model(m).mass; %[CO2] In Production
    driveDistTotal=1; %[km] %total driving distance
    
    model(m).Iu_energy = material.disposalEnergy * model(m).mass; %[J]
    model(m).Iu_km_kg= material.disposalCO2 * model(m).mass; %[co2/km/kg]
    model(m).Iu_CO_2=model(m).Iu_km_kg * driveDistTotal * model(m).mass; %[co2/km]
    model(m).Ie_energy= material.eolEnergy .* model(m).mass; %[J]
    model(m).Ie_CO_2 = material.eolCO2 * model(m).mass; %[J]
    model(m).prise = material.price * model(m).mass; %[SEK]
    
    end
  

% %% production phase
% Ip_energy=material.productionEnergy*mass; %[J] In Production
% Ip_CO_2 = material.productionCO2 * mass; %[CO2] In Production
% 
% 
% %% use phase
% driveDistTotal=1; %[km] %total driving distance
% 
% Iu_energy = material.disposalEnergy * mass; %[J]
% Iu_km_kg= material.disposalCO2 * mass; %[co2/km/kg] 
% Iu_CO_2=Iu_km_kg * driveDistTotal * mass; %[co2/km] 
% 
% 
% %% end-of-life phase
% Ie_energy= material.eolEnergy .* mass; %[J] 
% Ie_CO_2 = material.eolCO2 * mass; %[J]
% 
% prise = material.price * mass; %[SEK]


%% visualise life cycle results
figure(3), clf
subplot(3,1,1)
bargraph(struct('title','Energy','xlabel',{{'Prod','Use','EoL Pot'}},'ylabel','[Joules]','legend',NAMN,'values',[model(:).Ip_energy;model(:).Iu_energy;model(:).Ie_energy]),1);
subplot(3,1,2)
bargraph(struct('title','CO_2','xlabel',{{'Prod','Use','EoL Pot'}},'ylabel','[kg]','legend',NAMN,'values',[model(:).Ip_CO_2;model(:).Iu_CO_2;model(:).Ie_CO_2]),1);
subplot(3,1,3)
bargraph(struct('title','Cost','xlabel',{{''}},'ylabel','[SEK]','legend',NAMN,'values',[model(:).prise]));

%% Compare Two Pie Charts
labels = {'Production CO_{2}','Disposal CO_{2}','eolCO_{2}'};
figure(4)
t = tiledlayout (3,3);

for m = 1:i
    
    pie3(nexttile,[model(m).Ip_CO_2,model(m).Ie_CO_2,model(m).Iu_CO_2])
    title(table2cell(model(m).materialNamn))
    legend(labels)
end





%% FUNCTIONS



function bargraph(result,showtotal)
%this function produces a bar graph of the values in 'result' - see 
%examples for usage. it should not need to be edited.
  if nargin==2
    result.values=[result.values;sum([result.values],1)];
    result.xlabel=[result.xlabel,{'Total'}];
    bar(result.values,'grouped')
    a=axis;
    hold on
    plot(repmat(size(result.values,1)-0.5,1,2),[a(3) a(4)]*10,'-k')
    axis(a)
    hold off
  elseif size(result.values,1)==1 & size(result.values,2)>1
    bar([result.values;result.values],'grouped')
    a=axis;
    a(2)=a(2)-1;
    axis(a)
  else
    bar(result.values,'grouped')
  end
  legend(result.legend','location','bestoutside')
  title(result.title)
  set(gca,'xticklabel',result.xlabel)
  ylabel(result.ylabel)
end

function MODEL = initModelStructArray(material,m,model)

model(m).F=1e5;
model(m).delta=1e-3;
model(m).L=1;
model(m).b=1;
model(m).h = sqrt((3/2) .* (model(m).F .* model(m).L)/(material.tensileStrength .*model(m).b));
model(m).E= material.youngsModulus;
model(m).rho=material.density;
model(m).mass = model(m).rho .* (model(m).b .* model(m).L .* model(m).h); 

model(m).Ip_energy=material.productionEnergy*model(m).mass; %[J] In Production
model(m).Ip_CO_2 = material.productionCO2 * model(m).mass; %[CO2] In Production
driveDistTotal=1; %[km] %total driving distance

model(m).Iu_energy = material.disposalEnergy * model(m).mass; %[J]
model(m).Iu_km_kg= material.disposalCO2 * model(m).mass; %[co2/km/kg] 
model(m).Iu_CO_2=model(m).Iu_km_kg * driveDistTotal * model(m).mass; %[co2/km]
model(m).Ie_energy= material.eolEnergy .* model(m).mass; %[J] 
model(m).Ie_CO_2 = material.eolCO2 * model(m).mass; %[J]
model(m).prise = material.price * model(m).mass; %[SEK]
MODEL = model(m); 
end 



function model=initModel(material)
%this function creates intialises/restores the model values to default
%values. The model variables are case specfic and they could be edited.
  model.F=1e5;
  model.delta=1e-3;
  model.L=1;
  model.b=1;
  model.h=1;
  model.E= material.youngsModulus;
  model.rho=material.density;
end

function mass=computeMass(model)
%this function should determine the mass of the model based on its
%geometric dimensions and material density.
  names=fieldnames(model);
  for i=1:numel(names)
    eval([names{i},'=model.',names{i},';'])
  end
  mass = model.rho .* (model.b .* model.L .* model.h); 
end

function model=computeOptimalVariable(model,material)
%this function should determine the optimal values for the geometric
%variable.
  names=fieldnames(model);
  for i=1:numel(names)
    eval([names{i},'=model.',names{i},';'])
  end
  model.h = sqrt((3/2) .* (model.F .* model.L)/(material.tensileStrength .*model.b));
  
end

function  [energy_consumption, CO2_footprint] = Product_Environmental_Footprint(material,vikt)
    energy_consumption = (material.disposalEnergy + material.eolEnergy + material.productionEnergy).*vikt;
    CO2_footprint = (material.disposalCO2 + material.eolCO2 + material.productionCO2).*vikt; 
end

function [fitresult, gof] = createFit(cost, Ip, prod_energy)
[xData, yData, zData] = prepareSurfaceData( cost, Ip, prod_energy );
% Set up fittype and options.
ft = 'biharmonicinterp';

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );

% Plot fit with data.
figure(1);
h = plot( fitresult, [xData, yData], zData );
legend( h, 'untitled fit 1', 'prod_energy vs. cost, Ip', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Prise', 'Interpreter', 'none' );
ylabel( 'Total CO_2', 'Interpreter', 'none' );
zlabel( 'Total Energy', 'Interpreter', 'none' );
grid off
view( -26.4, 55.5 );

end


