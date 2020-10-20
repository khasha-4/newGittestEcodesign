% This file is a part of ecoDesign, a example code to design and 
% assess the environmental performance of simply supported panel.
% 
% Copyright (C) 2020 Ciarán O'Reilly <ciaran@kth.se>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.% 
% 
% This provides an example design a mono-material simly-supported panel
% where for a selected material 1) the thickness is picked to meet a
% deflection constraint, and 2) the life cycle impact of the correspoinding 
% mass of material is evaluated.

%% prelims
clear all
clc
close all
addpath('.') %path to material database [you could pick a different material database]

%% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 14);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A3:N11";

% Specify column names and types
opts.VariableNames = ["VarName1", "price", "density", "youngsModulus", "poissonsRatio", "yieldStrength", "tensileStrength", "productionEnergy", "productionCO2", "disposalEnergy", "disposalCO2", "eolEnergy", "eolCO2","info"];
opts.VariableTypes = ["char", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double","char"];

% Specify variable properties
opts = setvaropts(opts, "VarName1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "VarName1", "EmptyFieldRule", "auto");

% Import the data
% materialsData=readtable('materialData.xlsx','ReadRowNames',true);
% materialsData=addvars(materialsData,materialsData.Row,'NewVariableName','info');
% materialsData=table2struct(materialsData);


% Import the data
materialsData=readtable('materialData.xlsx',opts, "UseExcel", false);
materialsData=table2struct(materialsData);

%% Clear temporary variables
clear opts
  
%% iterate material alternatives
i=0;
for m=[2 4]
  i=i+1;
  
  %% initiate model of the panel
  model=initModel; %define the model in function
  
  %% select a material
  model=setModelMaterial(model,materialsData(m));
  names(i)={materialsData(m).VarName1};
  
  %% determine height to fit constraint and the resulting mass
  model=computeOptimalVariable(model);
  h(i)=model.H;
  mass(i)=computeMass(model);
  
  %% compute LCE and LCCO2
  LCE(i,:)=computeLCE(model);
  LCCO2(i,:)=computeLCCO2(model);
  
%% end loop
end

%% visualise design results
figure(1), clf
subplot(2,1,1)
bargraph(struct('title','Thickness','xlabel',{{''}},'ylabel','[m]','legend',{names},'values',[h],'showtotal',0))
subplot(2,1,2)
bargraph(struct('title','Mass','xlabel',{{''}},'ylabel','[kg]','legend',{names},'values',[mass],'showtotal',0))

figure(2), clf
subplot(2,1,1)
bargraph(struct('title','Energy','xlabel',{{'Prod','Use','EoL Disp','EoL Pot'}},'ylabel','[GJ]','legend',{names},'values',LCE','showtotal',1))
subplot(2,1,2)
bargraph(struct('title','CO_2','xlabel',{{'Prod','Use','EoL Disp','EoL Pot'}},'ylabel','[tonne]','legend',{names},'values',LCCO2','showtotal',1))







%% FUNCTIONS

%%
function model=initModel
%this function creates intialises/restores the model values to default
%values. The model variables are case specfic and they could be edited.
  model.dmax=[5e-3];
  model.driveDistTotal=1e5;
  model.P=-1e4;
  model.L=1;
  model.B=1;
  model.H=0.05;
end

%%
function model=setModelMaterial(model,materialsData)
  model.E=materialsData.youngsModulus;
  model.rho=materialsData.density;
  model.EProd=materialsData.productionEnergy;
  model.EDisp=materialsData.disposalEnergy;
  model.EEoL=materialsData.eolEnergy;
  model.CO2Prod=materialsData.productionCO2;
  model.CO2Disp=materialsData.disposalCO2;
  model.CO2EoL=materialsData.eolCO2;
  model.Cost=materialsData.price;
end

%%
function model=computeOptimalVariable(model)
%this function should determine the optimal values for the geometric
%variable.
  model.H=abs(model.P*model.L^3/(4*model.E*model.B*model.dmax)).^(1/3);
end

%%
function mass=computeMass(model)
%this function should determine the mass of the model based on its
%geometric dimensions and material density.
  mass=model.L.*model.B.*model.H.*model.rho;
end

%%
function LCE=computeLCE(model)
  mass=computeMass(model);
  % production phase
  Ep_kg=model.EProd; %[J/kg]
  Ep=sum(Ep_kg.*mass);
  % use phase      %simple model based on fuel efficiency
  Eu_km_kg=11.5/100*33.7e6/1500; %[J/km/kg]
  Eu=sum(Eu_km_kg.*model.driveDistTotal.*mass); %[J]
  % disposal phase
  Ed_kg=model.EDisp; %[J/kg]
  Ed=sum(Ed_kg.*mass); %[J]
  % end-of-life phase
  Ee_kg=model.EEoL; %[J/kg]
  Ee=sum(Ee_kg.*mass); %[J]
  % total
  LCE=[Ep Eu Ed Ee]*1e-9;
end

%%
function LCCO2=computeLCCO2(model)
  mass=computeMass(model);
  % production phase
  CO2p_kg=model.CO2Prod; %[kg/kg]
  CO2p=sum(CO2p_kg.*mass); %[kg]
  % use phase      %simple model based on fuel efficiency
  CO2u_km_kg=11.5/100*2.31/1500; %[kg/km/kg]
  CO2u=sum(CO2u_km_kg.*model.driveDistTotal.*mass); %[kg]
  % disposal phase
  CO2d_kg=model.CO2Disp; %[kg/kg]
  CO2d=sum(CO2d_kg.*mass); %[kg]
  % end-of-life phase
  CO2e_kg=model.CO2EoL; %[kg/kg]
  CO2e=sum(CO2e_kg.*mass); %[kg]
  % total
  LCCO2=[CO2p CO2u CO2d CO2e]*1e-3;
end


%%
function bargraph(result)
%this function produces a bar graph of the values in 'result' - see 
%examples for usage. it should not need to be edited.
  if result.showtotal
    result.values=[result.values;sum([result.values],1)];
    result.xlabel=[result.xlabel,{'Total'}];
    bar(result.values,'grouped')
    a=axis;
    hold on
    plot(repmat(size(result.values,1)-0.5,1,2),[a(3) a(4)]*10,'-k')
    axis(a)
    hold off
  elseif size(result.values,1)==1&size(result.values,2)>1
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