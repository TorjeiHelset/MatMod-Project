close all

mrstModule add ad-core mrst-gui 

jsonfile = fileread('diffusionMainPartScaled.json');
jsonstruct = jsondecode(jsonfile);


paramobj = ReactionDiffusionInputParamsMainPart(jsonstruct);

paramobj.k_1 = paramobj.k_1*(mol/meter)*(1/second);
paramobj.k_1
