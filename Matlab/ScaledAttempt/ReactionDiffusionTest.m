classdef ReactionDiffusionTest < BaseModel

    properties
        
        N % Diffision model for component N
        R % Diffision model for component R
        C % Diffision model for result
        
        k1 % reaction constant of N + R -> C
        k_1 % reaction constant of C -> N + R
        
    end
    
    methods
        
        function model = ReactionDiffusionTest(paramobj)
            
            model = model@BaseModel();
            
            % All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            %% Setup the model using the input parameters
            fdnames = {'G'            , ...
                       'k1', 'k_1'};
            
            model = dispatchParams(model, paramobj, fdnames);

            model.N = DiffusionComponent(paramobj.N);
            model.R = DiffusionComponent(paramobj.R);
            model.C = DiffusionComponent(paramobj.C);
            
        end

        
        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            %% Temperature dispatch functions
            fn = @ReactionDiffusion.updateSourceTerm;
            
            inputnames = {{'N', 'c'}, ...
                          {'R', 'c'}, ...
                          {'C', 'c'}};
            model = model.registerPropFunction({{'N', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'R', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'C', 'source'} , fn, inputnames});
            
        end

        function forces = getValidDrivingForces(model);
            forces = getValidDrivingForces@PhysicalModel(model);
            forces.none = [];
        end
        
        function state = updateSourceTerm(model, state)

            k1 = model.k1;
            k_1 = model.k_1;

            vols = model.G.cells.volumes;
            
            cN = state.N.c;
            cR = state.R.c;
            cC = state.C.c;

            R1 = k1.*vols.*cN.*cR;
            R_1 = k_1.*vols.*cC;
            %R1 = k1.*cN.*cR;
            %R_1 = k1.*cN.*cR;

            state.N.source = 0;
            state.R.source = 0;
            state.C.source = 0;
            
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            
            state = model.initStateAD(state);

            state.N = model.N.updateFlux(state.N);
            state.R = model.R.updateFlux(state.R);
            state.C = model.C.updateFlux(state.C);
            
            state = model.updateSourceTerm(state);

            state.N = model.N.updateMassAccum(state.N, state0.N, dt);
            state.R = model.R.updateMassAccum(state.R, state0.R, dt);
            state.C = model.C.updateMassAccum(state.C, state0.C, dt);
            
            state.N = model.N.updateMassConservation(state.N);
            state.R = model.R.updateMassConservation(state.R);
            state.C = model.C.updateMassConservation(state.C);
            
            eqs = {}; types = {}; names = {};
            
            eqs{end + 1}   = state.N.massCons;
            names{end + 1} = 'massCons N';
            types{end + 1} = 'cell';
            
            eqs{end + 1}   = state.R.massCons;
            names{end + 1} = 'massCons R';
            types{end + 1} = 'cell';

            eqs{end + 1}   = state.C.massCons;
            names{end + 1} = 'massCons C';
            types{end + 1} = 'cell';
                        
            primaryVars = model.getPrimaryVariables();

            %% Setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end
        
        
        function state = initStateAD(model, state)
        % initialize a new cleaned-up state with AD variables

            % initStateAD in BaseModel erase all fields
            newstate = initStateAD@BaseModel(model, state);
            newstate.time = state.time;
            state = newstate;
            
        end 

        
        function primaryvarnames = getPrimaryVariables(model)

            primaryvarnames = {{'N', 'c'}, ...
                               {'R', 'c'}, ...
                               {'C', 'c'}};
            
        end
        

        function model = validateModel(model, varargin)
        % nothing special to do
        end


    end

end



%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
