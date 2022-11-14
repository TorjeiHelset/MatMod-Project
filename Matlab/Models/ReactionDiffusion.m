classdef ReactionDiffusion < BaseModel

    properties
        
        A % Diffision model for component A
        B % Diffision model for component B
        C % Diffision model for result
        
        k % reaction constant
        
    end
    
    methods
        
        function model = ReactionDiffusion(paramobj)
            
            model = model@BaseModel();
            
            % All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            %% Setup the model using the input parameters
            fdnames = {'G'            , ...
                       'k'};
            
            model = dispatchParams(model, paramobj, fdnames);

            model.A = DiffusionComponent(paramobj.A);
            model.B = DiffusionComponent(paramobj.B);
            model.C = DiffusionComponent(paramobj.C);
            
        end

        
        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            %% Temperature dispatch functions
            fn = @ReactionDiffusion.updateSourceTerm;
            
            inputnames = {{'A', 'c'}, ...
                          {'B', 'c'}};
            model = model.registerPropFunction({{'A', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'B', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'C', 'source'} , fn, inputnames});
            
        end

        function forces = getValidDrivingForces(model);
            forces = getValidDrivingForces@PhysicalModel(model);
            forces.none = [];
        end
        
        function state = updateSourceTerm(model, state)

            k = model.k;
            vols = model.G.cells.volumes;
            
            cA = state.A.c;
            cB = state.B.c;

            R = k.*vols.*cA.*cB;

            state.A.source = -R;
            state.B.source = -R;
            state.C.source = R;
            
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            
            state = model.initStateAD(state);

            state.A = model.A.updateFlux(state.A);
            state.B = model.B.updateFlux(state.B);
            state.C = model.C.updateFlux(state.C);
            
            state = model.updateSourceTerm(state);

            state.A = model.A.updateMassAccum(state.A, state0.A, dt);
            state.B = model.B.updateMassAccum(state.B, state0.B, dt);
            state.C = model.C.updateMassAccum(state.C, state0.C, dt);
            
            state.A = model.A.updateMassConservation(state.A);
            state.B = model.B.updateMassConservation(state.B);
            state.C = model.C.updateMassConservation(state.C);
            
            eqs = {}; types = {}; names = {};
            
            eqs{end + 1}   = state.A.massCons;
            names{end + 1} = 'massCons A';
            types{end + 1} = 'cell';
            
            eqs{end + 1}   = state.B.massCons;
            names{end + 1} = 'massCons B';
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

            primaryvarnames = {{'A', 'c'}, ...
                               {'B', 'c'}, ...
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
