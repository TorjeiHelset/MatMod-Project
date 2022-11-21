classdef ReactionDiffusionGliaCells < BaseModel

    properties
        
        N % Diffision model for component N
        R % Diffision model for component R
        C % Diffision model for result
        T % Diffusion model for component T
        Tb % Diffusion model for result Tb
        N_inac % Diffusion model for N_inac

        k1 % reaction constant of N + R -> C
        k_1 % reaction constant of C -> N + R

        k_t1 % reaction constant of T + C -> Tb
        k_t_1 % reaction constant of Tb -> T + C
        k_t2  % reaction constant of Tb -> T + N_inac
        
    end
    
    methods
        
        function model = ReactionDiffusionGliaCells(paramobj)
            
            model = model@BaseModel();
            
            % All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            %% Setup the model using the input parameters
            fdnames = {'G'            , ...
                       'k1', 'k_1', 'k_t1', 'k_t_1', 'k_t2'};
            
            model = dispatchParams(model, paramobj, fdnames);

            model.N = DiffusionComponent(paramobj.N);
            model.R = DiffusionComponent(paramobj.R);
            model.C = DiffusionComponent(paramobj.C);
            model.T = DiffusionComponent(paramobj.T);
            model.Tb = DiffusionComponent(paramobj.Tb);
            model.N_inac = DiffusionComponent(paramobj.N_inac);
            
        end

        
        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            %% Temperature dispatch functions
            fn = @ReactionDiffusion.updateSourceTerm;
            
            inputnames = {{'N', 'c'}, ...
                          {'R', 'c'}, ...
                          {'C', 'c'}, ...
                          {'T', 'c'}, ...
                          {'Tb', 'c'}, ...
                          {'N_inac', 'c'}};

            model = model.registerPropFunction({{'N', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'R', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'C', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'T', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'Tb', 'source'} , fn, inputnames});
            model = model.registerPropFunction({{'N_inac', 'source'} , fn, inputnames});
            
        end

        function forces = getValidDrivingForces(model);
            forces = getValidDrivingForces@PhysicalModel(model);
            forces.none = [];
        end
        
        function state = updateSourceTerm(model, state)

            k1 = model.k1;
            k_1 = model.k_1;
            k_t1 = model.k_t1;
            k_t_1 = model.k_t_1;
            k_t2 = model.k_t2;

            vols = model.G.cells.volumes;
            
            cN = state.N.c;
            cR = state.R.c;
            cC = state.C.c;
            cT = state.T.c;
            cTb = state.Tb.c;

            R1 = k1.*vols.*cN.*cR;
            R_1 = k_1.*vols.*cC;
            R_t1 = k_t1.*vols.*cT.*cN;
            R_t_1 = k_t_1.*vols.*cTb;
            R_t2 = k_t2.*vols.*cTb;

            %k_t1 % reaction constant of T + N -> Tb
            %k_t_1 % reaction constant of Tb -> T + N
            %k_t2  % reaction constant of Tb -> T + N_inac

            state.N.source = -R1 + R_1 -R_t1 + R_t_1;
            state.R.source = -R1 + R_1;
            state.C.source = R1 - R_1;
            state.T.source = -R_t1 + R_t_1 + R_t2;
            state.Tb.source = R_t1 - R_t_1 - R_t2;
            state.N_inac.source = R_t2;
            
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            
            state = model.initStateAD(state);

            state.N = model.N.updateFlux(state.N);
            state.R = model.R.updateFlux(state.R);
            state.C = model.C.updateFlux(state.C);
            state.T = model.T.updateFlux(state.T);
            state.Tb = model.Tb.updateFlux(state.Tb);
            state.N_inac = model.N_inac.updateFlux(state.N_inac);

            state = model.updateSourceTerm(state);

            state.N = model.N.updateMassAccum(state.N, state0.N, dt);
            state.R = model.R.updateMassAccum(state.R, state0.R, dt);
            state.C = model.C.updateMassAccum(state.C, state0.C, dt);
            state.T = model.T.updateMassAccum(state.T, state0.T, dt);
            state.Tb = model.Tb.updateMassAccum(state.Tb, state0.Tb, dt);
            state.N_inac = model.N_inac.updateMassAccum(state.N_inac, state0.N_inac, dt);

            state.N = model.N.updateMassConservation(state.N);
            state.R = model.R.updateMassConservation(state.R);
            state.C = model.C.updateMassConservation(state.C);
            state.T = model.T.updateMassConservation(state.T);
            state.Tb = model.Tb.updateMassConservation(state.Tb);
            state.N_inac = model.N_inac.updateMassConservation(state.N_inac);
            
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

            eqs{end + 1}   = state.T.massCons;
            names{end + 1} = 'massCons T';
            types{end + 1} = 'cell';

            eqs{end + 1}   = state.Tb.massCons;
            names{end + 1} = 'massCons Tb';
            types{end + 1} = 'cell';

            eqs{end + 1}   = state.N_inac.massCons;
            names{end + 1} = 'massCons N_inac';
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
                               {'C', 'c'}, ...
                               {'T', 'c'}, ...
                               {'Tb', 'c'}, ...
                               {'N_inac', 'c'}};
            
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
