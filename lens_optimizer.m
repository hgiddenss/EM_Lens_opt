classdef lens_optimizer < handle
    %LENS_OPTIMIZATION optimize a dielectric lens for generation of a flat top beam radiation pattern
    %   Detailed explanation goes here

    properties
        %Lens properties
        lens_width     % lens width in mm (default 250)
        t              % lens thickness in mm (default 50)
        nX             %number of pixels in the X direction (across the width of the lens) (default 20)
        nY             %number of pixels in the Y direction (across the thickness of the lens) (default 4)
        pX             %width of a single pixel in X (i.e lens_wdtih./nX) default = 5;
        pY             %width of a single pixel in Y (i.e t./nY) default = 5;
        focal_distance %distance of focal point behind the lens (default = 25)
        frequency      %frequency in GHz (default = 10)

        min_permittivity %minimum permittivity value (default 1)
        max_permittivity %maximum permittivity value (default 5)

        %Optimization properties
        options 
        data_direc = 'C:\Git 2\LensOptimisationML\surrogate_optimize_data_20230805142532'
    end
    properties (SetAccess = private, SetObservable = true) %can only be Set from within the class
        permittivity_matrix
        optimized = false;
        optimizationResults = struct('fval',[],'exitflag',[],'output',[],'population',[],'score',[],'trials',[])
    end
    properties (Dependent,SetAccess = private)
        wvlen % wavelength in mm, updates whenever frequency property is changed
    end
    properties
        net
        mdl

        inputdata
        outputdata
    end
    methods
        function obj = lens_optimizer(varargin)
            %LENS_OPTIMIZATION Construct an instance of the lens_optimizer
            
            %   Initiate the lens properties with input args / default parameters
            p = inputParser;
            p.addParameter('lens_width',250);
            p.addParameter('t',50);
            p.addParameter('nX',20);
            p.addParameter('nY',4);
            p.addParameter('frequency',10)
            p.addParameter('focal_distance',25);
            p.addParameter('min_permittivity',1);
            p.addParameter('max_permittivity',5);
            
            p.parse(varargin{:});
            params = p.Results;
            [obj.lens_width,obj.t,obj.nX,obj.nY,obj.frequency,obj.focal_distance,obj.min_permittivity,obj.max_permittivity]...
                = deal(params.lens_width,params.t,params.nX,params.nY,params.frequency,params.focal_distance,...
                params.min_permittivity,params.max_permittivity);
            
            %Calculate the other parameters
            obj.pX = obj.lens_width./obj.nX;
            obj.pY = obj.t./obj.nY;


            %intialize a matrix of permittivity values according to the lens properties. Ensure the lens is symmetrical.
            %Assume we are working with an even number of pixels for now.
            %assert(mod(obj.nX,2) == 0,'Error: Ensure the number of pixels is even')
            obj.permittivity_matrix = randi(5,[obj.nX,obj.nY]);

        end
        function value = get.wvlen(obj)
            value = 299792458./obj.frequency*1e-9*1e3;
        end
        function varargout = simulate(obj,varargin)
            %FBT_lens_optimizer.simulate() on its own runs the simulation with the permittivity_matrix properties 
            % and displays the E field in a new figure window
            %[E,H] = FTB_lens_optimizer.simulate() returns the E and H fields from the simulation
            %[E,H] = FTB_lens_optimizer.simulate(permittivity_matrix) runs the simulation with a custom permittivity
            %matrix
            %[E,H] = FTB_lens_optimizer.simulate('plot',false) returns the E and H fields but does not plot the results
            

            p = inputParser;
            p.addOptional('permittivity_matrix',[])
            p.addParameter('plot',true);
            p.parse(varargin{:});
            
            %assign a matrix of permittivity values, either from the input or from the property value, which will be
            %used in simulation.
            if isempty(p.Results.permittivity_matrix)
                permittivity_values = obj.permittivity_matrix;
            else
                permittivity_values = reshape(p.Results.permittivity_matrix,obj.nX,obj.nY);
            end
            assert(numel(permittivity_values) == obj.nX*obj.nY,'Error: wrong number of permittivity values');
            
            
            %initiate colors for display
            nColor = obj.max_permittivity - obj.min_permittivity + 1;
            colorMatrix = cbrewer('qual','Set1',nColor);

            materials = cell(nColor,1);

            %create a new "material" for each discrete permittivity value possible. The discrete permittivity values are
            %just equal to the value of the counter.
            for i = 1:nColor
                Er = [i./2 1 1; 1 i./2 1; 1 1 i./2]+0.5; %Override permittivty values for a maximum of 2.5
                Mu = [1 1 1; 1 1 1; 1 1 1];
                materials{i} = Material(['M',num2str(i)],colorMatrix(i,:),Er,Mu,true);
            end

            %Create the lens for the simulation with nX x nY pixels. The permittivity value of each pixel is determined
            %by the permittivity_matrix value. The lens_parameters matrix is a cell array that contains the dimensions
            %of each pixel as a "box" followed by the corresponding "material" properties
            lens_parameters = cell(0,1);
            for iY = 1:obj.nY
                for iX = 1:obj.nX
                    %i = iX+(iY-1)*obj.nX; %counter

                    X = [-obj.lens_width/2 -obj.lens_width/2+obj.pX]+obj.pX*(iX-1);
                    Y = [obj.pY*(iY-1) obj.pY*iY];
                    C = Box([X; Y; 0 1]);
                    %select the material with permittivity value matching that in the permittivity_matrix - the value in permittivity_matrix is the index of the materials array
                    lens_parameters{end+1} = materials{permittivity_values(iX,iY)}; %#ok<AGROW> 
                    lens_parameters{end+1} = C; %#ok<AGROW>

                    %X = fliplr(X.*-1); %Lens is symetrical along x = 0 axis so add the opposite cell
                    %Y = [obj.pY*(iY-1) obj.pY*iY];
                    %C = Box([X; Y; 0 1]);
                    %lens_parameters{end+1} = materials{permittivity_values(iX,iY)}; %#ok<AGROW>
                    %lens_parameters{end+1} = C; %#ok<AGROW>
                end
            end

            %Simulation domain properties:
            inspect_only = false;

            if isunix %If running on server, iterative method provides increase in speed by around 5x compared to direct calculation
                solveropts.method = 'iterative';
            else
                solveropts.method = 'direct';
            end
            if p.Results.plot
                %This launches the simulation in maxwellFDFD
                [E, H, obj_array, src_array] = maxwell_run(...
                    'OSC', 1e-9, obj.wvlen, ... %units and frequency
                    'DOM', {'vacuum', 'none', 1.0}, [-obj.wvlen*10 obj.wvlen*10; -obj.wvlen*10, obj.wvlen*10; 0, 1],...
                    [1 1 1], [BC.p BC.p BC.p], [10 10 0], 4, exp(-16),true,... %simulation domain, boundary conditions
                    'OBJ', lens_parameters{:},... %materials in the simulation domain that we previously set
                    'SRCJ', LineSrc(Axis.z, [0 -obj.focal_distance], Axis.z,10),...   % LineSrc(axis, intercept)',...  %  source/excitation
                    solveropts,inspect_only);
            else 
                %Generally you should not run functions in evalc - but this has been tested and is working properly.
                %Using evalc supresses the output text in the command window when running the optimisation, which is
                %annoying and hides progress of the optimiser. This is set for when the plot option is false (we dont
                %want to plot every optimiser field map after all) 
                evalc(['[E, H, obj_array, src_array] = maxwell_run(',...
                    '''OSC'', 1e-9, obj.wvlen,', ... %units and frequency
                    '''DOM'', {''vacuum'', ''none'', 1.0}, [-obj.wvlen*10 obj.wvlen*10; -obj.wvlen*10, obj.wvlen*10; 0, 1],',...
                    '[1 1 1], [BC.p BC.p BC.p], [10 10 0], 4, exp(-16),true,',... %simulation domain, boundary conditions
                    '''OBJ'', lens_parameters{:},',... %materials in the simulation domain that we previously set
                    '''SRCJ'', LineSrc(Axis.z, [0 -obj.focal_distance], Axis.z,10),',...   % LineSrc(axis, intercept)',...  %  source/excitation
                    'solveropts,inspect_only)']);

            end

            %set output arguments:
            varargout{1} = E;
            varargout{2} = H;
            varargout{3} = obj_array;
            varargout{4} = src_array;

            %plot with own function:
            if p.Results.plot
                %plot E field results using inbuilt function... This is faster
                hAx = obj.plotSetup;
                X = E{Axis.z}.grid3d.lall{1,1}'; %X-grid
                Y = E{Axis.z}.grid3d.lall{2,1};  %Y-grid
                E1 = E{3}.array(:,:,2)';
                %E1 = E1./max(E1(:)); %Normalize for plotting purposes
                s = surf(hAx,X,Y,-1*ones(size(E1)),real(E1)); 
                hAx.DataAspectRatio = [1 1 1]; hAx.NextPlot = 'add'; hAx.XLimSpec = 'tight'; hAx.YLimSpec = 'tight';
                s.EdgeColor = 'none';
                s.FaceColor = 'interp';
                colormap(hAx,flipud(cbrewer('div','RdYlBu',128)));
                cbar = colorbar(hAx); %#ok<NASGU> 
                hAx.CLim = [-1 1]*max(real(E1(:)));
                drawnow
                view([0 90])
                
            end

        end
        function hAx = plotSetup(obj)
            hAx = axes('parent',figure);
            hold on;
            nColor = obj.max_permittivity - obj.min_permittivity + 1;
            colorMatrix = cbrewer('qual','Set1',nColor);

            for iY = 1:obj.nY
                for iX = 1:obj.nX
                    %i = iX+(iY-1)*obj.nX; %counter

                    XX = [-obj.lens_width/2 -obj.lens_width/2+obj.pX]+obj.pX*(iX-1);
                    YY = [obj.pY*(iY-1) obj.pY*iY];

                    xBox = [XX(1) XX(2) XX(2) XX(1) XX(1)];
                    yBox = [YY(1) YY(1) YY(2) YY(2) YY(1)];
                    plot(hAx,xBox,yBox,'color',colorMatrix(obj.permittivity_matrix(iX,iY),:),'linewidth',2)
                    %plot(hAx,xBox*-1,yBox,'color',colorMatrix(obj.permittivity_matrix(iX,iY),:),'linewidth',2)

                end
            end
            %plot location of line source
            plot(hAx,0,-obj.focal_distance,'o','markerFaceColor',[1 0 0],'MarkerEdgeColor','none');
            %Set axis limits - these should be equal to simulation limits
            axis equal;
            xlim([-1 1]*obj.wvlen*10);
            ylim([-obj.wvlen*10 obj.wvlen*10]);
            
        end
        function S = fitnessFunction(obj,x,target_pattern)
            
            [E,~,obj_array,src_array] = obj.simulate(x,'plot',false);

             %I wrote this function to approximate the far field radiation pattern from the near field data from MaxwellFDFD
            [E_far,~] = near2FarApprox(E,obj_array,src_array);
            
            %calculate fitness score based on absolute value of E_far
            %S = mean(abs(abs(E_far)-target_pattern)).^2;
            S = mean(abs(E_far(target_pattern == 1)-1)).^2;
            idx = (target_pattern == 0.1) & E_far < 0.1;
            
            S = mean(abs(E_far(~idx) - target_pattern(~idx))).^2;

            %Save raw data to disk
            fid = fopen(fullfile(obj.data_direc,[datestr(now,'yyyymmddHHMMSS'),'.txt']),'wt');
            fprintf(fid,'%d \t',x(:));
            fprintf(fid,'\n');
            fprintf(fid,'%.2f \t',abs(E_far(:))); %we are only concerned with magnitude of radiation pattern
            fclose(fid);
            %store data in object...
                
        end
        function generate_optimized_training_data(obj)
            
            nParams = obj.nX.*obj.nY;
            IntCon = 1:numel(obj.permittivity_matrix); %three discrete permittivity values
            
            lb = ones(1,nParams)*obj.min_permittivity;
            ub = ones(1,nParams)*obj.max_permittivity;

            % set some default optimisation options:
            obj.options = optimoptions('surrogateopt');
            obj.options = optimoptions(obj.options,'MaxFunctionEvaluations', 2000);
            obj.options = optimoptions(obj.options,'MinSurrogatePoints', 200);
            obj.options = optimoptions(obj.options,'BatchUpdateInterval', 1);
            obj.options = optimoptions(obj.options,'Display', 'iter');
            obj.options = optimoptions(obj.options,'PlotFcn', '');
            obj.options.UseParallel = true;
            
            
            %beam at 0 degrees with beamwdidth of 30 degrees
            phi = -90:3:90;
            beamwidth = 30;
            idx = (phi) >= 0-beamwidth/2 & (phi)  <= 0+beamwidth/2 ;
            E_ideal = ones(size(phi));
            E_ideal(~idx) = E_ideal(~idx) * 0.1;

            f = @(x)fitnessFunction(obj,x,E_ideal);
            ts = tic;
            [x,obj.optimizationResults.fval,obj.optimizationResults.exitflag,obj.optimizationResults.output,...
                obj.optimizationResults.trials] = surrogateopt(f,lb,ub,IntCon,[],[],[],[],obj.options);
            toc(ts)
            [E,~,obj_array,src_array] = obj.simulate(x);
            [E_far,phi] = near2FarApprox(E,obj_array,src_array);
            figure; plot(phi,E_far); hold on; plot(phi,E_ideal);
            %return

            
            %beam at 0 degrees with beamwdidth of 60 degrees
            beamwidth = 60;
            idx = (phi) >= 0-beamwidth/2 & (phi)  <= 0+beamwidth/2 ;
            E_ideal = ones(size(phi));
            E_ideal(~idx) = E_ideal(~idx) * 0.1;
            [x,obj.optimizationResults.fval,obj.optimizationResults.exitflag,obj.optimizationResults.output,...
                obj.optimizationResults.trials] = surrogateopt(f,lb,ub,IntCon,[],[],[],[],obj.options);
            [E,~,obj_array,src_array] = obj.simulate(x);
            [E_far,phi] = near2FarApprox(E,obj_array,src_array);
            figure; plot(phi,E_far); hold on; plot(E_ideal);
            toc(ts)

            %beam at 30 degrees with beamwdidth of 20 degrees
            beamwidth = 20;
            idx = (phi) >= 30-beamwidth/2 & (phi)  <= 30+beamwidth/2 ;
            E_ideal = ones(size(phi))*0.1;
            E_ideal(idx) = 1;
            [x,obj.optimizationResults.fval,obj.optimizationResults.exitflag,obj.optimizationResults.output,...
                obj.optimizationResults.trials] = surrogateopt(f,lb,ub,IntCon,[],[],[],[],obj.options);
            [E,~,obj_array,src_array] = obj.simulate(x);
            [E_far,phi] = near2FarApprox(E,obj_array,src_array);
            figure; plot(phi,E_far); hold on; plot(phi,E_ideal);
            toc(ts)

            %split beam at -30 and +20 degrees with beamwdidth of 10 degrees
            beamwidth = 10;
            idx = (phi) >= -30-beamwidth/2 & (phi)  <= -30+beamwidth/2; 
            idx3 = (phi) >= 20-beamwidth/2 & (phi)  <= 20+beamwidth/2;
            idx = idx | idx3;
            E_ideal = ones(size(phi))*0.1;
            E_ideal(idx) = 1;
            [x,obj.optimizationResults.fval,obj.optimizationResults.exitflag,obj.optimizationResults.output,...
                obj.optimizationResults.trials] = surrogateopt(f,lb,ub,IntCon,[],[],[],[],obj.options);
            toc(ts)

            [E,~,obj_array,src_array] = obj.simulate(x);
            [E_far,phi] = near2FarApprox(E,obj_array,src_array);
            figure; plot(phi,E_far); hold on; plot(phi,E_ideal);
            
        end
        function generate_random_training_data(obj,nRandom)
            if nargin == 1
                nRandom = 500;
            end

            obj.inputdata = zeros(numel(obj.permittivity_matrix(:)),nRandom);
            obj.outputdata = zeros(121,nRandom);
            for i = 1:nRandom
                permittivity_values = randi(10,[obj.nX,obj.nY]);
                [E,~,obj_array,src_array] = obj.simulate(permittivity_values,'plot',false);
                [E_far,~] = near2FarApprox(E,obj_array,src_array);
                E_far = abs(E_far); %magnitude only important
                %Write data to disk
                fid = fopen(fullfile(obj.data_direc,[datestr(now,'yyyymmddHHMMSS'),'.txt']),'wt');
                fprintf(fid,'%d \t',permittivity_values(:));
                fprintf(fid,'\n');
                fprintf(fid,'%.2f \t',E_far(:));
                fclose(fid);
                %store data in object...
                %obj.inputdata(:,i) = permittivity_values(:);
                %obj.outputdata(:,i) = E_far(:);
            end
        end
        function read_data_from_disk(obj)
            info = dir(obj.data_direc);
            info(1:2) = [];
            filenames = {info.name};
            obj.inputdata = zeros(numel(obj.permittivity_matrix(:))+1,0);
            obj.outputdata = zeros(61,0);
            
            phi = (-90:3:90)';
            for i = 1:numel(filenames)
                try %#ok<TRYNC> 
                fid = fopen(fullfile(obj.data_direc,filenames{i}),'rt');

                line = str2double(strsplit(fgetl(fid))); %first line, input data, permittivity values
                line(isnan(line)) = []; %remove NaN value that appears at the end
                input = line; 
                
                line = str2double(strsplit(fgetl(fid))); %second line, output data, far field pattern
                line(isnan(line)) = []; %remove NaN value that appears at the end
                output = line';
                fclose(fid);
                
                input = [phi,repmat(input,numel(phi),1)];
                if size(input,1) == size(output,1) % ensure they have the same number of rows
                    obj.inputdata = [obj.inputdata;input];
                    obj.outputdata = [obj.outputdata;output];
                end

                end
            end
        end
        function reduce_dataset(obj,factor)
            %lens_obj.reduce_dataset(factor) Reduces the dataset by randomly sampling data. The new dataset will have
            %length(outputdata)*factor number of samples

            assert(factor >= 0 & factor <= 1,'Error: factor must be between 0 and 1')

            %randomly partition data
            hpartition = cvpartition(length(obj.outputdata),'Holdout',factor);
            
            %resample the data at the new rate
            obj.inputdata = obj.inputdata(test(hpartition),:);
            obj.outputdata = obj.outputdata(test(hpartition),:);

        end
        function create_net(obj)
            obj.net = fitrnet(obj.inputdata,obj.outputdata,"LayerSizes",[30 10]); 

        end
        function train_net(obj)
            obj.net = train(obj.net,obj.inputdata',obj.outputdata');
        end
        function output = predict_net(obj,map)
            phi = (-90:3:90)';
            input = [phi,repmat(map(:)',numel(phi),1)];
            output = zeros(numel(phi),1);
            for i = 1:numel(phi)
                output(i) = predict(obj.net,input(i,:));
            end
        end
        function [L] = create_model(obj)

            %partition the data set for training and test
            hpartition = cvpartition(length(lens_obj.outputdata),'Holdout',0.3);
            idxTrain = training(hpartition);
            idxTest = test(hpartition);

            Xtrain = lens_obj.inputdata(idxTrain,:);
            Ytrain = lens_obj.outputdata(idxTrain,:);

            Xtest = lens_obj.inputdata(idxTest,:);
            Ytest = lens_obj.outputdata(idxTest,:);

            %train the model
            fprintf('Training model\n')
            ts = tic;
            obj.mdl = fsrnca(Xtrain,Ytrain,'Standardize',true);
            toc(ts)
            
            fprintf('Calculating Loss\n')
            tic;
            L = loss(obj.mdl,Xtest,Ytest);
            toc(ts)

            %Plot the computed and actual responses for the test values
            fprintf('Generating Testing Statistics');
            ts = tic;
            Ypred = predict(obj.mdl,Xtest);
            toc(ts)
            hAx = axes('parent',figure);
            scatter(hAx,Ypred,Ytest,'o','markerfacecolor',[0.4 0.4 0.87],'markerfacealpha',0.1,'MarkerEdgeColor','none');
            xlabel('Predicted response')
            ylabel('Actual response')
            axis equal
            xlim([0 1])
            ylim([0 1])

        end
        function output = predict_model(obj,map)
            phi = (-90:3:90)';
            input = [phi,repmat(map(:)',numel(phi),1)];
            output = zeros(numel(phi),1);
            for i = 1:numel(phi)
                output(i) = predict(obj.mdl,input(i,:));
            end
        end
        function [x] = optimize_with_ML(obj,target_pattern,predictor)
            nparams = numel(obj.permittivity_matrix);
            IntCon = 1:nparams;
            opts = optimoptions('ga');
            % options setting
            opts = optimoptions(opts,'MaxGenerations', 250);
            opts = optimoptions(opts,'MaxStallGenerations', 20);
            %options = optimoptions(options,'FunctionTolerance', 0);
            %options = optimoptions(options,'ConstraintTolerance', 0);
            opts = optimoptions(opts,'Display', 'iter');
            opts = optimoptions(opts,'PlotFcn', 'gaplotbestf');
            opts.PopulationSize = 100;
            opts.UseParallel = true;
            opts.UseVectorized = false;
            lb = ones(1,nparams)*1;
            ub = ones(1,nparams)*5;

            f = @(x)fitness_function_ML(obj,x,target_pattern,predictor);

            [x,obj.optimizationResults.fval,obj.optimizationResults.exitflag,obj.optimizationResults.output,...
                obj.optimizationResults.population,obj.optimizationResults.score] =...
                ga(f,nparams,[],[],[],[],lb,ub,[],IntCon,opts);

            obj.permittivity_matrix = reshape(x,obj.nX,obj.nY);
        end
        function S = fitness_function_ML(obj,x,target_pattern,predictor)
            
            switch lower(predictor)
                case {'net','nueralnet'}
                    E_far_nnet = obj.predict_net(x(:));
                    S = mean(abs(abs(E_far_nnet(:))-target_pattern(:))).^2;
                case {'mdl','model'}
                    E_far_mdl = obj.predict_model(x(:)');
                    S = mean(abs(abs(E_far_mdl(:))-target_pattern(:))).^2;
            end
        end
    end

end

