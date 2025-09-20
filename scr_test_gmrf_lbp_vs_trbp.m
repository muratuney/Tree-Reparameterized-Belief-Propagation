% This file and the library are licensed under the T&Cs of 
% Creative Commons BY-NC-SA
%
% If you benefit from this code, please also cite the paper that details the algorithm:
%  
% Murat Uney, Simo Sarkka, Angel Garcia-Fernandez, "Tree Reparameterized Belief Propagation for Gaussian Markov Random Fields"
% submitted to the 2026 IEEE International Conference on Acoustics, Speech, and Signal Processing 
%
% This script tests Tree-Reparameterized Gaussian Belief
% Propagation against Loopy Belief Propagation on an Undirected Network.
%

clear all;

global DEBUG_GMRF 
global DEBUG_GMRF_CARRAY 
global DEBUG_VERBOSE
DEBUG_GMRF = 1;
DEBUG_GMRF_CARRAY = {};
DEBUG_VERBOSE = 1;

% Add all the folders and subfolders under the current directory
addpath(genpath(cd));

% Create a Gaussian MRF
N = 10^2; % number of nodes
V = [1:N]; % The set of vertices
%[Lambda, E] = gaussgrid( N , 0.2, 1);
[Lambda, E] = attractrepulsegmrf( N , 0.15, 1); % Find the information/precision matrix of an attractive/repulse model

mu_x = randn(N,1); % Pick a random mean vector
C_x = Lambda^-1;   % The covariance of the model is the inverse of the information/precision matrix

p_x = cpdf( gk( C_x, mu_x ) ); % Create a Gaussian pdf object

% Find the marginals of the Gaussian pdf object as baselines
for i=1: length( V )
     marginals(i) = p_x.marginalise(i); % store the myopic posteriors
end

% Find the node and edge potentials
dims = ones(N,1); % dimensionality of the variables associated with V
[nodepots,edgepots] = gpot.findgmrfpotentials( mu_x, C_x, E, dims ); % Find edge potential functions given the whole model

% Configure a gmrf object for TRP-BP and loopy BP
mygmrfcfg = gmrfcfg;
mygmrfcfg.itermax = 20; % number of iterations
mygmrfcfg.V = V;        % The vertex set
mygmrfcfg.E = E;        % The edge set
mygmrfcfg.edgepots = edgepots; % Edgepotentials

nodecfgs = gmrfnodecfg([]); % Configure the node objects
for i=1:numel(V)
    mynodecfg = gmrfnodecfg;
    mynodecfg.nodepot = nodepots(i); % Node potentials
    nodecfgs(i) = mynodecfg;
end

mygmrfcfg.nodes = nodecfgs;

% Create the Gaussiam MRF object
mygraph = gmrf( mygmrfcfg );

% Perform Loopy Belief Propagation
mygraph = mygraph.lbp;

% Get the BP results, i.e., estimates for poterior marginals
margs_lbp = [ mygraph.nodes.state ]; % take the state fields of the nodes array in an array of its own 

% Find the LBP error given by the Kullback-Leibler divergences of node
% beliefs from the actual marginals
klderrors_lbp = [];
numberofmessages_lbp = [];
for icnt=1:length(DEBUG_GMRF_CARRAY)
    numberofmessages_lbp = [numberofmessages_lbp, DEBUG_GMRF_CARRAY{icnt}.numberofmessages ];
    for i=1:length(V)
        marg_lbpobj = DEBUG_GMRF_CARRAY{icnt}.nodes(i).state;
        if ~isempty( marg_lbpobj)
            klderrors_lbp(icnt,i) = marginals(i).kld( marg_lbpobj );
        else
            klderrors_lbp(icnt,i) = inf;
        end
    end
end
kldsum_lbp = sum(klderrors_lbp,2)/N;

%% TRMP
% Get two spanning trees
E_row = treeovergrid( N, 'row' );
E_col = treeovergrid( N, 'col' );
messagetrees = {E_row, E_col};

% Configure gmrf object for tree reparameterised bp
mygmrfcfg = gmrfcfg;
mygmrfcfg.itermax = 20;
mygmrfcfg.V = V;
mygmrfcfg.E = E;
mygmrfcfg.edgepots = edgepots;

nodecfgs = gmrfnodecfg([]);
for i=1:numel(V)
    mynodecfg = gmrfnodecfg;
    mynodecfg.nodepot = nodepots(i);
    nodecfgs(i) = mynodecfg;
end

mygmrfcfg.nodes = nodecfgs;


mygraph = gmrf( mygmrfcfg );
% Perform Tree reparameterised belief propagation
DEBUG_GMRF_CARRAY = {};
mygraph = mygraph.trbp( messagetrees );

% Get the TRBP results, i.e., estimates for poterior marginals
margs_trbp = [ mygraph.nodes.state ]; % take the state fields of the nodes array in an array of its own

% Find the TRP-BP error given by the Kullback-Leibler divergences of node
% beliefs from the actual marginals
klderrors_trbp = [];
numberofmessages_trbp = [];
for icnt=1:length(DEBUG_GMRF_CARRAY)
    numberofmessages_trbp = [numberofmessages_trbp, DEBUG_GMRF_CARRAY{icnt}.numberofmessages ];
    for i=1:length(V)
        marg_bpobj = DEBUG_GMRF_CARRAY{icnt}.nodes(i).state;
        if ~isempty( marg_bpobj)
            klderrors_trbp(icnt,i) = marginals(i).kld( marg_bpobj );
        else
            klderrors_trbp(icnt,i) = inf;
        end
    end
end
kldsum_trbp = sum(klderrors_trbp,2)/N;
disp(sprintf('For LBP, the final KLD sum over nodes is %g achived by %d messages', kldsum_lbp(end), numberofmessages_lbp(end) ));
disp(sprintf('For TRBP, the final KLD sum over nodes is %g achived by %d messages', kldsum_trbp(end), numberofmessages_trbp(end) ));


KLDfigure = newfigure;
hold on
grid on
l1 = plot(numberofmessages_lbp, kldsum_lbp, 'Linestyle','--','Color','k','Marker','x' );
xlabel('Number of messages')
ylabel('Average KLD')
l2 = plot(numberofmessages_trbp, kldsum_trbp, 'Linestyle',':','Color','b','Marker','+' );
legend('LBP','TRBP')

