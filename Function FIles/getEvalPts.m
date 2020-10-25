function [Xi_pos,Eta_pos] = getEvalPts(xiRes,etaRes,Xi,Eta)

    stepXi = 1 / xiRes;
    stepEta = 1 / etaRes;
    
    % Evaluation points for plots
    Xi_pos    = 0 : stepXi : 1;      % Evaluate basis functions at each xi in Xi_store
    Eta_pos   = 0 : stepEta : 1;  % Evaluate basis functions at each eta in Eta_store
    
    if(nargin == 4)
        Xi_pos    = unique(sort([Xi_pos Xi]));              % <- Makes sure that the knots are included (In order to plot the knots, otherwise we could ignore this)
        Eta_pos   = unique(sort([Eta_pos Eta]));            % <- Makes sure that the knots are included (In order to plot the knots)
    end
    

end