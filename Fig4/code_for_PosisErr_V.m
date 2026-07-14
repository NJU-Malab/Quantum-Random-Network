%%%%%%%%%%%% Fidelity calculation for target state (111+122+331+342)
%%%%%%%%%%%% FCAA: 17x16 experimental photon distribution matrix
%%%%%%%%%%%% Columns: basis ordering; Rows: different measurement settings
%%%%%%%%%%%% Rows: {122, 121, 132, 131, 142, 141, 112, 111, 322, 321, 332, 331, 342, 341, 312, 311}
%%%%%%%%%%%% Columns: 
%%%%%%%%%%%%    row1-{Computational basis}, 
%%%%%%%%%%%%    row2-{I^1 prod X^{1,2} prod X^{1,2}},
%%%%%%%%%%%%    row3-{I^1 prod Y^{1,2} prod Y^{1,2}},
%%%%%%%%%%%%    row4-{I^1 prod X^{3,4} prod X^{1,2}},
%%%%%%%%%%%%    row5-{I^1 prod Y^{3,4} prod Y^{1,2}},
%%%%%%%%%%%%    row6-{X^{1,3} prod X^{1,3} prod I^1},
%%%%%%%%%%%%    row7-{Y^{1,3} prod Y^{1,3} prod I^1},
%%%%%%%%%%%%    row8-{X^{1,3} prod X^{2,4} prod I^2},
%%%%%%%%%%%%    row9-{X^{1,3} prod X^{2,4} prod I^2},
%%%%%%%%%%%%    row10-{X^{1,3} prod X^{2,3} prod X^{2,1}},
%%%%%%%%%%%%    row11-{X^{1,3} prod Y^{2,3} prod Y^{2,1}},
%%%%%%%%%%%%    row12-{Y^{1,3} prod X^{2,3} prod Y^{2,1}},
%%%%%%%%%%%%    row13-{Y^{1,3} prod Y^{2,3} prod X^{2,1}},
%%%%%%%%%%%%    row14-{X^{1,3} prod X^{1,4} prod X^{1,2}},
%%%%%%%%%%%%    row15-{X^{1,3} prod Y^{1,4} prod Y^{1,2}},
%%%%%%%%%%%%    row16-{Y^{1,3} prod X^{1,4} prod Y^{1,2}},
%%%%%%%%%%%%    row17-{Y^{1,3} prod Y^{1,4} prod X^{1,2}},
%%%%%%%%%%%%    A1C1+A2C2+A3C3+A4C4，B1D1→A1，B1D2→A2，B2D1→A3，B2D2→A4
%%%%%%%%%%%%



close all; clear;
load('.\origin_data.mat');   % load FCAA

% Monte Carlo parameters
N_sim = 1000;                            % Number of Monte Carlo simulations (each simulation = one Poisson‑noised copy of FCAA)
Fidelity = zeros(N_sim, 1);              % Store fidelity for each simulation
ave = zeros(N_sim, 17);                  % Store the 17 measured ratios (one per row) for each simulation
FCAA_orig = FCAA;                        % Preserve original experimental data for generating noisy copies
% ==================== CONFIGURATION FOR EACH ROW ====================
% config(row) is a structure with fields:
%   .cols      : column indices used to compute the ratio (and coff)
%   .yita      : weighting vector for the ratio (empty for row 1)
%   .use_coff  : logical, true if this row contributes to the coff scaling factor
%
% Row 1 (special): ratio = sum(cols)/sum(all columns), no yita, no coff scaling.
% Rows 2-17: ratio = sum(yita .* data_sub) / sum(data_sub), and coff = sum(subset)/sum(row).

% Define configuration for each row (1..17)
% Fields: cols (indices for sub-array), yita (weight vector), use_coff (whether to compute coff)

% For row 1: special (no coff, sum of specific cols over whole row)
config(1) = struct('cols', [1,8,12,13], 'cofc',[1,8,12,13],'yita', [],'yitaz', [], 'use_coff', false);
% Rows 2-3: C1C2
for r = 2:3, config(r) = struct('cols', [1,2,7,8], 'cofc',[1,2,7,8],'yita', [1,-1,-1,1], 'yitaz', [1,-1,-1,1],'use_coff', true); end
% Rows 4-5: C3C4
for r = 4:5, config(r) = struct('cols', [11,12,13,14], 'cofc',[11,12,13,14],'yita', [1,-1,-1,1], 'yitaz', [-1,1,1,-1],'use_coff', true); end
% Rows 6-7: C1C3
for r = 6:7, config(r) = struct('cols', [14,6,10,2],  'cofc',[4,8,12,16],'yita', [1,-1,-1,1], 'yitaz', [-1,1,1,-1],'use_coff', true); end
% Rows 8-9: C2C4
for r = 8:9, config(r) = struct('cols', [1,9,5,13], 'cofc',[1,9,5,13],'yita', [1,-1,-1,1], 'yitaz', [1,-1,-1,1],'use_coff', true); end
% Rows 10-13: C2C3
for r = 10:13, config(r) = struct('cols', [1,2,5,6,9,10,13,14], 'cofc',[1,2,3,4,9,10,11,12],'yita', [1,-1,-1,1,-1,1,1,-1], 'yitaz', [1,0,-1,0,0,-1,0,1],'use_coff', true); end
% Rows 14-17: C1C4
for r = 14:17, config(r) = struct('cols', [1,2,5,6,9,10,13,14], 'cofc',[5,6,7,8,13,14,15,16],'yita', [1,-1,-1,1,-1,1,1,-1], 'yitaz', [0,-1,0,1,1,0,-1,0],'use_coff', true); end

% ==================== FIDELITY WEIGHTS ====================
% Base weights y0 for the 17 observables (before multiplying by coff), the weights for 2-17 observables is y0.*coff, where coff is 0.5 for ideal case.
% These correspond to the ideal target state
% The weights are derived from the quantum mechanical expectation values.

y0 = [1/4, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, ...
      1/16, 1/16, 1/16, 1/16, 1/16, 1/16, 1/16, 1/16];

for rot = 1:N_sim
    FCAA = possrand(FCAA_orig);    % add Poisson noise
    coff = zeros(1,17);
    coff(1) = 1;                   % row1 coff fixed to 1
    FCA=FCAA(1,:);
    for row = 1:17
        data_row = FCAA(row, :);
        cfg = config(row);
        % Compute the ratio stored in ave(rot,row)
        if row == 1
            % Row1: ratio = sum(cols) / sum(all columns)
            ub = sum(data_row(cfg.cols));
            vb = sum(data_row);
        else
            % Other rows: ratio = sum(yita .* data_sub) / sum(data_sub)
            sub = data_row(cfg.cols);
            ub = sum(cfg.yita .* sub);
            vb = sum(sub);
        end
        if vb == 0
            ave(rot, row) = 0;
        else
            ave(rot, row) = round(ub / vb, 4);
        end
        
        % Compute coff for this row if needed (used for fidelity weighting)
        if cfg.use_coff
            total_row = sum(data_row);
            if total_row ~= 0
                coff(row) = sum(data_row(cfg.cols)) / total_row;
            end
        end

    end
    yy = ave(rot, :)';
    for ind= 1: 6                 % index order for C1C2,C3C4,C1C3,C2C4,C2C3,C1C4
        search_cofc=[2 4 6 8 10 17];
        zind=config(search_cofc(ind)).cofc;
        Vzind(ind)=sum(FCA(zind).*config(search_cofc(ind)).yitaz)/sum(FCA(zind).*abs(config(search_cofc(ind)).yitaz));
    end





    %%% add variable coffe 
    coffe(1)=coff(1);
    for r=2:17
    coffe(r)=sum(FCA(config(r).cofc))/sum(FCA);
    end
    %%%
    
    % Fidelity for this simulation
    W(rot)=sum(abs(yy([2:9]))+abs(yy([10:17]))/2)+sum(Vzind);
end

% Compute final statistics

WW = sort(W, 'descend');
W_mean = mean(WW);
W_std = std(WW);
% Display (optional)
fprintf('Mean W = %.5f ± %.5f \n', W_mean, W_std);

%%%%%%%%%%%% Local function definition  %%%%%%%%%%%%
function FCAA_Pos = possrand(FCAA)
    % Generate Poisson random copies for each matrix element
    kd = zeros(17,16);
    for i = 1:17
        for j = 1:16
            kd(i,j) = poissrnd(FCAA(i,j));
        end
    end
    FCAA_Pos = kd;
end
