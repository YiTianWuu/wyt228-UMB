function [] = WriteFilterVar2File(VarMolList_filter, VarMol_filter, VarMolList, VarMol, input_ng, filename)
% Write filtered variant information to file
% 2 files: "libx_varlist.xls" including positions in EnR regions
% "libx_varlist_genomic.xls" including genomic positions


% Conversion yield of different plexes
CY_allplex = [0.659482759 0.407758621 0.379310345 0.318965517 0.509913793 0.240086207 0.353017241 0.337068966 0.320258621 0.824568966 0.460775862 0.377155172 0.292672414 0.102155172 0.234051724 0.579741379 0.359913793 0.331034483 0.164224138 0.556896552 0.238362069 0.363362069 0.032327586 0.378448276 0.174137931 0.634913793 0.444827586 0.275862069];


%% Make printable cell arrays
% columns: (mut name, count, vaf) * # of mut in this plex
Varfilter_print = cell(size(VarMolList_filter,1),3);
Var_print = cell(size(VarMolList_filter,1),3);

for tarnum = 1:size(VarMolList_filter,1) % each target
    curCY = CY_allplex(tarnum);
    for mutnum = 1:size(VarMolList_filter,2)
        if ~isempty(VarMolList_filter{tarnum,mutnum})
            Varfilter_print{tarnum,3*mutnum-2} = VarMolList_filter{tarnum,mutnum}; % mut name
            Varfilter_print{tarnum,3*mutnum-1} = VarMol_filter(tarnum,mutnum); % count
            Varfilter_print{tarnum,3*mutnum} = VarMol_filter(tarnum,mutnum) / (input_ng*300*2*curCY); % vaf
        end
    end
end

for tarnum = 1:size(VarMolList,1) % each target
    curCY = CY_allplex(tarnum);
    for mutnum = 1:size(VarMolList,2)
        if ~isempty(VarMolList{tarnum,mutnum})
            Var_print{tarnum,3*mutnum-2} = VarMolList{tarnum,mutnum}; % mut name
            Var_print{tarnum,3*mutnum-1} = VarMol(tarnum,mutnum); % count
            Var_print{tarnum,3*mutnum} = VarMol(tarnum,mutnum) / (input_ng*300*2*curCY); % vaf
        end
    end
end


%% Write genomic positions of mutations

%columns:Chr, EnR first base pos, direction(fwd=1, rvs=2)
EnRpos_dir = csvread('EnRpos_dir_AMLQBDA.csv');

Varfilter_print_genome = Varfilter_print;
Var_print_genome = Var_print;

for tarnum = 1:size(VarMolList,1) % each target
    for mutnum = 1:floor(size(Varfilter_print,2)/3)
        if ~isempty(Varfilter_print{tarnum,mutnum*3-2})
            Varfilter_print_genome{tarnum,mutnum*3-2} = EnR2GenomePos(Varfilter_print{tarnum,mutnum*3-2}, EnRpos_dir(tarnum,1), EnRpos_dir(tarnum,2), EnRpos_dir(tarnum,3));
        end
    end
    
    
    for mutnum = 1:floor(size(Var_print,2)/3)
        if ~isempty(Var_print{tarnum,mutnum*3-2})
            Var_print_genome{tarnum,mutnum*3-2} = EnR2GenomePos(Var_print{tarnum,mutnum*3-2}, EnRpos_dir(tarnum,1), EnRpos_dir(tarnum,2), EnRpos_dir(tarnum,3));
        end
    end
end


filestr = [filename(6:end-22) '_varlist_genomic.xls'];
writecell(Varfilter_print_genome,filestr,'Sheet','Filtered');
writecell(Var_print_genome,filestr,'Sheet','BeforeFilter');


filestr = [filename(6:end-22) '_varlist.xls'];
writecell(Varfilter_print,filestr,'Sheet','Filtered');
writecell(Var_print,filestr,'Sheet','BeforeFilter');
