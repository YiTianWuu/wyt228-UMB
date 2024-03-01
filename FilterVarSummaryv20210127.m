function [VarMolList_filter VarMol_filter] = FilterVarSummaryv20210127(VarMolList, VarMol, input_ng, fix_thresh, vaf_thresh)
% Remove error based on different threshold for different mutations
% different threshold for A>G, G>A mutations, homopolymer indels, and others
% Adjust threshold based on input ng

% Standard mutations need to be >= fixed count threshold, and >= VAF threshold
% Complex mutations only need to be >= complex threshold


complex_thresh = 6; % for complex mutations, save those with count>=fix_thresh
complex_vaf_thresh = vaf_thresh / 5; % for complex mutations, vaf threshold = standard vaf thresh/5
homopolymer_vaf_thresh = vaf_thresh * 2; % for homopolymer indels, vaf threshold = 2*standard vaf thresh

% VAF threshold for different single-base substitutions
vaf_thresh_hash = vaf_thresh * ones(4,4); % vaf_thresh_hash(wtbase,mutbase)
vaf_thresh_hash(1,4) = vaf_thresh * 2; % ATCG = 1234
vaf_thresh_hash(2,3) = vaf_thresh * 2; 
vaf_thresh_hash(4,1) = vaf_thresh * 2; 
vaf_thresh_hash(3,2) = vaf_thresh * 2; 

% Conversion yield of different plexes
CY_allplex = [0.659482759 0.407758621 0.379310345 0.318965517 0.509913793 0.240086207 0.353017241 0.337068966 0.320258621 0.824568966 0.460775862 0.377155172 0.292672414 0.102155172 0.234051724 0.579741379 0.359913793 0.331034483 0.164224138 0.556896552 0.238362069 0.363362069 0.032327586 0.378448276 0.174137931 0.634913793 0.444827586 0.275862069];

VarMolList_filter = VarMolList;
VarMol_filter = VarMol;

homopolymerlist = csvread('homopolymer_in_EnR_AMLQBDA.csv'); % 3 columns: tarnum, homopolymer startpos in EnR, homopolymer endpos in EnR
cursor_HP = 1; % homopolymerlist cursor

for curtar = 1:size(VarMolList,1)
    
    if curtar > homopolymerlist(cursor_HP,1) && cursor_HP < size(homopolymerlist,1)
        cursor_HP = cursor_HP + 1;
    end
    
    for curmut = 1:size(VarMolList,2)
        
        if ~isempty(VarMolList{curtar,curmut})
            
            [mutpos wtbase mutbase] = ConvertMutInfo(VarMolList{curtar,curmut});
            
            cur_mol_thresh = fix_thresh;
            
            if mutpos == -1 % not single base substitution
                [mutpos indel_seq homoindel] = ConvertIndelInfo(VarMolList{curtar,curmut});
                
                % complex indel and mutations
                cur_vaf_thresh = complex_vaf_thresh;
                cur_mol_thresh = complex_thresh;
                
                % homopolymer indel
                if curtar == homopolymerlist(cursor_HP,1) && homoindel == 1 && mutpos >= homopolymerlist(cursor_HP,2) && mutpos <= homopolymerlist(cursor_HP,3)
                    cur_vaf_thresh = homopolymer_vaf_thresh;
                    cur_mol_thresh = fix_thresh;
                end
            else % single base substitution
                cur_vaf_thresh = vaf_thresh_hash(wtbase,mutbase);
                cur_mol_thresh = fix_thresh;
            end
            
            curCY = CY_allplex(curtar);
            cur_VAF = VarMol(curtar,curmut) / (input_ng*290*2*curCY);
            
            % if lower than VAF threshold or count threshold, remove mutation
            if cur_VAF < cur_vaf_thresh || VarMol(curtar,curmut) < cur_mol_thresh
                VarMol_filter(curtar,curmut) = 0;
                VarMolList_filter{curtar,curmut} = [];
            end
            
        end
    end
end

% compress list
for curtar = 1:size(VarMolList_filter,1)
    allmutcount = 0;
    for curmut = 1:size(VarMolList_filter,2)
        if ~isempty(VarMolList_filter{curtar,curmut})
            allmutcount = allmutcount + 1;
            
            curmutinfo = VarMolList_filter{curtar,curmut};
            VarMolList_filter{curtar,curmut} = [];
            
            curmutcount = VarMol_filter(curtar,curmut);
            VarMol_filter(curtar,curmut) = 0;
            
            VarMolList_filter{curtar,allmutcount} = curmutinfo;
            VarMol_filter(curtar,allmutcount) = curmutcount;
            
        end
        
    end
end


end


function [mutpos wtbase mutbase] = ConvertMutInfo(mutstr)

mutpos = -1;
wtbase = -1;
mutbase = -1;

if ~isempty(strfind(mutstr, 'Ins')) || ~isempty(strfind(mutstr, 'Del')) || ~isempty(strfind(mutstr, ','))
    return % indel or multiple mutations
else
    temp_info = strsplit(mutstr,'>');
    mutpos = str2num(temp_info{1}(1:end-1));
    
    if isempty(mutpos) || length(temp_info{2})~=1 % not single base mutation
        mutpos = -1;
        wtbase = -1;
        mutbase = -1;
        return
    end
    
    curwtbase = temp_info{1}(end);
    curmutbase = temp_info{2}(1);
    wtbase = GCAT2num(curwtbase);
    mutbase = GCAT2num(curmutbase);
end


end


function [mutpos indel_seq homoindel] = ConvertIndelInfo(mutstr)

mutpos = -1;
indel_seq = [];
homoindel = 1; % default "True", whether indel sequence only consists of the same base

ins_ind = strfind(mutstr, 'Ins');
del_ind = strfind(mutstr, 'Del');

if ~isempty(strfind(mutstr, ','))
    return % complex mutation
elseif ~isempty(ins_ind)
    mutpos = str2num(mutstr(1:ins_ind-1));
    indel_seq = mutstr(ins_ind+3:end);
elseif ~isempty(del_ind)
    mutpos = str2num(mutstr(1:del_ind-1));
    indel_seq = mutstr(del_ind+3:end);
end

if ~isempty(indel_seq)
    if length(unique(indel_seq)) > 1
        homoindel = 0;
    end
    
end


end


function curnum = GCAT2num(curbase)

curnum = 0;

if curbase == 'A' || curbase == 'a'
    curnum = 1;
elseif curbase == 'T' || curbase == 't'
    curnum = 2;
elseif curbase == 'C' || curbase == 'c'
    curnum = 3;
elseif curbase == 'G' || curbase == 'g'
    curnum = 4;
end

end