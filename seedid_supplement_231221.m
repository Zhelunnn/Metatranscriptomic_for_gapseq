function [output_model,core_reactions,seedrxn_notincorrected_formatted,meta_idof_matchedseed_id_instep_sum,output_unfound_metaid,notexist_in_model_list]=seedid_supplement_231221(varargin)
''''''' trying to add those found and expressed (default top 100 % expression level) reactions based on biocyc id, but were not included in the model due to metacyc db inconsistencies with modelseed db, these reactions are listed in binname_expressedRxn_notinModel.tsv''''';
'''''''''''instructions: step1:get the top 90 % good_blast expressed reactions(metacyc id)''''';
'''''''step2:link the metacyc id as many as possible with seed id using 1) metacyc id to kegg to seed id 2) metacyc id to seed id using mnxref_seed-other.tsv. ';
'''''''''''''3) metacyc id to EC number to seed id(which is not accurate), use modelseed_db_reactions.tsv to correct it. 4)metacyc id to seed id using modelseed_db_reactions.tsv';
'''''''''''''5) and 6) concatentate results obtained above and retrieve reactions information from all-Reactions.tbl file';
'''''''step3:add matched seed reactions into the model, using seed_reactions_corrected.formatted_reat_add.tsv and modelseed_db_github_formatted4curation.tsv';
%   add those reactions that are found by gapfind (have metacyc id) and expressed according to reads count file but are not added in the model 
%   try to find the matched seed rxns using gapseq database, those failed reactions are stored in optput_unfound_metaid

% input
%       reaction_reads_count='bin_14_061016_reaction_reads_count.tsv';
%       seed_reactions_corrected='seed_reactions_corrected.tsv';
%       mnxref_seed_other='mnxref_seed-other.tsv';
%       meta_rea='meta_rea.tbl';
%       seed_Enzyme_Class='seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv';
%       model='bin.14_connected_edges-draft.xml';
%       all_reactions_file='bin.14_connected_edges-all-Reactions.tbl'
%       seed_reactions_corrected_formatted="seed_reactions_corrected.formatted_reat_add.tsv";
%       expressed_level='0.9';

% output
%       outputmodel:    model with as many as possible reactions added (reactions can not directly get seed id and reactions that are not included in the model)
%       output_unfoundd_metaid: a table containing the reactions can not find the seed id by any means
% only focus on the top 90% highly expressed reactions(selectedRows_good_blast_dereplicate)

%{ 
reaction_reads_count='bin_14_061016_reaction_reads_count.tsv';
seed_reactions_corrected='seed_reactions_corrected.tsv';
mnxref_seed_other='mnxref_seed-other.tsv';
meta_rea='meta_rea.tbl';
seed_Enzyme_Class='seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv';
model='bin.14_connected_edges-draft.xml';
all_reactions_file='bin.14_connected_edges-all-Reactions.tbl'
seed_reactions_corrected_formatted="seed_reactions_corrected.formatted_reat_add.tsv";
expressed_level='0.9'; 
%}


% 1. match the metacyc id with kegg id (1 or 2) using meta_rea file
% 2. match the kegg id (1 or 2) with seed id using both seed_reaction_corrected and mnxref_seed_other file
function [output,meta_id_notin_metarea,no_kegg_ids]=extract_seedid(no_seedrxn,mnxref_seed_db,seedid_meta_rea,rxn_equation_db,meta_id_notin_metarea,no_kegg_ids)
    row=seedid_meta_rea(ismember(strrep(seedid_meta_rea.id,'|',''),no_seedrxn),:);% no_seedrxn is metacyc rxn id
    output='';
    % if the metacyc id is not in meta_rea.tbl file, return '' and add id in meta_id_notin_metarea
    if isempty(row)
        meta_id_notin_metarea{end+1}=no_seedrxn;
        return;
    end
    % search the matched seed id by kegg id
    if isempty(row.kegg{1})
        no_kegg_ids{end+1}=no_seedrxn;
        return;
    else
        kegg_id=strsplit(row.kegg{1},',');
        if length(kegg_id)==1
            % find the matched row in seed reaction corrected file using kegg_id
            % Find the indices of the matching rows
            matchingIndices = ismember(rxn_equation_db.abbreviation, kegg_id{1});
            matchingIndices_mnxref=ismember(mnxref_seed_db.other,kegg_id{1});

            % Check if there are any matches (multiple matched rows) and get only the first one
            if any(matchingIndices)
                firstMatchIndex = find(matchingIndices, 1, 'first');
                matched_row = rxn_equation_db(firstMatchIndex, :);
                output=matched_row.id{1};
            % check the kegg in mnxref_seed_db other column
            elseif any(matchingIndices_mnxref)
                firstMatchIndex_mnxref = find(matchingIndices_mnxref, 1, 'first');
                matched_row_mnxref = mnxref_seed_db(firstMatchIndex_mnxref, :);
                output=matched_row_mnxref.seed{1};
            else
                return
            end

        % sometimes one metacyc id can link to two kegg id based on meta_rea file
        elseif length(kegg_id)==2
            % find the matched index using rxn_equation_db.abbreviation
            matchedIndices_1=ismember(rxn_equation_db.abbreviation,kegg_id{1});
            matchedIndices_2=ismember(rxn_equation_db.abbreviation,kegg_id{2});
            % combined the two logical array of matched result in seed_reaction_corrected.tsv file
            combinedMatchedIndices = matchedIndices_1 | matchedIndices_2;

            % find the matched index using mnxref_seed_db.other
            matchedIndices_3=ismember(mnxref_seed_db.other,kegg_id{1});
            matchedIndices_4=ismember(mnxref_seed_db.other,kegg_id{2});
            % combined the two logical array of matched result in mnxref_seed_other.tsv file
            combinedMatchedIndices_mnxref = matchedIndices_3 | matchedIndices_4;

            if all(~combinedMatchedIndices) & all(~combinedMatchedIndices_mnxref)
                return

            elseif any(combinedMatchedIndices_mnxref) % 1.for mnxref_seed_other, only use the first matched rxn
                firstMatchIndex_mnxref = find(combinedMatchedIndices_mnxref, 1, 'first');
                matched_row_mnxref = mnxref_seed_db(firstMatchIndex_mnxref, :);
                output=matched_row_mnxref.seed{1};

            else % 2.mete_rea, only use the first matched rxn
                firstMatchIndex_combined = find(combinedMatchedIndices, 1, 'first');
                matched_row_combined = rxn_equation_db(firstMatchIndex_combined, :);
                output=matched_row_combined.id{1};
            end
        else
            disp('There are ',length(kegg_id),' kegg ids for this metacyc rxn, adjust it manually')
            %output=rxn_equation_db(ismember(rxn_equation_db.abbreviation{1},kegg_id{1}),:).id{1};
        end
    end
end



%%%% parameter without default value
p=inputParser;
addParameter(p, 'reaction_reads_count',NaN);
addParameter(p, 'model',NaN);
addParameter(p, 'all_reactions_file',NaN);
addParameter(p, 'table_output',false); % if true a table contain all the unfound metacyc id will be written in a table
addParameter(p, 'output_file_name','output.tsv');

%%%% parameter with default value 
addParameter(p, 'expressed_level',1,@isnumeric);
if isfolder('/srv/scratch/z5245780/software/gapseq/gapseq_1.3') % if the script is running on katana
    addParameter(p, 'mnxref_seed_other','/srv/scratch/z5245780/software/gapseq/gapseq_1.3/dat/mnxref_seed-other.tsv');
    addParameter(p, 'seed_reactions_corrected','/srv/scratch/z5245780/software/gapseq/gapseq_1.3/dat/seed_reactions_corrected.tsv');
    addParameter(p, 'meta_rea','/srv/scratch/z5245780/software/gapseq/gapseq_1.3/dat/meta_rea.tbl');
    addParameter(p, 'seed_Enzyme_Class','/srv/scratch/z5245780/software/gapseq/gapseq_1.3/dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv');
    addParameter(p, 'seed_reactions_corrected_formatted','/srv/scratch/z5245780/DB/custom/seed_reactions_corrected.formatted_reat_add.tsv');
    addParameter(p, 'seed_reactions_github_formatted','/srv/scratch/z5245780/DB/custom/modelseed_db_github_formatted4curation.tsv');
    seed_rxn_github=readtable('/srv/scratch/z5245780/DB/custom/modelseed_db_reactions.tsv',"FileType","text","Delimiter",'\t','ReadVariableNames',true,'VariableNamingRule','preserve');

else % if the scriot is running locally
    addParameter(p, 'mnxref_seed_other','mnxref_seed-other.tsv');
    addParameter(p, 'seed_reactions_corrected','seed_reactions_corrected.tsv');
    addParameter(p, 'meta_rea','meta_rea.tbl');
    addParameter(p, 'seed_Enzyme_Class','seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv');
    addParameter(p, 'seed_reactions_corrected_formatted','seed_reactions_corrected.formatted_reat_add.tsv');
    addParameter(p, 'seed_reactions_github_formatted','modelseed_db_github_formatted4curation.tsv');
    seed_rxn_github=readtable('modelseed_db_reactions.tsv',"FileType","text","Delimiter",'\t','ReadVariableNames',true,'VariableNamingRule','preserve');

end
    addParameter(p, 'expressed_rxn_notinmodel_file_generated',NaN);% if true a table contain all the expressed seed id but not included in the model will be written in a table

parse(p, varargin{:});

reaction_reads_count = p.Results.reaction_reads_count;
model = p.Results.model;
all_reactions_file=p.Results.all_reactions_file;
table_output_logical=p.Results.table_output;
output_file_name=p.Results.output_file_name;

expressed_level=p.Results.expressed_level;
mnxref_seed_other=p.Results.mnxref_seed_other;
seed_reactions_corrected=p.Results.seed_reactions_corrected;
meta_rea=p.Results.meta_rea;
seed_Enzyme_Class=p.Results.seed_Enzyme_Class;
seed_reactions_corrected_formatted=p.Results.seed_reactions_corrected_formatted;
seed_reactions_github_formatted=p.Results.seed_reactions_github_formatted;
% read the summary table of each model and preserve the column name 
file_reads_count=readtable(reaction_reads_count, 'FileType', 'text', 'Delimiter', '\t','ReadVariableNames', true,'VariableNamingRule','preserve');
expressed_rxn_notinmodel_file_generated=p.Results.expressed_rxn_notinmodel_file_generated;
%file_reads_count(1:4,:)=[]; % remove the first 4 rows 

% ##################### prepare the expressed table ###################
% ########## only consider the rxn annotated with good blast from reads_count file/ summary file ###########
selectedRows_good_blast = file_reads_count(contains(file_reads_count.Biocyc_ID, 'good_blast'), :);

% ####### disregard of the enzyme subunit (dereplicate the table based on rxn id (characters before the first underscore in first column) ############
parsedIDs_metacyc= cellfun(@(x) strtok(x, '_'), selectedRows_good_blast.Biocyc_ID, 'UniformOutput', false);
[~, indices] = unique(parsedIDs_metacyc); % unique returns the unique entries from this list, along with the indices of their first occurrences the result will be reordered
%selectedRows_good_blast_dereplicate =sortrows(selectedRows_good_blast(indices, :),12,'descend'); % Create a new table with only the unique entries for reads count file
selectedRows_good_blast_dereplicate = sortrows(selectedRows_good_blast(indices, :),'TPM_mean','descend'); % Create a new table with only the unique entries for summary file

% get the top 90 % /100 % expressed reactions    
if isempty(selectedRows_good_blast_dereplicate.TPM_mean)|ismember(selectedRows_good_blast_dereplicate.TPM_mean,'NA')
    disp('TPM_mean is empty');
    return;
else
    TPM_cumsum=cumsum(selectedRows_good_blast_dereplicate.TPM_mean); %cumulative sum
    TPM_sum=sum(selectedRows_good_blast_dereplicate.TPM_mean);
    index_top90=find(TPM_cumsum>=TPM_sum*(expressed_level-1e-10),0.9,'first'); % 1e-10 a small tolerance to avoid float error
end

selectedRows_good_blast_dereplicate=selectedRows_good_blast_dereplicate(1:index_top90,:);

% ######################################################################
% operation for expressed table: get the biocyc id from expressed reactions file by extract the character before the first underscore
biocyc_ids_expressed=cellfun(@(x) strtok(x,'_'),selectedRows_good_blast_dereplicate{:,1},'UniformOutput',false);

% the metacyc rxn id that have (seedrxn_expressed_table) or don't have (no_seedrxn_expressed_table) matched modelseed id
no_seedrxn_expressed_table=selectedRows_good_blast_dereplicate(strcmp(selectedRows_good_blast_dereplicate.modelseed_id,''),:);

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''';
'''''''get seed ids of those metacyc reactions that dont have the matched seed id in gapseq database''''''''''''';
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''';
% 1.find the seed id by meta_rea file and mnxref_seed_other (found_metacycid_meta_rea)
no_seedrxns=cellfun(@(x) strtok(x,'_'),no_seedrxn_expressed_table{:,1},'UniformOutput',false);

mnxref_seed_db=readtable(mnxref_seed_other,'FileType','text','Delimiter','\t','ReadVariableNames', true,'VariableNamingRule','preserve');
rxn_equation_db=readtable(seed_reactions_corrected,'FileType','text','Delimiter','\t','ReadVariableNames', true,'VariableNamingRule','preserve');
seedid_meta_rea=readtable(meta_rea,'FileType','text','Delimiter','\t','ReadVariableNames', true,'VariableNamingRule','preserve');

meta_id_notin_metarea={}; % metacyc id not included in meta_rea file
no_kegg_ids={};
[found_seedid_meta_rea, meta_id_notin_metarea, no_kegg_ids]=cellfun(@(x) extract_seedid(x, mnxref_seed_db, seedid_meta_rea, rxn_equation_db, meta_id_notin_metarea, no_kegg_ids), no_seedrxns,'UniformOutput',false );

'''''''''''''''make the matched metacyc id and seed id cellarray''''''''''''''';
no_empty_indice=find(cellfun(@(x) ~isempty(x),found_seedid_meta_rea)); % indice of seed id 
step1_seed_id=found_seedid_meta_rea(no_empty_indice);
step1_metacyc_id=no_seedrxns(no_empty_indice);% Retrieve the MetaCyc IDs of reactions that correspond to the SEED IDs identified in Step 1.

meta_idof_matchedseed_id_instep1=cell(length(step1_metacyc_id),2); % store the metacyc and matched seed id
meta_idof_matchedseed_id_instep1(1:length(step1_metacyc_id),1)=step1_metacyc_id;
meta_idof_matchedseed_id_instep1(1:length(step1_seed_id),2)=step1_seed_id;

%disp(meta_idof_matchedseed_id_instep1)
% add the TPM_mean to the second column of found_seedid_meta_rea
for i=1:length(found_seedid_meta_rea)
    ele=found_seedid_meta_rea(i);
    if ~isempty(ele{1})
        found_seedid_meta_rea{i,2}=no_seedrxn_expressed_table{i,'TPM_mean'};
    end
end

% metacyc rxn id that can not find the matched seed id using meta_rea file
% get the unfound metacyc id first and then remove the empty element
unfound_metacycid_metarea=no_seedrxns(cellfun(@isempty, found_seedid_meta_rea(:,1)));
% remove the empty element 
meta_id_notin_metarea=meta_id_notin_metarea(~cellfun(@isempty, meta_id_notin_metarea));
found_seedid_meta_rea_nonempty=found_seedid_meta_rea(~cellfun(@isempty,found_seedid_meta_rea(:,1)),:);

% 2.find the seed id of those unfound seed id rxn in last step by mnxref_seed-other.tsv file 
found_seed_id_mnxref_seed={};
meta_idof_matchedseed_id_instep2={}; % the MetaCyc IDs of reactions that correspond to the SEED IDs identified in Step 2
for meta_id_index=1:length(unfound_metacycid_metarea)
    matched_index=ismember(mnxref_seed_db.other,unfound_metacycid_metarea{meta_id_index});
    seed_id=mnxref_seed_db.seed(matched_index);

    if ~isempty(seed_id)
        % find the index of meta id that has the seed id in mnxref_seed file and add the tpm_mean
        tpm_index=ismember(no_seedrxns,unfound_metacycid_metarea{meta_id_index});
        tpm_mean=no_seedrxn_expressed_table{tpm_index,'TPM_mean'};

        nextRow = size(found_seed_id_mnxref_seed, 1) + 1;
        found_seed_id_mnxref_seed{nextRow,1} = seed_id;
        found_seed_id_mnxref_seed{nextRow,2} = tpm_mean;

        '''''''''''''''make the matched metacyc id and seed id cellarray''''''''''''''';       
        meta_idof_matchedseed_id_instep2{nextRow,1}=unfound_metacycid_metarea{meta_id_index,1};% Retrieve the MetaCyc IDs of reactions that correspond to the SEED IDs identified in Step 2 
        meta_idof_matchedseed_id_instep2{nextRow,2}=seed_id{1};
    end
end
%disp(meta_idof_matchedseed_id_instep2)
% 3. find the seed id from last step using meta_rea and seed_Enzyme_Class_Reactions_Aliases_unique_edited
%First, retrieve the EC (Enzyme Commission) number using the 'meta_rea' dataset. Then, locate the corresponding SEED identifier by searching in the 'seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv' file.

EC_seed_db=readtable(seed_Enzyme_Class,'FileType','text','Delimiter','\t','ReadVariableNames', true,'VariableNamingRule','preserve');
no_EC_metaid={}; % store the meta rxn id that match with EC number
meta_id_withEC_butnoseed={}; % metacyc id that can find the eEC number but can not find seed id using seed_Enzyme_Class_Reactions_Aliases_unique_edited
seed_ids_found={};
meta_idof_matchedseed_id_instep3={}; % the MetaCyc IDs of reactions that correspond to the SEED IDs identified in Step 3

for i = 1:length(unfound_metacycid_metarea)
    meta_id_ec = unfound_metacycid_metarea{i};  % Extract meta_cyc_id from cell
    
    % remove the | in seedid_meta_rea and find the match index
    matched_index = ismember(strrep(seedid_meta_rea.id, '|', ''), meta_id_ec); 
    ECs = seedid_meta_rea.ec(matched_index);

    if ~isempty(ECs)
        combined_EC_index = false(size(EC_seed_db.("External ID")));
        EC_split = strsplit(ECs{1}, ',');
        
        % one EC may contain two EC numbers
        for j=1:length(EC_split)            
            EC_index=ismember(EC_seed_db.("External ID"),strrep(EC_split{j},'EC-',''));
            if any(EC_index)
                combined_EC_index = combined_EC_index | EC_index;
                seed_ids=EC_seed_db.("MS ID")(combined_EC_index);

                next_row=size(seed_ids_found,1)+1;
                seed_ids_found{next_row,1}=strsplit(seed_ids{1},'|');
                seed_ids_found{next_row,2}=no_seedrxn_expressed_table{ismember(no_seedrxns,meta_id_ec),'TPM_mean'};

                '''''''''''''''make the matched metacyc id and seed id cellarray''''''''''''''';       
                meta_idof_matchedseed_id_instep3{next_row,1}=meta_id_ec;
                meta_idof_matchedseed_id_instep3{next_row,2}=seed_ids{1};
            else
                meta_id_withEC_butnoseed{end+1}=meta_id_ec;
            end
        end
    else
        no_EC_metaid{end+1}=meta_id_ec;
    end
end
%disp(meta_idof_matchedseed_id_instep3)
unfound_metaid=[meta_id_withEC_butnoseed,no_EC_metaid];

''''''''''''''''' there are some issues in step3 because of database version or gapseq''''''''''''';
''''''''''''''''' correct the step 3 matching result ''''''''''''''''';
% using aliases in modelseed_db_reactions.tsv to get the metacyc rxn id

for i=1:size(meta_idof_matchedseed_id_instep3,1)
    metacyc_id=meta_idof_matchedseed_id_instep3{i,1};
    
    indice=cellfun(@(x) contains(x,metacyc_id),seed_rxn_github.aliases,'UniformOutput', true);
    if any(indice)
        meta_idof_matchedseed_id_instep3{i,2}=seed_rxn_github.id{indice,1};% replace the old seed id with corrected seed id in matching table
        seed_ids_found{i,1}=seed_rxn_github.id(indice);% replace the old seed id with corrected seed id in id found cellarray (seed id will be added into models)
        
    end
end

% 4. find the correspongding modelseed id against unfound metacyc id using modelseed db reactions file from modelseed github
%seed_rxn_github=readtable('/srv/scratch/z5245780/DB/custom/modelseed_db_reactions.tsv',"FileType","text","Delimiter",'\t','ReadVariableNames',true,'VariableNamingRule','preserve');

% some seed rxn may in this form metaid.c/metaid.p so include all these rxns if metacyc id is found
rxnid_otherDB = cellfun(@(x) strsplit(x, '.'), seed_rxn_github.abbreviation, 'UniformOutput', false);
metacycid_beforedot = cellfun(@(x) x{1}, rxnid_otherDB, 'UniformOutput', false); % get all rxn if metacyc id is found
metacyc_id_aliases=cellfun(@(x) regexp(x,'MetaCyc:\s*(RXN-\d+)','tokens'),seed_rxn_github.aliases,'UniformOutput', false);
% retrieve the string in the nested cell
function out = extractContentOrEmptyString(cellArray)
    if iscell(cellArray) && ~isempty(cellArray)
        out = cellArray{1}{1};  % Extract the first element if it's a non-empty cell array
    else
        out = '';           % Return an empty string for empty cells
    end
end
meta_id_aliases_unnested=cellfun(@(x) extractContentOrEmptyString(x), metacyc_id_aliases,'UniformOutput', false);

meta_idof_matchedseed_id_instep4={}; % the MetaCyc IDs of reactions that correspond to the SEED IDs identified in Step 4
seed_id_modelseedgithub={};% store the mapped seed id 
unfound_metaid_final = {};
for i=1:length(unfound_metaid)
    meta_id=unfound_metaid{i};
    indices=ismember(metacycid_beforedot,meta_id);
    if any(indices)
        seed_id_github=seed_rxn_github.id{indices};
        seed_id_modelseedgithub{end+1,1}=seed_id_github;
        
        '''''''''''''''make the matched metacyc id and seed id cellarray'''''''''''''''; 
        next_row=size(meta_idof_matchedseed_id_instep4,1)+1;
        meta_idof_matchedseed_id_instep4{next_row,1}=meta_id;
        meta_idof_matchedseed_id_instep4{next_row,2}=seed_id_github;
    else
        % check metacyc id in aliases column
        indices_aliases = false(size(meta_id_aliases_unnested));

        % Iterate over each element in meta_id_aliases_unnested,unfinished!!!
        for j = 1:length(meta_id_aliases_unnested)
            % Split the string by comma to get individual IDs
            splitIds = strsplit(meta_id_aliases_unnested{j}, ',');

            % Check if meta_id is one of the split IDs
            indices_aliases(j) = any(ismember(splitIds, meta_id));
            
        end        
        if any(indices_aliases)
            seed_id_github=seed_rxn_github.id{indices_aliases};
            seed_id_modelseedgithub{end+1,1}=seed_id_github;
       
            '''''''''''''''make the matched metacyc id and seed id cellarray'''''''''''''''; 
            next_row=size(meta_idof_matchedseed_id_instep4,1)+1;
            meta_idof_matchedseed_id_instep4{next_row,1}=meta_id;
            meta_idof_matchedseed_id_instep4{next_row,2}=seed_id_github;
        else
            unfound_metaid_final{end+1}=meta_id;
        end
    end
end
%disp(meta_idof_matchedseed_id_instep4)



% 5.merge found seed id in above steps
if isempty(found_seed_id_mnxref_seed)
    found_seed_id_mnxref_seed_IDs=seed_id_modelseedgithub;
else
    found_seed_id_mnxref_seed_IDs=found_seed_id_mnxref_seed(:,1);% extract found id from step 2
    found_seed_id_mnxref_seed_IDs=[found_seed_id_mnxref_seed_IDs;seed_id_modelseedgithub];
end
meta_idof_matchedseed_id_instep_sum=[meta_idof_matchedseed_id_instep1;meta_idof_matchedseed_id_instep2;meta_idof_matchedseed_id_instep3;meta_idof_matchedseed_id_instep4]
% extract and flatten the found id from step 3 and step 1
if ~isempty(seed_ids_found)
    seed_ids_found_IDs = cellfun(@(c) c(:), seed_ids_found(:,1), 'UniformOutput', false);
else
    seed_ids_found_IDs = {}; % Assign an empty cell array or handle it accordingly
end
% Ensure intermediate variables are non-empty and vertically concatenable
part1 = found_seedid_meta_rea(~cellfun(@isempty, found_seedid_meta_rea(:, 1)));
part2 = vertcat(found_seed_id_mnxref_seed_IDs{:});
part3 = vertcat(seed_ids_found_IDs{:});

% Handle empty parts
if isempty(part1), part1 = {}; end
if isempty(part2), part2 = {}; end
if isempty(part3), part3 = {}; end
% Concatenate parts
final_seed_ids_founds = [part1; part2; part3];



% 6.retrieve the reactions information from all-Reactions.tbl file generated by gapseq for unfound meta id rxns
all_reactions_tbl=readtable(all_reactions_file,"FileType","text","Delimiter",'\t','ReadVariableNames',true,'VariableNamingRule','preserve');
all_reactions_tbl=all_reactions_tbl(ismember(all_reactions_tbl.status,'good_blast'),:); % only keep the good_blast reactions

unfound_reactions_information=cellfun(@(x) all_reactions_tbl(find(ismember(all_reactions_tbl.rxn,x),1,'first'),:),unfound_metaid_final,UniformOutput=false);
output_unfound_metaid=vertcat(unfound_reactions_information{:});

if table_output_logical
    writetable(output_unfound_metaid, output_file_name, 'Delimiter', '\t', 'FileType', 'text');
end


'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''';
'''''''''''''''''''''''''add all the missed seed reactions in the model''''''''''''''''''''''';


addpath('/srv/scratch/z5245780/software/cobratoolbox');
initCobraToolbox;
model=readCbModel(model);
% operation for model:remove the META: from model.rxnbiocycid of each cell
modified_biocyc=cellfun(@(x) strrep(x, 'META:',''), model.rxnBioCycID, 'UniformOutput',false);

% get the serial number of matched rxn from model and switch result from cell to matrix
rxn_serial=cell2mat(cellfun(@(x) find(strcmp(x,modified_biocyc),1),biocyc_ids_expressed,'UniformOutput',false));
% note: the entry in expressed reaction file also included the subunit difference.. 
% so there must be rxn with more than one subunit, which means there are same biocyc Id in variable biocyc_ids_expressed...
% I use setdiff to get the matched id and unmatched id from two cell, the total count of both will be less than the count of biocyc_ids_expressed

% the reactions that can not find the matched seedid using model.rxnbiocycID
rxn_added_manu=setdiff(biocyc_ids_expressed,modified_biocyc(rxn_serial));

%%%%%%%%%% add the expressed reactions directly to the model that missed in the link %%%%%%%%%
rxn_db_file = readtable(seed_reactions_corrected_formatted, 'Delimiter', '\t', 'FileType', 'text');
seed_github_formatted4curation_db=readtable(seed_reactions_github_formatted, 'Delimiter', '\t', 'FileType', 'text');

% the metacyc rxn id that have matched modelseed id (seedrxn_expressed_table)
seedrxn_expressed_table=selectedRows_good_blast_dereplicate(~strcmp(selectedRows_good_blast_dereplicate.modelseed_id,''),:)
seedrxn_tobeadded=rowfun(@(x) strsplit(char(x),' '),seedrxn_expressed_table,'InputVariables', 'modelseed_id','OutputFormat', 'cell');
% Flatten each cell's content
flattenedContent = cellfun(@(x) x(:), seedrxn_tobeadded, 'UniformOutput', false);
% Concatenate to get a single list
flattened_seedrxn_tobeadded = unique(vertcat(flattenedContent{:}));

% get the rxn that are not present in model.rxns
notexist_in_model_list=flattened_seedrxn_tobeadded(~ismember(strcat(flattened_seedrxn_tobeadded, '_c0'),model.rxns));

notexist_in_model_list=[notexist_in_model_list;final_seed_ids_founds];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% output the reactions that are expressed but not in the model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

''''' before ouput seed rxns that are expressed but not in model, further check if it has counterparts included in the model''''';
 % to do this, 1. get its biocyc id in all-Reactions.tbl file. !! the all-REactions.tbl file need to be dereplicate based on biocycid 2. get seed id in the dbhit under this metacyc id
 % 3. check if any seed id except this seed rxn is already in the model, if not append this seed id in the list
 
 % dereplicate the all-Reactions.tbl based on rxn column
 rxnColumn = all_reactions_tbl.rxn;
 [~, uniqueIndices] = unique(rxnColumn, 'stable'); % 'stable' keeps the first occurrence
    
 % Keep only the rows corresponding to the unique indices
 all_reactions_tbl_deduplicated = all_reactions_tbl(uniqueIndices, :);

 list_refined = {}; % the list contain expressed seed id that are not included in the model, this time also excluded the counterparts of reactions in the models
 disp(model.rxns)
 for i = 1:length(notexist_in_model_list)
    rxn = notexist_in_model_list{i};  % Get the seed reaction
    disp(rxn)
    for j = 1:length(all_reactions_tbl_deduplicated.rxn)  % Loop through all reactions in the all_reactions_tbl
        rxns = all_reactions_tbl_deduplicated.dbhit{j};  % Get the dbhit for this reaction
        rxns_split = strsplit(rxns, ' ');

        % Check if the seed reaction is in the dbhit
        if ismember(rxn, rxns_split)
              % Split the dbhit by space
            
            % Create a flag that will check if any of the other reactions are in the model
            all_other_reactions_not_in_model = true;  % Assume that all others are not in the model
            
            for k = 1:length(rxns_split)  % Loop through the split reactions
                if ~strcmp(rxns_split{k}, rxn)  % If it's not the current rxn

                    % a self check to avoid including counterparts in list_refined but always keeping the first reaction in dbhit
                    if any(ismember(rxns_split{k},list_refined))
                        all_other_reactions_not_in_model = false; 
                        break
                    end
                    % Check if this reaction is in the model
                    if any(ismember([rxns_split{k},'_c0'], model.rxns))
                        all_other_reactions_not_in_model = false;  % If any other reaction is in the model, set flag to false
                        break;  % No need to check further, exit the loop
                    end
                end
            end

            % If all other reactions are not in the model, add the current rxn to the list
            if all_other_reactions_not_in_model
                list_refined{end+1} = rxn;  % Append to the refined list
            end
        end
    end
end

% generate the file containg rxns that are expressed but not in model ##### 26112024 added
if ~isnan(expressed_rxn_notinmodel_file_generated)
    list_refined=unique(list_refined);
    %writecell(notexist_in_model_list, expressed_rxn_notinmodel_file_generated, 'FileType', 'text', 'Delimiter', '\t');
    writecell(list_refined, expressed_rxn_notinmodel_file_generated, 'FileType', 'text', 'Delimiter', '\t'); % counterparts already have been excluded
    %writetable(seedrxn_expressed_table, 'expressed_rxn_inmodel.tsv', 'FileType', 'text', 'Delimiter', '\t'); % didn't exluded those counterparts
end



absent_rxn_ids_in_db = {}; % for the reactions that are not included in seed_reaction_corrected file
%length(notexistlist)/3
%{
for i=1:length(notexist_in_model_list)
    matched_row_fromrxndbfile=rxn_db_file(strcmp(rxn_db_file.Rid,strcat(notexist_in_model_list{i},'_zhelun')),:);

    % check if the reaction is included in seed_reaction_corrected file
    if isempty(matched_row_fromrxndbfile)
        try
            matched_seed_github_formatted=seed_github_formatted4curation_db(strcmp(seed_github_formatted4curation_db.Rid,strcat(notexist_in_model_list{i},'_zhelun')),:);
        catch exception
            matched_seed_github_formatted='';
            disp(notexist_in_model_list{i})
            fprintf('An error occurred: %s\n', exception.message);
        end
        if isempty(matched_seed_github_formatted)
            absent_rxn_ids_in_db{end+1} = notexist_in_model_list{i};
            continue;
        else
            matched_row_fromrxndbfile= matched_seed_github_formatted;
        end
    end

    % extract rxn information for addReaction
    %rxnid=matched_row_fromrxndbfile.Rid{1};
    rxnid=strrep(matched_row_fromrxndbfile.Rid{1},'_zhelun','_c0');
    reactionName=matched_row_fromrxndbfile.react_name{1};
    
    % Split met_id into a cell array
    metaboliteList= strsplit(matched_row_fromrxndbfile.met_id{1}, ',');
    stoichCoeffList=str2double(strsplit(matched_row_fromrxndbfile.met_Scoef{1}, ','));
    reversible=strcmp(matched_row_fromrxndbfile.react_rev{1},'true');
    lowerBound=matched_row_fromrxndbfile.lowbnd;
    upperBound=matched_row_fromrxndbfile.uppbnd;
    %fprintf('rxnid: %s\nreactionName: %s\nmetaboliteList: %s\nstoichCoeffList: %s\nreversible: %d\nlowerBound: %f\nupperBound: %f\n', ...
    %rxnid, reactionName, strjoin(metaboliteList, ', '), strjoin(arrayfun(@num2str, stoichCoeffList, 'UniformOutput', false), ', '), reversible, lowerBound, upperBound);
    model=addReaction(model,rxnid,'reactionName',reactionName,'metaboliteList',metaboliteList,'stoichCoeffList',stoichCoeffList,'reversible',reversible,'lowerBound',lowerBound,'upperBound',upperBound);
end 
%}
% note: new metabolites were not introduced into the model when 3 new reactions were added, 
% which means the form of met ids in added reactions are corresponding to id form in model


output_model=model;
%tissuemodel=fastcore(model,combined_indices)

%%%%%%%%%instruction: after as many as possible reactions supplemented into the model,we need to retrieve the core reactions of the model
%%%%%% retrieve the core reactions (top 90 % expressed reactions as many as possible %%%%%% 
% retrieve the expressed seed rxn id from variable 'selectedRows_good_blast_dereplicate', modelseed_id_unique can be the core set reactions
modelseed_id_cell=cellfun(@(x) strsplit(x,' '),selectedRows_good_blast_dereplicate.modelseed_id,'UniformOutput', false);
modelseed_id=cellfun(@(x) x(:), modelseed_id_cell, 'UniformOutput', false);

core_reactions_raw=unique(vertcat(modelseed_id{:})); % modelseed id of corerxns
core_reactions=core_reactions_raw(~cellfun(@(x) isempty(x), core_reactions_raw));

seedrxn_notincorrected_formatted=absent_rxn_ids_in_db; % reactions that are not included by the seed_reactions_corrected_formatted file

end

%{
% read and combine all the unfound table and dereplicate the combined table based on first column
bin_names = 'bin_3 bin_4 bin_5 bin_7 bin_11 bin_12 bin_13 bin_14 bin_15 bin_16 bin_17 bin_19 bin_25 bin_26 bin_32 bin_34 bin_35 bin_37 bin_49';
bin_names_list = strsplit(bin_names, ' ');
combine = [];

for i = 1:length(bin_names_list)
    dataTable = readtable(strcat(bin_names_list{i}, '_unfound_metaid.tsv'), 'FileType', 'text', 'Delimiter', '\t', 'VariableNamingRule', 'preserve');
    
    % Check if 'ec' column is not a cell, and convert it to a cell if needed
    if ~iscell(dataTable.ec)
        dataTable.ec = num2cell(dataTable.ec);
    end
    
    % Add a new column 'File' with the current file name
    dataTable.File = repmat({bin_names_list{i}}, size(dataTable, 1), 1);

    combine = [combine; dataTable];
end

[~, uniqueIdx] = unique(combine(:,1), 'rows', 'stable');
deduplicatedData = combine(uniqueIdx, :);
writetable(deduplicatedData, 'combined_deduplicated.tsv', 'FileType', 'text', 'Delimiter', '\t');
%}