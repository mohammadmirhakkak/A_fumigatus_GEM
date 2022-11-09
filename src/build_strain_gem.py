import pandas as pd
import cobra
from cobra import Model, Metabolite, Reaction
import re
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis import fastcc
import copy
import timeit
import pickle
import glob


start = timeit.default_timer()

#function -> gene ID unification (model genes vs predicted genes)
def unify_gene_id(matrix):
    new_index = list()
    for i in range(matrix.shape[0]):
        s = matrix.index[i]
        s = 'AFUA_'+s[3:].upper()
        new_index.append(s)
    matrix.index = new_index
    return matrix

#function -> make model based on gene presence data using a draft model containing exchanges and biomass related reactions
def identify_rxn_presence(model,model_info,gene_presence):

    for i in range(gene_presence.shape[0]):
        sub_df = model_info[model_info['GENE ASSOCIATION'].str.contains(gene_presence.index[i])]
        related_rxn_locs = list(sub_df.index)

        for j in related_rxn_locs:
            s = model_info.loc[j,'GENE ASSOCIATION']
            s = s.replace(gene_presence.index[i],str(gene_presence[gene_presence.index[i]]))
            model_info.loc[j,'GENE ASSOCIATION'] = s

    #rxn_presence is a binarized vector showing which reactions should be there and which not regarding the
    #threshold we take
    rxn_presence = list()
    no_info_genes= list()
    for i in range(model_info.shape[0]):
        s = model_info.iloc[i,9]

        #take care if there is no similarity information for a gene in the model. in this case I just
        #keep them out from all the isolates (replace them with zero)
        while 'AFUA' in s:
            loc = s.find('AFUA')
            no_info_gene = s[loc:loc + 12]
            s = s.replace(no_info_gene,'0')
            no_info_genes.append(no_info_gene)
            
        
        rxn_presence.append(eval(s))

    no_info_genes = list(set(no_info_genes))

    rxn_presence = pd.DataFrame({'present':rxn_presence},index = model_info['ID'])

    return model_info,rxn_presence,no_info_genes

#for minimal cut sets 
def del_rxn_give_rest(del_rxns,del_ids,permanent_deletion):
    del_prior_rxns_new = []
    del_prior_ids_new = []
    del_rxns_new = []
    del_ids_new = []
        
    for i in range(len(del_rxns)):
        if i<permanent_deletion:
            del_rxns[i].bounds = (0,0)
            del_prior_rxns_new.append(del_rxns[i])
            del_prior_ids_new.append(del_ids[i])
        else:
            rest_rxns_new.append(del_rxns[i])
            rest_ids_new.append(del_ids[i])
                
    return del_prior_rxns_new,del_prior_ids_new,bounds,del_rxns_new,del_ids_new

def del_one_rest_fva(del_rxns,next_one):
    bound = del_rxns[next_one].bounds
    del_rxns[next_one].bounds = (0,0)
    single_del_id = del_rxns[next_one].id
    return del_rxns,bound,single_del_id

def random_deletion(del_rxns,del_ids,rand_selections,index):
    del_prior_rxns_new = []
    del_prior_ids_new = []
    del_rxns_new = []
    del_ids_new = []
    
    for i in range(len(del_rxns)):
        if i in rand_selections[index]:
            del_rxns[i].bounds = (0,0)
            del_prior_rxns_new.append(del_rxns[i])
            del_prior_ids_new.append(del_ids[i])
        else:
            del_rxns_new.append(del_rxns[i])
            del_ids_new.append(del_ids[i])
        
    return del_prior_rxns_new,del_prior_ids_new,del_rxns_new,del_ids_new

def max_min_grRule(a):
    sign_locs = []
    before_sign_close_pr = 0
    before_sign_open_pr = 0
    sign = 'none'
    i=len(a)-1
    while i>=0:

        if sign=='none' and a[i-3:i]=='and':# and (i-3,i): not in sign_locs:
            sign = 'and'
            sign_loc = i
        if sign=='none' and a[i-2:i]=='or':# and (i-2,i) not in sign_locs:
            sign = 'or'
            sign_loc = i

        if sign in ['and','or'] and a[i]==')':
            before_sign_close_pr += 1
            pr_occured = True
        if sign in ['and','or'] and a[i]=='(':
            before_sign_open_pr += 1
            pr_occured = True

        if before_sign_open_pr == before_sign_close_pr + 1:
            if a[i-3:i] not in ['max','min']:
                if sign=='and':
                    a = a[:i]+'min'+a[i:]
                    #a = a[:i+3] + a[i+3:].replace('and',',')
                if sign=='or':
                    a = a[:i]+'max'+a[i:]
                    #a = a[:i+3] + a[i+3:].replace('or',',')

            #renew initial values
            i=sign_loc
            sign = 'none'
            before_sign_open_pr = 0
            before_sign_close_pr = 0

        i-=1

    return a


def score_gene_replacement(s,iso_name,score_mat):

    while 'AFUA' in s:
        loc = s.find('AFUA')
        gene = s[loc:loc + 12]
        score = str(score_mat.loc[gene,iso_name])
        s = s.replace(gene,score)
    return s

def set_minimal_media(m):
    for i in m.reactions:
        if i.id[:2]=='EX':
            if i.id in ['EX_C00001[e]','EX_C00007[e]','EX_C00009[e]','EX_C00059[e]','EX_C00031[e]','EX_C01342[e]','EX_CHEBI29033[e]']:
                i.lower_bound = -1000
            else:
                i.lower_bound = 0
    return m

#import identity matrix which was completed in tongta_data_connection folder, Aspergilli model and model_info
#I must keep in mind that every time I'm using the last version of identity matrix, model and model_info
model = cobra.io.read_sbml_model('mohammadmirhakkak/A_fumigatus_GEM/GEMs/Pan_aspergillus_fumigatus.xml')
model_info = pd.read_excel('mohammadmirhakkak/A_fumigatus_GEM/dat/MM_Af_strainGEMs_Supplementary_TableS7.xlsx',sheet_name = 'Reactions',index_col=0)
identity_score_48 = pd.read_csv('mohammadmirhakkak/A_fumigatus_GEM/dat/identity_score_48.csv',index_col=0)
identity_score_252 = pd.read_csv('mohammadmirhakkak/A_fumigatus_GEM/dat/identity_score_252.csv',index_col=0)

identity_score = pd.merge(identity_score_48,identity_score_252,how='outer', left_index=True, right_index=True)
identity_score = identity_score.fillna(0)

#ID unification to have the same IDs as Model IDs
identity_score = unify_gene_id(identity_score)

#unique list of genes in the model
grRules = list(model_info['GENE ASSOCIATION'])
model_genes = list()
for i in range(len(grRules)):
    if not pd.isna(grRules[i]):
        for m in re.finditer('AFUA', grRules[i]):
            model_genes.append(grRules[i][m.start():m.start()+12])

#filter identity score based on the genes present in the model
identity_score = identity_score[identity_score.index.isin(model_genes)]


#binarize identity_score based on the defined threshold (95%)
identity_score_bin = identity_score.copy()
identity_score_bin[identity_score >= 95] = 1
identity_score_bin[identity_score < 95] = 0


#create a vector of grRules with 0 and 1 to replacements of genes
#The aim is to calculate reaction presence and absence
#also make a dataframe as many rows as model_info and as many columns as number of isolates
#this will indicate which reactions should be preliminary present in either of the isolates
#0 -> absent, 1 -> present, nan -> not associated to genes. So will not be picked up only if essential for growth of models
rxn_presence = pd.DataFrame(index = model_info['ID'],columns=identity_score.columns)
for isolate in identity_score.columns:
    grRules_0_1 = grRules.copy()
    model_genes = list()
    for i in range(len(grRules_0_1)):
        if not pd.isna(grRules_0_1[i]):
            while 'AFUA' in grRules_0_1[i]:
                AFUA_loc = grRules_0_1[i].find('AFUA')
                selected_gene = grRules_0_1[i][AFUA_loc:AFUA_loc+12]
                on_off = str(identity_score_bin.loc[selected_gene,isolate])
                grRules_0_1[i] = grRules_0_1[i].replace(selected_gene,on_off)

            rxn_presence.loc[rxn_presence.index[i],isolate] = eval(grRules_0_1[i])

        else:

            rxn_presence.loc[rxn_presence.index[i],isolate] = grRules_0_1[i]



#do FVA to detect the essential reactions. This time on minimal media. so do not relax the bounds
#for i in model.reactions:
#    if i.id[:2]=='EX':
#        i.lower_bound = -1000
model = set_minimal_media(model)
r_bio = model.reactions.get_by_id('R3629')
r_bio.bounds = (0.01,0.01)
fva = flux_variability_analysis(model,model.reactions)
single_essential = fva.query('minimum*maximum>0.0000000001').index
r_bio.bounds = (0,1000)

#active part of network
active_net = set(fva.query('abs(minimum)>0.00001 or abs(maximum)>0.00001').index)

#Keep the essential rxns in for all the models
for i in single_essential:
    ind = list(rxn_presence.index).index(i)
    rxn_presence.iloc[ind,] = 1
    
#ind = list(rxn_presence.index).index('r0154YCM606R00509[m]')
#rxn_presence.iloc[ind,] = 1
#ind = list(rxn_presence.index).index('r0436YCM606R00310[m]')
#rxn_presence.iloc[ind,] = 1
#ind = list(rxn_presence.index).index('r0271YCM606R02333[c]')
#rxn_presence.iloc[ind,] = 1

#set the minimal media with glucose for the base model. I did it above
#model = set_minimal_media(model)

###requirments for the following model generator LOOP (Add all non-genome associated rxns)
non_genome_rxns_ids = model_info[model_info['GENE ASSOCIATION'].isna()].ID.values

non_genome_rxns = []
for i in non_genome_rxns_ids:
    r = model.reactions.get_by_id(i)
    non_genome_rxns.append(r.copy())
###requirments for the following model generator LOOP (Add all non-genome associated rxns)

rxn_presence = rxn_presence.drop(non_genome_rxns_ids)

print('Draft GEM generation started\n')

#LOOP for generating draft models for all the isolates
minimal_cut_sets = []
minimal_cut_sets_concat = []
GRs = []
counter = 0
for isolate in rxn_presence.columns:

    counter+=1

    #draft model containing only exchange and biomass-related reactions
    draft_model = Model(isolate)

    draft_model.add_reactions(non_genome_rxns)

    draft_model.objective = 'R3629'


    #add reactions to each draft isolate model regarding the rxn_presence matrix
    for i in rxn_presence.index:
        if rxn_presence.loc[i,isolate]==1:
            if i in model.reactions:
                r = model.reactions.get_by_id(i)
            else:
                print(i+' is not in the model')
                break
            draft_model.add_reaction(r.copy())

    #set the minimal media with glucose
    draft_model = set_minimal_media(draft_model)
    
    #test if the model can grow. If not, find the minimal cut sets
    sol = draft_model.optimize()
    GRs.append(sol.objective_value)
    if sol.objective_value<0.01:

        model_copied = copy.deepcopy(model)
        r_bio = model_copied.reactions.get_by_id('R3629')
        r_bio.bounds = (0.01,0.01)
        
        del_rxn_ids = rxn_presence.loc[rxn_presence.index.isin(active_net),isolate]
        del_rxn_ids = list(del_rxn_ids[del_rxn_ids==0].index)
            

        del_rxns = []
        bounds = []
        for r_del_id in del_rxn_ids:
            r_del = model_copied.reactions.get_by_id(r_del_id)
            del_rxns.append(r_del)
            bounds.append(r_del.bounds)
        
        '''
        #THIS is to check if there is repeated mcs in the new model. So I put them in minimal_cut_sets beforehand without calculations
        all_minimal_cut_sets = []
        for i in minimal_cut_sets:
            all_minimal_cut_sets = all_minimal_cut_sets+i
        mcs_s = []
        exclude_ids = []
        for i in all_minimal_cut_sets:
            if len(i - set(del_rxn_ids))==0:
                mcs_s.append(i)
                exclude_ids = exclude_ids+list(i)
        minimal_cut_sets_concat.append(mcs_s)
        exclude_ids_locs = []
        for i in del_rxn_ids:
            if i in exclude_ids:
                exclude_ids_locs.append(del_rxn_ids.index(i))
        exclude_ids_locs.reverse()        
        for i in exclude_ids_locs:
            del bounds[i]
            del del_rxn_ids[i]
            del del_rxns[i]
        '''

        still_mcs = True
        MCSs = []
        while still_mcs:
            
            j = len(del_rxns)
            dynamic_index = [0]*len(del_rxns)
            from_beggining = False
            i = 0
            while i < j:
                
                #delete rxns one-by-one
                if dynamic_index[i] == 0:
                    del_rxns[i].bounds = (0,0)
                    dynamic_index[i] = 1

                #bring back the reaction and calculate growth. if still zero ->  reaction is not in MCS    
                else:
                    del_rxns[i].bounds = bounds[i]
                    sol = model_copied.optimize()
                    if sol.objective_value < 0.01 or sol.status=='infeasible':
                        dynamic_index[i] = 0
                    else:
                        del_rxns[i].bounds = (0,0)

                if from_beggining:
                    i += 1
                    continue
                    
                fva = flux_variability_analysis(model_copied,del_rxns)
    
                if fva.query('minimum*maximum>0.0000000001').shape != (0,2):

                    #Maybe more than one reaction identified by fva. so I loop through it
                    for each in fva.query('minimum*maximum>0.0000000001').index:
                    
                        #take the index of fva mandatory rxn
                        ind = del_rxn_ids.index(each)

                        #take that reaction
                        r_del = model_copied.reactions.get_by_id(each)

                        #knock out the reaction
                        r_del.bounds = (0,0)

                        #indicate that the reaction is deleted (1) in the dynamic_index 
                        dynamic_index[ind] = 1
                    
                    #the one obtained via FVA is the last reaction in MCS. and j is one before that. as certainly these two
                    #are in MCS, we keep track of bringing back the reactions until 'j'.
                    j = i

                    #start from beginning
                    i = 0
                    from_beggining = True
                    continue

                i += 1
    
                if i==len(del_rxns):
                    still_mcs = False

            #the last one is not mcs. So should be ignored
            if len(dynamic_index) != sum(dynamic_index):
                mcs = []
                locs = [i for i,x in enumerate(dynamic_index) if x==1]
                locs.reverse()
                for i in locs:
                    mcs.append(del_rxn_ids[i])
                    del_rxns[i].bounds = bounds[i]
                    del del_rxns[i]
                    del del_rxn_ids[i]
                    del bounds[i]

                MCSs.append(mcs)

    else:
        MCSs = [[]]
                
            
    minimal_cut_sets.append(MCSs)











    #go through MCSs and bring back the reactions having highest score in each MCS (one per MCS)
    must_kept = []
    for mcs in MCSs:

        if len(mcs)==0:
            continue
        
        rxn_scores = []
        for rxn in mcs:
                
            gr_rule = model_info.loc[model_info.ID==rxn,'GENE ASSOCIATION'].values[0]
            gr_rule = '(' + gr_rule + ')'
                
            command_string = max_min_grRule(gr_rule)
            command_string = command_string.replace('and',',')
            command_string = command_string.replace('or',',')

            command_string = score_gene_replacement(command_string, isolate, identity_score)

            rxn_scores.append(eval(command_string))

        df_score = pd.DataFrame(index = mcs,data = {'score':rxn_scores})
        df_score.sort_values(by = 'score',ascending=False)
            
        must_kept.append(df_score.index[1])

    
    for i in must_kept:
        
        r_kept = model_copied.reactions.get_by_id(i)
        
        r_kept_to_add = Reaction(r_kept.id)
        draft_model.add_reaction(r_kept_to_add)
        r_kept_to_add.reaction = r_kept.reaction
        

    #when I add missing reactions, I need to make sure that there is no new metabolite added to the model. So compartment should be assigned to save the model. It is so unlikely. But just to make sure
    for i in draft_model.metabolites:
        if i.compartment==None:
            met_added = model.metabolites.get_by_id(i.id)
            i.compartment = met_added.compartment

    draft_model.name = draft_model.id
    
    #print if there is problem with one model
    sol = draft_model.optimize()
    if sol.objective_value<0.01:
        print('problem with ' + draft_model.id)
    
    cobra.io.write_sbml_model(draft_model,'mohammadmirhakkak/A_fumigatus_GEM/GEMs/draft_strain_GEMs/'+draft_model.id+'.xml')
    f'Draft GEM #{counter} was generated and saved.'

pickle.dump(minimal_cut_sets,open('mohammadmirhakkak/A_fumigatus_GEM/res/minimal_cut_sets.pkl','wb'))

stop = timeit.default_timer()

print('Time: ', stop - start)

print('Draft GEM generation finished.\n')  


#FASTCC following
start = timeit.default_timer()

print('FASTCC started.')

model_dir = glob.glob('mohammadmirhakkak/A_fumigatus_GEM/GEMs/draft_strain_GEMs/*.xml')

counter = 0
for dir_m in model_dir:
    counter+=1
    strain = cobra.io.read_sbml_model(dir_m)
    sol = strain.optimize()
    if sol.objective_value<0.01:
        print('original' + strain.name + 'is not viable')
    r_bio = strain.reactions.get_by_id('R3629')
    r_bio.lower_bound = 0.1
    strain_fastcc = fastcc(strain)
    sol = strain_fastcc.optimize()
    if sol.objective_value<0.01:
        print('fastcc-derived' + strain.name + 'is not viable')
    strain_fastcc.name = strain.name
    strain_fastcc.id = strain.name
    r_bio.lower_bound = 0
    cobra.io.write_sbml_model(strain_fastcc,'mohammadmirhakkak/A_fumigatus_GEM/GEMs/strain_GEMs/'+strain.name+'.xml')
    f'Curated GEM #{counter} was generated and saved.'

stop = timeit.default_timer()

print('FASTCC finished.')

print('Time for FASTCC: ', stop - start)
