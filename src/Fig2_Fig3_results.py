import cobra
import pandas as pd
import glob


model_dir = glob.glob('mohammadmirhakkak/A_fumigatus_GEM/GEMs/strain_GEMs/*.xml')
model_dir_252 = []
model_dir_48 = []
for i in model_dir:
    if 'NRZ' in i:
        model_dir_252.append(i)
    else:
        s = i.split('/')[-1][0]
        if s.isdigit():
            model_dir_252.append(i)
        else:
            model_dir_48.append(i)



rxns_300_1_0 = pd.DataFrame()
for i in model_dir:

    model = cobra.io.read_sbml_model(i)
    rxns_in_model = [r.id for r in model.reactions]
    df = pd.DataFrame({model.name:[1]*len(rxns_in_model)},index = rxns_in_model)
    rxns_300_1_0 = pd.merge(rxns_300_1_0,df,right_index = True,left_index = True,how = 'outer')

rxns_300_1_0 = rxns_300_1_0.fillna(0)

rxns_300_1_0.to_csv("mohammadmirhakkak/A_fumigatus_GEM/res/rxns_300.csv")




# take the 252 GEMs of German isolates
cols_300 = rxns_300_1_0.columns
cols_252 = []
for i in cols_300:
    if 'NRZ' in i:
        cols_252.append(i)
    else:
        s = i.split('/')[-1][0]
        if s.isdigit():
            cols_252.append(i)
        else:
            pass
rxns_252_1_0 = rxns_300_1_0.loc[:,rxns_300_1_0.columns.isin(cols_252)]

rxns_252_1_0.to_csv("mohammadmirhakkak/A_fumigatus_GEM/res/rxns_252.csv")


#################
#### Fig 2.a ####
#################
sgem = pd.read_csv('mohammadmirhakkak/A_fumigatus_GEM/dat/metadata_reformat_supplS3.csv')
sgem['number'] = [0] * sgem.shape[0]

for mdir in model_dir_252:

    model = cobra.io.read_sbml_model(mdir)

    mid = model.name.replace('AB01-','')
    subsgem = sgem[sgem['Strain information']==mid]


    # num genes
    number.append(len(model.genes))

    # num total rxns
    number.append(len(model.reactions))
    ind = subsgem[subsgem['GEM information']=='total reactions'].index[0]
    sgem.iloc[ind,5] = len(model.reactions)

    # num unique rxns
    old2new = pd.read_csv("mohammadmirhakkak/A_fumigatus_GEM/dat/old2new.csv",index_col=0)
    #change the reaction IDs to the old ones.
    for i in model.reactions:
        i.id = old2new.loc[old2new.new==i.id,'old'].values[0]
    rxn_ids = [i.id for i in model.reactions]
    unq_rxn_ids = set([i.split('[')[0] for i in rxn_ids])
    number.append(len(unq_rxn_ids))
    ind = subsgem[subsgem['GEM information']=='unique reactions'].index[0]
    sgem.iloc[ind,5] = len(unq_rxn_ids)

    # num gene-associated rxns
    c=0
    for i in model.reactions:
        if len(i.genes)>0:
            c+=1
    number.append(c)
    ind = subsgem[subsgem['GEM information']=='gene-associated reactions'].index[0]
    sgem.iloc[ind,5] = c

    # num orphan rxns
    # relax all the exchange reactions
    for i in model.reactions:
        if i.id[:2]=='EX':
            i.bounds = (-1000,1000)
    orphan_rxns = cobra.flux_analysis.find_blocked_reactions(model)
    number.append(len(orphan_rxns))
    ind = subsgem[subsgem['GEM information']=='orphan reactions'].index[0]
    sgem.iloc[ind,5] = len(orphan_rxns)

    # num exchange rxns
    c = 0
    for i in model.reactions:
        if i.id[:2]=='EX':
            c+=1
    number.append(c)
    ind = subsgem[subsgem['GEM information']=='exchange reactions'].index[0]
    sgem.iloc[ind,5] = c

    # num demand rxns
    c = 0
    for i in model.reactions:
        if i.id[:2]=='DM':
            c+=1
    number.append(c)
    ind = subsgem[subsgem['GEM information']=='demand reactions'].index[0]
    sgem.iloc[ind,5] = c

    # num total mets
    number.append(len(model.metabolites))
    ind = subsgem[subsgem['GEM information']=='total metabolites'].index[0]
    sgem.iloc[ind,5] = len(model.metabolites)

    # num unique metabolites
    met_ids = [i.id for i in model.metabolites]
    unq_met_ids = set([i.split('[')[0] for i in met_ids])
    number.append(len(unq_met_ids))
    ind = subsgem[subsgem['GEM information']=='unique metabolites'].index[0]
    sgem.iloc[ind,5] = len(unq_met_ids)

sgem.to_csv('mohammadmirhakkak/A_fumigatus_GEM/res/metadata_reformat_supplS3_v2.csv')



###############
### Fig 2.b ###
###############

rxns_bin = rxns_252_1_0.copy()

rxns_info = pd.read_excel('mohammadmirhakkak/A_fumigatus_GEM/dat/MM_Af_strainGEMs_Supplementary_TableS7.xlsx',sheet_name = 'Reactions',index_col=0)
rxns_info.loc[rxns_info.SUBSYSTEM=='Demand','SUBSYSTEM'] = 'Other'

rxns_info.SUBSYSTEM = rxns_info.SUBSYSTEM.fillna('Other')

core = rxns_bin[rxns_bin.sum(axis=1)==252].index
path_core = rxns_info.loc[rxns_info.ID.isin(core),'SUBSYSTEM'].values

acc = rxns_bin[rxns_bin.sum(axis=1)<252].index
path_acc = rxns_info.loc[rxns_info.ID.isin(acc),'SUBSYSTEM'].values

# categorize pathways
path_dict = dict()
path_dict['2-Oxocarboxylic acid metabolism'] = 'Carbohydrate m.'
path_dict['Alanine, aspartate, asparagine, glutamate and glutamine metabolism'] = 'Amino acid m.'
path_dict['Amino sugar and nucleotide sugar metabolism'] = 'Carbohydrate m.'
path_dict['Arachidonic acid metabolism'] = 'Lipid m.'
path_dict['Arginine and proline metabolism'] = 'Amino acid m.'
path_dict['Ascorbate and aldarate metabolism'] = 'Carbohydrate m.'
path_dict['Biosynthesis of secondary metabolites'] = 'Secondary metabolites m.'
path_dict['Biotin metabolism'] = 'Cofactors and vitamins m.'
path_dict['Butanoate metabolism'] = 'Carbohydrate m.'
path_dict['C5-Branched dibasic acid metabolism'] = 'Carbohydrate m.'
path_dict['Citrate cycle (TCA cycle)'] = 'Carbohydrate m.'
path_dict['Conversion'] = 'Carbohydrate m.'
path_dict['Cyanoamino acid metabolism'] = 'Amino acid m.'
path_dict['Cysteine and methionine metabolism'] = 'Amino acid m.'
path_dict['Energy metabolism'] = 'Energy m.'
path_dict['Exchange'] = 'Exchange'
path_dict['Fatty acid metabolism'] = 'Lipid m.'
path_dict['Folate biosynthesis'] = 'Cofactors and vitamins m.'
path_dict['Glutathione metabolism'] = 'Amino acid m.'
path_dict['Glycerolipid metabolism'] = 'Lipid m.'
path_dict['Glycerophospholipid metabolism'] = 'Lipid m.'
path_dict['Glycine, serine and threonine metabolism'] = 'Amino acid m.'
path_dict['Glycolysis / Gluconeogenesis'] = 'Carbohydrate m.'
path_dict['Glyoxylate and dicarboxylate metabolism'] = 'Carbohydrate m.'
path_dict['Histidine metabolism'] = 'Amino acid m.'
path_dict['Inositol phosphate metabolism'] = 'Carbohydrate m.'
path_dict['L-Arabinose/Arabitol and D-Xylose/D,L-Xylulose/Xylitol metabolism'] = 'Carbohydrate m.'
path_dict['Lysine metabolism'] = 'Amino acid m.'
path_dict['Mannose/Mannitol, Fructose and Sorbose/Sorbitol metabolism'] = 'Carbohydrate m.'
path_dict['Metabolism of terpenoids and polyketides'] = 'Terpenoids and polyketides m.'
path_dict['Methane metabolism'] = 'Energy m.'
path_dict['Nicotinate and nicotinamide metabolism'] = 'Cofactors and vitamins m.'
path_dict['Nitrogen metabolism'] = 'Energy m.'
path_dict['Nucleotide salvage pathway'] = 'Nucleotide m.'
path_dict['One carbon pool by folate'] = 'Cofactors and vitamins m.'
path_dict['Oxidative phosphorylation'] = 'Energy m.'
path_dict['Pantothenate and CoA biosynthesis'] = 'Cofactors and vitamins m.'
path_dict['Pentose and glucuronate interconversions'] = 'Carbohydrate m.'
path_dict['Pentose phosphate pathway'] = 'Carbohydrate m.'
path_dict['Phenylalanine, tyrosine and tryptophan metabolism'] = 'Amino acid m.'
path_dict['Phenylpropanoid biosynthesis'] = 'Secondary metabolites m.'
path_dict['Phosphonate and phosphinate metabolism'] = 'Amino acid m.'
path_dict['Polysaccharide metabolism (Starch, Cellulose, Chitin, and Xylan)'] = 'Carbohydrate m.'
path_dict['Porphyrin and chlorophyll metabolism'] = 'Cofactors and vitamins m.'
path_dict['Propanoate metabolism'] = 'Carbohydrate m.'
path_dict['Purine metabolism'] = 'Nucleotide m.'
path_dict['Pyrimidine metabolism'] = 'Nucleotide m.'
path_dict['Pyruvate metabolism'] = 'Carbohydrate m.'
path_dict['Riboflavin metabolism'] = 'Cofactors and vitamins m.'
path_dict['Selenocompound metabolism'] = 'Amino acid m.'
path_dict['Sphingolipid metabolism'] = 'Lipid m.'
path_dict['Steroid biosynthesis'] = 'Lipid m.'
path_dict['Sulfur metabolism'] = 'Energy m.'
path_dict['Taurine and hypotaurine metabolism'] = 'Amino acid m.'
path_dict['Thiamine metabolism'] = 'Cofactors and vitamins m.'
path_dict['Transport'] = 'Transport'
path_dict['Ubiquinone and other terpenoid-quinone biosynthesis'] = 'Cofactors and vitamins m.'
path_dict['Valine, leucine and isoleucine metabolism'] = 'Amino acid m.'
path_dict['Vitamin B6 metabolism'] = 'Cofactors and vitamins m.'
path_dict['Xenobiotics biodegradation and metabolism'] = 'Xenobiotics biodegradation m.'
path_dict['alpha-Linolenic acid metabolism'] = 'Lipid m.'
path_dict['beta-Alanine metabolism'] = 'Amino acid m.'
path_dict['Other'] = 'Other'


path_core = [path_dict[i] for i in path_core]
path_acc = [path_dict[i] for i in path_acc]

parent_path = list(set(path_core+path_acc))

path_count_core = pd.DataFrame(pd.DataFrame(path_core).value_counts())
path_count_core.columns = ['Core']
path_count_acc = pd.DataFrame(pd.DataFrame(path_acc).value_counts())
path_count_acc.columns = ['Accessory']

path_count = pd.merge(path_count_core,path_count_acc,left_index = True, right_index = True, how = 'outer')

path_count = path_count.fillna(0)

path_count['sum'] = path_count['Core'] + path_count['Accessory']

path_count = path_count.sort_values(by = 'sum',ascending=False)

path_count.loc['Other','Core'] = path_count.loc['Other','Core'] + path_count.iloc[-3:,0].sum()
path_count.loc['Other','Accessory'] = path_count.loc['Other','Accessory'] + path_count.iloc[-3:,1].sum()
path_count.loc['Other','sum'] = path_count.loc['Other','sum'] + path_count.iloc[-3:,2].sum()

path_count = path_count.iloc[:-3,:]
path_count = path_count.sort_values(by = 'sum',ascending=False)

path_count_v2 = pd.DataFrame({'num':list(path_count.Core)+list(path_count.Accessory),'Reactome':['Core']*path_count.shape[0]+['Accessory']*path_count.shape[0],'Pathway':list(path_count.index)*2})

path_modified = [i[0] for i in path_count_v2.Pathway]
path_count_v2['Pathway'] = path_modified

path_count_v2.to_csv('mohammadmirhakkak/A_fumigatus_GEM/res/path_count_v2.csv')






###############
### Fig 2.c ###
###############

identity_score_bin = rxns_252_1_0.copy()
acc_5 = pd.DataFrame()

boolean = (0<identity_score_bin.sum(axis=1)).values & (identity_score_bin.sum(axis=1)<=0.2*252).values
sub_rxns_bin = identity_score_bin[boolean]
sub_rxns_bin = sub_rxns_bin.sum()
sub_rxns_bin = pd.DataFrame(sub_rxns_bin,columns = ['num'])
sub_rxns_bin['Accessory'] = ['(0, 20%]']*252
sub_rxns_bin['isolates'] = sub_rxns_bin.index.values
acc_5 = pd.concat([acc_5,sub_rxns_bin])

boolean = (0.2*252<identity_score_bin.sum(axis=1)).values & (identity_score_bin.sum(axis=1)<=0.4*252).values
sub_rxns_bin = identity_score_bin[boolean]
sub_rxns_bin = sub_rxns_bin.sum()
sub_rxns_bin = pd.DataFrame(sub_rxns_bin,columns = ['num'])
sub_rxns_bin['Accessory'] = ['(20%, 40%]']*252
sub_rxns_bin['isolates'] = sub_rxns_bin.index.values
acc_5 = pd.concat([acc_5,sub_rxns_bin])

boolean = (0.4*252<identity_score_bin.sum(axis=1)).values & (identity_score_bin.sum(axis=1)<=0.6*252).values
sub_rxns_bin = identity_score_bin[boolean]
sub_rxns_bin = sub_rxns_bin.sum()
sub_rxns_bin = pd.DataFrame(sub_rxns_bin,columns = ['num'])
sub_rxns_bin['Accessory'] = ['(40%, 60%]']*252
sub_rxns_bin['isolates'] = sub_rxns_bin.index.values
acc_5 = pd.concat([acc_5,sub_rxns_bin])

boolean = (0.6*252<identity_score_bin.sum(axis=1)).values & (identity_score_bin.sum(axis=1)<=0.8*252).values
sub_rxns_bin = identity_score_bin[boolean]
sub_rxns_bin = sub_rxns_bin.sum()
sub_rxns_bin = pd.DataFrame(sub_rxns_bin,columns = ['num'])
sub_rxns_bin['Accessory'] = ['(60%, 80%]']*252
sub_rxns_bin['isolates'] = sub_rxns_bin.index.values
acc_5 = pd.concat([acc_5,sub_rxns_bin])

boolean = (0.8*252<identity_score_bin.sum(axis=1)).values & (identity_score_bin.sum(axis=1)<1*252).values
sub_rxns_bin = identity_score_bin[boolean]
sub_rxns_bin = sub_rxns_bin.sum()
sub_rxns_bin = pd.DataFrame(sub_rxns_bin,columns = ['num'])
sub_rxns_bin['Accessory'] = ['(80%, 100%)']*252
sub_rxns_bin['isolates'] = sub_rxns_bin.index.values
acc_5 = pd.concat([acc_5,sub_rxns_bin])

acc_5.index = range(1260)

acc_5.to_csv("mohammadmirhakkak/A_fumigatus_GEM/res/acc_5_rxns.csv")






###############
### Fig 3.a ###
###############


# take only the network reactions

old2new = pd.read_csv("mohammadmirhakkak/A_fumigatus_GEM/dat/old2new.csv",index_col=0)

rxns_bin = rxns_252_1_0.copy()

old_rxns_index = [old2new.loc[old2new.new==i,'old'].values[0] for i in rxns_bin.index]
rxns_bin.index = old_rxns_index
excluded_rxns = []
for i in old_rxns_index:
	if i.startswith('t_') or i.startswith('pot_t_') or i.startswith('t_') or i.startswith('EX_') or i.startswith('DM_') or '_formation' in i or 'Growth' in i:
		excluded_rxns.append(i)

rxns_bin = rxns_bin[~rxns_bin.index.isin(excluded_rxns)]


names = []
group = []
for i in rxns_bin.columns:

    if '-NRZ-' in i:
        group.append('Clinical')
    else:
        group.append('Environmental')

    names.append(i)


jac = pd.DataFrame(columns = names,index = names)
for i in range(rxns_bin.shape[1]):
    for j in range(i,rxns_bin.shape[1]):
        A = set(rxns_bin.index[rxns_bin.loc[:,names[i]]==1])
        B = set(rxns_bin.index[rxns_bin.loc[:,names[j]]==1])
        common = len(A & B)
        union = len(A | B)
        jac.loc[names[i],names[j]] = common/union
        jac.loc[names[j],names[i]] = common/union

df_groups = pd.DataFrame({'origin':group},index = names)

jac.to_csv("mohammadmirhakkak/A_fumigatus_GEM/res/jac_252.csv")
df_groups.to_csv("mohammadmirhakkak/A_fumigatus_GEM/res/group_jac_252.csv")




###############
### Fig 3.b ###
###############

rxns_bin = rxns_252_1_0.copy()

env = list(df_groups.query("origin == 'Environmental'").index)
cli = list(df_groups.query("origin == 'Clinical'").index)

#just to see the percentage of those eight reactions in each group
env_prop = rxns_bin.loc[:,rxns_bin.columns.isin(env)].sum(axis=1)/len(env)
cli_prop = rxns_bin.loc[:,rxns_bin.columns.isin(cli)].sum(axis=1)/len(cli)

diff = env_prop - cli_prop
substantial_diff = diff[abs(diff) >= 0.2]
substantial_diff_ind = list(substantial_diff.index)

#rxns_bin.to_csv("Documents/Aspergillus_fumigatus/works_ncom_revision/8_rxns_300.csv")

# fisher's test
import numpy as np
from scipy.stats import fisher_exact
a = rxns_bin.loc[rxns_bin.index.isin(substantial_diff_ind),rxns_bin.columns.isin(env)]
b = rxns_bin.loc[rxns_bin.index.isin(substantial_diff_ind),rxns_bin.columns.isin(cli)]
#loop for each reaction of 8
pval = []
for i in substantial_diff_ind:
    table = np.array([[a.loc[i,:].sum(), a.shape[1] - a.loc[i,:].sum()], [b.loc[i,:].sum(), b.shape[1] - b.loc[i,:].sum()]])
    oddsr, p = fisher_exact(table, alternative='two-sided')
    pval.append(p)
substantial_diff = pd.DataFrame(substantial_diff)
substantial_diff['pvalue'] = pval

sig_diff = substantial_diff.query("pvalue < 0.05")

rxns_bin_sig = rxns_bin.loc[rxns_bin.index.isin(sig_diff.index),:]

# drop exchange, transport, and demand reactions
a = old2new.loc[old2new.new.isin(list(rxns_bin_sig.index)),]
a.index = a.new
rxns_bin_sig = pd.merge(rxns_bin_sig,a,left_index = True, right_index = True, how = 'outer')
rxns_bin_sig = rxns_bin_sig[~rxns_bin_sig.old.str.startswith('pot_t_')]
rxns_bin_sig = rxns_bin_sig[~rxns_bin_sig.old.str.startswith('t_')]
rxns_bin_sig = rxns_bin_sig[~rxns_bin_sig.old.str.startswith('EX_')]
rxns_bin_sig = rxns_bin_sig[~rxns_bin_sig.old.str.startswith('DM_')]
rxns_bin_sig = rxns_bin_sig.drop('old',axis = 1)
rxns_bin_sig = rxns_bin_sig.drop('new',axis = 1)

rxns_bin_sig.to_csv("mohammadmirhakkak/A_fumigatus_GEM/res/sig_rxns_252.csv")



###############
### Fig 3.c ### catabolic activity prediction. The code for the decision tree is implemented in R.
###############

#import biolog to model metabolite IDs mapping file
df_sources = pd.read_csv('mohammadmirhakkak/A_fumigatus_GEM/dat/df_sources.csv',index_col = 0)


# First of all, catabolic capability of the models should be calculated
#first let's make a dataframe to record the results
catabolic = pd.DataFrame(columns = list(df_sources.Substrates.values),index = range(len(model_dir_252)))  
ind = list(catabolic.index)
catabolic.index = ind

sources = df_sources['Source'].values

#rename the columns and label the with C, N, P, S
col_names = list(catabolic.columns)
for i in range(len(sources)):
    if sources[i]=='carbon':
        col_names[i] = col_names[i] + '_c'
    elif sources[i]=='nitrogen':
        col_names[i] = col_names[i] + '_n'
    elif sources[i]=='phosphorus':
        col_names[i] = col_names[i] + '_p'
    elif sources[i]=='sulfur':
        col_names[i] = col_names[i] + '_s'
catabolic.columns = col_names


ind = []
#import strain models
counter = -1
group = []
for x,i in enumerate(model_dir_252):

    counter+=1
    
    isolate_model = cobra.io.read_sbml_model(i)
    isolate_model.id = isolate_model.name

    #add the row name and category for the related strain model.
    ind = list(catabolic.index)
    ind[counter] = isolate_model.id
    catabolic.index = ind

    if isolate_model.id.startswith('AB01'):
        group.append('clinical')
        st_name = isolate_model.id.replace('AB01-','')
    else:
        group.append('environmental')
        st_name = isolate_model.id
    
    try:
        water = isolate_model.reactions.get_by_id('EX_C00001[e]')
        oxygen = isolate_model.reactions.get_by_id('EX_C00007[e]')
        iron = isolate_model.reactions.get_by_id('EX_CHEBI29033[e]')
        phosphate = isolate_model.reactions.get_by_id('EX_C00009[e]')
        sulfate = isolate_model.reactions.get_by_id('EX_C00059[e]')
        ammonia = isolate_model.reactions.get_by_id('EX_C01342[e]')
        glucose = isolate_model.reactions.get_by_id('EX_C00031[e]')
    except:
        catabolic.iloc[counter,:] = 0
        continue
    

    #first inhibit all the influxes
    for i in isolate_model.reactions:
        if i.id.startswith('EX_'):
            i.lower_bound = 0

    #water and oxygen should be available for all conditions
    water.lower_bound = -1000
    oxygen.lower_bound = -1000
    iron.lower_bound = -1000
    

    for subs in range(df_sources.shape[0]):

        met = df_sources.ID.iloc[subs]

        if df_sources.Source.iloc[subs] == 'carbon':
            #test carbon sources
            ammonia.lower_bound = -10
            phosphate.lower_bound = -10
            sulfate.lower_bound = -10
            #if corresponding exchange reaction is in the model take it, otherwise make it (but first check if the metabolite is available)
            if 'EX_'+met+'[e]' in isolate_model.reactions:
            
                r = isolate_model.reactions.get_by_id('EX_'+met+'[e]')

                r.lower_bound = -10
                sol = isolate_model.optimize()
                obj = sol.objective_value
                r.lower_bound = 0

                if sol.status=='optimal':
                    catabolic.iloc[counter,subs] = obj
                else:
                    catabolic.iloc[counter,subs] = 0
                
                #except KeyError:
                #    r = Reaction('Ex_'+met)
                #    isolate_model.add_reaction(r)
                #    r.reaction = met + ' -->'
                
            else:
        
                catabolic.iloc[counter,subs] = 0
        

        

        
        if df_sources.Source.iloc[subs] == 'nitrogen':
            #test nitrogen sources
            glucose.lower_bound = -10
            ammonia.lower_bound = 0
            phosphate.lower_bound = -10
            sulfate.lower_bound = -10
            #if corresponding exchange reaction is in the model take it, otherwise make it (but first check if the metabolite is available)
            if 'EX_'+met+'[e]' in isolate_model.reactions:
            
                r = isolate_model.reactions.get_by_id('EX_'+met+'[e]')

                r.lower_bound = -10
                sol = isolate_model.optimize()
                obj = sol.objective_value
                r.lower_bound = 0
            
                if sol.status=='optimal':
                    catabolic.iloc[counter,subs] = obj
                else:
                    catabolic.iloc[counter,subs] = 0
                
                #except KeyError:
                #    r = Reaction('Ex_'+met)
                #    isolate_model.add_reaction(r)
                #    r.reaction = met + ' -->'
                
            else:
        
                catabolic.iloc[counter,subs] = 0




        if df_sources.Source.iloc[subs] == 'phosphorus':
            #test phosphorus sources
            glucose.lower_bound = -10
            ammonia.lower_bound = -10
            phosphate.lower_bound = 0
            sulfate.lower_bound = -10
            #if corresponding exchange reaction is in the model take it, otherwise make it (but first check if the metabolite is available)
            if 'EX_'+met+'[e]' in isolate_model.reactions:
            
                r = isolate_model.reactions.get_by_id('EX_'+met+'[e]')

                r.lower_bound = -10
                sol = isolate_model.optimize()
                obj = sol.objective_value
                r.lower_bound = 0
        
        
                if sol.status=='optimal':
                    catabolic.iloc[counter,subs] = obj
                else:
                    catabolic.iloc[counter,subs] = 0
                
                #except KeyError:
                #    r = Reaction('Ex_'+met)
                #    isolate_model.add_reaction(r)
                #    r.reaction = met + ' -->'
                    
            else:
        
                catabolic.iloc[counter,subs] = 0



        if df_sources.Source.iloc[subs] == 'sulfur':
            #test sulfur sources
            glucose.lower_bound = -10
            ammonia.lower_bound = -10
            phosphate.lower_bound = -10
            sulfate.lower_bound = 0
            #if corresponding exchange reaction is in the model take it, otherwise make it (but first check if the metabolite is available)
            if 'EX_'+met+'[e]' in isolate_model.reactions:
            
                r = isolate_model.reactions.get_by_id('EX_'+met+'[e]')

                r.lower_bound = -10
                sol = isolate_model.optimize()
                obj = sol.objective_value
                r.lower_bound = 0
            
                ##put the value in the table
                col = df_sources.loc[df_sources.ID==met,'Substrates'].iloc[0]
                row = isolate_model.id
        
                if sol.status=='optimal':
                    catabolic.iloc[counter,subs] = obj
                else:
                    catabolic.iloc[counter,subs] = 0
                
                #except KeyError:
                #    r = Reaction('Ex_'+met)
                #    isolate_model.add_reaction(r)
                #    r.reaction = met + ' -->'
                
            else:
        
                catabolic.iloc[counter,subs] = 0


#dataframe of groups according to the rows
df_groups = pd.DataFrame({'niche':group})    
df_groups.to_csv('mohammadmirhakkak/A_fumigatus_GEM/res/df_groups_catabolic_activity.csv') 


#save the results. the version must be changed accordingly
catabolic.to_csv("mohammadmirhakkak/A_fumigatus_GEM/res/catabolic_activity.csv")
#since we have duplicates in column, I read the file again to separate the duplicated ones automatically.
#catabolic = pd.read_csv("Documents/Aspergillus_fumigatus/stage2_strain_model/results/catabolic_activity_models_v3.csv",index_col=0)




###############
### Fig 3.d ### flux variability analysis in minimal media. The code for the ML is implemented in R 
###############

#import isolate info
isolate_info = pd.read_csv('mohammadmirhakkak/A_fumigatus_GEM/dat/tree_cluster_metadata_20200826.csv',index_col = 0)


#strain model directions
model_dir = glob.glob('mohammadmirhakkak/A_fumigatus_GEM/GEMs/strain_GEMs/*.xml')
model_dir_252 = []
model_dir_48 = []
for i in model_dir:
    if 'NRZ' in i:
        model_dir_252.append(i)
    else:
        s = i.split('/')[-1][0]
        if s.isdigit():
            model_dir_252.append(i)
        else:
            model_dir_48.append(i)




diet_mm = {'EX_CHEBI29033[e]':-1000,
           'EX_C00001[e]':-1000,
           'EX_C00007[e]':-1000,
           'EX_C00009[e]':-1000,
           'EX_C00059[e]':-1000,
           'EX_C01342[e]':-1000,
           'EX_C00031[e]':-10}


def set_diet(model,diet):

    for i in model.reactions:
        if i.id[:2]=='EX':
            if i.id in diet.keys():
                i.lower_bound = diet[i.id]
            else:
                i.lower_bound = 0

    return model






#################################################
### set up minimal media diet and perform FVA ###
#################################################

min_fva = pd.DataFrame()
max_fva = pd.DataFrame()
group = []
names = []
#import strain models
for i in model_dir_252:

    isolate_model = cobra.io.read_sbml_model(i)

    isolate_model = set_diet(isolate_model,diet_mm)

    if '-NRZ-' in isolate_model.name:
        group.append('Clinical')
    else:
        group.append('Environmental')

    names.append(isolate_model.name)

    fva = flux_variability_analysis(isolate_model,fraction_of_optimum = 0.9, loopless = True, pfba_factor = 1.1)


    min_ = pd.DataFrame(fva['minimum'])
    max_ = pd.DataFrame(fva['maximum'])

    min_.columns = [isolate_model.name]
    max_.columns = [isolate_model.name]
    
    min_fva = pd.merge(min_fva,min_,right_index=True,left_index=True,how='outer')
    max_fva = pd.merge(max_fva,max_,right_index=True,left_index=True,how='outer')


min_fva.to_csv('mohammadmirhakkak/A_fumigatus_GEM/res/min_fva_mm.csv')
max_fva.to_csv('mohammadmirhakkak/A_fumigatus_GEM/res/max_fva_mm.csv')

#dataframe of groups according to the rows
df_groups = pd.DataFrame({'niche':group},index = names)    
df_groups.to_csv('mohammadmirhakkak/A_fumigatus_GEM/res/df_groups_fit_fva_mm.csv')
