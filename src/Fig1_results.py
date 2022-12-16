import pandas as pd
import cobra
import glob
import numpy as np
import pickle

#################
#### Fig. 1b ####
#################

model = cobra.io.read_sbml_model('mohammadmirhakkak/A_fumigatus_GEM/GEMs/Pan_aspergillus_fumigatus.xml')

number = []
gem_info = ["total genes","total rxns","unique rxns","gene-associated rxns","orphan rxns","exchange rxns","demand rxns","total mets","unique mets"]
information = ["genes"]*1 + ["reactions"]*6 + ["metabolites"]*2

# num genes (total and unique)
number.append(len(model.genes))

# num total rxns
number.append(len(model.reactions))

# num unique rxns
old2new = pd.read_csv("mohammadmirhakkak/A_fumigatus_GEM/dat/old2new.csv",index_col=0)
#change the reaction IDs to the old ones.
for i in model.reactions:
    i.id = old2new.loc[old2new.new==i.id,'old'].values[0]
rxn_ids = [i.id for i in model.reactions]
unq_rxn_ids = set([i.split('[')[0] for i in rxn_ids])
number.append(len(unq_rxn_ids))

# num gene-associated rxns
c=0
for i in model.reactions:
    if len(i.genes)>0:
        c+=1
number.append(c)

# num orphan rxns
# relax all the exchange reactions
for i in model.reactions:
    if i.id[:2]=='EX':
        i.bounds = (-1000,1000)
orphan_rxns = cobra.flux_analysis.find_blocked_reactions(model)
number.append(len(orphan_rxns))

# num exchange rxns
c = 0
for i in model.reactions:
    if i.id[:2]=='EX':
        c+=1
number.append(c)

# num demand rxns
c = 0
for i in model.reactions:
    if i.id[:2]=='DM':
        c+=1
number.append(c)

# num total mets
number.append(len(model.metabolites))

# num unique metabolites
met_ids = [i.id for i in model.metabolites]
unq_met_ids = set([i.split('[')[0] for i in met_ids])
number.append(len(unq_met_ids))


gem_info = pd.DataFrame({'number':number, 'gem_info':gem_info, 'information':information})


gem_info.to_csv("mohammadmirhakkak/A_fumigatus_GEM/res/pan_gem_statistics.csv")


#################
#### Fig. 1d ####
#################

rxns_info = pd.read_excel('mohammadmirhakkak/A_fumigatus_GEM/dat/MM_Af_strainGEMs_Supplementary_TableS7.xlsx',sheet_name = 'Reactions',index_col=0)

rxns_info.SUBSYSTEM = rxns_info.SUBSYSTEM.fillna('Other')

# categorize pathways
path_dict = dict()
path_dict['2-Oxocarboxylic acid metabolism'] = 'Carbohydrate metabolism'
path_dict['Alanine, aspartate, asparagine, glutamate and glutamine metabolism'] = 'Amino acid metabolism'
path_dict['Amino sugar and nucleotide sugar metabolism'] = 'Carbohydrate metabolism'
path_dict['Arachidonic acid metabolism'] = 'Lipid metabolism'
path_dict['Arginine and proline metabolism'] = 'Amino acid metabolism'
path_dict['Ascorbate and aldarate metabolism'] = 'Carbohydrate metabolism'
path_dict['Biosynthesis of secondary metabolites'] = 'Biosynthesis of secondary metabolites'
path_dict['Biotin metabolism'] = 'Metabolism of cofactors and vitamins'
path_dict['Butanoate metabolism'] = 'Carbohydrate metabolism'
path_dict['C5-Branched dibasic acid metabolism'] = 'Carbohydrate metabolism'
path_dict['Citrate cycle (TCA cycle)'] = 'Carbohydrate metabolism'
path_dict['Conversion'] = 'Carbohydrate metabolism'
path_dict['Cyanoamino acid metabolism'] = 'Amino acid metabolism'
path_dict['Cysteine and methionine metabolism'] = 'Amino acid metabolism'
path_dict['Energy metabolism'] = 'Energy metabolism'
path_dict['Exchange'] = 'Exchange'
path_dict['Fatty acid metabolism'] = 'Lipid metabolism'
path_dict['Folate biosynthesis'] = 'Metabolism of cofactors and vitamins'
path_dict['Glutathione metabolism'] = 'Amino acid metabolism'
path_dict['Glycerolipid metabolism'] = 'Lipid metabolism'
path_dict['Glycerophospholipid metabolism'] = 'Lipid metabolism'
path_dict['Glycine, serine and threonine metabolism'] = 'Amino acid metabolism'
path_dict['Glycolysis / Gluconeogenesis'] = 'Carbohydrate metabolism'
path_dict['Glyoxylate and dicarboxylate metabolism'] = 'Carbohydrate metabolism'
path_dict['Histidine metabolism'] = 'Amino acid metabolism'
path_dict['Inositol phosphate metabolism'] = 'Carbohydrate metabolism'
path_dict['L-Arabinose/Arabitol and D-Xylose/D,L-Xylulose/Xylitol metabolism'] = 'Carbohydrate metabolism'
path_dict['Lysine metabolism'] = 'Amino acid metabolism'
path_dict['Mannose/Mannitol, Fructose and Sorbose/Sorbitol metabolism'] = 'Carbohydrate metabolism'
path_dict['Metabolism of terpenoids and polyketides'] = 'Metabolism of terpenoids and polyketides'
path_dict['Methane metabolism'] = 'Energy metabolism'
path_dict['Nicotinate and nicotinamide metabolism'] = 'Metabolism of cofactors and vitamins'
path_dict['Nitrogen metabolism'] = 'Energy metabolism'
path_dict['Nucleotide salvage pathway'] = 'Nucleotide metabolism'
path_dict['One carbon pool by folate'] = 'Metabolism of cofactors and vitamins'
path_dict['Oxidative phosphorylation'] = 'Energy metabolism'
path_dict['Pantothenate and CoA biosynthesis'] = 'Metabolism of cofactors and vitamins'
path_dict['Pentose and glucuronate interconversions'] = 'Carbohydrate metabolism'
path_dict['Pentose phosphate pathway'] = 'Carbohydrate metabolism'
path_dict['Phenylalanine, tyrosine and tryptophan metabolism'] = 'Amino acid metabolism'
path_dict['Phenylpropanoid biosynthesis'] = 'Biosynthesis of secondary metabolites'
path_dict['Phosphonate and phosphinate metabolism'] = 'Amino acid metabolism'
path_dict['Polysaccharide metabolism (Starch, Cellulose, Chitin, and Xylan)'] = 'Carbohydrate metabolism'
path_dict['Porphyrin and chlorophyll metabolism'] = 'Metabolism of cofactors and vitamins'
path_dict['Propanoate metabolism'] = 'Carbohydrate metabolism'
path_dict['Purine metabolism'] = 'Nucleotide metabolism'
path_dict['Pyrimidine metabolism'] = 'Nucleotide metabolism'
path_dict['Pyruvate metabolism'] = 'Carbohydrate metabolism'
path_dict['Riboflavin metabolism'] = 'Metabolism of cofactors and vitamins'
path_dict['Selenocompound metabolism'] = 'Amino acid metabolism'
path_dict['Sphingolipid metabolism'] = 'Lipid metabolism'
path_dict['Steroid biosynthesis'] = 'Lipid metabolism'
path_dict['Sulfur metabolism'] = 'Energy metabolism'
path_dict['Taurine and hypotaurine metabolism'] = 'Amino acid metabolism'
path_dict['Thiamine metabolism'] = 'Metabolism of cofactors and vitamins'
path_dict['Transport'] = 'Transport'
path_dict['Ubiquinone and other terpenoid-quinone biosynthesis'] = 'Metabolism of cofactors and vitamins'
path_dict['Valine, leucine and isoleucine metabolism'] = 'Amino acid metabolism'
path_dict['Vitamin B6 metabolism'] = 'Metabolism of cofactors and vitamins'
path_dict['Xenobiotics biodegradation and metabolism'] = 'Xenobiotics biodegradation and metabolism'
path_dict['alpha-Linolenic acid metabolism'] = 'Lipid metabolism'
path_dict['beta-Alanine metabolism'] = 'Amino acid metabolism'
path_dict['Other'] = 'Other'
path_dict['NAD and NADP Conversion'] = 'Nucleotide metabolism'
path_dict['Galactose/Galactitol metabolism'] = 'Carbohydrate metabolism'
path_dict['Glycan biosynthesis'] = 'Glycan biosynthesis and metabolism'
path_dict['D-Alanine metabolism'] = 'Amino acid metabolism'
path_dict['Alkaloid metabolism'] = 'Biosynthesis of secondary metabolites'
path_dict['Retinol metabolism'] = 'Metabolism of cofactors and vitamins'
path_dict['Hormone metabolism'] = 'Lipid metabolism'
path_dict['D-Arginine and D-ornithine metabolism'] = 'Amino acid metabolism'
path_dict['Flavonoid biosynthesis'] = 'Biosynthesis of secondary metabolites'
path_dict['Ether lipid metabolism'] = 'Lipid metabolism'
path_dict['Primary bile acid biosynthesis'] = 'Lipid metabolism'
path_dict['Aminoacyl-tRNA biosynthesis'] = 'Other'
path_dict['Anthocyanin biosynthesis'] = 'Biosynthesis of secondary metabolites'
path_dict['Linoleic acid metabolism'] = 'Lipid metabolism'
path_dict['Penicillin and cephalosporin biosynthesis'] = 'Biosynthesis of secondary metabolites'
path_dict['Caffeine metabolism'] = 'Biosynthesis of secondary metabolites'
path_dict['Staurosporine biosynthesis'] = 'Biosynthesis of secondary metabolites'
path_dict['Biosynthesis of fumiquinazolines'] = 'Biosynthesis of secondary metabolites'
path_dict['Non protein amino acid biosynthesis'] = 'Amino acid metabolism'
path_dict['Cholesterol metabolism'] = 'Lipid metabolism'
path_dict['Lipoic acid metabolism'] = 'Metabolism of cofactors and vitamins'
path_dict['Hypoglycin biosynthesis'] = 'Biosynthesis of secondary metabolites'
path_dict['Aflatoxin biosynthesis'] = 'Biosynthesis of secondary metabolites'
path_dict['Glucosinolate biosynthesis'] = 'Biosynthesis of secondary metabolites'
path_dict['Biosynthesis of antibiotics'] = 'Metabolism of terpenoids and polyketides'
path_dict['Melanogenesis'] = 'Other'
path_dict['Demand'] = 'Other'

subsystem = rxns_info.SUBSYSTEM
path = [path_dict[i] for i in subsystem]

path_count = pd.DataFrame(pd.DataFrame(path).value_counts())

path_count['pathway'] = path_count.index

path_modified = [i[0] for i in path_count.pathway]
path_count['pathway'] = path_modified

path_count.columns = ['num','pathway']

new_index = list(path_count.pathway)
path_count.index = new_index
del new_index[new_index.index('Other')]
new_index.append('Other')
path_count = path_count.reindex(new_index)

path_count.to_csv('mohammadmirhakkak/A_fumigatus_GEM/res/path_count_pan_model.csv')



#################
#### Fig. 1e ####
#################

#count number of compartments
model = cobra.io.read_sbml_model('mohammadmirhakkak/A_fumigatus_GEM/GEMs/Pan_aspergillus_fumigatus.xml')
old2new = pd.read_csv("mohammadmirhakkak/A_fumigatus_GEM/dat/old2new.csv",index_col=0)

#change the reaction IDs to the old ones.
for i in model.reactions:
    i.id = old2new.loc[old2new.new==i.id,'old'].values[0]

compartments = []
for r in model.reactions:
    if '[' in r.id:
        if r.id.startswith('EX_'):
            compartments.append('Exchange')
        else:
            comp = [met.compartment for met in r.metabolites]
            unq_comp = list(set(comp))
            num_unq_comp = [unq_comp.count(i) for i in unq_comp]
            compartments.append(unq_comp[num_unq_comp.index(max(num_unq_comp))])
    else:
        compartments.append('Transport')
compartments = pd.DataFrame(compartments)
compartments.iloc[compartments[0]=='cytoplasm',0] = 'Cytoplasm'
compartments.iloc[compartments[0]=='mitochondrion',0] = 'Mitochondrion'
compartments.iloc[compartments[0]=='extracellular',0] = 'Extracellular space'
compartments.iloc[compartments[0]=='nucleus',0] = 'Nucleus'
compartments.iloc[compartments[0]=='peroxisome',0] = 'Peroxisome'
compartments.iloc[compartments[0]=='endoplasmic_reticulum',0] = 'Endoplasmic reticulum'
compartments.iloc[compartments[0]=='lipid_particle',0] = 'Lipid particle, Vacuole, and Golgi'
compartments.iloc[compartments[0]=='vacuole',0] = 'Lipid particle, Vacuole, and Golgi'
compartments.iloc[compartments[0]=='golgi',0] = 'Lipid particle, Vacuole, and Golgi'
compartments = pd.DataFrame(compartments[0].value_counts())
compartments.to_csv('mohammadmirhakkak/A_fumigatus_GEM/res/num_compartments.csv')



#################
#### Fig. 1f ####
#################

model_before = cobra.io.read_sbml_model('mohammadmirhakkak/A_fumigatus_GEM/GEMs/draft_aspergillus_fumigatus.xml')
model_after = cobra.io.read_sbml_model('mohammadmirhakkak/A_fumigatus_GEM/GEMs/Pan_aspergillus_fumigatus.xml')

#change the reaction IDs to the old ones.
for i in model_after.reactions:
    i.id = old2new.loc[old2new.new==i.id,'old'].values[0]

df_growth_whole = pd.read_csv('mohammadmirhakkak/A_fumigatus_GEM/dat/biolog_afu.csv',index_col=0)

r = model_after.reactions.get_by_id('maintenance[c]')
r.bounds = (0,1000)

#inhibit influxes
for i in model_after.reactions:
    if i.id[:2]=="EX":
        if i.lower_bound<0:
            i.lower_bound=0
            
iron_ex = model_after.reactions.get_by_id("EX_CHEBI29033[e]")
water_ex = model_after.reactions.get_by_id("EX_C00001[e]")
oxygen_ex = model_after.reactions.get_by_id("EX_C00007[e]")
ammonia_ex = model_after.reactions.get_by_id("EX_C01342[e]")
phosphate_ex = model_after.reactions.get_by_id("EX_C00009[e]")
sulfate_ex = model_after.reactions.get_by_id("EX_C00059[e]")
glucose_ex = model_after.reactions.get_by_id("EX_C00031[e]")

iron_ex.lower_bound = -1000
water_ex.lower_bound = -1000
oxygen_ex.lower_bound = -1000

tested_strain = list(set(df_growth_whole['Strain']))

'''
mutants

niaD: AFUA_1G12830
r527AORYZAER00794[c]: AFUA_5G10420

pyrG: AFUA_2G08360
r0821YCM606R00965[c] : AFUA_2G08360

MET2: AFUA_5G07210
r0549YCM606R01776[c]: AFUA_1G15350 or AFUA_5G07210

LYS4: AFUA_5G08890
r0027YCM606R03444[m]: AFUA_5G08890
r0542YCM606R04371[m]: AFUA_5G08890
'''

#first of all, define the transporter reactions if they are essential
#For instance, if a compound is inside the cell it should be braught to the extracellular space to be tested.
#However, I do it only if the data shows growth. If the data does not show the growth the logic is that there should not be any transporter
#for that and the model correctly does not have that. But if the data shows growth, I add that to the ex space and add the corresponding transport reaction
#Same logic applies for the compounds that are not in the model, but are in the gap-filling reactions.
'''
df_growth = df_growth_whole[df_growth_whole['Strain']=='Af293']

for i in range(df_growth.shape[0]):
    
    met_id = df_growth['ID'].values[i]
    is_out = df_growth['in_gap_filled'].values[i]
    gr_status = df_growth['Growth_status_stat'].values[i]

    if type(met_id)!=float and gr_status==1 and not met_id+'[e]' in model.metabolites and is_out==0:

        #skip if the compound is repeated, because compounds can be present in more than one source, so the algorithm has already made the reactions 
        if 'EX_'+met_id+'[e]' in model.reactions:
            continue

        #make one proper transporter for the model to simulate if the model is able to growth on substrate or not.
        r = Reaction('EX_'+met_id+'[e]')
        model.add_reactions([r])
        r.reaction = met_id + '[e]' + ' -->'
        m = model.metabolites.get_by_id(met_id+'[e]')
        m.compartment = 'extracellular'
        df_growth['in_gap_filled'].values[i] = 0
            
        r = Reaction('t_'+met_id+'_c_e')
        model.add_reactions([r])
        model.add_reaction = met_id + '[e]' + ' <=> ' + met_id + '[c]'
'''
#regarding this logic, I calculate the accuracy by conunting the following conditions like below:

#data:yes, model:yes -> basically, all data-yes coditions should be present in ex space and have proper transporters regardless of if the model can predict growth or not.

#data:yes, model:no -> same as above

#data:no, model:yes -> It can be only the ones with transporter reactions because of secretion of the compound and model showed the growth

#data:no, model:no -> EITHER the ones with transporter reactions and model did not show the growth (for this case, the transporter should be removed (they might have been added from the yeast), unless that is a secretable compound like formate)
                    # OR the ones not present in extracellular space (which is right)  



df_growth_whole['prediction'] = np.nan
df_growth_whole['solution_status'] = np.nan

#filter the ones present in gap-filling area
df_growth_whole = df_growth_whole.query('in_gap_filled!=1')

for st in tested_strain:

    df_growth = df_growth_whole[df_growth_whole['Strain']==st]
    
    #gene KO
    if st=='niaD':
        del_rxns = ['r527AORYZAER00794[c]']
    elif st=='pyrG':
        del_rxns = ['r0821YCM606R00965[c]']
    elif st=='MET2':
        del_rxns = ['r0549YCM606R01776[c]','r0549YCM606R01776[m]']
    elif st=='LYS4':
        del_rxns = ['r0027YCM606R03444[m]','r0542YCM606R04371[m]']
    elif st=='Af293':
        del_rxns = []

    del_bounds = []
    for r_id in del_rxns:
        r = model_after.reactions.get_by_id(r_id)
        del_bounds.append(r.bounds)
        r.bounds = (0,0)



    gr = []
    status = []
    for i in range(df_growth.shape[0]):

        source = df_growth['Source'].values[i]
        plate = df_growth['Plate'].values[i]
        is_out = df_growth['in_gap_filled'].values[i]
        met_id = df_growth['ID'].values[i]
        gr_status = df_growth['Growth_status_stat'].values[i]

        if not pd.isna(met_id) and 'EX_'+met_id+'[e]' in model_after.reactions:

            r = model_after.reactions.get_by_id("EX_"+met_id+'[e]')
            r.lower_bound=-1

            #set up the media with respect to the tested source
            if source=='carbon':

                phosphate_ex.lower_bound=-1
                sulfate_ex.lower_bound=-1
                ammonia_ex.lower_bound=-1

            elif source=='nitrogen':

                phosphate_ex.lower_bound=-1
                sulfate_ex.lower_bound=-1
                glucose_ex.lower_bound=-1

            elif source=='phosphorus':

                glucose_ex.lower_bound=-1
                sulfate_ex.lower_bound=-1
                ammonia_ex.lower_bound=-1

            elif source=='sulfur':

                glucose_ex.lower_bound=-1
                phosphate_ex.lower_bound=-1
                ammonia_ex.lower_bound=-1

            sol = model_after.optimize()
            gr.append(sol.objective_value)
            status.append(sol.status)
            r.lower_bound=0

            #inhibit influxes
            glucose_ex.lower_bound = 0
            phosphate_ex.lower_bound = 0
            ammonia_ex.lower_bound = 0
            sulfate_ex.lower_bound = 0

        else:

            #print warnings
            if st==tested_strain[0] and type(met_id)!=float and 'EX_'+met_id+'[e]' not in model_after.reactions and is_out==1:
                
                print(met_id + ' is out of the compartmentalized model, so was not tested')
                
            elif st==tested_strain[0] and type(met_id)!=float and 'EX_'+met_id+'[e]' not in model_after.reactions and gr_status==1:

                print(met_id + ' is not present in the extracellular space, so was not tested')

            if type(met_id)!=float and gr_status==0 and 'EX_'+met_id+'[e]' not in model_after.reactions:
                gr.append(0)
                status.append('was_not_simulated')
            else:
                gr.append(np.nan)
                status.append(np.nan)

        
    df_growth_whole.loc[df_growth_whole['Strain']==st,'prediction'] = gr
    df_growth_whole.loc[df_growth_whole['Strain']==st,'solution_status'] = status

    #gene KO, bring back
    for ind,r_id in enumerate(del_rxns):
        r = model_after.reactions.get_by_id(r_id)
        r.bounds = del_bounds[ind]


#calculate the accuracy
#first one for data, second one for prediction. For instance 'yn' means experimental yes and prediction no

#Af293
df_growth = df_growth_whole[df_growth_whole['Strain']=='Af293']
df_growth = df_growth[df_growth['Source']=='carbon']
df_growth = df_growth.dropna()
Af293_c_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
Af293_c_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
Af293_c_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
Af293_c_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
Af293_c_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='Af293']
df_growth = df_growth[df_growth['Source']=='nitrogen']
df_growth = df_growth.dropna()
Af293_n_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
Af293_n_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
Af293_n_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
Af293_n_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
Af293_n_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='Af293']
df_growth = df_growth[df_growth['Source']=='phosphorus']
df_growth = df_growth.dropna()
Af293_p_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
Af293_p_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
Af293_p_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
Af293_p_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
Af293_p_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='Af293']
df_growth = df_growth[df_growth['Source']=='sulfur']
df_growth = df_growth.dropna()
Af293_s_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
Af293_s_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
Af293_s_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
Af293_s_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
Af293_s_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

#niaD
df_growth = df_growth_whole[df_growth_whole['Strain']=='niaD']
df_growth = df_growth[df_growth['Source']=='carbon']
df_growth = df_growth.dropna()
niaD_c_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
niaD_c_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
niaD_c_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
niaD_c_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
niaD_c_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='niaD']
df_growth = df_growth[df_growth['Source']=='nitrogen']
df_growth = df_growth.dropna()
niaD_n_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
niaD_n_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
niaD_n_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
niaD_n_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
niaD_n_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='niaD']
df_growth = df_growth[df_growth['Source']=='phosphorus']
df_growth = df_growth.dropna()
niaD_p_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
niaD_p_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
niaD_p_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
niaD_p_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
niaD_p_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='niaD']
df_growth = df_growth[df_growth['Source']=='sulfur']
df_growth = df_growth.dropna()
niaD_s_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
niaD_s_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
niaD_s_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
niaD_s_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
niaD_s_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]


#pyrG
df_growth = df_growth_whole[df_growth_whole['Strain']=='pyrG']
df_growth = df_growth[df_growth['Source']=='carbon']
df_growth = df_growth.dropna()
pyrG_c_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
pyrG_c_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
pyrG_c_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
pyrG_c_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
pyrG_c_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='pyrG']
df_growth = df_growth[df_growth['Source']=='nitrogen']
df_growth = df_growth.dropna()
pyrG_n_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
pyrG_n_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
pyrG_n_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
pyrG_n_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
pyrG_n_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='pyrG']
df_growth = df_growth[df_growth['Source']=='phosphorus']
df_growth = df_growth.dropna()
pyrG_p_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
pyrG_p_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
pyrG_p_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
pyrG_p_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
pyrG_p_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='pyrG']
df_growth = df_growth[df_growth['Source']=='sulfur']
df_growth = df_growth.dropna()
pyrG_s_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
pyrG_s_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
pyrG_s_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
pyrG_s_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
pyrG_s_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]


#MET2
df_growth = df_growth_whole[df_growth_whole['Strain']=='MET2']
df_growth = df_growth[df_growth['Source']=='carbon']
df_growth = df_growth.dropna()
MET2_c_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
MET2_c_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
MET2_c_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
MET2_c_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
MET2_c_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='MET2']
df_growth = df_growth[df_growth['Source']=='nitrogen']
df_growth = df_growth.dropna()
MET2_n_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
MET2_n_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
MET2_n_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
MET2_n_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
MET2_n_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='MET2']
df_growth = df_growth[df_growth['Source']=='phosphorus']
df_growth = df_growth.dropna()
MET2_p_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
MET2_p_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
MET2_p_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
MET2_p_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
MET2_p_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='MET2']
df_growth = df_growth[df_growth['Source']=='sulfur']
df_growth = df_growth.dropna()
MET2_s_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
MET2_s_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
MET2_s_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
MET2_s_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
MET2_s_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]


#LYS4
df_growth = df_growth_whole[df_growth_whole['Strain']=='LYS4']
df_growth = df_growth[df_growth['Source']=='carbon']
df_growth = df_growth.dropna()
LYS4_c_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
LYS4_c_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
LYS4_c_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
LYS4_c_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
LYS4_c_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='LYS4']
df_growth = df_growth[df_growth['Source']=='nitrogen']
df_growth = df_growth.dropna()
LYS4_n_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
LYS4_n_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
LYS4_n_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
LYS4_n_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
LYS4_n_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='LYS4']
df_growth = df_growth[df_growth['Source']=='phosphorus']
df_growth = df_growth.dropna()
LYS4_p_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
LYS4_p_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
LYS4_p_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
LYS4_p_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
LYS4_p_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='LYS4']
df_growth = df_growth[df_growth['Source']=='sulfur']
df_growth = df_growth.dropna()
LYS4_s_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
LYS4_s_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
LYS4_s_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
LYS4_s_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
LYS4_s_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]


gr_pred_after = pd.DataFrame({'Strain':['Af293']*4 + ['niaD']*4 + ['pyrG']*4 + ['MET2']*4 + ['LYS4']*4,
                             'Source':['C','N','P','S']*5,
                             'TP':[Af293_c_acc.iloc[0,0],Af293_n_acc.iloc[0,0],Af293_p_acc.iloc[0,0],Af293_s_acc.iloc[0,0],
                                  niaD_c_acc.iloc[0,0],niaD_n_acc.iloc[0,0],niaD_p_acc.iloc[0,0],niaD_s_acc.iloc[0,0],
                                  pyrG_c_acc.iloc[0,0],pyrG_n_acc.iloc[0,0],pyrG_p_acc.iloc[0,0],pyrG_s_acc.iloc[0,0],
                                  MET2_c_acc.iloc[0,0],MET2_n_acc.iloc[0,0],MET2_p_acc.iloc[0,0],MET2_s_acc.iloc[0,0],
                                  LYS4_c_acc.iloc[0,0],LYS4_n_acc.iloc[0,0],LYS4_p_acc.iloc[0,0],LYS4_s_acc.iloc[0,0]],
                             'FN':[Af293_c_acc.iloc[0,1],Af293_n_acc.iloc[0,1],Af293_p_acc.iloc[0,1],Af293_s_acc.iloc[0,1],
                                  niaD_c_acc.iloc[0,1],niaD_n_acc.iloc[0,1],niaD_p_acc.iloc[0,1],niaD_s_acc.iloc[0,1],
                                  pyrG_c_acc.iloc[0,1],pyrG_n_acc.iloc[0,1],pyrG_p_acc.iloc[0,1],pyrG_s_acc.iloc[0,1],
                                  MET2_c_acc.iloc[0,1],MET2_n_acc.iloc[0,1],MET2_p_acc.iloc[0,1],MET2_s_acc.iloc[0,1],
                                  LYS4_c_acc.iloc[0,1],LYS4_n_acc.iloc[0,1],LYS4_p_acc.iloc[0,1],LYS4_s_acc.iloc[0,1]],
                             'FP':[Af293_c_acc.iloc[1,0],Af293_n_acc.iloc[1,0],Af293_p_acc.iloc[1,0],Af293_s_acc.iloc[1,0],
                                  niaD_c_acc.iloc[1,0],niaD_n_acc.iloc[1,0],niaD_p_acc.iloc[1,0],niaD_s_acc.iloc[1,0],
                                  pyrG_c_acc.iloc[1,0],pyrG_n_acc.iloc[1,0],pyrG_p_acc.iloc[1,0],pyrG_s_acc.iloc[1,0],
                                  MET2_c_acc.iloc[1,0],MET2_n_acc.iloc[1,0],MET2_p_acc.iloc[1,0],MET2_s_acc.iloc[1,0],
                                  LYS4_c_acc.iloc[1,0],LYS4_n_acc.iloc[1,0],LYS4_p_acc.iloc[1,0],LYS4_s_acc.iloc[1,0]],
                             'TN':[Af293_c_acc.iloc[1,1],Af293_n_acc.iloc[1,1],Af293_p_acc.iloc[1,1],Af293_s_acc.iloc[1,1],
                                  niaD_c_acc.iloc[1,1],niaD_n_acc.iloc[1,1],niaD_p_acc.iloc[1,1],niaD_s_acc.iloc[1,1],
                                  pyrG_c_acc.iloc[1,1],pyrG_n_acc.iloc[1,1],pyrG_p_acc.iloc[1,1],pyrG_s_acc.iloc[1,1],
                                  MET2_c_acc.iloc[1,1],MET2_n_acc.iloc[1,1],MET2_p_acc.iloc[1,1],MET2_s_acc.iloc[1,1],
                                  LYS4_c_acc.iloc[1,1],LYS4_n_acc.iloc[1,1],LYS4_p_acc.iloc[1,1],LYS4_s_acc.iloc[1,1]]})



#FOR BEFORE REFINEMENT

#inhibit influxes
for i in model_before.reactions:
    if i.id[:2]=="Ex":
        if i.lower_bound<0:
            i.lower_bound=0
            
iron_ex = model_before.reactions.get_by_id("Ex_CHEBI29033")
water_ex = model_before.reactions.get_by_id("Ex_C00001")
oxygen_ex = model_before.reactions.get_by_id("Ex_C00007")
ammonia_ex = model_before.reactions.get_by_id("Ex_C01342")
phosphate_ex = model_before.reactions.get_by_id("Ex_C00009")
sulfate_ex = model_before.reactions.get_by_id("Ex_C00059")
glucose_ex = model_before.reactions.get_by_id("Ex_C00031")

iron_ex.lower_bound = -1000
water_ex.lower_bound = -1000
oxygen_ex.lower_bound = -1000

tested_strain = list(set(df_growth_whole['Strain']))

'''
mutants

niaD: AFUA_1G12830
r527AORYZAER00794[c]: AFUA_1G12830 or AFUA_5G10420

pyrG: AFUA_2G08360
r0821YCM606R00965[c] : AFUA_2G08360

MET2: AFUA_5G07210
r0549YCM606R01776[c]: AFUA_1G15350 or AFUA_5G07210

LYS4: AFUA_5G08890
r0027YCM606R03444[m]: AFUA_5G08890
r0542YCM606R04371[m]: AFUA_5G08890
'''

#first of all, define the transporter reactions if they are essential
#For instance, if a compound is inside the cell it should be braught to the extracellular space to be tested.
#However, I do it only if the data shows growth. If the data does not show the growth the logic is that there should not be any transporter
#for that and the model correctly does not have that. But if the data shows growth, I add that to the ex space and add the corresponding transport reaction
#Same logic applies for the compounds that are not in the model, but are in the gap-filling reactions.
'''
df_growth = df_growth_whole[df_growth_whole['Strain']=='Af293']

for i in range(df_growth.shape[0]):
    
    met_id = df_growth['ID'].values[i]
    is_out = df_growth['in_gap_filled'].values[i]
    gr_status = df_growth['Growth_status_stat'].values[i]

    if type(met_id)!=float and gr_status==1 and not met_id+'[e]' in model.metabolites and is_out==0:

        #skip if the compound is repeated, because compounds can be present in more than one source, so the algorithm has already made the reactions 
        if 'EX_'+met_id+'[e]' in model.reactions:
            continue

        #make one proper transporter for the model to simulate if the model is able to growth on substrate or not.
        r = Reaction('EX_'+met_id+'[e]')
        model.add_reactions([r])
        r.reaction = met_id + '[e]' + ' -->'
        m = model.metabolites.get_by_id(met_id+'[e]')
        m.compartment = 'extracellular'
        df_growth['in_gap_filled'].values[i] = 0
            
        r = Reaction('t_'+met_id+'_c_e')
        model.add_reactions([r])
        model.add_reaction = met_id + '[e]' + ' <=> ' + met_id + '[c]'
'''
#regarding this logic, I calculate the accuracy by conunting the following conditions like below:

#data:yes, model:yes -> basically, all data-yes coditions should be present in ex space and have proper transporters regardless of if the model can predict growth or not.

#data:yes, model:no -> same as above

#data:no, model:yes -> It can be only the ones with transporter reactions because of secretion of the compound and model showed the growth

#data:no, model:no -> EITHER the ones with transporter reactions and model did not show the growth (for this case, the transporter should be removed (they might have been added from the yeast), unless that is a secretable compound like formate)
                    # OR the ones not present in extracellular space (which is right)  



df_growth_whole['prediction'] = np.nan
df_growth_whole['solution_status'] = np.nan

#filter the ones present in gap-filling area
df_growth_whole = df_growth_whole.query('in_gap_filled!=1')

for st in tested_strain:

    df_growth = df_growth_whole[df_growth_whole['Strain']==st]
    
    #gene KO
    if st=='niaD':
        del_rxns = ['r527AORYZAER00794']
    elif st=='pyrG':
        del_rxns = ['r0821YCM606R00965']
    elif st=='MET2':
        del_rxns = ['r0549YCM606R01776','r0549YCM606R01776']
    elif st=='LYS4':
        del_rxns = ['r0027YCM606R03444','r0542YCM606R04371']
    elif st=='Af293':
        del_rxns = []

    del_bounds = []
    for r_id in del_rxns:
        r = model_before.reactions.get_by_id(r_id)
        del_bounds.append(r.bounds)
        r.bounds = (0,0)



    gr = []
    status = []
    for i in range(df_growth.shape[0]):

        source = df_growth['Source'].values[i]
        plate = df_growth['Plate'].values[i]
        is_out = df_growth['in_gap_filled'].values[i]
        met_id = df_growth['ID'].values[i]
        gr_status = df_growth['Growth_status_stat'].values[i]

        if not pd.isna(met_id) and 'Ex_'+met_id in model_before.reactions:

            r = model_before.reactions.get_by_id("Ex_"+met_id)
            r.lower_bound=-1

            #set up the media with respect to the tested source
            if source=='carbon':

                phosphate_ex.lower_bound=-1
                sulfate_ex.lower_bound=-1
                ammonia_ex.lower_bound=-1

            elif source=='nitrogen':

                phosphate_ex.lower_bound=-1
                sulfate_ex.lower_bound=-1
                glucose_ex.lower_bound=-1

            elif source=='phosphorus':

                glucose_ex.lower_bound=-1
                sulfate_ex.lower_bound=-1
                ammonia_ex.lower_bound=-1

            elif source=='sulfur':

                glucose_ex.lower_bound=-1
                phosphate_ex.lower_bound=-1
                ammonia_ex.lower_bound=-1

            sol = model_before.optimize()
            gr.append(sol.objective_value)
            status.append(sol.status)
            r.lower_bound=0

            #inhibit influxes
            glucose_ex.lower_bound = 0
            phosphate_ex.lower_bound = 0
            ammonia_ex.lower_bound = 0
            sulfate_ex.lower_bound = 0

        else:

            #print warnings
            if st==tested_strain[0] and type(met_id)!=float and 'Ex_'+met_id not in model_before.reactions and is_out==1:
                
                print(met_id + ' is out of the compartmentalized model, so was not tested')
                
            elif st==tested_strain[0] and type(met_id)!=float and 'Ex_'+met_id not in model_before.reactions and gr_status==1:

                print(met_id + ' is not present in the extracellular space, so was not tested')

            if type(met_id)!=float and gr_status==0 and 'Ex_'+met_id not in model_before.reactions:
                gr.append(0)
                status.append('was_not_simulated')
            else:
                gr.append(np.nan)
                status.append(np.nan)

        
    df_growth_whole.loc[df_growth_whole['Strain']==st,'prediction'] = gr
    df_growth_whole.loc[df_growth_whole['Strain']==st,'solution_status'] = status

    #gene KO, bring back
    for ind,r_id in enumerate(del_rxns):
        r = model_before.reactions.get_by_id(r_id)
        r.bounds = del_bounds[ind]


#calculate the accuracy
#first one for data, second one for prediction. For instance 'yn' means experimental yes and prediction no

#Af293
df_growth = df_growth_whole[df_growth_whole['Strain']=='Af293']
df_growth = df_growth[df_growth['Source']=='carbon']
df_growth = df_growth.dropna()
Af293_c_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
Af293_c_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
Af293_c_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
Af293_c_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
Af293_c_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='Af293']
df_growth = df_growth[df_growth['Source']=='nitrogen']
df_growth = df_growth.dropna()
Af293_n_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
Af293_n_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
Af293_n_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
Af293_n_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
Af293_n_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='Af293']
df_growth = df_growth[df_growth['Source']=='phosphorus']
df_growth = df_growth.dropna()
Af293_p_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
Af293_p_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
Af293_p_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
Af293_p_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
Af293_p_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='Af293']
df_growth = df_growth[df_growth['Source']=='sulfur']
df_growth = df_growth.dropna()
Af293_s_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
Af293_s_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
Af293_s_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
Af293_s_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
Af293_s_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

#niaD
df_growth = df_growth_whole[df_growth_whole['Strain']=='niaD']
df_growth = df_growth[df_growth['Source']=='carbon']
df_growth = df_growth.dropna()
niaD_c_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
niaD_c_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
niaD_c_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
niaD_c_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
niaD_c_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='niaD']
df_growth = df_growth[df_growth['Source']=='nitrogen']
df_growth = df_growth.dropna()
niaD_n_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
niaD_n_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
niaD_n_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
niaD_n_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
niaD_n_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='niaD']
df_growth = df_growth[df_growth['Source']=='phosphorus']
df_growth = df_growth.dropna()
niaD_p_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
niaD_p_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
niaD_p_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
niaD_p_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
niaD_p_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='niaD']
df_growth = df_growth[df_growth['Source']=='sulfur']
df_growth = df_growth.dropna()
niaD_s_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
niaD_s_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
niaD_s_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
niaD_s_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
niaD_s_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]


#pyrG
df_growth = df_growth_whole[df_growth_whole['Strain']=='pyrG']
df_growth = df_growth[df_growth['Source']=='carbon']
df_growth = df_growth.dropna()
pyrG_c_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
pyrG_c_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
pyrG_c_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
pyrG_c_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
pyrG_c_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='pyrG']
df_growth = df_growth[df_growth['Source']=='nitrogen']
df_growth = df_growth.dropna()
pyrG_n_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
pyrG_n_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
pyrG_n_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
pyrG_n_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
pyrG_n_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='pyrG']
df_growth = df_growth[df_growth['Source']=='phosphorus']
df_growth = df_growth.dropna()
pyrG_p_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
pyrG_p_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
pyrG_p_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
pyrG_p_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
pyrG_p_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='pyrG']
df_growth = df_growth[df_growth['Source']=='sulfur']
df_growth = df_growth.dropna()
pyrG_s_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
pyrG_s_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
pyrG_s_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
pyrG_s_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
pyrG_s_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]


#MET2
df_growth = df_growth_whole[df_growth_whole['Strain']=='MET2']
df_growth = df_growth[df_growth['Source']=='carbon']
df_growth = df_growth.dropna()
MET2_c_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
MET2_c_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
MET2_c_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
MET2_c_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
MET2_c_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='MET2']
df_growth = df_growth[df_growth['Source']=='nitrogen']
df_growth = df_growth.dropna()
MET2_n_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
MET2_n_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
MET2_n_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
MET2_n_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
MET2_n_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='MET2']
df_growth = df_growth[df_growth['Source']=='phosphorus']
df_growth = df_growth.dropna()
MET2_p_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
MET2_p_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
MET2_p_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
MET2_p_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
MET2_p_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='MET2']
df_growth = df_growth[df_growth['Source']=='sulfur']
df_growth = df_growth.dropna()
MET2_s_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
MET2_s_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
MET2_s_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
MET2_s_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
MET2_s_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]


#LYS4
df_growth = df_growth_whole[df_growth_whole['Strain']=='LYS4']
df_growth = df_growth[df_growth['Source']=='carbon']
df_growth = df_growth.dropna()
LYS4_c_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
LYS4_c_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
LYS4_c_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
LYS4_c_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
LYS4_c_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='LYS4']
df_growth = df_growth[df_growth['Source']=='nitrogen']
df_growth = df_growth.dropna()
LYS4_n_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
LYS4_n_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
LYS4_n_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
LYS4_n_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
LYS4_n_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='LYS4']
df_growth = df_growth[df_growth['Source']=='phosphorus']
df_growth = df_growth.dropna()
LYS4_p_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
LYS4_p_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
LYS4_p_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
LYS4_p_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
LYS4_p_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]

df_growth = df_growth_whole[df_growth_whole['Strain']=='LYS4']
df_growth = df_growth[df_growth['Source']=='sulfur']
df_growth = df_growth.dropna()
LYS4_s_acc = pd.DataFrame(data = {'pred_yes':[0,0],'pred_no':[0,0]},index = ['data_yes','data_no'])
LYS4_s_acc.iloc[0,0] = df_growth.query('Growth_status_stat==1 and prediction>=0.001').shape[0]
LYS4_s_acc.iloc[0,1] = df_growth.query('Growth_status_stat==1 and prediction<0.001').shape[0]
LYS4_s_acc.iloc[1,0] = df_growth.query('Growth_status_stat==0 and prediction>=0.001').shape[0]
LYS4_s_acc.iloc[1,1] = df_growth.query('Growth_status_stat==0 and prediction<0.001').shape[0]


gr_pred_before = pd.DataFrame({'Strain':['Af293']*4 + ['niaD']*4 + ['pyrG']*4 + ['MET2']*4 + ['LYS4']*4,
                             'Source':['C','N','P','S']*5,
                             'TP':[Af293_c_acc.iloc[0,0],Af293_n_acc.iloc[0,0],Af293_p_acc.iloc[0,0],Af293_s_acc.iloc[0,0],
                                  niaD_c_acc.iloc[0,0],niaD_n_acc.iloc[0,0],niaD_p_acc.iloc[0,0],niaD_s_acc.iloc[0,0],
                                  pyrG_c_acc.iloc[0,0],pyrG_n_acc.iloc[0,0],pyrG_p_acc.iloc[0,0],pyrG_s_acc.iloc[0,0],
                                  MET2_c_acc.iloc[0,0],MET2_n_acc.iloc[0,0],MET2_p_acc.iloc[0,0],MET2_s_acc.iloc[0,0],
                                  LYS4_c_acc.iloc[0,0],LYS4_n_acc.iloc[0,0],LYS4_p_acc.iloc[0,0],LYS4_s_acc.iloc[0,0]],
                             'FN':[Af293_c_acc.iloc[0,1],Af293_n_acc.iloc[0,1],Af293_p_acc.iloc[0,1],Af293_s_acc.iloc[0,1],
                                  niaD_c_acc.iloc[0,1],niaD_n_acc.iloc[0,1],niaD_p_acc.iloc[0,1],niaD_s_acc.iloc[0,1],
                                  pyrG_c_acc.iloc[0,1],pyrG_n_acc.iloc[0,1],pyrG_p_acc.iloc[0,1],pyrG_s_acc.iloc[0,1],
                                  MET2_c_acc.iloc[0,1],MET2_n_acc.iloc[0,1],MET2_p_acc.iloc[0,1],MET2_s_acc.iloc[0,1],
                                  LYS4_c_acc.iloc[0,1],LYS4_n_acc.iloc[0,1],LYS4_p_acc.iloc[0,1],LYS4_s_acc.iloc[0,1]],
                             'FP':[Af293_c_acc.iloc[1,0],Af293_n_acc.iloc[1,0],Af293_p_acc.iloc[1,0],Af293_s_acc.iloc[1,0],
                                  niaD_c_acc.iloc[1,0],niaD_n_acc.iloc[1,0],niaD_p_acc.iloc[1,0],niaD_s_acc.iloc[1,0],
                                  pyrG_c_acc.iloc[1,0],pyrG_n_acc.iloc[1,0],pyrG_p_acc.iloc[1,0],pyrG_s_acc.iloc[1,0],
                                  MET2_c_acc.iloc[1,0],MET2_n_acc.iloc[1,0],MET2_p_acc.iloc[1,0],MET2_s_acc.iloc[1,0],
                                  LYS4_c_acc.iloc[1,0],LYS4_n_acc.iloc[1,0],LYS4_p_acc.iloc[1,0],LYS4_s_acc.iloc[1,0]],
                             'TN':[Af293_c_acc.iloc[1,1],Af293_n_acc.iloc[1,1],Af293_p_acc.iloc[1,1],Af293_s_acc.iloc[1,1],
                                  niaD_c_acc.iloc[1,1],niaD_n_acc.iloc[1,1],niaD_p_acc.iloc[1,1],niaD_s_acc.iloc[1,1],
                                  pyrG_c_acc.iloc[1,1],pyrG_n_acc.iloc[1,1],pyrG_p_acc.iloc[1,1],pyrG_s_acc.iloc[1,1],
                                  MET2_c_acc.iloc[1,1],MET2_n_acc.iloc[1,1],MET2_p_acc.iloc[1,1],MET2_s_acc.iloc[1,1],
                                  LYS4_c_acc.iloc[1,1],LYS4_n_acc.iloc[1,1],LYS4_p_acc.iloc[1,1],LYS4_s_acc.iloc[1,1]]})

gr_pred_after.to_csv('mohammadmirhakkak/A_fumigatus_GEM/res/gr_pred_curated.csv')
gr_pred_before.to_csv('mohammadmirhakkak/A_fumigatus_GEM/res/gr_pred_draft.csv')




#################
#### Fig. 1g ####
#################
#Gene essentiality

#after refinement
essential_genes = pd.read_csv("mohammadmirhakkak/A_fumigatus_GEM/dat/essentialAndNon_genes_hu_et_al.csv")
#aaa_essential_genes = pd.read_csv("Documents/Aspergillus_fumigatus/works_ncom_revision/aromatic_essential_genes.csv",index_col=0)
model_after = cobra.io.read_sbml_model('mohammadmirhakkak/A_fumigatus_GEM/GEMs/Pan_aspergillus_fumigatus.xml')

#change the reaction IDs to the old ones.
for i in model_after.reactions:
    i.id = old2new.loc[old2new.new==i.id,'old'].values[0]

#inhibit influxes
for i in model_after.reactions:
    if i.id[:2]=="EX":
        if i.lower_bound<0:
            i.lower_bound=0
            
iron_ex = model_after.reactions.get_by_id("EX_CHEBI29033[e]")
water_ex = model_after.reactions.get_by_id("EX_C00001[e]")
oxygen_ex = model_after.reactions.get_by_id("EX_C00007[e]")
ammonia_ex = model_after.reactions.get_by_id("EX_C01342[e]")
phosphate_ex = model_after.reactions.get_by_id("EX_C00009[e]")
sulfate_ex = model_after.reactions.get_by_id("EX_C00059[e]")
glucose_ex = model_after.reactions.get_by_id("EX_C00031[e]")

iron_ex.lower_bound = -1000
water_ex.lower_bound = -1000
oxygen_ex.lower_bound = -1000
glucose_ex.lower_bound = -1
sulfate_ex.lower_bound = -1000
phosphate_ex.lower_bound = -1000
ammonia_ex.lower_bound = -1000

df_essentiality = pd.DataFrame({'gene':essential_genes['gene'],'data_essential':essential_genes['essential']})
df_essentiality = df_essentiality.drop_duplicates()

unq_genes = list(df_essentiality['gene'])
predicted_essential = list()
objective_value = list()

for i in unq_genes:
    
    sub_df = essential_genes[essential_genes['gene']==i]
    rxns_bounds = list()
    rxns = list()

    #CHS2: no specific Afu gene relates to that. They are all in Or relationship
    if i=='CHS2':
        objective_value.append(model_after.slim_optimize())
        if sol.objective_value < 0.001:
            predicted_essential.append(1)
        else:
            predicted_essential.append(0)
        continue
    
    for j in range(sub_df.shape[0]):
        r = model_after.reactions.get_by_id(sub_df['reaction'].iloc[j])
        rxns_bounds.append(r.bounds)
        rxns.append(r)

        r.bounds = (0,0)

    sol = model_after.optimize()
    objective_value.append(sol.objective_value)

    if sol.objective_value<0.001:
        predicted_essential.append(1)
    else:
        predicted_essential.append(0)

    for i in range(len(rxns)):
        rxns[i].bounds = rxns_bounds[i]



df_essentiality['predicted_essential'] = predicted_essential

"""
df_essentiality_aaa = pd.DataFrame({'gene':aaa_essential_genes['gene'],'data_essential':aaa_essential_genes['essential']})
predicted_essential = []

for i in aaa_essential_genes.afum_gene.values:
    model_copied = model_after.copy()
    model_copied.genes.get_by_id(i).knock_out()
    gr = model_copied.slim_optimize()
    if gr < 0.001:
        predicted_essential.append(1)
    else:
        predicted_essential.append(0)

df_essentiality_aaa['predicted_essential'] = predicted_essential

df_essentiality = pd.concat([df_essentiality,df_essentiality_aaa])
"""

acc_after = pd.DataFrame(data = {'pred_essential':[0,0],'pred_non_essential':[0,0]},index = ['data_essential','data_non_essential'])
acc_after.iloc[0,0] = df_essentiality.query('data_essential==1 and predicted_essential==1').shape[0]
acc_after.iloc[0,1] = df_essentiality.query('data_essential==1 and predicted_essential==0').shape[0]
acc_after.iloc[1,0] = df_essentiality.query('data_essential==0 and predicted_essential==1').shape[0]
acc_after.iloc[1,1] = df_essentiality.query('data_essential==0 and predicted_essential==0').shape[0]

acc_after.to_csv("mohammadmirhakkak/A_fumigatus_GEM/res/Fig1g_gene_essentiality.csv")



#################
#### Fig. 1g ####
#################
import numpy as np
import seaborn as sns
import matplotlib.pylab as plt

model = cobra.io.read_sbml_model('mohammadmirhakkak/A_fumigatus_GEM/GEMs/Pan_aspergillus_fumigatus.xml')

iron_ex = model.reactions.get_by_id("EX_CHEBI29033[e]")
water_ex = model.reactions.get_by_id("EX_C00001[e]")
oxygen_ex = model.reactions.get_by_id("EX_C00007[e]")
ammonia_ex = model.reactions.get_by_id("EX_C01342[e]")
phosphate_ex = model.reactions.get_by_id("EX_C00009[e]")
sulfate_ex = model.reactions.get_by_id("EX_C00059[e]")
glucose_ex = model.reactions.get_by_id("EX_C00031[e]")

iron_ex.lower_bound = -1000
water_ex.lower_bound = -1000
ammonia_ex.lower_bound = -1000
phosphate_ex.lower_bound = -1000
sulfate_ex.lower_bound = -1000


#calibrate the O2 uptake value for normoxic
glucose_ex.lower_bound = - 0.25
o2_upt = 0
o2_uptakes = [o2_upt]
oxygen_ex.lower_bound = o2_upt
gr = [model.slim_optimize()]
for i in range(70):
    o2_upt-=0.01
    o2_uptakes.append(o2_upt)
    oxygen_ex.lower_bound = o2_upt
    gr.append(model.slim_optimize())

df_growth = pd.DataFrame({'o2':abs(np.array(o2_uptakes)),'growth':gr})

sns.lineplot(data = df_growth, x='o2', y='growth')
plt.xlabel('O2 uptake [mmol/grDW/hr]')
plt.ylabel('growth [1/hr]')
plt.grid()
plt.show()

# Normoxic constrains
oxygen_ex.lower_bound = - 14
glucose_ex.lower_bound = - 0.25
norm_exper = 0.013
norm_fba = model.slim_optimize()

# Hypoxic constrains
oxygen_ex.lower_bound = - 0.14
glucose_ex.lower_bound = - 0.206
hypo_exper = 0.011
hypo_fba = model.slim_optimize()

df_growth = pd.DataFrame({"growth_rate":[round(norm_fba,3),norm_exper,round(hypo_fba,3),hypo_exper],
    "condition":['Normoxic','Normoxic','Hypoxic','Hypoxic'],
    "type":['Simulation','Experiment','Simulation','Experiment']})

df_growth.to_csv("mohammadmirhakkak/A_fumigatus_GEM/res/norm_hypo_growth.csv")


# acetate, ethanol, lactate predictions
from numpy import arange
import seaborn as sns

acet_ex = model.reactions.get_by_id('EX_C00033[e]')
ethan_ex = model.reactions.get_by_id('EX_C00469[e]')
slact_ex = model.reactions.get_by_id('EX_C00186[e]')
rlact_ex = model.reactions.get_by_id('EX_C00256[e]')

df = pd.DataFrame()

# Normoxic constrains
oxygen_ex.lower_bound = - 14
glucose_ex.lower_bound = - 0.25
acet_exper_n = 0.004
ethan_exper_n = 0.003
lact_exper_n = 0

for i in arange(0.9,1.01,0.01):
    fva = cobra.flux_analysis.flux_variability_analysis(model,[acet_ex,ethan_ex,slact_ex,rlact_ex],fraction_of_optimum = i)
    sub_df = pd.DataFrame(fva['maximum'])
    sub_df['metabolite'] = ['acetate','ethanol','(S)-lactate','(R)-lactate']
    sub_df['condition'] = ['normoxic']*4
    sub_df['fraction_of_optimum'] = [i]*4
    df = pd.concat([df,sub_df])

# Hypoxic constrains
oxygen_ex.lower_bound = - 0.14
glucose_ex.lower_bound = - 0.206
acet_exper_h = 0.015
ethan_exper_h = 0
lact_exper_h = 0.002
hypo_fba = model.slim_optimize()

for i in arange(0.9,1.01,0.01):
    fva = cobra.flux_analysis.flux_variability_analysis(model,[acet_ex,ethan_ex,slact_ex,rlact_ex],fraction_of_optimum = i)
    sub_df = pd.DataFrame(fva['maximum'])
    sub_df['metabolite'] = ['acetate','ethanol','(S)-lactate','(R)-lactate']
    sub_df['condition'] = ['hypoxic']*4
    sub_df['fraction_of_optimum'] = [i]*4
    df = pd.concat([df,sub_df])

cols = list(df.columns)
cols[0] = 'secretion'
df.columns = cols

df.fraction_of_optimum = round(df.fraction_of_optimum,2)

sns.barplot(data = df.query("condition == 'normoxic'"), x = 'fraction_of_optimum', y = 'secretion', hue = 'metabolite')
plt.xlabel('Fraction of optimum')
plt.ylabel('secretion [mmol/grDW/hr]')
plt.title('byproduct secretions in normoxic condition')
plt.grid()
plt.show()

sns.barplot(data = df.query("condition == 'hypoxic'"), x = 'fraction_of_optimum', y = 'secretion', hue = 'metabolite')
plt.xlabel('Fraction of optimum')
plt.ylabel('secretion [mmol/grDW/hr]')
plt.title('byproduct secretions in hypoxic condition')
plt.grid()
plt.show()