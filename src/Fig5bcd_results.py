##################
##### Fig 5a #####
##################

import pandas as pd
import cobra
import glob
from cobra.flux_analysis import flux_variability_analysis


#import MAMBO predictions
met_profile = pd.read_csv("mohammadmirhakkak/A_fumigatus_GEM/dat/mambo_media.csv",index_col=0)


#import id_mapping
id_mapping = pd.read_csv("mohammadmirhakkak/A_fumigatus_GEM/dat/id_mapping.csv",index_col=0)


#import sample info
abun_data = pd.read_csv("mohammadmirhakkak/A_fumigatus_GEM/dat/clinical_data.csv",index_col=0)



#convert row names of met_profile to aspergillus exchange IDs
new_ids = []
for i in met_profile.index:
    s = i[3:-3]
    if s in id_mapping.agora.values:
        s = id_mapping.loc[id_mapping.agora==s,'coreco'].iloc[0]
        new_ids.append('EX_'+s+'[e]')
    else:
        new_ids.append(s)
met_profile.index = new_ids




#make dictionary of diet out of met_profile
met_profile = -1 * met_profile
diets = met_profile.to_dict()


#strain model directions
isolates_dir = glob.glob('mohammadmirhakkak/A_fumigatus_GEM/GEMs/strain_GEMs/AB01-*.xml')


#import isolate info
isolate_info = pd.read_csv('mohammadmirhakkak/A_fumigatus_GEM/dat/tree_cluster_metadata_20200826.csv',index_col = 0)


drug = []
group = []
df_gr = pd.DataFrame(index=range(49),columns=range(2))

df_gr.columns = ['before','pos']

for x,isolate_dir in enumerate(isolates_dir):

    
    isolate = cobra.io.read_sbml_model(isolate_dir)
    isolate.id = isolate.name
    #do not use defective models
    if len(isolate.reactions)==0:
        continue

    #counter_isolates+=1
    if isolate.name.startswith('AB01'):
        group.append('clinical')
        st_name = isolate.id.replace('AB01-','')
        drug.append(isolate_info.loc[isolate_info.Sample==st_name,'overall_resistant'].iloc[0])
    else:
        group.append('environmental')
        st_name = isolate.id
        drug.append(isolate_info.loc[isolate_info.Sample==st_name,'overall_resistant'].iloc[0])


    for y,sample in enumerate(met_profile.columns):

        #diet setup
        for r in isolate.reactions:
            if r.id[:2]=='EX':
                if r.id in diets[sample].keys():
                    r.lower_bound = diets[sample][r.id]
                else:
                    r.lower_bound = 0
        sol = isolate.optimize()
        df_gr.iloc[x,y] = sol.objective_value

df_gr.columns = met_profile.columns

df_isolates = pd.DataFrame({'group':group,'drug':drug})

sub_df1 = pd.DataFrame(df_gr[df_isolates['group']=='clinical'].median())
sub_df2 = pd.DataFrame(df_gr[df_isolates['group']=='environmental'].median())

df_gr_groups = pd.concat([sub_df1,sub_df2],axis=1)
df_gr_groups = sub_df1

df_gr_groups = df_gr_groups.reindex(abun_data['Ifd.Nummer'].values)

df_gr_groups_reformat = pd.DataFrame({'growth_rate':list(df_gr.iloc[:,0].values)+list(df_gr.iloc[:,1].values),'group':group*2,'drug':drug*2,'cohort':['before']*49+['pos']*49})

isolate_names = [isolates_dir[i].split('/')[5].replace('.xml','') for i in range(49)]

df_gr_groups_reformat['isolate_name'] = isolate_names*2

df_gr_groups_reformat.to_csv("mohammadmirhakkak/A_fumigatus_GEM/res/fba_mambo.csv")


####################
##### Fig 5c,d #####
####################

#import isolate info
isolate_info = pd.read_csv('mohammadmirhakkak/A_fumigatus_GEM/dat/tree_cluster_metadata_20200826.csv',index_col = 0)


#strain model directions
model_dir = glob.glob('mohammadmirhakkak/A_fumigatus_GEM/GEMs/strain_GEMs/AB01-*.xml')


#import mambo diet
diet_median = pd.read_csv("mohammadmirhakkak/A_fumigatus_GEM/dat/mambo_media.csv",index_col=0)


#import id_mapping
id_mapping = pd.read_csv("mohammadmirhakkak/A_fumigatus_GEM/dat/id_mapping.csv",index_col=0)

agora_ids = ['EX_' + i + '(e)' for i in list(id_mapping.agora)]

diet_median = diet_median[diet_median.index.isin(agora_ids)]

agora = [i[3:-3] for i in list(diet_median.index)]

coreco = []
for i in agora:
    coreco.append(id_mapping.loc[id_mapping.agora==i,'coreco'].iloc[0])
    

coreco_ids = ['EX_' + i + '[e]' for i in coreco]

diet_median.index = coreco_ids


diet_pos = diet_median.pos.to_dict()
diet_before = diet_median.before.to_dict()



def set_diet(model,diet):

    for i in model.reactions:
        if i.id[:2]=='EX':
            if i.id in diet.keys():
                i.lower_bound = -1 * diet[i.id]
            else:
                i.lower_bound = 0

    return model




#########################################
### set up 'pos' diet and perform FVA ###
#########################################

min_fva = pd.DataFrame()
max_fva = pd.DataFrame()
group = []
drug = []
#import strain models
for dir_m in model_dir:
    
    isolate_model = cobra.io.read_sbml_model(dir_m)
    isolate_model.id = isolate_model.name

    isolate_model = set_diet(isolate_model,diet_pos)

    sol = isolate_model.optimize()
    if sol.status=='infeasible':
        continue
    
    if isolate_model.id.startswith('AB01'):
        group.append('clinical')
        st_name = isolate_model.id.replace('AB01-','')
        drug.append(isolate_info.loc[isolate_info.Sample==st_name,'overall_resistant'].iloc[0])
    else:
        group.append('environmental')
        st_name = isolate_model.id
        drug.append(isolate_info.loc[isolate_info.Sample==st_name,'overall_resistant'].iloc[0])

    fva = flux_variability_analysis(isolate_model,fraction_of_optimum = 0.9, loopless = True, pfba_factor = 1.1)

    

    min_ = pd.DataFrame(fva['minimum'])
    max_ = pd.DataFrame(fva['maximum'])

    min_.columns = [isolate_model.name]
    max_.columns = [isolate_model.name]
    
    min_fva = pd.merge(min_fva,min_,right_index=True,left_index=True,how='outer')
    max_fva = pd.merge(max_fva,max_,right_index=True,left_index=True,how='outer')

min_fva.to_csv('mohammadmirhakkak/A_fumigatus_GEM/res/min_fva_pos_mambo.csv')
max_fva.to_csv('mohammadmirhakkak/A_fumigatus_GEM/res/max_fva_pos_mambo.csv')

#dataframe of groups according to the rows
df_groups = pd.DataFrame({'niche':group,'drug':drug})    
df_groups.to_csv('mohammadmirhakkak/A_fumigatus_GEM/res/df_groups_fva_pos_mambo.csv') 





############################################
### set up 'before' diet and perform FVA ###
############################################

min_fva = pd.DataFrame()
max_fva = pd.DataFrame()
group = []
drug = []
#import strain models
for dir_m in model_dir:
    
    isolate_model = cobra.io.read_sbml_model(dir_m)
    isolate_model.id = isolate_model.name

    isolate_model = set_diet(isolate_model,diet_before)

    sol = isolate_model.optimize()
    if sol.status=='infeasible':
        continue
    
    if isolate_model.id.startswith('AB01'):
        group.append('clinical')
        st_name = isolate_model.id.replace('AB01-','')
        drug.append(isolate_info.loc[isolate_info.Sample==st_name,'overall_resistant'].iloc[0])
    else:
        group.append('environmental')
        st_name = isolate_model.id
        drug.append(isolate_info.loc[isolate_info.Sample==st_name,'overall_resistant'].iloc[0])

    fva = flux_variability_analysis(isolate_model,fraction_of_optimum = 0.9, loopless = True, pfba_factor = 1.1)

    

    min_ = pd.DataFrame(fva['minimum'])
    max_ = pd.DataFrame(fva['maximum'])

    min_.columns = [isolate_model.name]
    max_.columns = [isolate_model.name]
    
    min_fva = pd.merge(min_fva,min_,right_index=True,left_index=True,how='outer')
    max_fva = pd.merge(max_fva,max_,right_index=True,left_index=True,how='outer')

min_fva.to_csv('mohammadmirhakkak/A_fumigatus_GEM/res/min_fva_before_mambo.csv')
max_fva.to_csv('mohammadmirhakkak/A_fumigatus_GEM/res/max_fva_before_mambo.csv')

#dataframe of groups according to the rows
df_groups = pd.DataFrame({'niche':group,'drug':drug})    
df_groups.to_csv('mohammadmirhakkak/A_fumigatus_GEM/res/df_groups_fva_before_mambo.csv')