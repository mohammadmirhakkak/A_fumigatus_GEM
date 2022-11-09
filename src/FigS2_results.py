import pandas as pd
import re


def unify_gene_id(matrix):
    new_index = list()
    for i in range(matrix.shape[0]):
        s = matrix.index[i]
        s = 'AFUA_'+s[3:].upper()
        new_index.append(s)
    matrix.index = new_index
    return matrix

identity_score = pd.read_csv('mohammadmirhakkak/A_fumigatus_GEM/dat/identity_score_252.csv',index_col=0)
model_info = pd.read_excel('mohammadmirhakkak/A_fumigatus_GEM/dat/Suppl/MM_Af_strainGEMs_Supplementary_TableS7.xlsx',sheet_name = 'Reactions',index_col=0)

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


identity_score_bin = identity_score.copy()
identity_score_bin[identity_score >= 95] = 1
identity_score_bin[identity_score < 95] = 0

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

acc_5.to_csv("mohammadmirhakkak/A_fumigatus_GEM/res/acc_5_genes.csv")