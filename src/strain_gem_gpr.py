import cobra
import pandas as pd
import glob



#function -> gene ID unification (model genes vs predicted genes)
def unify_gene_id(matrix):
    new_index = list()
    for i in range(matrix.shape[0]):
        s = matrix.index[i]
        s = 'AFUA_'+s[3:].upper()
        new_index.append(s)
    matrix.index = new_index
    return matrix


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


def adjust_gpr(s,deleted):

    s = s.replace(deleted,'')

    inconsistent_gpr = True

    while inconsistent_gpr:

        inconsistent_gpr = False

        loc_rparen = ''
        loc_lparen = ''

        if '  and' in s:
            loc = s.find('  and')
            innerP = False
            for j in range(loc,len(s)):# go forward and find ')'
                if s[j]=='(':
                    innerP = True
                if s[j]==')' and innerP:
                    innerP = False
                    continue
                if s[j]==')' and not innerP:
                    loc_rparen = j
                    break
            if loc_rparen == '':
                s = ''
                break

            innerP = False
            for j in range(loc,-1,-1):# go backward and find '('
                if s[j]==')':
                    innerP = True
                if s[j]=='(' and innerP:
                    innerP = False
                    continue
                if s[j]=='(' and not innerP:
                    loc_lparen = j
                    break
            if loc_lparen == '':
                s = ''
                break

            s = s[:loc_lparen] + s[loc_rparen+1:]
            inconsistent_gpr = True
    
        if 'and  ' in s:
            loc = s.find('and  ')
            innerP = False
            for j in range(loc,len(s)):# go forward and find ')'
                if s[j]=='(':
                    innerP = True
                if s[j]==')' and innerP:
                    innerP = False
                    continue
                if s[j]==')' and not innerP:
                    loc_rparen = j
                    break
            if loc_rparen == '':
                s = ''
                break

            innerP = False
            for j in range(loc,-1,-1):# go backward and find '('
                if s[j]==')':
                    innerP = True
                if s[j]=='(' and innerP:
                    innerP = False
                    continue
                if s[j]=='(' and not innerP:
                    loc_lparen = j
                    break
            if loc_lparen == '':
                s = ''
                break

            s = s[:loc_lparen] + s[loc_rparen+1:]
            inconsistent_gpr = True
    
        if s.endswith(' or '):
            s = s[:-4]
            inconsistent_gpr = True

        if s.startswith(' or '):
            s = s[4:]
            inconsistent_gpr = True

        if ' or  or ' in s:
            s = s.replace(' or  or ',' or ')
            inconsistent_gpr = True
        
        if '  or ' in s:
            s = s.replace('  or ',' ')
            inconsistent_gpr = True

        if ' or  ' in s:
            s = s.replace(' or  ',' ')
            inconsistent_gpr = True

        if s.endswith('and '):
            s = ''
            break

        if s.startswith(' and'):
            s = ''
            break

    return s


model_dir = glob.glob('mohammadmirhakkak/A_fumigatus_GEM/GEMs/strain_GEMs/*.xml')


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


#binarize identity_score based on the user-defined threshold (95%)
identity_score_bin = identity_score.copy()
identity_score_bin[identity_score >= 95] = 1
identity_score_bin[identity_score < 95] = 0


for i in model_dir:

    isolate_model = cobra.io.read_sbml_model(i)

    mid = i.split('/')[-1][:-4]

    del_genes = list(identity_score_bin.index[identity_score_bin.loc[:,mid]==0])

    for gene in del_genes:
        if gene in isolate_model.genes:
            genes_still_in.append(gene) # store the genes that need to be deleted

            # Modify the gpr accordingly
            rxns = isolate_model.genes.get_by_id(gene).reactions
            for rxn in rxns:
                gpr = rxn.gene_reaction_rule
                gpr = adjust_gpr(gpr,gene)
                if gpr != '':
                    rxn.gene_reaction_rule = gpr

    # To remove the unused genes
    genes_still_in = list()
    for gene in isolate_model.genes:
        if len(gene.reactions)==0:
            genes_still_in.append(gene) # store the genes that need to be deleted
    cobra.manipulation.delete.remove_genes(isolate_model,genes_still_in) # remove permanently

    # adjust the model ID according to SBML valid ID and save the model
    adj_mid = mid.replace('-','_')
    adj_mid = adj_mid.replace('.','_')
    adj_mid = 'Afu_' + adj_mid
    isolate_model.id = adj_mid

    cobra.io.write_sbml_model(isolate_model,"mohammadmirhakkak/A_fumigatus_GEM/GEMs/strain_GEMs/" + mid + ".xml")








        


