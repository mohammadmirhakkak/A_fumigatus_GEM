import pandas as pd


model_dir = glob.glob('mohammadmirhakkak/A_fumigatus_GEM/GEMs/strain_GEMs/*.xml')
model_dir_252 = []
for i in model_dir:
    if 'NRZ' in i:
        model_dir_252.append(i)
    else:
        s = i.split('/')[-1][0]
        if s.isdigit():
            model_dir_252.append(i)


genes_252 = pd.DataFrame()
for i in model_dir_252:

    model = cobra.io.read_sbml_model(i)
    genes_in_model = [g.id for g in model.genes]
    df = pd.DataFrame({model.name:[1]*len(genes_in_model)},index = genes_in_model)
    genes_252 = pd.merge(genes_252,df,right_index = True,left_index = True,how = 'outer')

genes_252 = genes_252.fillna(0)


acc_5 = pd.DataFrame()

boolean = (0<genes_252.sum(axis=1)).values & (genes_252.sum(axis=1)<=0.2*252).values
sub_genes_bin = genes_252[boolean]
sub_genes_bin = sub_genes_bin.sum()
sub_genes_bin = pd.DataFrame(sub_genes_bin,columns = ['num'])
sub_genes_bin['Accessory'] = ['(0, 20%]']*252
sub_genes_bin['isolates'] = sub_genes_bin.index.values
acc_5 = pd.concat([acc_5,sub_genes_bin])

boolean = (0.2*252<genes_252.sum(axis=1)).values & (genes_252.sum(axis=1)<=0.4*252).values
sub_genes_bin = genes_252[boolean]
sub_genes_bin = sub_genes_bin.sum()
sub_genes_bin = pd.DataFrame(sub_genes_bin,columns = ['num'])
sub_genes_bin['Accessory'] = ['(20%, 40%]']*252
sub_genes_bin['isolates'] = sub_genes_bin.index.values
acc_5 = pd.concat([acc_5,sub_genes_bin])

boolean = (0.4*252<genes_252.sum(axis=1)).values & (genes_252.sum(axis=1)<=0.6*252).values
sub_genes_bin = genes_252[boolean]
sub_genes_bin = sub_genes_bin.sum()
sub_genes_bin = pd.DataFrame(sub_genes_bin,columns = ['num'])
sub_genes_bin['Accessory'] = ['(40%, 60%]']*252
sub_genes_bin['isolates'] = sub_genes_bin.index.values
acc_5 = pd.concat([acc_5,sub_genes_bin])

boolean = (0.6*252<genes_252.sum(axis=1)).values & (genes_252.sum(axis=1)<=0.8*252).values
sub_genes_bin = genes_252[boolean]
sub_genes_bin = sub_genes_bin.sum()
sub_genes_bin = pd.DataFrame(sub_genes_bin,columns = ['num'])
sub_genes_bin['Accessory'] = ['(60%, 80%]']*252
sub_genes_bin['isolates'] = sub_genes_bin.index.values
acc_5 = pd.concat([acc_5,sub_genes_bin])

boolean = (0.8*252<genes_252.sum(axis=1)).values & (genes_252.sum(axis=1)<1*252).values
sub_genes_bin = genes_252[boolean]
sub_genes_bin = sub_genes_bin.sum()
sub_genes_bin = pd.DataFrame(sub_genes_bin,columns = ['num'])
sub_genes_bin['Accessory'] = ['(80%, 100%)']*252
sub_genes_bin['isolates'] = sub_genes_bin.index.values
acc_5 = pd.concat([acc_5,sub_genes_bin])

acc_5.index = range(1260)

acc_5.to_csv("mohammadmirhakkak/A_fumigatus_GEM/res/acc_5_genes.csv")
