import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import seaborn as sns; sns.set_theme(color_codes=True)



###########################################################################################################
###########################################################################################################
############################################ FUNCTIONS ####################################################
###########################################################################################################
###########################################################################################################

# generating F input for the model using output frequencies of phyloWGS
## when we have phyloWGS then we should change the code according to the tree
def phylowgs_F(F_epsilon,useFraction,tree_json,which_tree):
    F_temp = []
    beta_gamma = F_epsilon[0][1]
    l = 200
    if useFraction==True:
        with open(tree_json) as json_file:
            tree_freq = json.load(json_file)
            for population in range(K):
                F_temp.append([l*tree_freq.get('trees').get(str(which_tree)).get('populations').get(str(population)).get('cellular_prevalence')[0],beta_gamma])
        F_temp = np.array(F_temp)

        F_temp[0][0] = F_temp[0][0] - F_temp[1][0]
        F_temp[1][0] = F_temp[1][0] - F_temp[2][0]
        F_temp[2][0] = F_temp[2][0] - F_temp[3][0] - F_temp[8][0]
        F_temp[3][0] = F_temp[3][0] - F_temp[4][0] - F_temp[6][0]
        F_temp[4][0] = F_temp[4][0] - F_temp[5][0]
        F_temp[6][0] = F_temp[6][0] - F_temp[7][0]
    else:
        F_temp = [[2,beta_gamma],[2,beta_gamma],[2,beta_gamma],[2,beta_gamma],[2,beta_gamma],[2,beta_gamma]]
        F_temp = np.array(F_temp)

    F = []
    for ff in range(len(F_temp)):
        if(F_temp[ff][0]/l<0.2):
            F.append([20,1])
        elif (F_temp[ff][0]/l < 0.7):
            F.append([30,1])
        elif (F_temp[ff][0]/l < 1):
            F.append([40,1])
        else:
            print('F has problem')
    F = np.array(F)
    print(F)
    print(F_epsilon)
    return(F)
    # F could be None, In that case, it will be generated using dirichlet distribution

#### this depend on the sections that we have and the arrangement we want for plotting
def adjusting_sections(n_s_data):
    x_max = np.max([np.max(n_s_data[n_s_data.section == 'P1.2'].x)])
    y_max = np.max([np.max(n_s_data[n_s_data.section == 'P1.2'].y), np.max(n_s_data[n_s_data.section == 'P3.3'].y)])

    # n_s_data.loc[n_s_data.section == 'P1.3', 'y'] = n_s_data[n_s_data.section=='P1.3'].y + y_max+1
    n_s_data.loc[n_s_data.section == 'P2.4', 'x'] = n_s_data[n_s_data.section == 'P2.4'].x + x_max + 1
    n_s_data.loc[n_s_data.section == 'P2.4', 'y'] = n_s_data[n_s_data.section == 'P2.4'].y + y_max + 1
    n_s_data.loc[n_s_data.section == 'P3.3', 'x'] = n_s_data[n_s_data.section == 'P3.3'].x + x_max + 1
    return n_s_data

## when we have phyloWGS then we should change the code according to the tree
def generate_C(method, mutations_locations, location, selected_tree,selected_canopy_tree):
    # generating C_ik out of input mutations of phyloWGS and output tree of phyloWGS
    if method == 'phylowgs':
        I = len(mutations_locations)
        with open(location + '.mutass/' + str(selected_tree) + '.json') as json_file:
            tree = json.load(json_file)
            K = len(tree.get('mut_assignments')) + 1
            C_ik = pd.DataFrame(np.zeros(shape=(I, K)))
            for i in range(1, K):
                one = tree.get('mut_assignments').get(str(i)).get('ssms')
                one_mut = [int(element.split('s', 1)[1]) for element in one]
                C_ik[i][one_mut] = 1
        if WES_pos_diff_ST == 0:
            C_ik = C_ik.set_index(mutations_locations)
        else:
            C_ik = C_ik.set_index(mutations_locations_minus1)
        C_ik[2][C_ik[1] == 1] = 1
        C_ik[3][C_ik[2] == 1] = 1
        C_ik[8][C_ik[2] == 1] = 1

        C_ik[4][C_ik[3] == 1] = 1
        C_ik[5][C_ik[4] == 1] = 1

        C_ik[6][C_ik[3] == 1] = 1
        C_ik[7][C_ik[6] == 1] = 1
        print(C_ik)
        print("C_ik")
        print(len(C_ik))
    elif method == 'canopy':
        ## keep in mind that the position here in canopy  is the position in phyloWGS minus 1
        C_ik = pd.read_csv(selected_canopy_tree, index_col=0, sep='\t')
        K = len(C_ik.columns)
    return K, C_ik

###### it would work only for canopy
def generate_F(F_eps,F_file,K):
    F_epsilon = np.tile(F_eps, (K, 1))
    #F = phylowgs_F(F_epsilon,data['Gamma']['F_fraction'],data['C_variation']['phyloWGS']+'.summ.json',data['C_variation']['tree'])
    F_temp = np.array(pd.read_csv(F_file).x)
    print(F_temp)
    F = []
    for ff in range(len(F_temp)):
        #if(F_temp[ff]<0.2):
        #    F.append([20,1])
        #elif (F_temp[ff] < 0.7):
        #    F.append([30,1])
        #elif (F_temp[ff] < 1):
        #    F.append([40,1])
        #else:
        #    print('F has problem')
        for i in range(1,21):
            if (F_temp[ff] < 0.05*i):
                F.append([5*i, 1])
                break
    #F = scale*F_temp
    #F[F<2*F_eps]=2*F_eps
    F = np.array(F)
    return F_epsilon,F
##########################################################################################################

#
def generate_st_wes(frames,wes_file,offset):
    st =  pd.concat(frames)
    ### making wes data structure ( the important thing is that it should have a gene column associiated with location on genome )
    wes = pd.read_csv(wes_file, delimiter="\t")
    wes['gene_minus1'] = wes.gene.str.split(expand=True)[0] + ' '+  (wes.gene.str.split(expand=True)[1].astype(int) - 1).astype(str)

    ### merge st and wes using gene2 column in st and gene column in wes
    gene2 = st['refContig'].astype(str)+' '+(st['refPos']+1).astype(str)
    gene3 = st['refContig'].astype(str) + ' ' + (st['refPos']).astype(str)
    st = pd.concat([st,gene2 ], axis=1)
    st = pd.concat([st, gene3], axis=1)
    st.columns = list(st.columns)[:-2]+['gene2','gene3']
    if(offset==1):
        st_wes = pd.merge(left = st,right = wes,validate = 'many_to_one', how="left",left_on=['gene3'],right_on=['gene']).dropna()
    else:
        st_wes = pd.merge(left=st, right=wes, validate='many_to_one', how="left", left_on=['gene2'],
                          right_on=['gene']).dropna()
    #sections = np.unique(st['source'])
    return st_wes.dropna()

# generating n_s out of the n file from Igor
# add x,y,barcode columns to Igor files
def n_barcode_merge(n_file,sections,barcode):
    n_s_data = pd.read_csv(n_file, delimiter = ",")
    n_s_data.nuclei = n_s_data.nuclei.astype(str).str.split(" ", expand=True, )[0]
    n_s_data = n_s_data[n_s_data.section.isin(sections)]

    coordinate = n_s_data.coordinates.str.split("x", expand=True)
    coordinate = coordinate.rename(columns={0: 'x', 1: 'y'})
    n_s_data = pd.concat([n_s_data, coordinate], axis=1)

    n_s_data['x'] = n_s_data['x'].astype(int)
    n_s_data['y'] = n_s_data['y'].astype(int)

    n_s_data = pd.merge(left=n_s_data, right=barcode, validate='many_to_one', how="left", left_on=['x', 'y'],right_on=['x', 'y'])
    n_s_data['barcode'] = n_s_data.section + '_' + n_s_data['barcode']
    return(n_s_data)


# generating A and D

# Input:
# st_wes:
# offset:
# spots_order:

# output:
# A_df and D_df are dataframes with columns i, s and value.
# A and D are matrices which have been made by pivote function.
def generate_A_D(st_wes,offset,spots_order):
    A_df = []
    D_df = []
    if(offset==1):
        mut_id = 'gene3'
    else:
        mut_id = 'gene2'
    for mutation in st_wes.groupby([mut_id,'spot']).groups:
         #print(mut_id)
        reads = st_wes[np.logical_and(st_wes[mut_id] == mutation[0],st_wes['spot'] == mutation[1])]
        ref = (reads[reads['refAllele'] == reads['base']].cnt.sum())
        alt = (reads[reads['refAllele'] != reads['base']].cnt.sum())
        A_df.append([mutation[0], mutation[1], alt])
        D_df.append([mutation[0], mutation[1], ref + alt])

    A_df = pd.DataFrame(A_df, columns=['i','s','value'])
    D_df = pd.DataFrame(D_df, columns=['i','s','value'])
    A_df = A_df[A_df.s.isin(spots_order)]
    D_df = D_df[D_df.s.isin(spots_order)]

    A_df = A_df.reindex(range(len(A_df['s'])))
    D_df = D_df.reindex(range(len(D_df['s'])))

    A = A_df.pivot(index='i', columns='s', values='value').fillna(0)
    D = D_df.pivot(index='i', columns='s', values='value').fillna(0)
    return A,D,A_df,D_df

def get_proper_spots(section, cells_file):
    cells = pd.read_csv(cells_file,  header=0, names=['section', 'spot', 'n_cells'])
    cells = cells.loc[cells['section']=='P'+section]
    return list(cells['spot'])


##### Analysing the results



def plot_prior_inferred_n(n_lambda, inferred_n, visualization_dir):
    plt.hist(np.array(n_lambda), bins=50)
    plt.savefig(visualization_dir + '/' + 'n_lambda.png', dpi=500)
    plt.close()
    plt.hist(np.array(inferred_n), bins=50)
    plt.savefig(visualization_dir + '/' + '_n_inferred.png', dpi=500)
    plt.close()

    plt.style.use('seaborn-deep')
    n_lists = [list(n_lambda), np.array(inferred_n.transpose()).tolist()]
    bins = np.linspace(0, np.max(n_lists), 30)
    plt.hist(n_lists, label=['Lambda', 'Inferred n'])
    plt.legend(loc='upper right')
    plt.savefig(visualization_dir + '/' +  'lambda_and_inferred.png', dpi=500)
    plt.close()



def plot_spots_each_clone(K, n_s_data, n_lambda, tum_1, section,visualization_dir):
    for clone in range(K):
        plt.scatter(n_s_data['x'][n_lambda.index.tolist()], n_s_data['y'][n_lambda.index.tolist()],
                    c=tum_1.inferred_H[:, clone], cmap='tab20b', vmin=0, vmax=1)
        plt.colorbar()
        plt.savefig(visualization_dir + '/' + section + '_clone_' + str(clone) + '_H.png', dpi=500)
        plt.close()

def plot_pie_chart_H(n_s_data, tum_1, col,section,tree_clones,):
    fontP = FontProperties()
    fontP.set_size('xx-small')
    fig, ax = plt.subplots()
    counter = 0
    for j in range(len(n_s_data['x'])):
        xi = n_s_data['x'][j]
        yi = n_s_data['y'][j]
        values = np.divide(tum_1.inferred_H[counter, :],np.sum(tum_1.inferred_H[counter, :])) #[float(i) / sum(tum_1.inferred_H[counter, :]) for i in tum_1.inferred_H[counter, :]]
        if (counter == 0):
            ax.pie(values, radius=0.5, colors=col, center=(xi, yi))
            ax.set(ylabel='', title='', aspect='equal')
            counter = counter + 1
        else:
            ax.pie(values, radius=0.5, colors=col, center=(xi, yi))
            ax.set(ylabel='', title='', aspect='equal')
            counter = counter + 1
    ax.autoscale()
    ax.legend(loc='upper left', bbox_to_anchor=(1.001, 1), labels=tree_clones, prop=fontP)

    plt.savefig(visualization_dir + '/' + section + '_pie_chart_H.png', dpi=500)
    plt.close()


def plot_pie_chart_HZ(n_s_data, tum_1, col,reads):
    fontP = FontProperties()
    fontP.set_size('xx-small')
    fig, ax = plt.subplots()
    counter = 0
    print(len(n_s_data['x']))
    print(len(tum_1.inferred_H[:,1]))
    for j in range(len(n_s_data['x'])):
        xi = n_s_data['x'][j]
        yi = n_s_data['y'][j]
        values = [float(i) / sum(tum_1.inferred_H[counter, :] * tum_1.inferred_Z[counter, :]) for i in
                  tum_1.inferred_H[counter, :] * tum_1.inferred_Z[counter, :]]
        if np.sum(np.isnan(values)) > 0:
            values = [float(i) / sum(tum_1.inferred_H[counter, :]) for i in
                      tum_1.inferred_H[counter, :]]
        if (counter == 0):
            ax.pie(values, radius=0.5, colors=col, center=(xi, yi))#,wedgeprops={'alpha':reads[j]/np.max(reads[j])})
            ax.set(ylabel='', title='', aspect='equal')
            counter = counter + 1
        else:
            ax.pie(values, radius=0.5, colors=col, center=(xi, yi),wedgeprops={'alpha':reads[j]/np.max(reads[j])})
            ax.set(ylabel='', title='', aspect='equal')
            counter = counter + 1

    ax.autoscale()
    ax.legend(loc='upper left', bbox_to_anchor=(1.001, 1), labels=tree_clones, prop=fontP)

    plt.savefig(visualization_dir + '/' + section + '_pie_chart_HZ.png', dpi=500)
    plt.close()


# merge number of cells from Igor and the barcode file (with x and y)
# file: file of the cell numbers
# sections_n_file:
# barcode: barcoode file which has barcode x y
# n_sampling: True or False, if we are considering n Fix or not.

# output:
# n_s_data: dataframe with barcode, number of cells, x, y
# spots_order: It is actually n_s_data['barcode'] --- maybe redundant
# S: number of spots in the merged dataframe
# n_lambda: again n_s_data.nuclei ---- can be redundant
# n = None or n_lamda based on n_sampling

def generate_n_lambda(file,sections_n_file,barcode,n_sampling):
    if file is not None:
        n_s_data = n_barcode_merge(file,sections_n_file,barcode)
        print("generate_n_lambda - number of Nans/Infs:")
        print(np.sum(n_s_data.isin([np.nan, np.inf, -np.inf])))
        n_s_data = n_s_data[~n_s_data.isin([np.nan, np.inf, -np.inf]).any(1)]
        #n_s_data = np.unique(n_s_data['barcode'])
        print(n_s_data[n_s_data['barcode'].duplicated()])
        spots_order = n_s_data['barcode']
        S = len(spots_order)
        n_lambda = n_s_data.nuclei.astype(float).reindex(n_s_data.index)
        if n_sampling is True:
            n = None
        else:
            n = n_lambda
    else:
        raise('You didn\'t give any prior file for number of cells in the spots')
    return n_s_data,spots_order,S,n_lambda,n

def update_spots_lambda(n_s_data):
    print(n_s_data[n_s_data['barcode'].duplicated()])
    spots_order = n_s_data['barcode'].reset_index(drop=True)
    S = len(spots_order)
    n_lambda = n_s_data['nuclei'].reset_index(drop=True)
    return spots_order,S,n_lambda


###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################


