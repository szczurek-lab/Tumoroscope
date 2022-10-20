
from tumoroscope.simulation import simulation as sim
from tumoroscope.visualization import visualization as vis
from tumoroscope.tumoroscope import tumoroscope as tum
from tumoroscope import constants
import pickle
import numpy as np
import random
from multiprocessing import Pool
import time
import json
import re
import pandas as pd
import glob, os
import sys
import matplotlib.pyplot as plt
import seaborn as sns

def run_in_parallel(args):
    return(args.gibbs_sampling(seed = random.randint(1,100) ,min_iter=min_iter,max_iter = max_iter,burn_in=burn_in, batch = batch,simulated_data= sample_1,n_sampling=True,F_fraction=F_fraction,theta_variable = theta_variable, pi_2D=pi_2D,th=th,every_n_sample=every_n_sample,changes_batch=changes_batch,var_calculation=var_calculation))


def run_in_parallel_n(args):
    return(args.gibbs_sampling(seed=random.randint(1, 100),min_iter=min_iter, max_iter=max_iter, burn_in=burn_in, batch=batch,
                                   simulated_data=sample_1, n_sampling=False,F_fraction=F_fraction,theta_variable = theta_variable, pi_2D=pi_2D,th=th,every_n_sample=every_n_sample,changes_batch=changes_batch,var_calculation=var_calculation))


pi_2D = True
th = 0.8 # threshhold for Z
constants.VISUALIZATION = 'visualization'

constants.CHAINS = 1
constants.CORES = 35
every_n_sample = 5
l=1
changes_batch= 500


onLaptop = 'False'
if onLaptop=='False':
    constants.RESULTS = sys.argv[3] #'Resutls_low_reads_corrected_'
    constants.RESULTS = constants.RESULTS + str(sys.argv[2])
    number = str(sys.argv[1])
    config_file =  sys.argv[4]
    noise = str(sys.argv[5])
else:
    constants.RESULTS = 'Resutls_simulation_REALTEMP_'
    constants.RESULTS = constants.RESULTS + '1'
    number = '1'
    config_file =  'configs/configs_500_real_low'
    noise = '1'





result_txt = constants.RESULTS+'/results_'
result_txt_n = constants.RESULTS+'/results_n_'
optimal_rate = 0.40
test = constants.RESULTS
begin= time.time()

if not os.path.exists(constants.RESULTS):
    os.makedirs(constants.RESULTS)
    print("Directory ", constants.RESULTS, " Created ")
else:
    print("Directory ", constants.RESULTS, " already exists")


# intialise
result = {'Name': [],
          'N status': [],
          'H_SEE': [],
          'phi_SEE': [],
          'pi_SEE': [],
          'PZ_SEE': [],
          'n_SEE': [],
          'h_theta_SEE': []
          }
result_df = pd.DataFrame(result)
#os.chdir("configs")
for file in glob.glob(config_file +"/*.json"):
    file_name = re.split("/|.json",file)[2]
    vis_1 = vis(constants.RESULTS + '/n_var_' + file_name + '_' + constants.VISUALIZATION)
    with open(file) as json_data_file:
        data = json.load(json_data_file)

    K = data['structure']['K']
    S = data['structure']['S']
    I = data['structure']['I']
    theta = data['structure']['theta']

    p_c_binom = data['C_variation']['p_c_binom']
    C_temp = data['C_variation']['C']
    repeat_temp = data['C_variation']['repeat']
    if C_temp is None:
        C = None
    else:
        C_t = []
        for ct in range(len(C_temp)):
            C_t.append(np.tile(C_temp[ct], (repeat_temp[ct], 1)))
        C = np.concatenate(C_t)
        #C = np.concatenate((np.tile(C_temp[0], (repeat_temp[0], 1)),np.tile(C_temp[1], (repeat_temp[1], 1)),np.tile(C_temp[2], (repeat_temp[2], 1)),np.tile(C_temp[3], (repeat_temp[3], 1)),np.tile(C_temp[4], (repeat_temp[4], 1))), axis=0)
    vis_1.heatmap_seaborn(C, 'C_seaborn', 'clones', 'mutations', False, 0.5)


    sns.set_theme()
    ax = sns.heatmap(C, annot=False, linewidths=0.5)
    ax.set(xlabel='clones', ylabel='mutations')
    fig = ax.get_figure()
    fig.savefig(constants.RESULTS+'/C_seaborn_new' + '.png')
    plt.close()

    n_sampling = data['n_variation']['n_sampling']

    if data['n_variation']['n'] is not None:
        n = np.array(data['n_variation']['n'])
    else:
        n = None

    Z = data['Z_variation']['Z']
    avarage_clone_in_spot = data['Z_variation']['avarage_clone_in_spot']

    if onLaptop == 'True':
        max_iter = np.int(data['smapling']['max_iter']/10)
        min_iter = np.int(data['smapling']['min_iter']/10)
        burn_in = np.int(data['smapling']['burn_in']/10)
        batch = np.int(data['smapling']['batch']/10)
        var_calculation = np.int(data['smapling']['min_iter'] * 0.09)
    else:
        max_iter = np.int(data['smapling']['max_iter'])
        min_iter = np.int(data['smapling']['min_iter'])
        burn_in = np.int(data['smapling']['burn_in'])
        batch = np.int(data['smapling']['batch'])
        var_calculation = np.int(data['smapling']['min_iter']*0.9)

    phi_gamma = np.array(data['Gamma']['phi_gamma'])
        # F could be None, In that case, it will be generated using dirichlet distribution
    F_epsilon = np.tile(data['Gamma']['F_epsilon'], (K, 1))
    F_fraction =  data['Gamma']['F_fraction']

    F = np.tile(data['Gamma']['F'], (K, 1))  # np.array([[9,2],[9,2],[9,2],[9,2],[9,2],[9,2]])
    phi_gamma_selected = np.array(data['Gamma']['phi_gamma_selected'])
        # F could be None, In that case, it will be generated using dirichlet distribution
    F_epsilon_selected = np.tile(data['Gamma']['F_epsilon_selected'], (K, 1))
    F_selected = np.tile(data['Gamma']['F_selected'], (K, 1))  # np.array([[9,2],[9,2],[9,2],[9,2],[9,2],[9,2]])
    gamma = data['theta']['gamma']
    gamma_sampling = data['theta']['gamma_sampling']
    theta_variable = data['theta']['theta_variable']

    n_lambda = np.tile(data['n_variation']['n_lambda'], (S))
    while(True):
        sample_1 = sim(K=K,S=S,r=phi_gamma[0],p=phi_gamma[1],I=I,F=F,D=None,A=None,C=C,avarage_clone_in_spot=avarage_clone_in_spot,random_seed=random.randint(1,100),F_epsilon= F_epsilon,n=n,p_c_binom=p_c_binom,theta=theta,Z = Z,n_lambda=n_lambda,F_fraction=F_fraction,theta_variable=theta_variable,gamma=gamma, pi_2D=pi_2D)
        if  np.mean(np.sum(sample_1.D, axis=0)) > int(data['Gamma']['mean_read']) - 2 and np.mean(np.sum(sample_1.D, axis=0)) < int(data['Gamma']['mean_read']) + 2:
            break

    if noise=='0':
        n_lambda_tum = sample_1.n
    else:
        b = np.random.binomial(n=1,p=0.5,size=S)
        noise_pois = np.random.poisson(lam=int(noise),size=S)
        n_lambda_tum = sample_1.n+noise_pois*b+noise_pois*(b-1)
        n_lambda_tum[n_lambda_tum<1]=1



    vis_1.visualizing_simulated_data(sample_1)
    pickle.dump(sample_1, open(constants.RESULTS+'/sample_1_'+file_name, 'wb'))

    tum_objs = []
    for cc in range(constants.CHAINS):
        tum_objs.append(tum(name = constants.RESULTS +'/'+ file_name+'_chain_'+str(constants.CHAINS),K=sample_1.K,S=sample_1.S,r=phi_gamma_selected[0],p=phi_gamma_selected[1],I=sample_1.I,avarage_clone_in_spot=sample_1.avarage_clone_in_spot,F = F_selected,C = sample_1.C,A = sample_1.A,D = sample_1.D,F_epsilon=F_epsilon_selected,optimal_rate=optimal_rate,n_lambda = n_lambda_tum,gamma = gamma_sampling, pi_2D=pi_2D,result_txt=result_txt+file_name+'.txt'))

    args_map = tum_objs

    with Pool(constants.CORES) as pool:
        tum_all = pool.map(run_in_parallel, args_map)
        #result_txt = result_txt+str(pool.name)


    pickle.dump(tum_all, open(constants.RESULTS+'/tum_all_'+file_name, 'wb'))

    tum_average = tum(name = constants.RESULTS +'/'+file_name+'avarage', K=sample_1.K,S=sample_1.S,r=sample_1.r,p=sample_1.p,I=sample_1.I,avarage_clone_in_spot=sample_1.avarage_clone_in_spot,F = sample_1.F,C = sample_1.C,A = sample_1.A,D = sample_1.D,F_epsilon=sample_1.F_epsilon,optimal_rate=optimal_rate,n_lambda = n_lambda_tum,gamma=gamma_sampling, pi_2D=pi_2D,result_txt=result_txt+file_name+'.txt')

    tum_average.max_iter = max_iter
    tum_average.burn_in = burn_in
    tum_average.time =  time.time() - begin

    tum_average.inferred_Z = np.sum([c.inferred_Z for c in tum_all],axis=0)/constants.CHAINS
    tum_average.inferred_P_Z = np.sum([c.inferred_P_Z for c in tum_all],axis=0)/constants.CHAINS
    tum_average.inferred_H =np.sum([c.inferred_H for c in tum_all],axis=0)/constants.CHAINS
    tum_average.inferred_phi =np.sum([c.inferred_phi for c in tum_all],axis=0)/constants.CHAINS
    tum_average.inferred_pi =np.sum([c.inferred_pi for c in tum_all],axis=0)/constants.CHAINS
    tum_average.inferred_n =np.sum([c.inferred_n for c in tum_all],axis=0)/constants.CHAINS
    tum_average.inferred_h_theta = np.sum([c.inferred_h_theta for c in tum_all], axis=0) / constants.CHAINS
    tum_average.loglik =np.sum([c.loglik for c in tum_all],axis=0)/constants.CHAINS
    # TODO: these two should be updated in a way that the mean calculated after sampling not here and we do not need max iter here
    tum_average.decision_matrix_G =np.sum([c.decision_matrix_G for c in tum_all],axis=0)/constants.CHAINS
    tum_average.decision_matrix_phi =np.sum([c.decision_matrix_phi for c in tum_all],axis=0)/constants.CHAINS
    tum_average.sigma_phi =np.sum([c.sigma_phi for c in tum_all],axis=0)/constants.CHAINS
    tum_average.sigma_G =np.sum([c.sigma_G for c in tum_all],axis=0)/constants.CHAINS
    tum_average.H_SEE = np.sum([c.H_SEE for c in tum_all],axis=0)/constants.CHAINS
    tum_average.phi_SEE = np.sum([c.phi_SEE for c in tum_all],axis=0)/constants.CHAINS
    tum_average.pi_SEE = np.sum([c.pi_SEE for c in tum_all],axis=0)/constants.CHAINS
    tum_average.PZ_SEE = np.sum([c.PZ_SEE for c in tum_all],axis=0)/constants.CHAINS
    tum_average.n_SEE = np.sum([c.n_SEE for c in tum_all],axis=0)/constants.CHAINS
    tum_average.h_theta_SEE = np.sum([c.h_theta_SEE for c in tum_all], axis=0) / constants.CHAINS
    tum_average.theta_variable = theta_variable
    tum_average.last_convergence = np.sum([c.last_convergence for c in tum_all],axis=0)/constants.CHAINS

    average_dict = {'Name': file_name,
          'N status': 'Variable',
          'H_SEE': tum_average.H_SEE,
          'phi_SEE': tum_average.phi_SEE,
          'pi_SEE': tum_average.pi_SEE,
          'PZ_SEE': tum_average.PZ_SEE,
          'n_SEE': tum_average.n_SEE,
          'h_theta_SEE': tum_average.h_theta_SEE
          }
    result_df = result_df.append(average_dict,ignore_index=True)

    #vis_1.visualizing_inferred_variables(tum_average)
    pickle.dump(tum_average, open(constants.RESULTS+'/tum_average', 'wb'))
    #tum_average.save_in_txt(sample_1,constants.RESULTS+'/average_'+file_name+'.txt')

    #for c in range(constants.CHAINS):
    #    vis_1.vector_plot(np.arange(0, 1, 0.05),tum_all[c].Z_SEE,'Z_SEE'+str(c))
    #vis_1.vector_plot(np.arange(0, 1, 0.05),tum_average.Z_SEE,'Z_SEE_all')
    #vis_1.error_boxplot(tum_all)


    vis_1.likelihood_all(tum_all,'likelihood_all')

    ###### n fix

    vis_2 = vis(constants.RESULTS+'/n_fix_'+file_name+'_'+constants.VISUALIZATION)


    tum_objs_n = []
    for cc in range(constants.CHAINS):
        tum_objs_n.append(tum(name = constants.RESULTS +'/'+file_name+'_n_chain_'+str(constants.CHAINS),K=sample_1.K,S=sample_1.S,r=sample_1.r,p=sample_1.p,I=sample_1.I,avarage_clone_in_spot=sample_1.avarage_clone_in_spot,F = sample_1.F,C = sample_1.C,A = sample_1.A,D = sample_1.D,F_epsilon=sample_1.F_epsilon,optimal_rate=optimal_rate,n_lambda = n_lambda_tum,gamma=gamma_sampling, pi_2D=pi_2D, result_txt= result_txt+file_name+'fixed_n.txt'))

    args_map = tum_objs_n

    with Pool(constants.CORES) as pool:
        tum_all_n = pool.map(run_in_parallel_n, args_map)
        #result_txt = result_txt+str(pool.name)


    pickle.dump(tum_all_n, open(constants.RESULTS+'/tum_all_n_'+file_name, 'wb'))

    tum_average_n = tum(name= constants.RESULTS +'/'+file_name+'_n_average',K=sample_1.K,S=sample_1.S,r=sample_1.r,p=sample_1.p,I=sample_1.I,avarage_clone_in_spot=sample_1.avarage_clone_in_spot,F = sample_1.F,C = sample_1.C,A = sample_1.A,D = sample_1.D,F_epsilon=sample_1.F_epsilon,optimal_rate=optimal_rate,n_lambda = n_lambda_tum,gamma=gamma_sampling, pi_2D=pi_2D,result_txt=result_txt+file_name+'fixed_n.txt')

    tum_average_n.max_iter = max_iter
    tum_average_n.burn_in = burn_in
    tum_average_n.time =  time.time() - begin

    tum_average_n.inferred_Z = np.sum([c.inferred_Z for c in tum_all_n],axis=0)/constants.CHAINS
    tum_average_n.inferred_P_Z = np.sum([c.inferred_P_Z for c in tum_all_n],axis=0)/constants.CHAINS
    tum_average_n.inferred_H =np.sum([c.inferred_H for c in tum_all_n],axis=0)/constants.CHAINS
    tum_average_n.inferred_phi =np.sum([c.inferred_phi for c in tum_all_n],axis=0)/constants.CHAINS
    tum_average_n.inferred_pi =np.sum([c.inferred_pi for c in tum_all_n],axis=0)/constants.CHAINS
    tum_average_n.inferred_n =np.sum([c.inferred_n for c in tum_all_n],axis=0)/constants.CHAINS
    tum_average_n.inferred_h_theta = np.sum([c.inferred_h_theta for c in tum_all_n], axis=0) / constants.CHAINS
    tum_average_n.loglik =np.sum([c.loglik for c in tum_all_n],axis=0)/constants.CHAINS
    # TODO: these two should be updated in a way that the mean calculated after sampling not here and we do not need max iter here
    tum_average_n.decision_matrix_G =np.sum([c.decision_matrix_G for c in tum_all_n],axis=0)/constants.CHAINS
    tum_average_n.decision_matrix_phi =np.sum([c.decision_matrix_phi for c in tum_all_n],axis=0)/constants.CHAINS
    tum_average_n.sigma_phi =np.sum([c.sigma_phi for c in tum_all_n],axis=0)/constants.CHAINS
    tum_average_n.sigma_G =np.sum([c.sigma_G for c in tum_all_n],axis=0)/constants.CHAINS
    tum_average_n.H_SEE = np.sum([c.H_SEE for c in tum_all_n],axis=0)/constants.CHAINS
    tum_average_n.phi_SEE = np.sum([c.phi_SEE for c in tum_all_n],axis=0)/constants.CHAINS
    tum_average_n.pi_SEE = np.sum([c.pi_SEE for c in tum_all_n],axis=0)/constants.CHAINS
    tum_average_n.PZ_SEE = np.sum([c.PZ_SEE for c in tum_all_n],axis=0)/constants.CHAINS
    tum_average_n.n_SEE = np.sum([c.n_SEE for c in tum_all_n],axis=0)/constants.CHAINS
    tum_average_n.h_theta_SEE = np.sum([c.h_theta_SEE for c in tum_all_n], axis=0) / constants.CHAINS
    tum_average_n.theta_variable = theta_variable
    tum_average_n.last_convergence = np.sum([c.last_convergence for c in tum_all_n],axis=0)/constants.CHAINS

    average_n_dict = {'Name': file_name,
          'N status': 'Fix',
          'H_SEE': tum_average_n.H_SEE,
          'phi_SEE': tum_average_n.phi_SEE,
          'pi_SEE': tum_average_n.pi_SEE,
          'PZ_SEE': tum_average_n.PZ_SEE,
          'n_SEE': tum_average_n.n_SEE,
          'h_theta_SEE': tum_average_n.h_theta_SEE
          }
    result_df = result_df.append(average_n_dict, ignore_index=True)

    #vis_2.visualizing_inferred_variables(tum_average_n)
    #tum_average_n.save_in_txt(sample_1,constants.RESULTS+'/n_average_'+file_name+'.txt')
    pickle.dump(tum_average_n, open(constants.RESULTS+'/tum_average_n', 'wb'))
    #for c in range(constants.CHAINS):
    #    vis_2.vector_plot(np.arange(0, 1, 0.05),tum_all_n[c].Z_SEE,'Z_SEE'+str(c))
    #vis_2.vector_plot(np.arange(0, 1, 0.05),tum_average_n.Z_SEE,'Z_SEE_all')
    #vis_2.error_boxplot(tum_all_n)
    #vis_2.likelihood_all(tum_all_n, 'likelihood_all')



    #tum_average.H_SEE = np.sum([c.H_SEE for c in tum_all],axis=0)/constants.CHAINS
    #tum_average.phi_SEE = np.sum([c.phi_SEE for c in tum_all],axis=0)/constants.CHAINS
    #tum_average.pi_SEE = np.sum([c.pi_SEE for c in tum_all],axis=0)/constants.CHAINS
    #tum_average.PZ_SEE = np.sum([c.PZ_SEE for c in tum_all],axis=0)/constants.CHAINS
    #tum_average.n_SEE = np.sum([c.n_SEE for c in tum_all],axis=0)/constants.CHAINS

print(result_df)
pickle.dump(result_df, open(constants.RESULTS + '/result_df_'+number, 'wb'))
vis_3 = vis.visualization(constants.RESULTS + '/' + constants.VISUALIZATION)
#vis_3.error_bar(result_df, 'H_SEE')
#vis_3.error_plot_connected( result_df, 'H_SEE')

#vis_3.error_bar(result_df, 'phi_SEE')
#vis_3.error_plot_connected( result_df, 'phi_SEE')

#vis_3.error_bar(result_df, 'PZ_SEE')
#vis_3.error_plot_connected( result_df, 'PZ_SEE')

#vis_3.error_bar(result_df, 'pi_SEE')
#vis_3.error_plot_connected( result_df, 'pi_SEE')

#vis_3.error_bar(result_df, 'n_SEE')
#vis_3.error_plot_connected( result_df, 'n_SEE')

#vis_3.error_boxplot(tum_all)
#vis_3.error_boxplot(tum_all_n)

