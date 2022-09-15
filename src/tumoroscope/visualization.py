import matplotlib.pyplot as plt
import os
from collections import Counter
import numpy as np
import scipy.stats as stats
import seaborn as sns
import pandas as pd
from tumoroscope import constants
from matplotlib.font_manager import FontProperties
import pickle

'''
    def plot_pie_chart_H(self,n_s_data, tum_1, col,section,tree_clones):
        fontP = FontProperties()
        fontP.set_size('xx-small')
        fig, ax = plt.subplots()
        #for xi, yi in zip(list(n_s_data['x'][n_lambda.index.tolist()]), list(n_s_data['y'][n_lambda.index.tolist()])):
        for j in range(len(n_s_data['x'])):
            values = np.divide(tum_1.inferred_H[j, :],sum(tum_1.inferred_H[j, :]))#[float(i) / sum(tum_1.inferred_H[counter, :]) for i in tum_1.inferred_H[counter, :]]
            ax.pie(values, radius=0.5, colors=col, center=(n_s_data['x'][j], n_s_data['y'][j]))
            ax.set(ylabel='', title='', aspect='equal')
        ax.autoscale()
        ax.legend(loc='upper left', bbox_to_anchor=(1.001, 1), labels=tree_clones, prop=fontP)
        plt.savefig(self.dir + '/' + section + '_pie_chart_H.png', dpi=500)
        plt.close()

    def plot_pie_chart_HZ(self,n_s_data, tum_1, col,section,tree_clones):
        fontP = FontProperties()
        fontP.set_size('xx-small')
        fig, ax = plt.subplots()
        counter = 0
        #for xi, yi in zip(list(n_s_data['x'][n_lambda.index.tolist()]), list(n_s_data['y'][n_lambda.index.tolist()])):
            #values = [float(i) / sum(tum_1.inferred_H[counter, :] * tum_1.inferred_Z[counter, :]) for i in
            #          tum_1.inferred_H[counter, :] * tum_1.inferred_Z[counter, :]]
        for j in range(len(n_s_data['x'])):
            values = np.divide(tum_1.inferred_H[counter, :] * tum_1.inferred_Z[counter, :],sum(tum_1.inferred_H[counter, :] * tum_1.inferred_Z[counter, :]))#[float(i) / sum(tum_1.inferred_H[counter, :]) for i in tum_1.inferred_H[counter, :]]
            xi =  n_s_data['x'][j]
            yi = n_s_data['y'][j]
            if np.sum(np.isnan(values)) > 0:
                values = [float(i) / sum(tum_1.inferred_H[counter, :]) for i in
                          tum_1.inferred_H[counter, :]]
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

        plt.savefig(self.dir + '/' + section + '_pie_chart_HZ.png', dpi=500)
        plt.close()
'''


class visualization:
    def __init__(self,dir):
        self.dir = dir
        # Create target directory & all intermediate directories if don't exists
        if not os.path.exists(self.dir):
            os.makedirs(self.dir)
            print("Directory ", self.dir, " Created ")
        else:
            print("Directory ", self.dir, " already exists")
        self.plot_colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

    def barplot_matrix(self,A,name):
        temp = (Counter(A.flatten("C")))
        vals = list(temp.values())
        keys = list(temp.keys())
        fig = plt.figure()
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        ax.bar(keys,vals)
        plt.savefig(self.dir + '/'+name+'.png')
        plt.close()

    def hist_matrix(self,X,num_bins,name):
        fig = plt.figure()
        fig.add_axes([0.1,0.1,0.8,0.8])
        n, bins, patches = plt.hist(X.flatten("C"), num_bins, facecolor='blue', alpha=0.5,edgecolor="black",linewidth=1.2)
        plt.savefig(self.dir + '/'+name+'.png')
        plt.close()

    def vector_plot(self,X,Y,name):
        plt.plot(X,Y)
        plt.savefig(self.dir + '/' + name + '.png')
        plt.close()

    def binary_bw(self,Z,name):
        #plt.figure()
        #plt.imshow(Z, cmap='Greys', extent=[1, len(Z[0]), 1, len(Z)])
        plt.matshow(Z, cmap='Greys')
        plt.savefig(self.dir + '/' + name + '.png')
        plt.close()

    def gamma(self,alpha,beta,name):
        #plt.figure()
        x = np.linspace (0, 40, 200)
        y1 = stats.gamma.pdf(x=x, a=alpha, scale=beta) #a is alpha, loc is beta???
        plt.plot(x, y1, "b-", label=(r'$\alpha=$'+str(alpha)+r'$, scale=$'+str(beta)))
        plt.legend()
        plt.savefig(self.dir + '/' + name + '.png')
        plt.close()
    def beta(self,alpha,beta,name):
        rv1 = stats.beta(alpha, beta)
        x = np.linspace(0, 1, 100)
        plt.plot(x, rv1.pdf(x))
        plt.savefig(self.dir + '/' + name + '.png')
        plt.close()

    def likelihood(self, loglik,name):
        #plt.figure()
        plt.plot(range(len(loglik)),loglik)
        plt.savefig(self.dir + '/' + name + '.png')
        plt.close()

    def likelihood_all(self, obj,name):
        #plt.figure()
        for c in obj:
            plt.plot(range(len(c.loglik)),c.loglik)
        plt.savefig(self.dir + '/' + name + '.png')
        plt.close()

    def heatmap_seaborn(self,data,name,x_label,y_label,annot,linewidths):
        sns.set_theme()
        ax = sns.heatmap(data, annot=annot, linewidths=linewidths)
        ax.set(xlabel=x_label, ylabel=y_label)
        fig = ax.get_figure()
        fig.savefig(self.dir + '/' + name + '.png')
        plt.close()

    def visualizing_simulated_data(self,sample_1):
        self.barplot_matrix(sample_1.A, 'A')
        self.barplot_matrix(sample_1.D, 'D')

        self.hist_matrix(sample_1.phi, 20, 'phi')
        self.hist_matrix(sample_1.G, 20, 'G')
        self.hist_matrix(sample_1.H, 20, 'H')

        self.binary_bw(sample_1.A, 'A_bw')
        self.binary_bw(sample_1.D, 'D_bw')
        self.binary_bw(sample_1.C, 'C')
        self.binary_bw(sample_1.Z, 'Z_bw')
        self.binary_bw(sample_1.H, 'H_bw')
        self.binary_bw(sample_1.phi, 'phi_bw')

        self.heatmap_seaborn(sample_1.Z, 'Z_seaborn','clones','spots',False,0)
        self.heatmap_seaborn(sample_1.H, 'H_seaborn','clones','spots',False,0)
        self.heatmap_seaborn(sample_1.phi, 'phi_seaborn','clones','mutations',False,0)
        self.heatmap_seaborn(sample_1.A, 'A_seaborn','spots','mutations',False,0)
        self.heatmap_seaborn(sample_1.D, 'D_seaborn','spots','mutations',False,0)
        self.heatmap_seaborn(sample_1.C, 'C_seaborn', 'clones', 'mutations', False, 0.5)

        self.plot_F_gamma(sample_1.F_epsilon, sample_1.F,self.plot_colors)
        self.gamma(sample_1.r, sample_1.p, 'phi_gamma')
        self.gamma(sample_1.F_epsilon[0, 0], sample_1.F_epsilon[0, 1], 'H_gamma_F_epsilon')
        for k in range(len(sample_1.F)):
            self.gamma(alpha=sample_1.F[k, 0], beta=sample_1.F[k, 1], name=('H_gamma_' + str(k)))
        if sample_1.theta_variable == True:
            self.beta(sample_1.theta_alpha[0][0], sample_1.theta_beta[0][0], 'beta_0_0')


    def visualizing_real_data(self,A,D,C,K,F_epsilon,F,r,p,theta_variable,theta_alpha,theta_beta):
        self.barplot_matrix(A, 'A')
        self.barplot_matrix(D, 'D')
        self.binary_bw(A, 'A_bw')
        self.binary_bw(D, 'D_bw')
        self.binary_bw(C, 'C')
        self.heatmap_seaborn(A, 'A_seaborn','spots','mutations',False,0)
        self.heatmap_seaborn(D, 'D_seaborn','spots','mutations',False,0)
        self.heatmap_seaborn(C, 'C_seaborn', 'clones', 'mutations', False, 0.5)

        self.plot_F_gamma(F_epsilon, F)
        self.gamma(r, p, 'phi_gamma')
        self.gamma(F_epsilon[0, 0], F_epsilon[0, 1], 'H_gamma_F_epsilon')
        for k in range(K):
            self.gamma(alpha=F[k, 0], beta=F[k, 1], name=('H_gamma_' + str(k)))
        if theta_variable == True:
            self.beta(theta_alpha[0][0], theta_beta[0][0], 'beta_0_0')



    def visualizing_inferred_variables(self,tum_1):
        self.binary_bw(tum_1.inferred_Z, 'Z_bw_inferred')
        self.binary_bw(tum_1.inferred_P_Z, 'P_Z_bw_inferred')
        self.binary_bw(tum_1.inferred_H, 'H_bw_inferred')
        self.binary_bw(tum_1.inferred_phi, 'phi_bw_inferred')
        self.likelihood(tum_1.loglik, 'loglik')
        self.binary_bw(tum_1.decision_matrix_G / tum_1.max_iter, 'acceptance_G')
        self.binary_bw(tum_1.decision_matrix_phi / tum_1.max_iter, 'acceptance_phi')

        self.heatmap_seaborn(tum_1.inferred_P_Z, 'P_Z_bw_inferred_seaborn','clones','spots',False,0)
        self.heatmap_seaborn(tum_1.inferred_H, 'H_bw_inferred_seaborn','clones','spots',False,0)
        self.heatmap_seaborn(tum_1.inferred_phi, 'phi_bw_inferred_seaborn','clones','mutations',False,0)
        self.heatmap_seaborn(tum_1.decision_matrix_G / tum_1.max_iter, 'acceptance_G_seaborn','clones','spots',False,0)
        self.heatmap_seaborn(tum_1.decision_matrix_phi / tum_1.max_iter, 'acceptance_phi_seaborn','clones','mutations',False,0)
        #for rate in range(6):
        #    self.vector_plot(range(len(tum_1.convergence_rate[0][rate])),tum_1.convergence_rate[0][rate],'convergence_rate_'+str(rate*10))

    def error_bar(self,result_df,error_name):
        df = pd.DataFrame({'Fix': result_df[error_name][result_df['N status'] == 'Fix'].to_list(),
                      'Variable': result_df[error_name][result_df['N status'] == 'Variable'].to_list()},
                     index=result_df['Name'][result_df['N status'] == 'Fix']
                     )
        fig = df.plot.bar(rot=0).get_figure()
        fig.savefig(self.dir + '/error_barplot_'+error_name+'.png')
        plt.close()

    def error_plot_connected(self,result_df,error_name):
        plot1, = plt.plot(result_df['Name'][result_df['N status'] == 'Fix'],
                          result_df[error_name][result_df['N status'] == 'Fix'],
                          'xb-')
        plot2, = plt.plot(result_df['Name'][result_df['N status'] == 'Variable'],
                          result_df[error_name][result_df['N status'] == 'Variable'], 'xr-')
        plt.legend([plot1, plot2], ["Fix n", "Variable n"])

        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(self.dir + '/error_connected_plot_'+error_name+'.png')
        plt.close()

    def error_boxplot(self,tum_all):
        H_chains = []
        for c in range(constants.CHAINS):
            H_chains.append(tum_all[c].H_SEE)
        ax = sns.boxplot(data=H_chains)
        ax = sns.swarmplot(data=H_chains)
        ax.set(xticklabels=['H SEE'])
        fig = ax.get_figure()
        fig.savefig(self.dir + '/' + 'H_error_barplot.png')
        plt.close()

        phi_chains = []
        for c in range(constants.CHAINS):
            phi_chains.append(tum_all[c].phi_SEE)
        ax = sns.boxplot(data=phi_chains)
        ax = sns.swarmplot(data=phi_chains)
        ax.set(xticklabels=['phi SEE'])
        fig = ax.get_figure()
        fig.savefig(self.dir + '/' + 'phi_error_barplot.png')
        plt.close()

        pi_chains = []
        for c in range(constants.CHAINS):
            pi_chains.append(tum_all[c].pi_SEE)
        ax = sns.boxplot(data=pi_chains)
        ax = sns.swarmplot(data=pi_chains)
        ax.set(xticklabels=['pi SEE'])
        fig = ax.get_figure()
        fig.savefig(self.dir + '/' + 'pi_error_barplot.png')
        plt.close()

        n_chains = []
        for c in range(constants.CHAINS):
            n_chains.append(tum_all[c].n_SEE)
        ax = sns.boxplot(data=n_chains)
        ax = sns.swarmplot(data=n_chains)
        ax.set(xticklabels=['n SEE'])
        fig = ax.get_figure()
        fig.savefig(self.dir + '/' +  'n_error_barplot.png')
        plt.close()

    def plot_F_gamma(self,F_epsilon, F,plot_colors):
        x = np.linspace(0, 40, 1000)
        y1 = stats.gamma.pdf(x=x, a=F_epsilon[0][0], scale=1 / F_epsilon[0][1])  # a is alpha, loc is beta???
        plt.plot(x, y1, plot_colors[0],
                 label=(r'$F_\epsilon: \alpha=$' + str(F_epsilon[0][0]) + r'$, \beta=$' + str(F_epsilon[0][1])))
        for f_gamma in range(len(F)):
            y1 = stats.gamma.pdf(x=x, a=F[f_gamma][0], scale=1 / F[f_gamma][1])  # a is alpha, loc is beta???
            plt.plot(x, y1, plot_colors[f_gamma + 1], label=(
                        'F_C' + str(f_gamma) + r': $\alpha=$' + str(np.round(F[f_gamma][0], 2)) + r'$, \beta=$' + str(
                    F[f_gamma][1])))
        plt.legend()
        plt.savefig(self.dir + '/' + 'F_distributions.png', dpi=500)
        plt.close()

    def plot_C_ik(self,csv_file):
        df = pd.read_csv(csv_file, sep='\t')
        print(np.sum(df>0))
        sns.set(font_scale=1.5)
        g = sns.clustermap(df, cmap="Blues",  yticklabels=False)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=20, fontweight='bold')
        g.ax_row_dendrogram.set_visible(False)
        g.ax_col_dendrogram.set_visible(False)
        g.ax_row_dendrogram.set_xlim([0, 0])
        g.fig.suptitle('', fontsize=22)
        g.ax_heatmap.set_title('Genotype', fontsize=22, fontweight='bold')
        g.ax_heatmap.set_xlabel("")
        g.ax_heatmap.set_ylabel("Mutations", fontsize=20, fontweight='bold')
        g.ax_heatmap.yaxis.set_label_position('left')
        g.fig.subplots_adjust(right=0.75)
        g.ax_cbar.set_position((0.8, .2, .03, .4))
        #g.ax_heatmap.legend( title="Mean copy number", fontsize=18, fancybox=True, fontweight='bold')
        g.savefig(self.dir+"/C_ik.png", dpi=500)

    def plot_adj_far_corr(self,adj_corrs, far_corrs, tree_clones, sample_count, visualization_dir,title,labels):


        bar_width = 0.35
        col = (0.02, 0.20, 0.40)  # sns.color_palette()[0]
        x = np.arange(len(labels))
        plt.style.use('classic')
        fig, ax = plt.subplots()

        # index = np.arange(len(adj_corrs))
        #a = fig.gca()
        #a.set_frame_on(False)
        #a.spines["top"].set_visible(False)
        #a.spines["right"].set_visible(False)
        #a.spines["left"].set_visible(False)

        adj = ax.bar(x - bar_width / 2, adj_corrs, bar_width, label="adjacent", edgecolor=col, color=col)
        far = ax.bar(x + bar_width / 2, far_corrs, bar_width, label="far", edgecolor=col, color='white')

        ax.legend()
        #ax.set_frame_on(True)
        #ax.set_facecolor('white')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=14, fontweight='bold')
        ax.set_yticklabels(np.round(ax.get_yticks(), 2), fontsize=14, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlabel("", fontsize=14, fontweight='bold')
        ax.set_ylabel("Correlation", fontsize=12, fontweight='bold')
        # ax.set_xlabel('Category')
        # ax.set_ylabel('Incidence')
        # ax.set_title('Crime incidence by season, type')
        # ax.set_xticks(index + bar_width / 2)
        # ax.set_xticklabels(["ASB", "Violence", "Theft", "Public Order", "Drugs"])
        ax.legend()
        bar_color = adj[0].get_facecolor()
        for bar in adj:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.005,
                round(bar.get_height(), 2),
                color=bar_color,
                horizontalalignment='center'
            )
        # bar_color = far[0].get_facecolor()
        for bar in far:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.005,
                round(bar.get_height(), 2),
                color=bar_color,
                horizontalalignment='center'
            )
        fig.tight_layout()
        plt.savefig(visualization_dir + '/corr_bar_' + str(sample_count) + '.png', dpi=500)


    def sim_error_box(self, obj, var_name, file_name,title):
        props = dict(boxes="DarkGreen", whiskers="DarkOrange", medians="DarkBlue", caps="Gray")
        fig = plt.figure()
        #ax = obj.boxplot( notch=True,patch_artist=True)
        ax = obj.boxplot()

        #ax.grid(False)
        ax.grid(axis='x')
        ax.set_title(title, fontsize=18)
        ax.set_xlabel("", fontsize=16)
        ax.set_ylabel("mean absolute error (MAE)", fontsize=16)

        a = fig.gca()
        #a.set_frame_on(False)
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
        a.spines["left"].set_visible(False)

        #plt.rcParams["axes.grid"] = False
        plt.savefig(self.dir + '/'+ var_name + file_name, dpi=500)
        plt.close()


    def sim_error_box_2(self,df, var_name, file_name,title):
        fig, ax = plt.subplots()

        col = ['yellow', 'blue', 'blue', 'green', 'green']
        bplot = sns.boxplot(data=df)
        ax = fig.gca()
        #a.set_frame_on(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        for i in range(len(col)):
            mybox = bplot.artists[i]
            mybox.set_facecolor(col[i])
        ax.yaxis.grid(True)  # Hide the horizontal gridlines
        ax.xaxis.grid(False)
        col = (0.02, 0.20, 0.40)
        #plt.rc('xtick', labelsize=16,weight='bold',color=col)
        #plt.rc('ytick', labelsize=16,weight='bold',color=col)
        x = list(range(1, 6, 1))
        #y = list(range(0, 1, 0.02))
        ax.set_xticklabels(x, fontsize=18,fontweight='bold',color=col)
        ax.set_yticklabels(np.round(ax.get_yticks(),4), fontsize=14, fontweight='bold',color=col)
        ax.set_title(title, fontsize=18,color=col, fontweight='bold')
        ax.set_xlabel("", fontsize=18,color=col, fontweight='bold')
        ax.set_ylabel("mean absolute error (MAE)", fontsize=16,fontweight='bold',color=col)
        fig.tight_layout()
        plt.savefig(self.dir + '/' + var_name + file_name, dpi=500,figsize=(10, 6))
        plt.close()


    def sim_error_box_subplot(self,df,title,ax,col,y_axis,x_tick_labels):
        font_size = 6
        bplot = sns.boxplot(data=df)
        ##df.boxplot()
        ##ax = fig.gca()
        #ax.set_frame_on(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        for i in range(len(col)):
            mybox = bplot.artists[i]
            mybox.set_facecolor(col[i])
        ax.yaxis.grid(True)  # Hide the horizontal gridlines
        ax.xaxis.grid(False)
        col = (0.02, 0.20, 0.40)
        #plt.rc('xtick', labelsize=16,weight='bold',color=col)
        #plt.rc('ytick', labelsize=16,weight='bold',color=col)
        #x = list(range(1, 6, 1))
        #y = list(range(0, 1, 0.02))
        ax.set_xticklabels(x_tick_labels, fontsize=font_size,color=col,fontweight='bold') # removed fontweight='bold'
        ax.set_yticklabels(np.round(ax.get_yticks(),4), fontsize=font_size,color=col,fontweight='bold')
        ax.set_title(title, fontsize=font_size,color=col,fontweight='bold')
        ax.set_xlabel("", fontsize=font_size,color=col,fontweight='bold')
        ax.set_ylabel(y_axis, fontsize=font_size,color=col,fontweight='bold')
        #fig.tight_layout()
        return  ax


    def plot_pie_chart_H(self, n_s_data, inferred_H, col, section, tree_clones):
        fontP = FontProperties()
        fontP.set_size('xx-small')
        fig, ax = plt.subplots()
        # for xi, yi in zip(list(n_s_data['x'][n_lambda.index.tolist()]), list(n_s_data['y'][n_lambda.index.tolist()])):
        for j in range(len(n_s_data['x'])):
            values = np.divide(inferred_H.loc[j, :], sum(inferred_H.loc[j,:]))  # [float(i) / sum(tum_1.inferred_H[counter, :]) for i in tum_1.inferred_H[counter, :]]
            ax.pie(values.values, radius=0.5, colors=col, center=(n_s_data['x'][j], n_s_data['y'][j]))
            ax.set(ylabel='', title='', aspect='equal')
        ax.autoscale()
        ax.legend(loc='upper left', bbox_to_anchor=(1.001, 1), labels=tree_clones, prop=fontP)
        plt.savefig(self.dir  + '/' + section + '_pie_chart_H.png', dpi=500)
        plt.close()


    def plot_pie_chart_HZ(self, n_s_data, inferred_H,inferred_Z, col, section, tree_clones):
        fontP = FontProperties()
        fontP.set_size('xx-small')
        fig, ax = plt.subplots()
        counter = 0
        # for xi, yi in zip(list(n_s_data['x'][n_lambda.index.tolist()]), list(n_s_data['y'][n_lambda.index.tolist()])):
        # values = [float(i) / sum(tum_1.inferred_H[counter, :] * tum_1.inferred_Z[counter, :]) for i in
        #          tum_1.inferred_H[counter, :] * tum_1.inferred_Z[counter, :]]
        for j in range(len(n_s_data['x'])):
            values = np.divide(inferred_H.loc[counter, :] * inferred_Z.loc[counter, :],
                               sum(inferred_H.loc[counter, :] * inferred_Z.loc[counter,
                                                                  :]))  # [float(i) / sum(tum_1.inferred_H[counter, :]) for i in tum_1.inferred_H[counter, :]]
            xi = n_s_data['x'][j]
            yi = n_s_data['y'][j]
            if np.sum(np.isnan(values)) > 0:
                values = [float(i) / sum(inferred_H.loc[counter, :]) for i in
                          inferred_H.loc[counter, :]]
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

        plt.savefig(self.dir  + '/' + section + '_pie_chart_HZ.png', dpi=500)
        plt.close()

    def plot_variance_read_D(self,result_obj,number,type):
        fig, ax = plt.subplots()
        tum_1 = pickle.load(open(result_obj, 'rb'))
        batch = 1000
        every = 5
        iter = tum_1.iter

        last_convergence = tum_1.last_convergence
        converged_batch = [idx for idx, element in enumerate(last_convergence) if
                           last_convergence[idx] == max(last_convergence)]
        #print((iter+1)/every)
        samples = tum_1.H[converged_batch[0] * batch:int((iter+1)/every)]

        variances = np.var(samples, axis=0)
        DC = np.matmul(np.transpose(tum_1.D), tum_1.C)
        corr = np.corrcoef(DC, variances)[0][1]
        ax.set_title(type + '_D*C_' + str(number), fontsize=18, fontweight='bold')
        plt.scatter(DC, variances)
        plt.xlabel('Reads')
        plt.ylabel('Variance')
        plt.savefig(self.dir + '/variance_read_D_' +type+'_'+ str(number) + '.png', dpi=500)
        plt.close()
        return corr
    def plot_variance_read_A(self,result_obj,number,type):
        fig, ax = plt.subplots()
        tum_1 = pickle.load(open(result_obj, 'rb'))

        batch = 1000
        every = 5
        iter = tum_1.iter

        last_convergence = tum_1.last_convergence
        converged_batch = [idx for idx, element in enumerate(last_convergence) if
                           last_convergence[idx] == max(last_convergence)]
        samples = tum_1.H[converged_batch[0] * batch:int((iter+1)/every)]

        variances = np.var(samples, axis=0)
        AC = np.matmul(np.transpose(tum_1.A), tum_1.C)
        corr = np.corrcoef(AC, variances)[0][1]
        ax.set_title('corr: '+str(np.round(corr,2))+type + '_A*C_' + str(number), fontsize=18, fontweight='bold')
        plt.scatter(AC, variances)
        plt.xlabel('Reads')
        plt.ylabel('Variance')
        plt.savefig(self.dir + '/variance_read_A_' +type+'_'+ str(number) + '.png', dpi=500)
        plt.close()
        return corr

    def plot_variance_read_spots_A(self,result_obj,number,type,x_q,y_q,batch_size,every_size):
        fig, ax = plt.subplots()
        tum_1 = pickle.load(open(result_obj, 'rb'))

        batch = batch_size
        every = every_size
        iter = tum_1.iter

        last_convergence = tum_1.last_convergence
        converged_batch = [idx for idx, element in enumerate(last_convergence) if
                           last_convergence[idx] == max(last_convergence)]
        samples = tum_1.H[int((converged_batch[0] * batch)/every):int((iter+1)/every)]


        variances = np.var(samples, axis=0)
        spots_var = np.sum(variances, axis=1)
        corr = np.corrcoef(spots_var,np.sum(tum_1.A,axis=0))[0][1]

        if type == 'real':
            print('here')
        ax.set_title('corr: '+str(np.round(corr,2))+'                  '+type + '_A_sum_' + str(number), fontsize=18, fontweight='bold')
        plt.axhline(y=np.min(spots_var)+y_q * (np.max(spots_var)-np.min(spots_var)), color='gray', linestyle='--')
        plt.axvline(x=np.min(np.sum(tum_1.A,axis=0)) +x_q * (np.max(np.sum(tum_1.A,axis=0))-np.min(np.sum(tum_1.A,axis=0))), color='gray', linestyle='--')
        plt.scatter(np.sum(tum_1.A,axis=0), spots_var)
        plt.xlabel('Reads')
        plt.ylabel('Variance')
        plt.savefig(self.dir + '/variance_read_spots_A_' +type+'_'+ str(number) + '.png', dpi=500)
        plt.close()
        return corr

    def plot_variance_read_spots_D(self,result_obj,number,type,x_q,y_q,batch_size,every_size,save_plot):
        #fig, ax = plt.subplots()
        tum_1 = pickle.load(open(result_obj, 'rb'))[0]

        batch = batch_size
        every = every_size
        iter = tum_1.iter

        last_convergence = tum_1.last_convergence
        converged_batch = [idx for idx, element in enumerate(last_convergence) if
                           last_convergence[idx] == max(last_convergence)]
        samples = tum_1.H[int((converged_batch[0] * batch)/every):int((iter+1)/every)]

        variances = np.var(samples, axis=0)
        spots_var = np.sum(variances, axis=1)
        corr = np.corrcoef(spots_var,np.sum(tum_1.D, axis=0))[0][1]
        if save_plot==True:
            fig, ax = plt.subplots()
            ax.set_title('corr: '+str(np.round(corr,2))+'                  '+ type + '_D_sum_' + str(number), fontsize=18, fontweight='bold')
            plt.scatter(np.sum(tum_1.D,axis=0), spots_var)
            plt.axhline(y=np.min(spots_var)+y_q * (np.max(spots_var)-np.min(spots_var)), color='gray', linestyle='--')
            plt.axvline(x=np.min(np.sum(tum_1.D,axis=0)) + x_q * (np.max(np.sum(tum_1.D,axis=0))-np.min(np.sum(tum_1.D,axis=0))), color='gray', linestyle='--')
            plt.xlabel('Reads')
            plt.ylabel('Variance')
            plt.savefig(self.dir + '/variance_read_spots_D_' +type+'_'+ str(number) + '.png', dpi=500)
            plt.close()
        return corr

    def plot_variance_error_simulation(self,number, type,result_1,x_q,y_q,batch_size,every_size,fix,sim_res_all):
        fig, ax = plt.subplots()
        #if fix==True:
        #    result_1 = sim_res_all + str(number) + '/' + type + '_n_chain_1'
        #else:
        #    result_1 = sim_res_all + str(number) + '/' + type + '_chain_1'
        tum = pickle.load(open(result_1, 'rb'))[0]
        last_convergence = tum.last_convergence
        converged_batch = [idx for idx, element in enumerate(last_convergence) if
                           last_convergence[idx] == max(last_convergence)]
        batch = batch_size
        every = every_size
        iter = tum.iter

        samples_1 = tum.H[int((converged_batch[0] * batch)/every):int((iter+1)/every)]
        true_sample_1 = sim_res_all +str(number)+ '/sample_1_' + type
        true_sample = pickle.load(open(true_sample_1, 'rb'))
        variances_1 = np.abs(np.var(samples_1, axis=0))
        err = np.abs(true_sample.H - tum.inferred_H)


        corr = np.corrcoef(variances_1, err)[0][1]

        ax.set_title('corr: '+str(np.round(corr,2))+'                  '+'simulation_'+type+'_'+str(number), fontsize=18, fontweight='bold')
        plt.axhline(y=np.min(variances_1)+y_q*(np.max(variances_1)-np.min(variances_1)), color='gray', linestyle='--')
        plt.axvline(x=np.min(err)+x_q * (np.max(err)-np.min(err)), color='gray', linestyle='--')

        plt.scatter(err, variances_1)
        plt.xlabel('Error')
        plt.ylabel('Variance')
        plt.savefig(self.dir + '/variance_error_simulation_'+type+'_' + str(number) + '.png', dpi=500)
        plt.close()
        return corr

    def plot_read_error_simulation(self,number, type,result_1,x_q,y_q,batch_size,every_size,save_plot,true_sample):
        tum = pickle.load(open(result_1, 'rb'))[0]
        last_convergence = tum.last_convergence
        converged_batch = [idx for idx, element in enumerate(last_convergence) if
                           last_convergence[idx] == max(last_convergence)]
        batch = batch_size
        every = every_size
        iter = tum.iter

        #samples_1 = tum.H[int((converged_batch[0] * batch)/every):int((iter+1)/every)]
        #true_sample_1 = sim_res_all + str(number)  + '/sample_1_' + type
        #true_sample = pickle.load(open(true_sample_1, 'rb'))


        err = np.mean(np.abs(true_sample.H - tum.inferred_H),axis=1)
        D_reads = np.sum(tum.D,axis=0)
        A_reads = np.sum(tum.A, axis=0)
        #print(len(err))
        #print(len(D_reads))
        #print(len(A_reads))

        corr_A = np.corrcoef(A_reads, err)[0][1]
        corr_D = np.corrcoef(D_reads, err)[0][1]
        if save_plot==True:
            fig, ax = plt.subplots()
            ax.set_title('corr: '+str(np.round(corr_A,2))+'                  '+'simulation_'+type+'_'+str(number), fontsize=18, fontweight='bold')
            plt.axhline(y=np.min(A_reads)+y_q*(np.max(A_reads)-np.min(A_reads)), color='gray', linestyle='--')
            plt.axvline(x=np.min(err) + x_q * (np.max(err)-np.min(err)), color='gray', linestyle='--')

            plt.scatter(err, A_reads)
            plt.xlabel('Error')
            plt.ylabel('Reads')
            plt.savefig(self.dir + '/A_read_error_simulation_'+type+'_' + str(number) + '.png', dpi=500)
            plt.close()

            fig, ax = plt.subplots()
            ax.set_title('corr: '+str(np.round(corr_D,2))+'                  '+'simulation_'+type+'_'+str(number), fontsize=18, fontweight='bold')
            plt.axhline(y=np.min(D_reads)+y_q*(np.max(D_reads)-np.min(D_reads)), color='gray', linestyle='--')
            plt.axvline(x=np.min(err)+x_q * (np.max(err)-np.min(err)), color='gray', linestyle='--')

            plt.scatter(err, D_reads)
            plt.xlabel('Error')
            plt.ylabel('Reads')
            plt.savefig(self.dir + '/D_read_error_simulation_'+type+'_' + str(number) + '.png', dpi=500)
            plt.close()
        return corr_A,corr_D


    def plot_read_error_simulation_toghether(self,number, type,result_1,result_2,x_q,y_q,batch_size,every_size,save_plot,true_sample_1,true_sample_2):
        tum_1 = pickle.load(open(result_1, 'rb'))[0]
        err_1 = np.mean(np.abs(true_sample_1.H - tum_1.inferred_H),axis=1)
        D_reads_1 = np.sum(tum_1.D,axis=0)
        A_reads_1 = np.sum(tum_1.A, axis=0)

        tum_2 = pickle.load(open(result_2, 'rb'))[0]
        err_2 = np.mean(np.abs(true_sample_2.H - tum_2.inferred_H),axis=1)
        D_reads_2 = np.sum(tum_2.D,axis=0)
        A_reads_2 = np.sum(tum_2.A, axis=0)

        A_reads = [*A_reads_1 , *A_reads_2]
        D_reads = [*D_reads_1, *D_reads_2]
        err = [*err_1, *err_2]

        corr_A = np.corrcoef(A_reads, err)[0][1]
        corr_D = np.corrcoef(D_reads, err)[0][1]
        if save_plot==True:
            fig, ax = plt.subplots()
            ax.set_title('corr: '+str(np.round(corr_A,2))+'                  '+'simulation_'+type+'_'+str(number), fontsize=18, fontweight='bold')
            plt.axhline(y=np.min(A_reads)+y_q*(np.max(A_reads)-np.min(A_reads)), color='gray', linestyle='--')
            plt.axvline(x=np.min(err) + x_q * (np.max(err)-np.min(err)), color='gray', linestyle='--')

            plt.scatter(err, A_reads)
            plt.xlabel('Error')
            plt.ylabel('Reads')
            plt.savefig(self.dir + '/A_read_error_simulation_'+type+'_' + str(number) + '.png', dpi=500)
            plt.close()

            fig, ax = plt.subplots()
            ax.set_title('corr: '+str(np.round(corr_D,2))+'                  '+'simulation_'+type+'_'+str(number), fontsize=18, fontweight='bold')
            plt.axhline(y=np.min(D_reads)+y_q*(np.max(D_reads)-np.min(D_reads)), color='gray', linestyle='--')
            plt.axvline(x=np.min(err)+x_q * (np.max(err)-np.min(err)), color='gray', linestyle='--')

            plt.scatter(err, D_reads)
            plt.xlabel('Error')
            plt.ylabel('Reads')
            plt.savefig(self.dir + '/D_read_error_simulation_'+type+'_' + str(number) + '.png', dpi=500)
            plt.close()
        return corr_A,corr_D


    def return_plot_read_error(self,number, type,result_1,x_q,y_q,batch_size,every_size,ax,true_sample):
        tum = pickle.load(open(result_1, 'rb'))[0]
        last_convergence = tum.last_convergence
        converged_batch = [idx for idx, element in enumerate(last_convergence) if
                           last_convergence[idx] == max(last_convergence)]
        batch = batch_size
        every = every_size
        iter = tum.iter

        #samples_1 = tum.H[int((converged_batch[0] * batch)/every):int((iter+1)/every)]
        #true_sample_1 = sim_res_all + str(number)  + '/sample_1_' + type
        #true_sample = pickle.load(open(true_sample_1, 'rb'))


        err = np.mean(np.abs(true_sample.H - tum.inferred_H),axis=1)
        D_reads = np.sum(tum.D,axis=0)
        A_reads = np.sum(tum.A, axis=0)
        #print(len(err))
        #print(len(D_reads))
        #print(len(A_reads))

        corr_A = np.corrcoef(A_reads, err)[0][1]
        corr_D = np.corrcoef(D_reads, err)[0][1]

        ax.set_title('corr: '+str(np.round(corr_D,2))+'                  '+'simulation_'+type+'_'+str(number), fontsize=18, fontweight='bold')
        plt.axhline(y=np.min(D_reads)+y_q*(np.max(D_reads)-np.min(D_reads)), color='gray', linestyle='--')
        plt.axvline(x=np.min(err)+x_q * (np.max(err)-np.min(err)), color='gray', linestyle='--')

        plt.scatter(err, D_reads)
        plt.xlabel('Error')
        plt.ylabel('Reads')
        return ax


    def df_box(self,types,name,df):
        fig, ax = plt.subplots()
        boxplot = df.boxplot(column=types)
        plt.savefig(self.dir + name+ '.png', dpi=500)
        plt.close()

    def plot_spots_each_clone(self, n_s_data, inferred_H, section):
        import matplotlib as mpl
        fig = plt.figure(figsize=(50, 50))
        counter = 0
        color = 'magma'
        l = []
        for clone in inferred_H:
            ax = fig.add_subplot(3, 5, (counter + 1))
            l.append(ax.scatter(n_s_data['x'], n_s_data['y'], cmap=color, c=inferred_H[clone]))
            ax.set_title("clone" + clone, fontsize=9, fontweight='bold')
            plt.tick_params(
                axis='x',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                bottom=False,  # ticks along the bottom edge are off
                top=False,  # ticks along the top edge are off
                left=False,
                labelbottom=False, labelleft=False)
            plt.tick_params(
                axis='y',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                bottom=False,  # ticks along the bottom edge are off
                top=False,  # ticks along the top edge are off
                left=False,
                labelbottom=False, labelleft=False)
            counter = counter + 1
        # plt.tight_layout()
        # Create the sub-plots, assigning a different color for each line.
        # Also store the line objects created

        ax = fig.add_subplot(3, 40, 34)
        cmap = mpl.cm.magma
        norm = mpl.colors.Normalize(vmin=0, vmax=1)

        cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                        norm=norm,
                                        orientation='vertical')
        cb1.set_label('Fraction')

        plt.savefig(self.dir + '/Per_clone_' +section + '_H.png', dpi=100,figsize=(8, 6))
        plt.close()
