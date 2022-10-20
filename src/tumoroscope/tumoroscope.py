import pickle
import numpy as np
import operator as op
from functools import reduce
import scipy
from tumoroscope.simulation import simulation as sim
from scipy.stats import norm
import time
import math
import copy
import scipy.special as sc

class tumoroscope:

    def __init__(self,name, K, S, r, p, I, avarage_clone_in_spot,F,C,A,D,F_epsilon,optimal_rate,n_lambda,gamma,pi_2D,result_txt):
        self.name = name
        self.K=K#5
        self.S=S#10
        # S.alpha/(1+(alpha/K)) = 2*S
        #self.r=r#2
        #self.p=p#2

        adjust = np.tile(n_lambda, (I, 1))
        phi_samples = np.mean(D/(adjust.astype(float)),axis=1)#,(K,1)).transpose()
        phi_samples[phi_samples==0] = 10000
        phi_samples[phi_samples==10000] = np.min(phi_samples)
        N = I
        nomin = N*np.sum(phi_samples)
        denomin = ((N*np.sum(phi_samples*np.log(phi_samples)))-(np.sum(np.log(phi_samples))*np.sum(phi_samples)))
        self.r = nomin/denomin
        self.p = (denomin*(1/(N*N)))


        self.I=I#10
        self.avarage_clone_in_spot=avarage_clone_in_spot#2
        self.zeta=(avarage_clone_in_spot*K)/(K-avarage_clone_in_spot)

        #self.cons = 5
        self.F = F
        self.C = C
        self.A = A
        self.D = D
        self.converged = 0
        self.reject = 0
        self.acceptance = 0
        self.theta_0 = 0.95
        self.F_epsilon = F_epsilon
        self.sigma_G = None
        self.sigma_phi = None
        self.optimal_rate = optimal_rate
        self.n_lambda = n_lambda
        self.gamma = gamma
        self.result_txt = result_txt
        self.convergence_rate = np.array([])

        adjust = np.tile(n_lambda, (len(A), 1))
        self.sigma_phi_initial = (np.tile(np.mean(D/(adjust.astype(float)),axis=1),(K,1)).transpose())/10
        self.sigma_phi_initial[self.sigma_phi_initial==0] = 100000
        self.sigma_phi_initial[self.sigma_phi_initial==10000] = np.min(self.sigma_phi_initial)
        self.sigma_phi_changes = 0.1
        self.sigma_G_initial = np.tile(F[:,0]/10,(S,1))
        self.sigma_G_changes = 1
        self.sigma_n_initial = (n_lambda.astype(float))/10
        self.sigma_n_changes = 1


    def gibbs_sampling(self,seed,min_iter,max_iter,burn_in,batch,simulated_data,n_sampling,F_fraction,theta_variable,pi_2D,th,every_n_sample,changes_batch,var_calculation):

        self.changes_batch = changes_batch
        self.var_calculation = var_calculation
        convergence_counter = 0
        first = time.time()
        self.theta_variable = theta_variable
        self.max_iter = max_iter
        length = np.zeros(max_iter)
        length[0] = 0
        current_batch = 0
        iter = 0
        PZ=0
        n,H,G,pi,phi,Z,h_theta = self.initialization(K = self.K,S = self.S,r = self.r,p = self.p,I = self.I,
							F = self.F, avarage_clone_in_spot = self.avarage_clone_in_spot,
							seed = seed,max_iter = max_iter,batch = batch, D = self.D,
							 A = self.A, C = self.C,F_epsilon = self.F_epsilon,
							theta=self.theta_0,n_lambda = self.n_lambda,n_sampling = n_sampling,
							F_fraction=F_fraction,theta_variable = theta_variable,gamma=self.gamma,pi_2D=pi_2D,every_n_sample=every_n_sample)
        loglik = self.log_likelihood_model(A=self.A, D=self.D, H=H,phi=phi, C=self.C, n=n, I=self.I,theta=self.theta_0,theta_variable=theta_variable,h_theta=h_theta)
        n_pre,H_pre,G_pre,pi_pre,phi_pre,Z_pre,h_theta_pre,loglik_pre = self.save_and_update_variables(n,H,G,pi,phi,Z,PZ,h_theta,loglik,iter,every_n_sample,current_batch)


        for iter in range(1,(max_iter-1)):
            begin= time.time()
            self.iter = iter
            if n_sampling==False:
                Z,PZ = self.sample_Z_matrix_log( S = self.S, K = self.K, G = G_pre, Z = Z_pre,pi = pi_pre,F = self.F,F_epsilon = self.F_epsilon,burn_in = burn_in,pi_2D=pi_2D)
                phi = self.sample_phi( currentx= phi_pre, r = self.r, p = self.p, K = self.K, I=self.I,H=H_pre, C=self.C, A=self.A, D=self.D, n=self.n_lambda,theta=self.theta_0,optimal_rate = self.optimal_rate,theta_variable = theta_variable,h_theta = h_theta_pre,gamma=self.gamma,sigma_phi_initial=self.sigma_phi_initial,sigma_phi_changes=self.sigma_phi_changes,changes_batch=self.changes_batch,var_calculation=self.var_calculation)
                pi = self.sample_pi_matrix(alpha=self.zeta, K=self.K, Z=Z_pre,pi_2D=pi_2D)
                G = self.sample_G(currentx = G_pre, Z = Z_pre, F = self.F, K = self.K , S = self.S , F_epsilon = self.F_epsilon , p_ais = self.p_ais, p_dis = self.p_ais,phi=phi_pre, C=self.C, A=self.A, D=self.D, n=self.n_lambda, I=self.I,theta=self.theta_0,optimal_rate=self.optimal_rate,theta_variable = theta_variable,h_theta = h_theta_pre,gamma=self.gamma,sigma_G_initial=self.sigma_G_initial,sigma_G_changes=self.sigma_G_changes,changes_batch=self.changes_batch,var_calculation=self.var_calculation)
                H = self.sample_H(G = G_pre,S = self.S,K = self.K)

            else:
                #print("n")
                #n = self.n_lambda
                n = self.sample_n(currentx=n_pre,n_lambda=self.n_lambda,S=self.S,optimal_rate=self.optimal_rate,H=H_pre,phi=phi_pre,D=self.D,I=self.I,sigma_n_initial=self.sigma_n_initial,sigma_n_changes=self.sigma_n_changes,changes_batch=self.changes_batch,var_calculation=self.var_calculation)
                #print("Z")
                Z,PZ = self.sample_Z_matrix_log(S=self.S, K=self.K, G=G_pre,
                                                              Z=Z_pre, pi=pi_pre,
                                                              F=self.F, F_epsilon=self.F_epsilon, burn_in=burn_in,pi_2D=pi_2D)
                #print("phi")
                #phi = simulated_data.phi
                phi = self.sample_phi(currentx=phi_pre,r=self.r, p=self.p, K=self.K, I=self.I,
                                                                           H=H_pre, C=self.C,
                                                                           A=self.A, D=self.D, n=n_pre, theta=self.theta_0,
                                                                           optimal_rate=self.optimal_rate,theta_variable = theta_variable,h_theta = h_theta_pre,gamma=self.gamma,sigma_phi_initial=self.sigma_phi_initial,sigma_phi_changes=self.sigma_phi_changes,changes_batch=self.changes_batch,var_calculation=self.var_calculation)
                #print("pi")
                pi = self.sample_pi_matrix(alpha=self.zeta, K=self.K, Z=Z_pre,pi_2D=pi_2D)
                #print("G")
                G = self.sample_G(currentx=G_pre,Z=Z_pre, F=self.F, K=self.K,
                                                                       S=self.S, F_epsilon=self.F_epsilon,
                                                                       p_ais=self.p_ais, p_dis=self.p_ais,
                                                                       phi=phi_pre, C=self.C,
                                                                       A=self.A, D=self.D, n=n_pre, I=self.I,
                                                                       theta=self.theta_0,
                                                                       optimal_rate=(self.optimal_rate/2),theta_variable = theta_variable,h_theta = h_theta_pre,gamma=self.gamma,sigma_G_initial=self.sigma_G_initial,sigma_G_changes=self.sigma_G_changes,changes_batch=changes_batch,var_calculation=var_calculation)
                #print("H")
                H = self.sample_H(G=G_pre, S=self.S, K=self.K)
            loglik = self.log_likelihood_model(A=self.A, D=self.D, H=H,
                                                                  phi=phi, C=self.C, n=n, I=self.I,
                                                                  theta=self.theta_0,theta_variable=theta_variable,h_theta=h_theta)

            if theta_variable == True:
                h_theta = self.sample_theta( self.gamma, H_pre, phi_pre, self.C,self.A,self.D)

            n_pre,H_pre,G_pre,pi_pre,phi_pre,Z_pre,h_theta_pre,loglik_pre = self.save_and_update_variables(n,H,G,pi,phi,Z,PZ,h_theta,loglik,iter,every_n_sample,current_batch)

            if ((iter+1) % batch)== 0:
                #convergence = (np.sum(np.abs(self.geweke_z(self.H, int(np.round(iter/every_n_sample)), first=0.1, last=0.5))<2)/(self.S*self.K))*100
                convergences = self.test_convergence(self.H,iter,every_n_sample)
                #print("convergences[0]: " + str(convergences[0]))
                print("Batch " + str(int(iter / batch)) + " finished with "+str(convergences)+ "% convergence for [0,10,20,30,40,50]% burn-in.")
                if(convergence_counter == 0):
                    convergence_counter=1
                    self.convergence_rate = convergences
                else:
                    self.convergence_rate = np.dstack((self.convergence_rate,convergences))
                if(iter>min_iter and np.any(convergences>99.5)):
                    print("converged early")
                    print(convergences)
                    converged = [idx for idx, element in enumerate(convergences) if convergences[idx]>99.5]
                    #burnin = int(round((converged[0]/10)*(iter/every_n_sample)))
                    max_iter=iter
                    self.max_iter = iter
                    break
                current_batch = current_batch+1
            length[iter] = time.time() - begin


        print("************* Inference started ************")
        self.time = (time.time() - first)
        samples_count = int(np.round((iter+1)/every_n_sample))
        batch_count = int(np.round((iter+1)/batch))
        batch_n = int(round(batch/every_n_sample))
        print("The number of samples was "+str(samples_count))


        # TODO: I think there shouldnt be any inferred variable ( n_inferred is wrong --> I changed it to self.n in the following)
        #loglik = self.log_likelihood_model(A=self.A, D=self.D, H=H_pre,phi=phi_pre, C=self.C, n=n, I=self.I,theta=self.theta_0,theta_variable=theta_variable,h_theta=h_theta_pre)
        #self.inferred_D = np.mean(self.D[burn_in:][:][:][:],axis=0)
        last_convergence = self.test_convergence_batches(self.H,iter,every_n_sample,batch_count,samples_count,batch_n)
        self.last_convergence = last_convergence
        converged_batch = [idx for idx, element in enumerate(last_convergence) if last_convergence[idx]==max(last_convergence)]
        print(last_convergence)
        print(last_convergence[converged_batch[0]])
        batch_count = int(round((iter+1)/batch))

        print("batch_count: ",str(batch_count))
        print("batch_n: ",str(batch_n))
        print("samples_count: ",str(samples_count))
        num_of_samples_after_burnin = samples_count-(converged_batch[0])*batch_n
        print("Number of samples after discarding the burn-in ",str(num_of_samples_after_burnin))

        if n_sampling == True:
            self.inferred_n = np.sum(self.n_sum[converged_batch[0]:(batch_count)], axis=0)/num_of_samples_after_burnin
        else:
            self.inferred_n = self.n
        self.inferred_H = np.sum(self.H_sum[converged_batch[0]:(batch_count)], axis=0)/num_of_samples_after_burnin
        self.inferred_G = np.sum(self.G_sum[converged_batch[0]:(batch_count)], axis=0)/num_of_samples_after_burnin
        #self.inferred_Z = np.multiply(np.sum(self.Z_sum[converged_batch[0]:(batch_count)], axis=0)/(samples_count-(converged_batch[0])*batch_n)<th,1)
        #self.inferred_P_Z = self.p_z/self.counter_p_z
        self.inferred_P_Z = np.sum(self.PZ_sum[converged_batch[0]:(batch_count)], axis=0)/num_of_samples_after_burnin
        self.inferred_Z = np.multiply((np.sum(self.PZ_sum[converged_batch[0]:(batch_count)], axis=0)/num_of_samples_after_burnin) <th,1)
        self.inferred_phi = np.sum(self.phi_sum[converged_batch[0]:(batch_count)], axis=0)/num_of_samples_after_burnin
        self.inferred_pi = np.sum(self.pi_sum[converged_batch[0]:(batch_count)] ,axis=0)/num_of_samples_after_burnin
        self.inferred_h_theta =  np.sum(self.h_theta_sum[converged_batch[0]:(batch_count)], axis=0)/num_of_samples_after_burnin

        if (simulated_data is None):
            self.save_in_txt_only_results(self.result_txt)
        else:
            self.save_in_txt(simulated_data, self.result_txt)

        pickle.dump(self, open(self.name, 'wb'))
        return self



    ## TODO: when we are not sampling theta, we do not need to initialize it, memory consumption
    # ????? Here we are not using D and A for generating initial values (Not Complete)
    def initialization(self,K,S,r,p,I,F,avarage_clone_in_spot,seed,max_iter,batch,D,A,C,F_epsilon,theta,n_lambda,n_sampling,F_fraction,theta_variable,gamma,pi_2D,every_n_sample):
        self.iter_current =  0
        self.counter_p_z = 0
        if n_sampling==True:
            sample = sim(K=K, S=S, r=r, p=p, I=I, F=F, D=D, A=A, C=C,
                                    avarage_clone_in_spot=avarage_clone_in_spot, random_seed=seed, F_epsilon=F_epsilon,
                                    n=n_lambda, p_c_binom=None, theta=theta, Z=None, n_lambda=n_lambda,F_fraction=False,theta_variable = theta_variable,gamma=gamma,pi_2D=pi_2D)
            n = sample.n
            if np.sum(sample.n==0)>0:
                raise Exception("initial n has some zero")
            if pi_2D is True:
                self.calculating_zeta_s( n_lambda, avarage_clone_in_spot, K)
        else:
            sample = sim(K=K, S=S, r=r, p=p, I=I, F=F, D=D, A=A, C=C,
                                    avarage_clone_in_spot=avarage_clone_in_spot, random_seed=seed, F_epsilon=F_epsilon,
                                    n=n_lambda, p_c_binom=None, theta=theta, Z=None, n_lambda=n_lambda,F_fraction=F_fraction,theta_variable = theta_variable,gamma=gamma,pi_2D=pi_2D)
            n = sample.n


        keep = int(np.round(max_iter/every_n_sample))
        #self.Z = np.zeros((keep, S, K), dtype=np.float64)
        self.H = np.zeros((keep, S, K), dtype=np.float64)
        self.G = np.zeros((keep, S, K), dtype=np.float64)
        self.n = np.zeros((keep, S), dtype=np.float64)
        self.phi = np.zeros((keep, I, K), dtype=np.float64)

        self.loglik = np.zeros((keep-1), dtype=np.float64)
        self.decision_matrix_phi = np.zeros( (I, K), dtype=np.float64)
        self.decision_matrix_G = np.zeros( (S, K), dtype=np.float64)
        self.decision_matrix_n = np.zeros( (S), dtype=np.float64)
        self.p_z = np.zeros( (S, K), dtype=np.float64)
        self.geweke_z_all = 1000*np.ones(int(max_iter/batch), dtype=np.float64)

        #self.Z[0] = copy.deepcopy(sample.Z)
        #self.H[0] = copy.deepcopy(sample.H)

        n= n_lambda
        Z = copy.deepcopy(sample.Z)
        adjust = np.tile(n_lambda, (I, 1))
        phi = np.tile(np.mean(D/(adjust.astype(float)),axis=1),(K,1)).transpose()#copy.deepcopy(sample.phi)
        phi[phi==0]=10000
        phi[phi==10000]=np.min(phi)
        
        pi = copy.deepcopy(sample.pi)
        H = copy.deepcopy(sample.H)
        G = copy.deepcopy(sample.G)
        if theta_variable==True:
            h_theta  = copy.deepcopy(sample.h_theta)
        else:
            h_theta = 0

        batch_count = int(round(self.max_iter/batch))
        self.PZ_sum = np.zeros( (batch_count,S, K), dtype=np.float64)
        self.n_sum =  np.zeros((batch_count,S))
        self.H_sum = np.zeros((batch_count,S,K))
        self.G_sum = np.zeros((batch_count,S,K))
        self.pi_sum = np.zeros((batch_count,S,K))
        self.phi_sum = np.zeros((batch_count,I,K))
        self.Z_sum = np.zeros((batch_count,S,K))
        self.h_theta_sum = np.zeros((batch_count,I,S))


        #self.n_sum[0] = copy.deepcopy(sample.n)
        #self.H_sum[0] = copy.deepcopy(sample.H)
        #self.G_sum[0] = copy.deepcopy(sample.G)
        #self.pi_sum[0] = copy.deepcopy(sample.pi)
        #self.phi_sum[0] = copy.deepcopy(sample.phi)
        #self.Z_sum[0] = copy.deepcopy(sample.Z)
        #self.h_theta_sum[0] = copy.deepcopy(sample.h_theta)

        return n,H,G,pi,phi,Z,h_theta



    ##############################################################################################################
    ############################################  utils  ################################################
    ##############################################################################################################

    def finding_cutoff(self, sample_1):
        self.Z_SEE = []
        for cutoff in np.arange(0, 1, 0.05):
            self.Z_SEE.append(np.sum(np.abs(sample_1.Z - (self.inferred_P_Z<cutoff)*1)) / (len(sample_1.Z)*len(sample_1.Z[0])))


    def ncr(self,n, r):
        r = min(r, n - r)
        r = int(r)
        n = int(n)
        numer = reduce(op.mul, range(n, n - r, -1), 1)
        denom = reduce(op.mul, range(1, r + 1), 1)
        return numer // denom  # or / in Python 2

    def truncate(self,f, n):
        '''Truncates/pads a float f to n decimal places without rounding'''
        s = '%.12f' % f
        i, p, d = s.partition('.')
        return '.'.join([i, (d + '0' * n)[:n]])


    def geweke_z(self,samples, N, first=0.1, last=0.5):
        m = int(np.floor(first * N))
        n = int(np.ceil(last * N))
        A = samples[0:m]
        B = samples[n:(N-1)]

        A_mean = np.mean(A,0)
        B_mean = np.mean(B,0)
        A_var = np.var(A,0)
        B_var = np.var(B,0)

        min_var = np.power(10, np.float64(-50))
        Z = (A_mean - B_mean) / np.sqrt(A_var + B_var + min_var)
        if np.any(np.isnan(Z)):
            raise Exception("Geweke is NaN")
        return Z

    def my_floor(self,a, precision=0):
        return np.round(a - 0.5 * 10 ** (-precision), precision)

    def test_convergence(self,H,iter,every_n_sample):
        test_count = 6
        steps = 0.1
        N = int(np.round(iter/every_n_sample))
        convergences = np.zeros(test_count, dtype=np.float64)
        for i in range(test_count):
            #print("i is: ",i)
            #print(int(round(i*steps*N)))
            #print(N-1)
            convergences[i] = np.round((np.sum(np.abs(self.geweke_z(H[int(round(i*steps*N)):N], N- int(round(i*steps*N)), first=0.1, last=0.5))<2)/(self.S*self.K))*100,2)
        return convergences
    def test_convergence_batches(self,H,iter,every_n_sample,batch_count,N,batch_n):
        convergences = np.zeros(batch_count, dtype=np.float64)
        for i in range(batch_count):
            #print("i is: ",i)
            #print(int(round(i*steps*N)))
            #print(N-1)
            convergences[i] = np.round((np.sum(np.abs(self.geweke_z(H[int(round(i*batch_n)):N], (N-int(round(i*batch_n))), first=0.1, last=0.5))<2)/(self.S*self.K))*100,2)
        return convergences
    def save_and_update_variables(self,n,H,G,pi,phi,Z,PZ,h_theta,loglik,iter,every_n_sample,current_batch):
        n_pre = copy.deepcopy(n)
        H_pre = copy.deepcopy(H)
        G_pre = copy.deepcopy(G)
        pi_pre = copy.deepcopy(pi)
        phi_pre = copy.deepcopy(phi)
        Z_pre = copy.deepcopy(Z)
        h_theta_pre = copy.deepcopy(h_theta)
        loglik_pre = copy.deepcopy(loglik)
        if(((iter+1) % every_n_sample)== 0):
            self.H[self.iter_current] = copy.deepcopy(H)
            #self.Z[self.iter_current] = copy.deepcopy(Z)
            self.G[self.iter_current] = copy.deepcopy(G)
            self.n[self.iter_current] = copy.deepcopy(n)
            self.phi[self.iter_current] = copy.deepcopy(phi)

            self.loglik[self.iter_current] = copy.deepcopy(loglik)
            self.n_sum[current_batch] = self.n_sum[current_batch]+copy.deepcopy(n)
            self.H_sum[current_batch] = self.H_sum[current_batch]+copy.deepcopy(H)
            self.G_sum[current_batch] = self.G_sum[current_batch]+copy.deepcopy(G)
            self.pi_sum[current_batch] = self.pi_sum[current_batch]+copy.deepcopy(pi)
            self.phi_sum[current_batch] = self.phi_sum[current_batch]+copy.deepcopy(phi)
            self.Z_sum[current_batch] = self.Z_sum[current_batch]+copy.deepcopy(Z)
            self.PZ_sum[current_batch] = self.PZ_sum[current_batch]+copy.deepcopy(PZ)
            self.h_theta_sum[current_batch] = self.h_theta_sum[current_batch]+copy.deepcopy(h_theta)
            self.iter_current = self.iter_current + 1

        #self.n_sum = self.n_sum+copy.deepcopy(n)
        #self.H_sum = self.H_sum+copy.deepcopy(H)
        #self.G_sum = self.G_sum+copy.deepcopy(G)
        #self.pi_sum = self.pi_sum+copy.deepcopy(pi)
        #self.phi_sum = self.phi_sum+copy.deepcopy(phi)
        #self.Z_sum = self.Z_sum+copy.deepcopy(Z)
        #self.h_theta_sum = self.h_theta_sum+copy.deepcopy(h_theta)

        return n_pre,H_pre,G_pre,pi_pre,phi_pre,Z_pre,h_theta_pre,loglik_pre


    def calculating_zeta_s(self,n_lambda,avarage_clone_in_spot,K):
        density = (n_lambda) / (np.max(n_lambda))
        self.avarage_clone_in_spot_s = np.maximum(1,np.minimum(n_lambda,density * avarage_clone_in_spot))
        self.zeta = (self.avarage_clone_in_spot_s*K)/(K-self.avarage_clone_in_spot_s)

    def trunc_norm_sampling_vector(self,mu, sigma):
        n = len(mu)
        U = np.random.mtrand._rand.uniform(size=n)
        y = mu + sigma * sc.ndtri(U + sc.ndtr(-mu / sigma) * (1 - U))
        return y

    def trunc_norm_sampling_matrix(self,mu, sigma):
        n_col = mu.shape[1]
        n_row = mu.shape[0]

        U = np.random.mtrand._rand.uniform(size=(n_row, n_col))
        y = mu + sigma * sc.ndtri(U + sc.ndtr(-mu / sigma) * (1 - U))
        y[np.isinf(y)] = mu[np.isinf(y)]
        y[y<0] = mu[y<0]
        return y

    # TODO: take care of added epsilon
    def calculate_logp_ais(self, A, D, H, phi, C, theta):
        binom_p_nom = np.matmul(H, np.matrix.transpose(np.array(phi) * np.array(C)))
        binom_p_denom = np.matmul(H, np.matrix.transpose(np.array(phi)))
        binom_p = (binom_p_nom * theta) / binom_p_denom
        if np.sum(binom_p > 1) + np.sum(binom_p == 0) > 0:
            raise Exception("Why binom_p in G_target is more than one, or is zero")
        binom_p[np.isnan(binom_p)] = 0.0001
        p_ais = scipy.stats.binom.logpmf(k=A, n=D, p=np.matrix.transpose(binom_p))
        return p_ais

    def calculate_logp_dis(self, D, H, phi, I, n):
        binom_p_denom = np.matmul(H, np.matrix.transpose(np.array(phi)))
        p_dis = scipy.stats.poisson.logpmf(D, np.tile(n, (I, 1)) * np.matrix.transpose(binom_p_denom))
        return p_dis

    ##############################################################################################################
    ############################################  sampling theta  ################################################
    ##############################################################################################################

    def sample_theta(self, gamma, H, phi, C,A,D):
        alpha = np.matmul(gamma*phi*C,np.transpose(H))
        beta = np.matmul(gamma*phi*(1-C),np.transpose(H))
        alpha[alpha==0] = 0.01
        alpha[alpha == 1] = 0.99
        beta[beta==0] = 0.01
        beta[beta == 1] = 0.99
        h_theta = np.random.beta(alpha + A, beta + D - A)
        #if np.sum(h_theta[h_theta == 0])>0:
        #    raise Exception('There are some 0\'s in thetas which are generated in the simulation')
        counter = 0
        while np.any(h_theta == 0) or np.any(h_theta == 1):
            counter = counter+1
            print(counter)
            h_theta[h_theta == 0] = np.random.beta(alpha[h_theta == 0], beta[h_theta == 0])
            h_theta[h_theta == 1] = np.random.beta(alpha[h_theta == 1], beta[h_theta == 1])
        #h_theta = self.my_floor(np.random.beta(alpha+A,beta+D-A),5)
        #h_theta[h_theta==0] = 0.00001
        return h_theta

    ##############################################################################################################
    ################################################  sampling Z  ################################################
    ##############################################################################################################

    def sample_Z_matrix_log(self, S, K, G, Z,pi,F,F_epsilon,burn_in,pi_2D):
        if pi_2D is True:
            p_0 = np.log(1 - pi) + scipy.stats.gamma.logpdf(x=G, a=F_epsilon[:, 0],scale=F_epsilon[:, 1])
            p_1 = np.log(pi) + scipy.stats.gamma.logpdf(x=G, a=F[:, 0], scale=F[:, 1])
        else:
            p_0 = np.log(np.tile((1-pi), (S, 1))) + scipy.stats.gamma.logpdf(x=G, a=F_epsilon[:, 0], scale=F_epsilon[:, 1])
            p_1 = np.log(np.tile(pi, (S, 1))) + scipy.stats.gamma.logpdf(x=G, a=F[:, 0], scale=F[:, 1])

        p_frac = np.exp(p_0-p_1)
        if np.sum(np.isnan(p_frac))>0 :
            p_frac[np.isnan(p_frac)] = 1
            #raise Exception('p_frac is nan')
        if np.sum(np.isinf(p_frac)) > 0:
            p_frac[np.isinf(p_frac)] = 1
            #raise Exception('p_frac is inf')
        if np.sum(p_frac==0)>0:
            print("Z migth have problems")
            print(np.sum(p_frac==0))
        decision_matrix = np.random.uniform(low=0.0, high=1.0, size=(S,K)) < (p_frac / (1 + p_frac))
        Z[decision_matrix==1]=0
        Z[decision_matrix==0]=1
        PZ = (p_frac/(1+p_frac))
        if self.iter >= burn_in:
            self.p_z += PZ
            self.counter_p_z += 1
        return Z,PZ


    ##############################################################################################################
    ################################################  sampling N  ################################################
    ##############################################################################################################
    def next_sigma(self,sigma_current,optimal_acceptance,last_acceptance,last_iter):
        theta = np.log(sigma_current)

        cp_accepted=1/(last_iter*optimal_acceptance)
        theta_next_accepted = theta + cp_accepted
        cp_rejected=1/((1-optimal_acceptance)*last_iter)
        theta_next_rejected = theta - cp_rejected

        theta_next = theta_next_rejected
        theta_next[last_acceptance] = theta_next_accepted[last_acceptance]
        return(np.exp(theta_next))


    def sample_n(self,currentx,n_lambda,S,optimal_rate,H,phi,D,I,sigma_n_initial,sigma_n_changes,changes_batch,var_calculation):
        #check and update sigma_n

        if self.iter%changes_batch==0:
            #if self.iter%(5*changes_batch)==changes_batch and var_calculation:
            #    self.sigma_n=np.sqrt(np.var(self.n[0:(self.iter-1)],axis=0))
            #else:
                #self.sigma_n= self.sigma_n+((self.decision_matrix_n/self.iter)-optimal_rate)*sigma_n_changes
            self.sigma_n = self.sigma_n + ((self.decision_matrix_n / self.iter) - optimal_rate) * sigma_n_changes * self.sigma_n
        elif self.iter==1:
            self.sigma_n = sigma_n_initial
        '''if self.iter == 1:
            self.sigma_n = np.tile(sigma_n_initial, S)
        else:
            self.sigma_n = self.next_sigma(self.sigma_n,optimal_rate,self.last_acceptance_n,self.iter)'''

        if np.sum(self.sigma_n<=0)>0:
            self.sigma_n[self.sigma_n <= 0] = sigma_n_initial
            print("Number of self.sigma_n<=0 is "+str(np.sum(self.sigma_n <= 0)))


        proposedx = np.round(self.trunc_norm_sampling_vector(currentx, self.sigma_n))
        #while (True):
        #    eless0 = np.sum(proposedx <= 0.5)
        #    if (eless0==0):
        #        break
        #    else:
        #        proposedx[proposedx <= 0.5]  = np.round(self.trunc_norm_sampling_vector(currentx[proposedx <= 0.5] , self.sigma_n))
                #proposedx[proposedx <= 0.5] = currentx[proposedx <= 0.5] + np.random.normal(mu,self.sigma_n[proposedx <= 0.5],size=(eless0))

        current_p = self.n_target_matrix(x=currentx, H=H, phi=phi,D=D,I=I,n_lambda=n_lambda)
        proposed_p = self.n_target_matrix(x=proposedx, H=H, phi=phi,D=D,I=I,n_lambda=n_lambda)

        n_loglik = ((proposed_p + norm.logcdf(currentx,loc=0,scale=self.sigma_n)) - (current_p + norm.logcdf(proposedx,loc=0,scale=self.sigma_n)))
        decision_matrix_n = (np.log(np.random.uniform(low=0.0, high=1.0, size=(S))) < n_loglik)*1

        n = np.zeros((S), dtype=np.float64)
        n[decision_matrix_n==1] = proposedx[decision_matrix_n==1] # accept move with probabily min(1,A)
        n[decision_matrix_n==0] = currentx[decision_matrix_n==0] # otherwise "reject" move, and stay where we are
        self.last_acceptance_n = decision_matrix_n
        self.decision_matrix_n += decision_matrix_n
        if np.sum(n<=0)>0:
            print("We have sum n<=0")
        if np.any(np.isnan(n)):
            print("n has NaN")
        return np.round(n)

    def n_target_matrix(self,x, H, phi,D,I,n_lambda):
        n = np.round(x)
        p_dis = self.calculate_logp_dis(D,H,phi,I,n)
        p_n = scipy.stats.poisson.logpmf(n,n_lambda) + np.sum(p_dis, axis=0)
        return p_n

    ##############################################################################################################
    ################################################  sampling pi  ################################################
    ##############################################################################################################


    # the beta bernouli should be checked again
    def sample_pi_matrix(self,alpha,K,Z,pi_2D):
        if pi_2D is True:
            pi = np.random.beta(np.transpose(np.array([alpha / K] * K)) + Z, 2 - Z)
        else:
            pi = np.random.beta((alpha / K) + sum(Z), 1 + sum(1 - Z))
        if np.any(np.isinf(pi)) or np.any(np.isnan(pi)):
            print(pi[np.isnan(pi)])
            print(pi[np.isinf(pi)])
            raise Exception('pi has problem')
        return pi



    ##############################################################################################################
    ################################################  sampling H & G  ############################################
    ##############################################################################################################


    # we expect that all elements of H is more than 0
    def sample_H(self,G,S,K):
        H = G / np.transpose(np.tile(np.sum(G, axis=1), (K, 1)))
        return H

    # here, first we generate some random numbers using a normal distribution. The mean is 0 and variance is being learned during
    # the sampling steps. Then we accept or reject the sample based on the probability of the smaple.
    def sample_G(self,currentx,Z,F,K,S,F_epsilon,p_ais,p_dis, phi,C,A,D,n,I,theta,optimal_rate,theta_variable,h_theta,gamma,sigma_G_initial,sigma_G_changes,changes_batch,var_calculation):

        if self.iter%changes_batch==0:
            #if self.iter%(5*changes_batch)==changes_batch:
            #    self.sigma_G=np.sqrt(np.var(self.G[0:(self.iter-1)],axis=0))
            #else:
                #self.sigma_G= self.sigma_G+((self.decision_matrix_G/self.iter)-optimal_rate)*sigma_G_changes
            self.sigma_G = self.sigma_G + ((self.decision_matrix_G / self.iter) - optimal_rate) * sigma_G_changes * self.sigma_G
        elif self.iter==1:
            self.sigma_G = sigma_G_initial


        '''if self.iter == 1:
            self.sigma_G = np.tile(sigma_G_initial, (S, K))
        else:
            self.sigma_G = self.next_sigma(self.sigma_G,optimal_rate,self.last_acceptance_G,self.iter)'''

        if np.sum(self.sigma_G<=0)>0:
            print("Problem: sigma_G <=0 ")
            self.sigma_G[self.sigma_G <= 0] = sigma_G_initial

        proposedx =  self.trunc_norm_sampling_matrix(currentx, self.sigma_G)
        while (True):
            eless0 = np.sum(proposedx <= 0)
            if (eless0==0):
                break
            else:
                proposedx[proposedx <= 0] = self.trunc_norm_sampling_matrix(currentx[proposedx <= 0], self.sigma_G[proposedx <= 0])

        if theta_variable == False:
            current_p = self.G_target(x=currentx,Z = Z,F=F,F_epsilon = F_epsilon,K=K, phi=phi, C=C, A=A, D=D, n=n, I=I,theta=theta)
            proposed_p = self.G_target(x=proposedx,Z = Z,F=F,F_epsilon = F_epsilon,K=K,phi=phi, C=C, A=A, D=D, n=n, I=I,theta=theta)
        else:
            current_p = self.G_target_theta(x=currentx, Z=Z, F=F, F_epsilon=F_epsilon, K=K, phi=phi, C=C,
                                                        D=D, n=n, I=I ,gamma=gamma,h_theta=h_theta)
            proposed_p = self.G_target_theta(x=proposedx, Z=Z, F=F, F_epsilon=F_epsilon, K=K, phi=phi, C=C,
                                                         D=D, n=n, I=I,gamma=gamma,h_theta=h_theta)
        log_g_lik = (proposed_p + norm.logcdf(currentx,loc=0,scale=self.sigma_G)) - (current_p + norm.logcdf(proposedx,loc=0,scale=self.sigma_G))
        g_loglik = (log_g_lik.sum(axis=1))
        if np.any(np.isinf(g_loglik)):
            print("Problem: there are inf values in G: "+str(np.sum(np.isinf(g_loglik))))

        decision_matrix_G_1 = np.log(np.random.uniform(low=0.0, high=1.0, size=(S))) < g_loglik
        decision_matrix_G = np.transpose(np.tile(decision_matrix_G_1, (K, 1)))

        G = np.zeros((S, K), dtype=np.float64)
        G[decision_matrix_G==1] = proposedx[decision_matrix_G==1] # accept move with probabily min(1,A)
        G[decision_matrix_G==0] = currentx[decision_matrix_G==0]  # otherwise "reject" move, and stay where we are
        self.decision_matrix_G += decision_matrix_G
        self.last_acceptance_G = decision_matrix_G
        if np.sum(G==0):
            print("Problem: number of zero Gs: "+str(np.sum(G==0)))
            raise Exception("G got zero in sample_G_matrix")
        return G



    def G_target(self,x,Z,F,F_epsilon,K,phi,C,A,D,n,I,theta):
        H = x/np.transpose(np.tile(np.sum(x,axis=1),(K,1)))
        p_ais = self.calculate_logp_ais(A,D,H,phi,C,theta)
        p_dis = self.calculate_logp_dis( D, H, phi, I, n)
        p_gsk = scipy.stats.gamma.logpdf(a = np.power(F[:,0],Z) * np.power(F_epsilon[:,0],1-Z), x = x,scale=np.power(F[:,1],Z) * np.power(F_epsilon[:,1],1-Z))
        p_g = p_gsk + np.transpose(np.tile(np.sum(p_ais, axis=0)+np.sum(p_dis, axis=0), (K, 1)))
        return p_g

    #should be written
    def G_target_theta(self,x,Z,F,F_epsilon,K,phi,C,D,n,I,gamma,h_theta):
        H = x/np.transpose(np.tile(np.sum(x,axis=1),(K,1)))

        binom_p_denom = np.matmul(H, np.matrix.transpose(np.array(phi)))
        p_dis = scipy.stats.poisson.logpmf(D, np.tile(n, (I, 1)) * np.matrix.transpose(binom_p_denom))

        alpha = np.matmul(gamma*phi*C,np.transpose(H))
        beta = np.matmul(gamma*phi*(1-C),np.transpose(H))
        #TO DO: check the values, if they are reseanable and the place of alpha and beta are correct
        p_theta = scipy.stats.beta.logpdf(h_theta,alpha,beta)

        p_gsk = scipy.stats.gamma.logpdf(a = np.power(F[:,0],Z) * np.power(F_epsilon[:,0],1-Z), x = x,scale=np.power(F[:,1],Z) * np.power(F_epsilon[:,1],1-Z))
        p_g = p_gsk + np.transpose(np.tile(np.sum(p_theta, axis=0)+np.sum(p_dis, axis=0), (K, 1)))
        if (np.sum(np.isinf(p_g))>0):
            print(h_theta[np.isinf(p_theta)])
            print(alpha[np.isinf(p_theta)])
            print(beta[np.isinf(p_theta)])
            raise Exception('p_g is inf')
        return p_g


    ##############################################################################################################
    ################################################  sampling PHI  ##############################################
    ##############################################################################################################

    def sample_phi(self,currentx, r, p, K,I,H,C,A,D,n,theta,optimal_rate,theta_variable,h_theta,gamma,sigma_phi_initial,sigma_phi_changes,changes_batch,var_calculation):

        # this while is not efficient
        if self.iter%changes_batch==0:
            #self.sigma_phi= self.sigma_phi+((self.decision_matrix_phi/self.iter)-optimal_rate)*sigma_phi_changes
            #if self.iter%(5*changes_batch)==changes_batch and var_calculation:
            #    self.sigma_phi=np.sqrt(np.var(self.phi[0:(self.iter-1)],axis=0))
            #else:
            self.sigma_phi = self.sigma_phi + ((self.decision_matrix_phi / self.iter) - optimal_rate) * sigma_phi_changes *self.sigma_phi
        elif self.iter==1:
            self.sigma_phi = sigma_phi_initial

        '''if self.iter == 1:
            self.sigma_phi = np.tile(sigma_phi_initial, (I, K))
        else:
            self.sigma_phi = self.next_sigma(self.sigma_phi,optimal_rate,self.last_acceptance_phi,self.iter)'''

        if np.sum(self.sigma_phi<=0)>0:
            print("Problem: sigma_phi<=0 ")
            #print(self.sigma_phi)
            #print(sigma_phi_initial)
            np.putmask(self.sigma_phi,self.sigma_phi <= 0,sigma_phi_initial)
            #self.sigma_phi[self.sigma_phi <= 0] = sigma_phi_initial

        proposedx =  self.trunc_norm_sampling_matrix(currentx, self.sigma_phi)
        while (True):
            eless0 = np.sum(proposedx <= 0)
            if (eless0==0):
                break
            else:
                proposedx[proposedx <= 0] = self.trunc_norm_sampling_matrix(currentx[proposedx <= 0], self.sigma_phi[proposedx <= 0])

        if theta_variable==False:
            current_p = self.phi_target(phi=currentx,r=r,p=p,K=K, H=H, C=C, A=A, D=D, n=n, I=I,theta=theta)
            proposed_p = self.phi_target(phi=proposedx,r=r,p=p,K=K, H=H, C=C, A=A, D=D, n=n, I=I,theta=theta)
        else:
            current_p = self.phi_target_matrix_recal_log_theta(phi=currentx,r=r,p=p,K=K, H=H, C=C, D=D, n=n, I=I,h_theta=h_theta,gamma=gamma)
            proposed_p = self.phi_target_matrix_recal_log_theta(phi=proposedx,r=r,p=p,K=K, H=H, C=C, D=D, n=n, I=I,h_theta=h_theta,gamma=gamma)

        phi_loglik = ((proposed_p + norm.logcdf(currentx,loc=0,scale=self.sigma_phi)) - (current_p+ norm.logcdf(proposedx,loc=0,scale=self.sigma_phi)))
        if np.any(np.isinf(phi_loglik)):
            print("number of inf in phi_lik: "+str(np.sum(np.isinf(phi_loglik))))
        decision_matrix_phi = (np.log(np.random.uniform(low=0.0, high=1.0, size=(I,K))) < phi_loglik)*1

        phi = np.zeros((I, K), dtype=np.float64)
        phi[decision_matrix_phi==1] = proposedx[decision_matrix_phi==1] # accept move with probabily min(1,A)
        phi[decision_matrix_phi==0] = currentx[decision_matrix_phi==0]  # otherwise "reject" move, and stay where we are
        self.decision_matrix_phi += decision_matrix_phi
        self.last_acceptance_phi = decision_matrix_phi
        if np.sum(phi==0):
            print(np.sum(phi==0))
            print(self.iter)
            raise Exception("phi got zero in sample_phi_matrix")
        return phi

    def phi_target(self, phi, r, p, K, H, C, A, D, n, I,theta):

        p_ais = self.calculate_logp_ais(A,D,H,phi,C,theta)
        p_dis = self.calculate_logp_dis( D, H, phi, I, n)
        p_phi_ik = scipy.stats.gamma.logpdf(x=phi, a=r, scale=p)
        if np.any(p_phi_ik == 0):
            raise Exception("gamma of phi ik is zero in phi target matrix")
        p_phi = p_phi_ik + np.transpose(np.tile(np.sum(p_ais, axis=1) + np.sum(p_dis, axis=1), (K, 1)))
        return p_phi

    #should be written
    def phi_target_matrix_recal_log_theta(self, phi, r, p, K, H, C, D, n, I,h_theta,gamma):

        binom_p_denom = np.matmul(H, np.matrix.transpose(np.array(phi)))

        p_dis = scipy.stats.poisson.logpmf(D, np.tile(n, (I, 1)) * np.matrix.transpose(binom_p_denom))

        alpha = np.matmul(gamma*phi*C,np.transpose(H))
        beta = np.matmul(gamma*phi*(1-C),np.transpose(H))
        #TO DO: check the values, if they are reseanable and the place of alpha and beta are correct
        p_theta = scipy.stats.beta.logpdf(h_theta,alpha,beta)

        p_phi_ik = scipy.stats.gamma.logpdf(x=phi, a=r, scale=p)
        if np.sum(p_phi_ik == 0) > 0:
            raise Exception("gamma of phi ik is zero in phi target matrix")
        p_phi = p_phi_ik + np.transpose(np.tile(np.sum(p_theta, axis=1) + np.sum(p_dis, axis=1), (K, 1)))
        return p_phi



    def log_likelihood_model(self,A,D,H,phi,C,n,I,theta,theta_variable,h_theta):
        binom_p_nom = np.matmul(H, np.matrix.transpose(np.array(phi) * np.array(C)))
        binom_p_denom = np.matmul(H, np.matrix.transpose(np.array(phi)))
        binom_p = (binom_p_nom*theta)/binom_p_denom
        if np.sum(binom_p > 1)+ np.sum(binom_p == 0) > 0 :
            raise Exception("Why binom_p in G_target is more than one, or is zero")
        if theta_variable==False:
            self.p_ais = scipy.stats.binom.logpmf(k=A, n=D, p=np.matrix.transpose(binom_p))
        else:
            self.p_ais = scipy.stats.binom.logpmf(k=A, n=D, p=h_theta)
        #try:
        #    print(I)
        #    print(n)
        #    print(binom_p_denom)
        self.p_dis =  scipy.stats.poisson.logpmf(D,np.tile(n,(I,1))*np.matrix.transpose(binom_p_denom))
        #except Exception:
        #    print(" problem in likelihood and d")

        return np.sum(self.p_dis)+np.sum(self.p_ais)

    def save_in_txt_only_results(self,name):
        file1 = open(name, "w")
        file1.write("Time:"+str(self.time)+'\n')

        file1.write("\n\nInferred H:\n")
        np.savetxt(file1, self.inferred_H, fmt='%.3f')

        file1.write("\n\nInferred G:\n")
        np.savetxt(file1, self.inferred_G, fmt='%.2f')

        file1.write("\n\nInferred Z:\n")
        np.savetxt(file1, self.inferred_Z, fmt='%.2f')


        file1.write("\n\nInferred HZ:\n")
        inferred_HZ = self.inferred_H*self.inferred_Z
        inferred_HZ = inferred_HZ/np.transpose(np.tile(inferred_HZ.sum(axis=1),(len(self.inferred_H[0]),1)))
        np.savetxt(file1, inferred_HZ, fmt='%.3f')


        file1.write("\nInferred Probability of Z:\n")
        np.savetxt(file1, self.inferred_P_Z, fmt='%.6f')

        file1.write("\n\nInferred phi:\n")
        np.savetxt(file1, self.inferred_phi, fmt='%.6f')

        file1.write("\n\nacceptance rate, G:\n")
        np.savetxt(file1, self.decision_matrix_G / self.max_iter, fmt='%.3f')
        file1.write("\nacceptance rate, phi:\n")
        np.savetxt(file1, self.decision_matrix_phi / self.max_iter, fmt='%.3f')
        file1.write("\nacceptance rate, n:\n")
        np.savetxt(file1, self.decision_matrix_n / self.max_iter, fmt='%.3f')



        file1.write("\n\nsigma phi: \n")
        np.savetxt(file1, self.sigma_phi, fmt='%.6f')
        file1.write("\nsigma G:\n")
        np.savetxt(file1, self.sigma_G, fmt='%.6f')
        file1.write("\nsigma n:\n")
        np.savetxt(file1, self.sigma_n, fmt='%.6f')

        file1.write("\n\nInferred pi:\n")
        np.savetxt(file1, self.inferred_pi, fmt='%.2f')


        if self.theta_variable == True:
            file1.write("\n\nInferred H_theta:\n")
            np.savetxt(file1, self.inferred_h_theta, fmt='%.2f')


        file1.write("\n\nInferred n:\n")
        np.savetxt(file1, self.inferred_n, fmt='%.2f')

        file1.close()

    def save_in_txt(self,sample_1,name):
        file1 = open(name, "w")
        file1.write("Time:"+str(self.time)+'\n')

        file1.write("\nStandard Error of the Estimate H:\n")
        #self.H_SEE = np.sqrt(np.sum(np.power((sample_1.H - self.inferred_H),2))/(len(sample_1.H)*len(sample_1.H[0])))
        self.H_SEE = np.mean(np.abs(sample_1.H - self.inferred_H))
        file1.write(str(self.H_SEE))


        file1.write("\nStandard Error of the Estimate HZ:\n")
        simulated_HZ = sample_1.H*sample_1.Z
        simulated_HZ = simulated_HZ/np.transpose(np.tile(simulated_HZ.sum(axis=1),(len(sample_1.H[0]),1)))
        inferred_HZ = self.inferred_H*self.inferred_Z
        inferred_HZ = inferred_HZ/np.transpose(np.tile(inferred_HZ.sum(axis=1),(len(sample_1.H[0]),1)))

        #self.HZ_SEE = np.sqrt(np.sum(np.power((simulated_HZ-inferred_HZ),2))/(len(sample_1.H)*len(sample_1.H[0])))
        self.HZ_SEE = np.mean(np.abs(simulated_HZ - inferred_HZ))
        file1.write(str(self.HZ_SEE))


        file1.write("\nStandard Error of the Estimate phi:\n")
        #self.phi_SEE = np.sqrt(np.sum(np.power((sample_1.phi - self.inferred_phi),2))/(len(sample_1.phi)*len(sample_1.phi[0])))
        self.phi_SEE = np.mean(np.abs(sample_1.phi - self.inferred_phi))
        file1.write(str(self.phi_SEE))

        file1.write("\nStandard Error of the Estimate n:\n")
        try:
            #self.n_SEE = np.sqrt(np.sum(np.power((sample_1.n - self.inferred_n),2))/(len(sample_1.n)))
            self.n_SEE = np.mean(np.abs(sample_1.n - self.inferred_n))

        except Exception:
            print("n has problem in calculating error")
        file1.write(str(self.n_SEE))

        file1.write("\nStandard Error of the Estimate pi:\n")
        #self.pi_SEE = np.sqrt(np.sum(np.power((sample_1.pi - self.inferred_pi),2))/(len(sample_1.pi)))
        self.pi_SEE = np.mean(np.abs(sample_1.pi - self.inferred_pi))
        file1.write(str(self.pi_SEE))



        file1.write("\nStandard Error of the Estimate P_Z:\n")
        #self.PZ_SEE = np.sqrt(np.sum(np.power((sample_1.Z - (1-self.inferred_P_Z)),2))/(len(sample_1.Z)*len(sample_1.Z[0])))
        self.PZ_SEE = np.mean(np.abs(sample_1.Z - (1 - self.inferred_P_Z)))
        file1.write(str(self.PZ_SEE))

        if self.theta_variable == True:
            file1.write("\nStandard Error of the Estimate H_theta:\n")
            try:
                #self.h_theta_SEE = np.sqrt(np.sum(np.power((sample_1.h_theta - self.inferred_h_theta),2))/(len(sample_1.h_theta)*len(sample_1.h_theta[0])))
                self.h_theta_SEE = np.mean(np.abs(sample_1.h_theta - self.inferred_h_theta))

            except Exception:
                print("writing result problem")
            file1.write(str(self.h_theta_SEE))
        else:
            self.h_theta_SEE = 0


        self.finding_cutoff(sample_1)
        file1.write("\nPercent of wrong estimaiton of Z:\n")
        file1.write(str(np.arange(0, 1, 0.05)))
        file1.write(str(self.Z_SEE))

        file1.write("\n\nH:\n")
        np.savetxt(file1, sample_1.H, fmt='%.2f')
        file1.write("\nInferred H:\n")
        np.savetxt(file1, self.inferred_H, fmt='%.2f')

        file1.write("\n\nZ:\n")
        np.savetxt(file1, sample_1.Z, fmt='%.2f')
        file1.write("\nInferred Z:\n")
        np.savetxt(file1, self.inferred_Z, fmt='%.2f')

        file1.write("\n\nHZ:\n")
        np.savetxt(file1, simulated_HZ, fmt='%.2f')
        file1.write("\nInferred HZ:\n")
        np.savetxt(file1, inferred_HZ, fmt='%.2f')

        file1.write("\nInferred Probability of Z:\n")
        np.savetxt(file1, self.inferred_P_Z, fmt='%.2f')

        file1.write("\n\nphi:\n")
        np.savetxt(file1, sample_1.phi, fmt='%.2f')
        file1.write("\nInferred phi:\n")
        np.savetxt(file1, self.inferred_phi, fmt='%.2f')

        file1.write("\n\nacceptance rate, G:\n")
        np.savetxt(file1, self.decision_matrix_G / self.max_iter, fmt='%.2f')
        file1.write("\nacceptance rate, phi:\n")
        np.savetxt(file1, self.decision_matrix_phi / self.max_iter, fmt='%.2f')

        file1.write("\n\nsigma phi: \n")
        np.savetxt(file1, self.sigma_phi, fmt='%.2f')
        file1.write("\nsigma G:\n")
        np.savetxt(file1, self.sigma_G, fmt='%.2f')

        file1.write("\n\npi:\n")
        np.savetxt(file1, sample_1.pi, fmt='%.2f')
        file1.write("\nInferred pi:\n")
        np.savetxt(file1, self.inferred_pi, fmt='%.2f')

        file1.write("\n\nA:\n")
        np.savetxt(file1, sample_1.A, fmt='%.2f')
        file1.write("\nD:\n")
        np.savetxt(file1, sample_1.D, fmt='%.2f')

        if self.theta_variable == True:
            file1.write("\n\nInferred H_theta:\n")
            np.savetxt(file1, self.inferred_h_theta, fmt='%.2f')
            file1.write("\n\n H_theta:\n")
            np.savetxt(file1, sample_1.h_theta, fmt='%.2f')


        file1.write("\n\nn:\n")
        np.savetxt(file1, sample_1.n, fmt='%.2f')
        file1.write("\nInferred n:\n")
        np.savetxt(file1, self.inferred_n, fmt='%.2f')

        file1.close()

