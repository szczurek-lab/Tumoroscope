import numpy as np

class simulation:
    def __init__(self,K,S,r,p,I,F,D,A,C,avarage_clone_in_spot,random_seed,F_epsilon,n,p_c_binom,theta,Z,n_lambda,F_fraction,theta_variable,gamma,pi_2D):

        self.K=K
        self.S=S
        # S.alpha/(1+(alpha/K)) = 2*S
        self.r=r
        self.p=p
        self.I=I
        self.avarage_clone_in_spot=avarage_clone_in_spot#2
        self.alpha=(avarage_clone_in_spot*K)/(K-avarage_clone_in_spot)
        self.F_epsilon = F_epsilon
        self.beta = F_epsilon[0][1]
        self.p_c_binom = p_c_binom
        self.theta = theta
        self.gamma  = gamma
        self.theta_variable = theta_variable
        np.random.seed(random_seed)#(4)
        # we can define n for the simulated data or let it be generated usinig it's distribution
        if n is None:
            self.generate_n(n_lambda)
        else:
            self.n = n

        self.calculating_alpha_s( n_lambda, avarage_clone_in_spot, K)

        # we can give F as input for the simulated data or let it be random
        if F_fraction is True:
            self.generate_F(self.K,self.beta)
        else:
            self.F = np.array(F)

        # we can define C for the simulated data or let it be generated usinig it's distribution
        if C is None:
            self.generate_C(self.I, self.K,self.p_c_binom)
        else:
            self.C = C

        if pi_2D is True:
            self.generate_pi_2D(self.zeta_s, self.K)
        else:
            self.generate_pi(self.alpha, self.K)

        self.generate_phi(self.r,self.p,self.I,self.K)

        if Z is None:
            self.generate_Z(self.K, self.pi,pi_2D)
        else:
            self.Z = Z

        self.generate_G(self.F,self.S,self.K,self.F_epsilon,self.Z)
        self.generate_H( self.G,self.S,self.K,self.Z)


        # We can give D as input or we can generate it by it's distribution
        if D is None:
            self.generate_D(I = self.I, S = self.S, K = self.K, phi = self.phi, H = self.H,n = self.n)
        else:
            self.D = D

        # The model can have theta both variable and constant, also in each case,
        # we can give Alternated read counts as input or we can generate them using
        # the distribution which we assumed.
        #print("generate A")
        #print(theta_variable)
        if theta_variable == False:
            if A is None:

                self.generate_A(self.I, self.S, self.H, self.C, self.D,self.phi,self.theta)
            else:
                self.A = A
        else:
            if A is None:
                self.generate_theta(self.gamma, self.phi,self.C,self.H)
                #print("theta finished")
                self.generate_A_theta(self.h_theta,self.D)
            else:
                #print("A is not non")
                self.A = A
                self.calculate_theta(self.I,self.S,self.gamma,self.phi,self.C,self.H)
        print("simulation finished")


    def calculating_alpha_s(self,n_lambda,avarage_clone_in_spot,K):
        density = (n_lambda) / (np.max(n_lambda))
        self.avarage_clone_in_spot_s = np.maximum(1,np.minimum(n_lambda,density * avarage_clone_in_spot))
        self.zeta_s = (self.avarage_clone_in_spot_s*K)/(K-self.avarage_clone_in_spot_s)

    def dirichlet_sample(self,alphas):
        """
        Generate samples from an array of alpha distributions.
        """
        r = np.random.standard_gamma(alphas)
        return np.divide(r , r.sum(-1, keepdims=True))

    def my_floor(self,a, precision=0):
        return np.round(a - 0.5 * 10 ** (-precision), precision)

    # there is a normal clone, which does not have any mutaiton,
    # if in C, we only have normal clone in one spot, then
    # we would get very small alpha and then nearly 0 h_theta
    # actually beta could generate 0 values, but beta.pdf return inf if theta s 0
    def generate_theta(self,gamma, phi,C,H):
        self.theta_alpha = np.round(100*np.matmul(gamma*phi*C,np.transpose(H)),2)
        self.theta_beta = np.round(100*np.matmul(gamma*phi*(1-C),np.transpose(H)))
        self.theta_alpha[self.theta_alpha==0] = 0.01
        self.theta_alpha[self.theta_alpha == 1] = 0.99
        self.theta_beta[self.theta_beta==0] = 0.01
        self.theta_beta[self.theta_beta == 1] = 0.99
        self.h_theta = np.random.beta(self.theta_alpha, self.theta_beta)
        #alphas = np.transpose(np.array([self.theta_alpha.flatten(),self.theta_beta.flatten()]))
        #self.dirichlet_sample(alphas)
        counter = 0
        print(self.h_theta)
        ##while np.any(self.h_theta == 0) or np.any(self.h_theta == 1):
        ##    counter = counter+1
        ##    print(counter)
        ##    self.h_theta[self.h_theta == 0] = np.random.beta(self.theta_alpha[self.h_theta == 0], self.theta_beta[self.h_theta == 0])
        ##    self.h_theta[self.h_theta == 1] = np.random.beta(self.theta_alpha[self.h_theta == 1], self.theta_beta[self.h_theta == 1])

        #if np.sum(self.h_theta == 0) > 0:
         #   raise Exception('There are some 0\'s in thetas which are generated in the simulation')
        #self.h_theta = self.my_floor(np.random.beta(self.theta_alpha,self.theta_beta),5)
        #self.h_theta[self.h_theta == 0] = 0.00001

    # I am not sure if we should even support this.
    # This is regarding to the situation where we have A
    # but we dont have theta and we want to guess about it.
    # TODO: we should take care of zero D,A and D-A
    def calculate_theta(self,I,S,gamma,phi,C,H):
        self.h_theta = 0.5*np.ones((I,S),dtype=np.float128)
        alpha = np.matmul(gamma*phi*C,np.transpose(H))
        beta = np.matmul(gamma*phi*(1-C),np.transpose(H))
        alpha[alpha<=0]=1
        beta[beta<=0]=1
        self.h_theta = np.random.beta(gamma * alpha, gamma * beta)

        ##while np.any(self.h_theta == 0) or np.any(self.h_theta == 1):
        ##    self.h_theta[self.h_theta == 0] = np.random.beta(alpha[self.h_theta == 0], beta[self.h_theta == 0])
        ##    self.h_theta[self.h_theta == 1] = np.random.beta(alpha[self.h_theta == 1], beta[self.h_theta == 1])

        #if np.sum(self.h_theta == 0)>0:
         #   raise Exception('There are some 0\'s in thetas which are generated in the simulation (calculate part)')
        #self.h_theta = self.my_floor(np.random.beta(gamma * alpha, gamma * beta), 5)
        #self.h_theta[self.h_theta == 0] = 0.00001

    def generate_n(self,n_lambda):
        self.n = np.random.poisson(n_lambda)
        self.n[self.n==0] = (n_lambda[self.n==0])
        #while np.sum(self.n==0)>0:
        #    self.n[self.n==0] = np.random.poisson(n_lambda[self.n==0])

    # Here we use F for multiplying the fractions to a constant value
    # For example if we do not want to have F as input, then F_fraction would be true
    # and here we generate some fractions using Dirichlet distribution but gamma( 0.2,...)
    # would not generate big numbers and because we want the numbers to be very larger than
    # F_epsilon, then we need a constant like 40 (F[0]) to got multiplied by the 0.2
    # also for the second parameter of the gamma distribution, we need a number which is given
    # by F[0][1]
    # TODO: I need to re-code this part with introducing new variable, this is ugly
    def generate_F(self,K,beta):
        temp = np.random.dirichlet(np.ones(K))
        self.F = []
        for k in range(K):
            self.F.append([temp[k], beta])
        self.F = np.array(self.F)


    # I am here generating c_i in all k with arbitrary probabilities of each k would be zero.
    # based on the experiance I just put the more probability for K's of the leaf nodes
    def generate_C(self,I,K,p_c_binom):
        #self.C = np.random.binomial(size=(I,K), n=1, p=np.random.beta(0.5, 1, size=K))
        self.C = np.zeros((I,K),dtype=np.float64)
        for i in range(I):
            half = (np.random.binomial(size=(K-1), n=1, p=p_c_binom)+1)/ 2
            self.C[i][:(K-1)] = np.random.binomial(size=(K-1), n=1, p=p_c_binom) *half
            while self.C[i][:(K-1)].sum()==0:
                half = (np.random.binomial(size=(K-1), n=1, p=p_c_binom)+1)/ 2
                self.C[i][:(K-1)] = np.random.binomial(size=(K-1), n=1, p=p_c_binom) *half
        if(sum(self.C.sum(axis=1)==0)>0):
            print(self.C.sum(axis=1)==0)
            raise Exception("A row in C is zero(0)")
        #print(C)

    def generate_pi(self,zeta,K):
        self.pi = np.random.beta(zeta/K, 1, size=(K))
        while np.any(self.pi == 0):
            self.pi[self.pi == 0] = np.random.beta(zeta[self.pi == 0], 1, size=(K))
        #print(pi)

    def generate_pi_2D(self,zeta,K):
        alpha_sk = np.tile(zeta / K,(K,1))
        T_num = 1
        for t in range(T_num):
            temp = np.transpose(np.random.beta(alpha_sk, 1))
            while np.any(temp == 0):
                temp[temp == 0] = np.transpose(np.random.beta(alpha_sk[temp == 0], 1))
            if t==0:
                self.pi = temp
            else:
                self.pi = self.pi+temp
        self.pi = self.pi/T_num

    # Here we are using binomial to generate random z
    # because z is bernouli ( binomial with only one element )
    def generate_Z(self,K,pi,pi_2D):
        self.Z = np.random.binomial( n=1, p=pi)
        while any(np.sum(self.Z,axis=1)==0):
            if pi_2D is True:
                self.Z[np.sum(self.Z, axis=1) == 0,] =  np.random.binomial(size=(sum(np.sum(self.Z, axis=1) == 0),K), n=1, p=pi[np.sum(self.Z, axis=1) == 0,])
            else:
                self.Z[np.sum(self.Z, axis=1) == 0,] = np.random.binomial(size=(sum(np.sum(self.Z, axis=1) == 0), K), n=1,p=pi)



    #self.Z = H > 0.05
    #self.Z = self.Z.astype(np.int)
    #print(z)
    def generate_G(self,F,S,K,F_epsilon,Z):
        self.G = np.zeros((S,K),dtype=np.float64)
        for s in range(S):
            try:
                alpha = np.power(F[:,0],Z[s][:])*np.power(F_epsilon[:,0],(1-Z[s][:]))
            except Exception:
                print("F or F_fraction has a problem. The exception is happened in generating G in the simulation.")
            scale = np.power(F[:,1],Z[s][:])*np.power(F_epsilon[:,1],(1-Z[s][:]))
            self.G[s][:] = np.random.gamma(shape=alpha,scale=scale,size=K)
        #print(G)

    def generate_H(self,G,S,K,Z):
        #self.H = G / sum(sum(G))
        self.H=np.zeros((S,K),dtype=np.float64)
        #G_modified = G*Z
        #self.H = G_modified/np.transpose(np.tile(G_modified.sum(axis=1),(K,1)))
        self.H = G / np.transpose(np.tile(G.sum(axis=1), (K, 1)))
        #for s in range(S):
         #   sum_temp = sum(G[s])
         #   self.H[s][:] = (G[s]/sum_temp)
        if np.sum(self.H==0)>0:
            temp = int(np.sum(self.H==0))
            print(G)
            print(self.F)
            print(self.H)
            raise Exception(str(temp)+" elements in H are zero in simulation")


    def rand_bin_array(K, N):
        arr = np.zeros(N)
        arr[:K]  = 1
        np.random.shuffle(arr)
        return arr

    def generate_phi(self,r,p,I,K):
        self.phi = np.random.gamma(shape=r,scale = p,size=(I,K))
        if np.sum(self.phi ==0):
            print(np.sum(self.phi ==0))
            raise Exception("phi got zero in simulation of phi")
        #mask = np.random.choice([0, 1,2], size=(I,K), p=[3./10, 3./10,4./10])
        #phi_more = np.random.gamma(shape=r*4,scale = p,size=(I,K))
        #phi_less  = np.random.gamma(shape=r/4,scale = p,size=(I,K))
        #self.phi[mask==0] = phi_more[mask==0]
        #self.phi[mask==1] = phi_less[mask==1]
        #print(phi)

    def generate_D(self,I,S,K,phi,H,n):
        self.D = np.zeros((I,S),dtype=np.int64)
        for s in range(S):
            while(True):
                for i in range(I):
                    try:
                        self.D[i][s]=np.random.poisson(n[s]*sum(H[s][:]*phi[i][:]),size=1)
                    except Exception:
                    	print("problem in generating D")
                if np.sum(self.D,axis=0)[s]>0:
                    break
	


    def generate_A(self,I,S,H,C,D,phi,theta):
        self.A = np.zeros((I,S),dtype=np.int64)
        self.p_binom = np.matmul(phi*C,np.transpose(H))/np.matmul(phi,np.transpose(H))
        self.p_binom = theta*self.p_binom
        try:
            self.A = np.random.binomial(n=D,p=self.p_binom)
        except Exception:
            print("problem in generating A in simulation class")

    def generate_A_theta(self,h_theta,D):
        try:
            self.A = np.random.binomial(n=D,p=h_theta)
        except Exception:
            print("problem in generating A using theta in simulation class")



