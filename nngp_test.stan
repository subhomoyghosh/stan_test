 /* Latent NNGP model*/

  functions{
    
      matrix make_positive_definite(matrix m){
        
       int N = cols(m);
       matrix[N,N] es= eigenvectors_sym(m);
       vector[N] esv= eigenvalues_sym(m);
       real delta = 2 * 1e-5; 
       vector[N] tau; 
       matrix[N,N] dm;
       matrix[N,N] mdm;
       
         for(i in 1:N){
           tau[i] = fmax(0,delta-esv[i]);
         }
      dm = es * diag_matrix(tau) * es'; 
      
      return m + dm; 
      
     }

      
      real nngp_w1_lpdf(vector w1, real sigmasq1, real phi1, matrix NN_dist1,matrix NN_cor1,
                       matrix NN_distM1, matrix NN_corM1, int[,] NN_ind1, int N1, int M){

          vector[N1] V;
          vector[N1] I_Aw = w1;
          int dim;
          int h;

          for (i in 2:N1) {

              matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M]
              iNNdistM;
              matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M]
              iNNCholL;
              vector[ i < (M + 1)? (i - 1) : M] iNNcorr;
              vector[ i < (M + 1)? (i - 1) : M] v;
              row_vector[i < (M + 1)? (i - 1) : M] v2;

              dim = (i < (M + 1))? (i - 1) : M;

              if(dim == 1){iNNdistM[1, 1] = 1;}
              else{
                  h = 0;
                  for (j in 1:(dim - 1)){
                      for (k in (j + 1):dim){
                          h = h + 1;
                          iNNdistM[j, k] = exp(- phi1 * NN_distM1[(i - 1), h] * NN_corM1[(i - 1), h]);
                          iNNdistM[k, j] = iNNdistM[j, k];
                      }
                  }
                  for(j in 1:dim){
                      iNNdistM[j, j] = 1;
                  }
              }
              
              iNNdistM = 0.5*(iNNdistM + iNNdistM');
              iNNdistM = make_positive_definite(iNNdistM);

              iNNCholL = cholesky_decompose(iNNdistM);
              iNNcorr = to_vector(exp(- phi1 * (NN_dist1[(i - 1), 1:dim] .* NN_cor1[(i - 1), 1:dim]) ) );
              v = mdivide_left_tri_low(iNNCholL, iNNcorr);

              V[i] = 1 - dot_self(v);

              v2 = mdivide_right_tri_low(v', iNNCholL);

              I_Aw[i] = I_Aw[i] - v2 * w1[NN_ind1[(i - 1), 1:dim]];

          }
          V[1] = 1;
          return - 0.5 * ( 1 / sigmasq1 * dot_product(I_Aw, (I_Aw ./ V)) +
                          sum(log(V)) + N1 * log(sigmasq1));
      }
    
      
  }


  data {
      int<lower=1> N1;
      int<lower=1> N;
      int<lower=1> M;
      int<lower=1> P;
      vector[N] Y;
      matrix[N, P + 1] X;
      int NN_ind1[N1 - 1, M];
      matrix[N1 - 1, M] NN_dist1;
      matrix[N1 - 1, (M * (M - 1) / 2)] NN_distM1;
      matrix[N1 - 1, M] NN_cor1;
      matrix[N1 - 1, (M * (M - 1) / 2)] NN_corM1;
      matrix[N,N1] H1;
      vector[P + 1] uB;
      matrix[P + 1, P + 1] VB;
      real ss1;
      real st;
      real ap1;
      real bp1;
  }

  transformed data {
      cholesky_factor_cov[P + 1] L_VB;
      L_VB = cholesky_decompose(VB);
  }

  parameters{
      vector[P + 1] beta;
      real<lower = 0> sigma1;
      real<lower = 0> tau;
      real<lower = 0> phi1;
      vector<lower = 0>[N1] w1;
  }

  transformed parameters {
      real sigmasq1 = square(sigma1);
      real tausq = square(tau);
  }

  model{
      beta ~ multi_normal_cholesky(uB, L_VB);
      phi1 ~ gamma(ap1, bp1);
      sigma1 ~ normal(0, ss1);
      tau ~ normal(0, st);
      w1 ~ nngp_w1(sigmasq1, phi1, NN_dist1, NN_cor1, NN_distM1, NN_corM1, NN_ind1, N1, M);
      Y ~ normal(X * beta + H1 * w1, tau);
  }

