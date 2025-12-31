#define ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#define ARMA_DONT_USE_OPENMP
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
//#include <Rcpp.h>
#include <unordered_set>
#include <set>
#include <vector>
using namespace Rcpp;
using namespace arma;
#ifdef _OPENMP
#include <omp.h>
#endif

void R2cpp_index(std::vector<IntegerVector>& vec) {
  for(IntegerVector& v : vec) {
    for(int& x : v) {
      x -= 1;
    }
  }
}

void cpp2R_index(std::vector<IntegerVector>& vec) {
  for(IntegerVector& v : vec) {
    for(int& x : v) {
      x += 1;
    }
  }
}

// [[Rcpp::export]]
arma::mat calAv_cpp(const arma::vec & v,
                    const arma::mat & X){
  arma::mat Av = X*v;
  return(Av);
}

// [[Rcpp::export]]
arma::mat GetS(const arma::mat & X,
               const arma::sp_mat & K,
               const arma::vec & DVec,
               double lambda,
               size_t batch_size = 100000) {

  arma::sp_mat K_modified = K;

  for (arma::uword i = 0; i < DVec.n_elem; i++) {
    K_modified.diag()[i] -= lambda * DVec(i);
  }

  const bool use_batch = (K_modified.n_rows > batch_size);

  arma::mat result(X.n_rows, X.n_rows, arma::fill::zeros);

  if (use_batch) {
    size_t n_rows = K_modified.n_rows;
    size_t num_batches = (n_rows + batch_size - 1) / batch_size;

    for (size_t b = 0; b < num_batches; ++b) {
      size_t start_row = b * batch_size;
      size_t end_row = std::min(start_row + batch_size, n_rows);

      arma::sp_mat K_batch = K_modified.rows(start_row, end_row - 1);
      arma::mat X_batch = X.cols(start_row, end_row - 1);

      arma::mat tmp = K_batch * X.t();
      result += X_batch * tmp;
    }
  } else {
    arma::mat tmp = K_modified * X.t();
    result = X * tmp;
  }

  return result;
}


// [[Rcpp::export]]
arma::mat GetS_new(const arma::mat & X,
                   const arma::sp_mat & K,
                   const arma::sp_mat & B,
                   double & lambda,
                   size_t batch_size = 50000){
  arma::sp_mat K_modified = K -lambda * B;

  const bool use_batch = (K_modified.n_rows > batch_size);

  arma::mat result(X.n_rows, X.n_rows, arma::fill::zeros);

  if (use_batch) {
    size_t n_rows = K_modified.n_rows;
    size_t num_batches = (n_rows + batch_size - 1) / batch_size;

    for (size_t b = 0; b < num_batches; ++b) {
      size_t start_row = b * batch_size;
      size_t end_row = std::min(start_row + batch_size, n_rows);

      arma::sp_mat K_batch = K_modified.rows(start_row, end_row - 1);
      arma::mat X_batch = X.cols(start_row, end_row - 1);

      arma::mat tmp = K_batch * X.t();
      result += X_batch * tmp;
    }
  } else {
    arma::mat tmp = K_modified * X.t();
    result = X * tmp;
  }

  return result;
}

// [[Rcpp::export]]
arma::mat Getmat(const arma::mat& X,
                           const arma::sp_mat& K,
                           const arma::sp_mat& B,
                           size_t batch_size = 10000) {
  const size_t n_rows = X.n_rows;
  const size_t n_cols = X.n_cols;
  arma::mat result(n_rows, n_rows, fill::zeros);

  const size_t num_batches = (n_cols + batch_size - 1) / batch_size;

  const arma::mat Xt = X.t();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
// #pragma omp parallel for schedule(dynamic)
  for (size_t b = 0; b < num_batches; ++b) {
    const size_t start = b * batch_size;
    const size_t end = std::min((b+1)*batch_size, n_cols) - 1;

    const arma::sp_mat K_block = K.rows(start, end);
    const arma::sp_mat B_block = B.rows(start, end);

    arma::mat KX = K_block * Xt;
    arma::mat BX = B_block * Xt;
    arma::mat KBX = KX - BX;

    const arma::mat X_block = X.cols(start, end);

#ifdef _OPENMP
#pragma omp critical
#endif
    result += X_block * KBX;
  }

  return result;
}

arma::mat pairwiseEuclideanDistance(const arma::mat & A,
                                    const arma::mat & B) {
  int n = A.n_cols;
  int m = B.n_cols;

  arma::mat distances(n, m);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      double dist = norm(A.col(i) - B.col(j), 2);
      distances(i, j) = dist;
    }
  }

  return distances;
}

// [[Rcpp::export]]
arma::mat FindMedian_self(const arma::mat & X,
                     std::vector<IntegerVector> & idx){
  R2cpp_index(idx);

  int n = idx.size();
  arma::mat out(n,n);
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n; ++j){
      arma::mat tmp = X.submat(Rcpp::as<arma::uvec>(idx[i]),Rcpp::as<arma::uvec>(idx[j]));
      out.at(i,j) = arma::median(arma::vectorise(tmp));
    }
  }
  cpp2R_index(idx);
  return(out);
}

// [[Rcpp::export]]
arma::mat FindMedian(const arma::mat & X,
                     std::vector<IntegerVector> & idx1,
                     std::vector<IntegerVector> & idx2){
  R2cpp_index(idx1);
  R2cpp_index(idx2);
  int n1 = idx1.size();
  int n2 = idx2.size();
  arma::mat out(n1,n2);
  for(int i = 0; i < n1; ++i){
    for(int j = 0; j < n2; ++j){
      arma::mat tmp = X.submat(Rcpp::as<arma::uvec>(idx1[i]),Rcpp::as<arma::uvec>(idx2[j]));
      out.at(i,j) = arma::median(arma::vectorise(tmp));
    }
  }
  cpp2R_index(idx1);
  cpp2R_index(idx2);
  return(out);
}

// [[Rcpp::export]]
List PairwiseKernel_norm(const arma::mat & sim,
                                 const arma::mat & groupsim,
                                 std::vector<IntegerVector> & idx1,
                                 std::vector<IntegerVector> & idx2,
                                 const int & seed){

  arma::sp_mat out1(sim.n_rows, sim.n_cols);

  R2cpp_index(idx1);
  R2cpp_index(idx2);

  int nr = idx1.size();
  for(int i = 0; i < nr; ++i){

    arma::rowvec row = groupsim.row(i);
    if (arma::accu(row != 0) > 0) {
      int cross_group = arma::index_max(row);
      arma::uvec r_idx = Rcpp::as<arma::uvec>(idx1[i]);
      arma::uvec c_idx = Rcpp::as<arma::uvec>(idx2[cross_group]);

      for (arma::uword ii = 0; ii < r_idx.n_elem; ++ii) {
        arma::uword row_idx = r_idx(ii);
        arma::uvec col_indices = c_idx;
        arma::uvec col_idx = arma::shuffle(col_indices);
        out1(row_idx, col_idx(0)) += 1;
      }

      for (arma::uword jj = 0; jj < c_idx.n_elem; ++jj) {
        arma::uword col_idx = c_idx(jj);
        arma::uvec row_indices = r_idx;
        arma::uvec row_idx = arma::shuffle(row_indices);
        out1(row_idx(0), col_idx) += 1;
      }
    }
  }

  int f1 = arma::accu(out1);
  if(f1 != 0){
    out1 /= f1;
  }

  arma::sp_mat out2(sim.n_rows, sim.n_cols);
  int nc = idx2.size();
  for(int i = 0; i < nc; ++i){

    arma::colvec col = groupsim.col(i);
    if (arma::accu(col != 0) > 0) {
      int cross_group = arma::index_max(col);
      arma::uvec r_idx = Rcpp::as<arma::uvec>(idx1[cross_group]);
      arma::uvec c_idx = Rcpp::as<arma::uvec>(idx2[i]);

      for (arma::uword ii = 0; ii < r_idx.n_elem; ++ii) {
        arma::uword row_idx = r_idx(ii);
        arma::uvec col_indices = c_idx;
        arma::uvec col_idx = arma::shuffle(col_indices);
        out2(row_idx, col_idx(0)) += 1;
      }

      for (arma::uword jj = 0; jj < c_idx.n_elem; ++jj) {
        arma::uword col_idx = c_idx(jj);
        arma::uvec row_indices = r_idx;
        arma::uvec row_idx = arma::shuffle(row_indices);
        out2(row_idx(0), col_idx) += 1;
      }

    }

  }

  int f2 = arma::accu(out2);
  if(f2 != 0){
    out2 /= f2;
  }

  arma::sp_mat out = (out1+out2)*0.5*std::max(f1,f2);
  cpp2R_index(idx1);
  cpp2R_index(idx2);

  return List::create(Named("subK") = out,
                      Named("scaler") = std::max(f1,f2));

}

// [[Rcpp::export]]
arma::sp_mat Insert_submat(arma::sp_mat & X,
                           std::vector<arma::sp_mat> & submat,
                           IntegerVector & idx,
                           const bool & scale,
                           const double & scaler){

  int n = idx.size();
  IntegerVector dim(n);

  for(int i = 0;i < n; ++i){
    NumericVector p(i+1);
    for(int j = 0; j < i+1; ++j){
      p[j] = idx[j]-1;
    }
    dim[i] = std::accumulate(p.begin(),p.end(),0);
  }

  int count = 0;
  for(int i = 0; i < (n-1); ++i){
    for(int j = i+1; j < n; ++j){
      if(i == 0){
        X.submat(0,dim[j-1],dim[i],dim[j]) = submat[count];
        X.submat(dim[j-1],0,dim[j],dim[i]) = submat[count].t();
        count++;
      }else{
        X.submat(dim[i-1],dim[j-1],dim[i],dim[j]) = submat[count];
        X.submat(dim[j-1],dim[i-1],dim[j],dim[i]) = submat[count].t();
        count++;
      }
    }
  }

  if(scale == true){
    X *= scaler;
  }

  return(X);
}

// [[Rcpp::export]]
List SubKernel(const arma::vec& x,
               const arma::vec& p,
               const arma::vec& i,
               const arma::vec& q,
               int k) {

  std::vector<int> xx(x.size());
  int n = p.size() - 1;
  for (int c = 0; c < n; c++) {
    if (p[c+1] - p[c] <= k) {
      for (int j = p[c]; j < p[c + 1]; j++) {
        xx[j] = 1;
      }
    } else {
      arma::vec x_sub = arma::abs(x.subvec(p[c], p[c+1]-1) - q[c]);
      arma::uvec sorted_indices = arma::sort_index(x_sub);
      arma::uvec top_idx = sorted_indices.head(k) + p[c];
      for(unsigned int j = 0; j < top_idx.n_elem; ++j){
        xx[top_idx[j]] = 1;
      }
    }
  }
  return Rcpp::List::create(Named("i") = i + 1, Named("j") = p, Named("x") = xx);
}

// [[Rcpp::export]]
arma::mat FindGroupLink(const arma::mat & X,
                        const arma::mat & row_mat,
                        const arma::mat & col_mat,
                        const double Angle_var = 15,
                        const double max_Angle = 40){

  int nrow = X.n_rows;
  int ncol = X.n_cols;
  arma::mat out1(nrow,ncol,arma::fill::zeros);
  arma::mat out2 = out1;

  double threshold = std::cos(max_Angle);
  for(int i = 0; i < nrow; ++i){
    double v1 = std::cos(std::acos(row_mat(i,i) + Angle_var));
    double cutoff = std::max(threshold,v1);
    out1.row(i) = arma::conv_to<arma::rowvec>::from(X.row(i) >= cutoff);
  }

  for(int j = 0; j < ncol; ++j){
    double v1 = std::cos(std::acos(col_mat(j,j) + Angle_var));
    double cutoff = std::max(threshold,v1);
    out2.col(j) = arma::conv_to<arma::colvec>::from(X.col(j) >= cutoff);
  }

  arma::mat out = out1%out2;

  return(X%out);
}

// [[Rcpp::export]]
List sKernel_norm(const std::vector<int> & N_list,
                       const int & N,
                       const int & clust,
                       std::vector<std::vector<IntegerVector>> & idx){

  arma::sp_mat out(N,N);
  int n_batch = N_list.size();
  IntegerVector dim(n_batch);

  for(int i = 0;i < n_batch; ++i){
    int sum = 0;
    for(int j = 0; j <= i; ++j){
      sum += N_list[j];
    }
    dim[i] = sum;

  }

  IntegerVector scaler_global;

  for(int i = 0; i < clust; ++i){
    arma::sp_mat tmp_out(N,N);
    IntegerVector scaler;
    for(int b1 = 0; b1 < (n_batch-1); ++b1){

      std::vector<IntegerVector> idx1 = idx[b1];
      R2cpp_index(idx1);

      IntegerVector tmp = idx1[i];
      int tmpi = tmp[0];
      if(tmpi >= 0){
        for(int b2 = b1+1; b2 < n_batch; ++b2){

          std::vector<IntegerVector> idx2 = idx[b2];
          R2cpp_index(idx2);
          IntegerVector tmpp = idx2[i];
          int tmpj = tmpp[0];
          if(tmpj >= 0){
            arma::sp_mat sub_mat;
            if(b1 == 0){
              sub_mat = tmp_out.submat(0,dim[b2-1],dim[b1]-1,dim[b2]-1);
            }else{
              sub_mat = tmp_out.submat(dim[b1-1],dim[b2-1],dim[b1]-1,dim[b2]-1);
            }

             for (size_t rowmat = 0; rowmat < tmp.size(); ++rowmat) {
               int row = tmp[rowmat];
              for (size_t colmat = 0; colmat < tmpp.size(); ++colmat) {
                int col = tmpp[colmat];
                sub_mat(row, col) = 1;
              }
            }

            int tt = tmp.size() * tmpp.size();
            if(tt != 0){
              sub_mat /= tt;
            }
            scaler.push_back(tt);

            if(b1 == 0){
              tmp_out.submat(0,dim[b2-1],dim[b1]-1,dim[b2]-1) = sub_mat;
              tmp_out.submat(dim[b2-1],0,dim[b2]-1,dim[b1]-1) = sub_mat.t();
            }else{
              tmp_out.submat(dim[b1-1],dim[b2-1],dim[b1]-1,dim[b2]-1) = sub_mat;
              tmp_out.submat(dim[b2-1],dim[b1-1],dim[b2]-1,dim[b1]-1) = sub_mat.t();
            }

          }
          cpp2R_index(idx2);
        }
      }
      cpp2R_index(idx1);
    }

    int ss = max(scaler);
    tmp_out *= ss;
    int gg = arma::accu(tmp_out);
    if(gg != 0){
      tmp_out /= gg;
    }
    scaler_global.push_back(gg);
    out += tmp_out;
  }

  IntegerVector scaler_global2;
  arma::sp_mat out2(N,N);
  for(int i = 0; i < clust; ++i){
    arma::sp_mat tmp_out(N,N);
    IntegerVector scaler;
    for(int b1 = 0; b1 < n_batch; ++b1){

      std::vector<IntegerVector> idx1 = idx[b1];
      R2cpp_index(idx1);

      IntegerVector tmp = idx1[i];
      int tmpi = tmp[0];
      if(tmpi >= 0){

        arma::sp_mat sub_mat;
        if(b1 == 0){
          sub_mat = tmp_out.submat(0,0,dim[b1]-1,dim[b1]-1);
        }else{
          sub_mat = tmp_out.submat(dim[b1-1],dim[b1-1],dim[b1]-1,dim[b1]-1);
        }
        for (size_t rowmat = 0; rowmat < tmp.size(); ++rowmat) {
          int row = tmp[rowmat];
          for (size_t colmat = 0; colmat < tmp.size(); ++colmat) {
            int col = tmp[colmat];
            sub_mat(row, col) = 1;
          }
        }

        int tt = tmp.size() * tmp.size();
        if(tt != 0){
          sub_mat /= tt;
        }
        scaler.push_back(tt);

        if(b1 == 0){
          tmp_out.submat(0,0,dim[b1]-1,dim[b1]-1) = sub_mat;
        }else{
          tmp_out.submat(dim[b1-1],dim[b1-1],dim[b1]-1,dim[b1]-1) = sub_mat;
        }

      }
      cpp2R_index(idx1);
    }

    int ss = max(scaler);
    if(ss != 0){
      tmp_out *= ss;
    }
    int gg = arma::accu(tmp_out);
    if(gg != 0){
      tmp_out /= gg;
    }
    scaler_global2.push_back(gg);
    out2 += tmp_out;
  }
  int ts = max(scaler_global2);
  arma::vec KD = arma::vec(arma::sum(out,1));
  arma::vec BD = arma::vec(arma::sum(out2,1));
  uvec idxD = find(KD == 0);
  bool has_zeros = !idxD.empty();
  if(has_zeros){
    for (size_t i = 0; i < idxD.n_elem; ++i) {
      out(idxD[i], idxD[i]) = BD[idxD[i]];
    }
  }

  out *= ts;
  out2 *= ts;
  return List::create(Named("K") = out, Named("B") = out2);

}

arma::mat scaler_assign(const arma::mat & A,
                        const arma::mat & B,
                        const arma::vec & index) {

  if (A.n_cols == 0 || B.n_cols == 0) {
    Rcpp::stop("scaler_assign: ncol of A or B is 0.");
  }
  if (index.n_elem == 0) {
    Rcpp::stop("scaler_assign: index empty.");
  }

  arma::uvec idx_u = arma::conv_to<arma::uvec>::from(index);

  if (idx_u.max() >= B.n_cols) {
    Rcpp::Rcout << "scaler_assign: index.max() = " << idx_u.max()
                << ", B.n_cols = " << B.n_cols << std::endl;
    Rcpp::stop("scaler_assign: index 中存在超过 B 列数的索引。");
  }

  arma::rowvec la = arma::sqrt(arma::sum(arma::square(A), 0));
  arma::mat sub_B = B.cols(idx_u);
  arma::rowvec lb = arma::sqrt(arma::sum(arma::square(sub_B), 0));

  arma::uword n   = la.n_elem;
  arma::uword len = lb.n_elem;

  if (len == 0) {
    Rcpp::stop("scaler_assign: sub_B 没有列（len == 0），无法进行采样。");
  }

  arma::uvec sorted_indices = arma::sort_index(la, "descend");

  arma::rowvec sampled(n);
  for (arma::uword i = 0; i < n; ++i) {
    arma::uword idx = arma::randi<arma::uword>(arma::distr_param(0, len - 1));
    sampled(i) = lb(idx);
  }

  arma::rowvec sorted_sample = arma::sort(sampled, "descend");
  arma::rowvec scaler(n, arma::fill::zeros);
  for (arma::uword i = 0; i < n; ++i) {
    scaler(sorted_indices(i)) = sorted_sample(i);
  }

  arma::rowvec la_safe = la;
  la_safe.replace(0.0, 1.0);

  arma::mat out = A * arma::diagmat(1.0 / la_safe) * arma::diagmat(scaler);
  return out;
}

arma::vec get_index(const arma::mat & A,
                    const arma::mat & B,
                    const int & scaler_k){
  int n = A.n_cols;
  arma::mat dist = pairwiseEuclideanDistance(A,B);
  arma::vec scale_index(n*scaler_k);
  int pos = 0;
  for(arma::uword i = 0; i < n; ++i){
    arma::rowvec rowi = dist.row(i);
    arma::uvec new_idx = arma::sort_index(rowi,"ascend");
    scale_index.subvec(pos,pos+scaler_k-1) = arma::conv_to<arma::vec>::from(new_idx.subvec(0,scaler_k-1));
    pos += scaler_k;
  }
  return(scale_index);
}

arma::mat Get_var(const arma::mat & D,
                  const arma::mat & Z,
                  const arma::mat & O,
                  const int & k){

  int cells = Z.n_cols;
  int dim = Z.n_rows;
  arma::mat var(dim,cells);
  for(arma::uword i = 0; i < cells; ++i){

    arma::uvec sorted_indices = arma::sort_index(D.row(i), "ascend");
    arma::mat sub_O = O.cols(sorted_indices.head(k+1));
    arma::vec stddev(dim);

    for(arma::uword j = 0; j < dim; ++j){
      stddev(j) = arma::stddev(sub_O.row(j));
    }

    var.col(i) = stddev;
  }

  return(var);
}

// [[Rcpp::export]]
List modal_infer_one_one(const arma::mat & X,
                                 const arma::mat & Y,
                                 const int & k,
                                 const bool & L2,
                                 const arma::mat & Z,
                                 const arma::mat & T,
                                 const bool & do_scaler,
                                 const int & scaler_k,
                                 const bool & do_var,
                                 const int & var_k){

  arma::mat dist(X.n_cols,Y.n_cols);
  if(L2 == true){
    dist = 1-X.t()*Y;
  }else{
    dist = pairwiseEuclideanDistance(X,Y);
  }
  arma::sp_mat M(dist.n_rows,dist.n_cols);
  const arma::uword n = X.n_cols;

  arma::vec scale_index(X.n_cols*k);
  int pos = 0;

  arma::mat subdsit(n, k);
  for(arma::uword i = 0; i < n; ++i){
    arma::rowvec rowi = dist.row(i);
    std::partial_sort(rowi.begin(), rowi.begin() + k, rowi.end());
    subdsit.row(i) = rowi.subvec(0, k - 1);

    arma::rowvec rowii = dist.row(i);
    arma::uvec top_k_indices = arma::sort_index(rowii, "ascend");
    if(do_scaler == true){
      scale_index.subvec(pos,pos+k-1) = arma::conv_to<arma::vec>::from(top_k_indices.subvec(0,k- 1));//new
      pos += k;//new
    }
  }

  arma::mat subdsit1 = exp(-subdsit);
  arma::vec row_sum_w = sum(subdsit1, 1)+1e-10;

  arma::mat weight_mat = arma::diagmat(1 / row_sum_w)*subdsit1;

  for(arma::uword i = 0; i < n; ++i){
    arma::uvec sorted_indices = arma::sort_index(dist.row(i));
    for(int j = 0; j < k; ++j){
      M(i,sorted_indices[j]) = weight_mat(i,j);
    }
  }

  arma::mat exp = Z * M.t();
  arma::mat inferMat(T.n_rows,X.n_cols);

  if(do_var == true){
    arma::mat D = pairwiseEuclideanDistance(exp,T);
    arma::mat var = Get_var(D,exp,T,var_k);
    int cells = X.n_cols;
    int dim = T.n_rows;
    for(int i = 0; i < cells; ++i){
      arma::colvec tmp = exp.col(i);
      arma::colvec samples = arma::randn<arma::colvec>(dim);
      arma::colvec tmp_var = var.col(i);
      arma::colvec rand = samples % tmp_var;
      inferMat.col(i) = tmp + rand;
    }
  }else{
    inferMat += exp;
  }

  if(do_scaler){
    arma::mat input_Z = Z.cols(arma::conv_to<arma::uvec>::from(scale_index));
    arma::vec idx = get_index(input_Z,Z,scaler_k);
    Function unique_r("unique");
    arma::vec sub_idx = Rcpp::as<arma::vec>(unique_r(idx));
    arma::mat Mat = scaler_assign(inferMat,T,sub_idx);
    inferMat = Mat;
  }

  return List::create(Named("M") = M,
                      Named("infer") = inferMat);

}

// [[Rcpp::export]]
List modal_infer_one_multi(const std::vector<arma::mat> & X,
                                   const std::vector<arma::mat> & Y,
                                   const int & k,
                                   const bool & L2,
                                   const arma::mat & Z,
                                   const arma::mat & T,
                                   const bool & do_scaler,
                                   const int & scaler_k,
                                   const bool & do_var,
                                   const int & var_k){

  int l = X.size();
  arma::mat dist;
  if(L2 == true){
    for(int i = 0; i < l; ++i){
      arma::mat tmp = 1-X[i].t()*Y[i];
      if(i == 0){
        dist = tmp;
      }else{
        dist = arma::join_rows(dist,tmp);
      }
    }
  }else{
    Function cor_r("cor");
    for(int i = 0; i < l; ++i){
      arma::mat tmp = 1-Rcpp::as<arma::mat>(cor_r(X[i],Y[i]));
      if(i == 0){
        dist = tmp;
      }else{
        dist = arma::join_rows(dist,tmp);
      }
    }
  }

  const arma::uword n = dist.n_rows;
  arma::sp_mat M(n,dist.n_cols);
  arma::mat subdsit(n, k);

  arma::vec scale_index(X[0].n_cols*k);
  int pos = 0;

  for(arma::uword i = 0; i < n; ++i){
    arma::rowvec rowi = dist.row(i);
    std::partial_sort(rowi.begin(), rowi.begin() + k, rowi.end());
    subdsit.row(i) = rowi.subvec(0, k - 1);

    arma::rowvec rowii = dist.row(i);
    arma::uvec top_k_indices = arma::sort_index(rowii, "ascend");
    if(do_scaler == true){
      scale_index.subvec(pos,pos+k-1) = arma::conv_to<arma::vec>::from(top_k_indices.subvec(0,k- 1));//new
      pos += k;
    }
  }

  arma::mat subdsit1 = exp(-subdsit);
  arma::vec row_sum_w = sum(subdsit1, 1)+1e-10;
  arma::mat weight_mat = arma::diagmat(1 / row_sum_w)*subdsit1;

  for(arma::uword i = 0; i < n; ++i){
    arma::uvec sorted_indices = arma::sort_index(dist.row(i));
    for(int j = 0; j < k; ++j){
      M(i,sorted_indices[j]) = weight_mat(i,j);
    }
  }

  arma::mat exp = Z * M.t();
  arma::mat inferMat(T.n_rows,X[0].n_cols);
  if(do_var == true){
    arma::mat D = pairwiseEuclideanDistance(exp,T);
    arma::mat var = Get_var(D,exp,T,var_k);
    int cells = X[0].n_cols;
    int dim = T.n_rows;
    for(int i = 0; i < cells; ++i){
      arma::colvec tmp = exp.col(i);
      arma::colvec samples = arma::randn<arma::colvec>(dim);
      arma::colvec tmp_var = var.col(i);
      arma::colvec rand = samples % tmp_var;
      inferMat.col(i) = tmp + rand;
    }
  }else{
    inferMat += exp;
  }

  if (do_scaler) {
    if (scale_index.n_elem == 0) {
      Rcpp::stop("do_scaler = TRUE，但是 scale_index 为空。");
    }
    if (scale_index.max() >= (int)Z.n_cols) {
      Rcpp::Rcout << "scale_index.max() = " << scale_index.max()
                  << ", Z.n_cols = " << Z.n_cols << std::endl;
      Rcpp::stop("scale_index 中存在超过 Z 列数的索引，请检查 k / dist / Z 的一致性。");
    }

    arma::uvec nn_idx = arma::conv_to<arma::uvec>::from(scale_index);
    arma::mat input_Z = Z.cols(nn_idx);

    arma::vec idx = get_index(input_Z, T, scaler_k);

    Function unique_r("unique");
    arma::vec sub_idx = Rcpp::as<arma::vec>(unique_r(idx));

    arma::mat Mat = scaler_assign(inferMat, T, sub_idx);
    inferMat = Mat;
  }

  return List::create(Named("M") = M,
                      Named("infer") = inferMat);
}

void replace_with_vector(arma::mat& M, const arma::uvec& V) {

  for (arma::uword i = 0; i < M.n_rows; i++) {
    for (arma::uword j = 0; j < M.n_cols; j++) {
      int tmp = M(i, j);
      M(i, j) = V(tmp);
    }
  }
}

// [[Rcpp::export]]
List modal_infer_multi_one(const arma::mat & X,
                                   const std::vector<std::vector<arma::mat>> & Y,
                                   const int & k,
                                   const bool & L2,
                                   const arma::mat & Z,
                                   const arma::mat & T,
                                   const bool & do_scaler,
                                   const int & scaler_k,
                                   const bool & do_var,
                                   const int & var_k){

  const arma::uword n = X.n_cols;
  int l = Y.size();
  arma::mat subdsit(n, k);
  arma::mat dist(X.n_cols,Y[0][0].n_cols);
  if(L2 == true){
    dist= 1-X.t()*Y[0][0];
  }else{
    dist = pairwiseEuclideanDistance(X,Y[0][0]);
  }

  arma::sp_mat M(n,Y[l-1][0].n_cols);
  for(arma::uword i = 0; i < n; ++i){
    arma::rowvec rowi = dist.row(i);
    std::partial_sort(rowi.begin(), rowi.begin() + k, rowi.end());
    subdsit.row(i) = rowi.subvec(0, k - 1);
  }

  arma::mat subdsit1 = exp(-subdsit);
  arma::vec row_sum_w = sum(subdsit1, 1)+1e-10;
  arma::mat weight_mat = arma::diagmat(1 / row_sum_w)*subdsit1;

  arma::mat indices(n, k);
  for(arma::uword i = 0; i < n; ++i){
    arma::uvec tmp =  arma::sort_index(dist.row(i));
    tmp = tmp.head(k);
    indices.row(i) = arma::conv_to<arma::rowvec>::from(tmp.t());
  }

  int i = 0;
  while(i < (l-1)){
    arma::mat tmpdist(Y[i][1].n_cols,Y[i+1][0].n_cols);
    if(L2 == true){
      tmpdist = 1-Y[i][1].t()*Y[i+1][0];
    }else{
      tmpdist =  pairwiseEuclideanDistance(Y[i][1],Y[i+1][0]);
    }
    arma::uvec minIndices = arma::index_min(tmpdist, 1);
    replace_with_vector(indices,minIndices);
    i +=1;
  }


  for(arma::uword i = 0; i < n; ++i){
    for(int j = 0; j < k; ++j){
      M(i,indices(i,j)) = weight_mat(i,j);
    }
  }

  arma::mat exp = Z * M.t();
  arma::mat inferMat(T.n_rows,X.n_cols);
  if(do_var == true){
    arma::mat D = pairwiseEuclideanDistance(exp,T);
    arma::mat var = Get_var(D,exp,T,var_k);
    int dim = T.n_rows;
    for(int i = 0; i < n; ++i){
      arma::colvec tmp = exp.col(i);
      arma::colvec samples = arma::randn<arma::colvec>(dim);
      arma::colvec tmp_var = var.col(i);
      arma::colvec rand = samples % tmp_var;
      inferMat.col(i) = tmp + rand;
    }
  }else{
    inferMat += exp;
  }

  if (do_scaler) {
    arma::uvec scale_index = arma::conv_to<arma::uvec>::from( arma::vectorise(indices) );

    if (scale_index.max() >= Z.n_cols) {
      Rcpp::stop("scaler: scale_index 中存在 >= Z.n_cols 的索引，可能越界。");
    }

    arma::mat input_Z = Z.cols(scale_index);
    arma::vec idx_scale = get_index(input_Z, Z, scaler_k);
    Function unique_r("unique");
    arma::vec sub_idx = Rcpp::as<arma::vec>(unique_r(idx_scale));
    arma::uvec sub_idx_u = arma::conv_to<arma::uvec>::from(sub_idx);
    if (sub_idx_u.max() >= T.n_cols) {
      Rcpp::stop("scaler: sub_idx 中存在 >= T.n_cols 的索引，可能越界。");
    }

    arma::mat Mat = scaler_assign(inferMat, T, sub_idx);
    inferMat = Mat;
  }

  return List::create(Named("M") = M,
                      Named("infer") = inferMat);

}

// [[Rcpp::export]]
List modal_infer_multi_multi(const std::vector<arma::mat> & X,
                             const std::vector<std::vector<std::vector<arma::mat>>> & Y,
                             const int & k,
                             const bool & L2,
                             const arma::mat & Z,
                             const arma::mat & T,
                             const bool & do_scaler,
                             const int & scaler_k,
                             const bool & do_var,
                             const int & var_k){

  const arma::uword n = X[0].n_cols;
  int l = Y.size();
  Rcpp::Function cor_r("cor"), cumsum_r("cumsum");

  std::vector<int> cols_per_path(l);
  for (int i = 0; i < l; ++i) {
    cols_per_path[i] = Y[i][0][0].n_cols;
  }

  std::vector<int> col_start(l + 1);
  col_start[0] = 0;
  for (int i = 0; i < l; ++i) {
    col_start[i + 1] = col_start[i] + cols_per_path[i];
  }
  int total_start_cols = col_start[l];

  int ll = Y[0].size();
  Rcpp::IntegerVector na(l + 1);
  na[0] = 0;
  for (int i = 0; i < l; ++i) {
    na[i + 1] = Y[i][ll - 1][1].n_cols;
  }
  Rcpp::IntegerVector acc_na = cumsum_r(na);
  int total_arrive_cols = acc_na[l];

  if ((int)Z.n_cols != total_arrive_cols) {
    Rcpp::stop("modal_infer_multi_multi: Z.n_cols(%d) != total_arrive_cols(%d)，请检查 arrive.data 的构造。",
               Z.n_cols, total_arrive_cols);
  }

  arma::mat dist(n, total_start_cols, arma::fill::zeros);

  for (int i = 0; i < l; ++i) {
    arma::mat Xi = X[i];
    arma::mat Y0 = Y[i][0][0];
    arma::mat sub_dist;
    if (L2) {
      sub_dist = 1 - Rcpp::as<arma::mat>(cor_r(Xi, Y0));
    } else {
      sub_dist = 1 - Xi.t() * Y0;
    }

    if ((int)sub_dist.n_rows != (int)n) {
      Rcpp::stop("sub_dist.n_rows(%d) != n(%d) 在路径 %d", sub_dist.n_rows, n, i);
    }
    if (sub_dist.n_cols != (arma::uword)cols_per_path[i]) {
      Rcpp::stop("sub_dist.n_cols(%d) != cols_per_path[%d](%d)",
                 sub_dist.n_cols, i, cols_per_path[i]);
    }

    int cs = col_start[i];
    int ce = col_start[i + 1] - 1;
    if (ce >= (int)dist.n_cols) {
      Rcpp::stop("dist.submat 越界：路径 %d, ce=%d >= dist.n_cols=%d",
                 i, ce, dist.n_cols);
    }

    dist.submat(0, cs, n - 1, ce) = sub_dist;
  }

  arma::mat subdist(n, k);
  for (arma::uword i = 0; i < n; ++i) {
    arma::rowvec r = dist.row(i);
    std::partial_sort(r.begin(), r.begin() + k, r.end());
    subdist.row(i) = r.subvec(0, k - 1);
  }

  arma::mat subdist1 = arma::exp(-subdist);
  arma::vec row_sum_w = arma::sum(subdist1, 1) + 1e-10;
  arma::mat weight_mat = arma::diagmat(1.0 / row_sum_w) * subdist1;

  arma::mat indices(n, k, arma::fill::ones);
  for (arma::uword i = 0; i < n; ++i) {
    arma::uvec tmp = arma::sort_index(dist.row(i));
    tmp = tmp.head(k);
    indices.row(i) = arma::conv_to<arma::rowvec>::from(tmp.t());
  }

  std::vector<arma::uvec> IDX(l);
  int step = 0;
  while (step < (ll - 1)) {
    for (int i = 0; i < l; ++i) {
      arma::mat Y_curr = Y[i][step][1];
      arma::mat Y_next = Y[i][step + 1][0];

      arma::mat tmpdist(Y_curr.n_cols, Y_next.n_cols);
      if (L2) {
        tmpdist = 1 - Y_curr.t() * Y_next;
      } else {
        tmpdist = pairwiseEuclideanDistance(Y_curr, Y_next);
      }

      arma::uvec minIndices = arma::index_min(tmpdist, 1);
      if (step != 0) {
        arma::uvec tmp = IDX[i];
        IDX[i] = minIndices.elem(tmp);
      } else {
        IDX[i] = minIndices;
      }
    }
    ++step;
  }

  for (int i = 0; i < l; ++i) {
    IDX[i] += acc_na[i];
  }

  arma::uvec idx(total_start_cols);
  {
    arma::uword pos = 0;
    for (int i = 0; i < l; ++i) {
      arma::uvec &cur = IDX[i];
      arma::uword len_i = cur.n_elem;
      if (pos + len_i > (arma::uword)total_start_cols) {
        Rcpp::stop("modal_infer_multi_multi: idx 填充越界，路径 %d。", i);
      }
      idx.subvec(pos, pos + len_i - 1) = cur;
      pos += len_i;
    }
    if (pos != (arma::uword)total_start_cols) {
      Rcpp::stop("modal_infer_multi_multi: idx 填充结束 pos(%d) != total_start_cols(%d)。",
                 pos, total_start_cols);
    }
  }

  replace_with_vector(indices, idx);

  if (indices.max() >= total_arrive_cols) {
    Rcpp::stop("modal_infer_multi_multi: 映射后 indices.max()(%f) >= total_arrive_cols(%d)，索引越界。",
               indices.max(), total_arrive_cols);
  }

  arma::sp_mat M(n, total_arrive_cols);
  for (arma::uword i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      int col = (int)indices(i, j);
      M(i, col) = weight_mat(i, j);
    }
  }

  if ((int)Z.n_cols != (int)M.n_cols) {
    Rcpp::stop("modal_infer_multi_multi: Z.n_cols(%d) != M.n_cols(%d)。",
               Z.n_cols, M.n_cols);
  }

  arma::mat exp = Z * M.t();          // d × n
  arma::mat inferMat(T.n_rows, n, arma::fill::zeros);

  if (do_var) {
    arma::mat D   = pairwiseEuclideanDistance(exp, T);
    arma::mat var = Get_var(D, exp, T, var_k);
    int dim = T.n_rows;
    for (arma::uword i = 0; i < n; ++i) {
      arma::colvec base   = exp.col(i);
      arma::colvec noise  = arma::randn<arma::colvec>(dim) % var.col(i);
      inferMat.col(i)     = base + noise;
    }
  } else {
    inferMat = exp;
  }

  if (do_scaler) {
    arma::rowvec scale_index_row = arma::vectorise(indices).t();
    arma::vec    scale_index     = arma::conv_to<arma::vec>::from(scale_index_row);
    arma::uvec   scale_index_u   = arma::conv_to<arma::uvec>::from(scale_index);

    if (scale_index_u.max() >= Z.n_cols) {
      Rcpp::stop("modal_infer_multi_multi: scale_index.max()(%d) >= Z.n_cols(%d)。",
                 scale_index_u.max(), Z.n_cols);
    }

    arma::mat input_Z = Z.cols(scale_index_u);

    arma::vec idx_scale = get_index(input_Z, T, scaler_k);

    Rcpp::Function unique_r("unique");
    arma::vec sub_idx = Rcpp::as<arma::vec>(unique_r(idx_scale));

    arma::uvec sub_idx_u = arma::conv_to<arma::uvec>::from(sub_idx);
    if (sub_idx_u.max() >= T.n_cols) {
      Rcpp::stop("modal_infer_multi_multi: sub_idx.max()(%d) >= T.n_cols(%d)。",
                 sub_idx_u.max(), T.n_cols);
    }

    arma::mat Mat = scaler_assign(inferMat, T, sub_idx);
    inferMat = Mat;
  }

  return Rcpp::List::create(
    Rcpp::Named("M")     = M,
    Rcpp::Named("infer") = inferMat
  );
}

arma::sp_mat generateMatchMatrixArma(const CharacterVector& batch1,
                                     const CharacterVector& batch2) {
  int n1 = batch1.size();
  int n2 = batch2.size();
  if (n1 == 0 || n2 == 0) {
    return arma::sp_mat(n1, n2);
  }


  std::unordered_map<std::string, std::vector<int>> map1, map2;
  for (int i = 0; i < n1; ++i) {
    map1[as<std::string>(batch1[i])].push_back(i);
  }
  for (int j = 0; j < n2; ++j) {
    map2[as<std::string>(batch2[j])].push_back(j);
  }


  std::vector<arma::uword> sorted_rows, sorted_cols;
  std::vector<double> sorted_vals;
  for (const auto& pair : map1) {
    const std::string& label = pair.first;
    if (map2.count(label)) {
      for (int row : pair.second) {
        for (int col : map2.at(label)) {
          if (row >= n1 || col >= n2) {
            Rcpp::stop("索引越界：行 %d >= %d 或列 %d >= %d", row, n1, col, n2);
          }
          sorted_rows.push_back(static_cast<arma::uword>(row));
          sorted_cols.push_back(static_cast<arma::uword>(col));
          sorted_vals.push_back(1.0);
        }
      }
    }
  }

  std::vector<size_t> indices(sorted_rows.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
    if (sorted_cols[a] == sorted_cols[b]) {
      return sorted_rows[a] < sorted_rows[b];
    }
    return sorted_cols[a] < sorted_cols[b];
  });
  std::vector<arma::uword> final_rows, final_cols;
  std::vector<double> final_vals;
  for (auto idx : indices) {
    final_rows.push_back(sorted_rows[idx]);
    final_cols.push_back(sorted_cols[idx]);
    final_vals.push_back(sorted_vals[idx]);
  }
  arma::umat locations(2, final_rows.size());
  locations.row(0) = arma::uvec(final_rows).t();
  locations.row(1) = arma::uvec(final_cols).t();
  arma::vec values(final_vals);

  return arma::sp_mat(locations, values, n1, n2, true, true);
}

arma::sp_mat elementwiseMultiplySp(const arma::sp_mat& A,
                                   const arma::sp_mat& B) {
  if ((A.n_rows != B.n_rows) || (A.n_cols != B.n_cols)) {
    Rcpp::stop("两个矩阵的维度必须相同");
  }

  std::vector<arma::uword> rowIndices;
  std::vector<arma::uword> colIndices;
  std::vector<double> values;

  for (arma::sp_mat::const_iterator it = A.begin(); it != A.end(); ++it) {
    arma::uword i = it.row();
    arma::uword j = it.col();
    double prod = (*it) * B(i, j);
    if (prod != 0.0) {
      rowIndices.push_back(i);
      colIndices.push_back(j);
      values.push_back(prod);
    }
  }

  arma::uvec arma_rows(rowIndices);
  arma::uvec arma_cols(colIndices);
  arma::vec arma_vals(values);

  arma::umat locations(2, arma_rows.n_elem);
  locations.row(0) = arma_rows.t();
  locations.row(1) = arma_cols.t();

  arma::sp_mat result(locations, arma_vals, A.n_rows, A.n_cols, true, true);

  return result;
}


// [[Rcpp::export]]
arma::sp_mat filter_SNN(arma::sp_mat & SNN,
                        std::vector<IntegerVector> & idx,
                        const std::vector<CharacterVector> & meta_list){
  IntegerVector dim(idx.size());
  int sum = 0;
  for(int i = 0; i < idx.size(); ++i){
    IntegerVector tmp = idx[i];
    sum += tmp.size();
    dim[i] = sum;
  }

  int n = meta_list.size();
  for(int i = 0; i < n; ++i){
    CharacterVector batch_i = meta_list[i];
    int start_i = (i == 0) ? 0 : dim[i-1];
    int end_i = dim[i] - 1;

    for(int j = 0; j < n; ++j){
      CharacterVector batch_j = meta_list[j];
      int start_j = (j == 0) ? 0 : dim[j-1];
      int end_j = dim[j] - 1;

      if(i != j){
        arma::sp_mat tmp = SNN.submat(start_i,start_j,end_i,end_j);
        arma::sp_mat filter_mat = generateMatchMatrixArma(batch_i,batch_j);
        arma::sp_mat res = elementwiseMultiplySp(tmp,filter_mat);
        SNN.submat(start_i,start_j,end_i,end_j) = res;
      }
    }
  }
  return SNN;
}


// [[Rcpp::export]]
List K_SNN(arma::sp_mat & SNN,
           std::vector<IntegerVector> & idx,
           const int & k,
           const double & lambda) {

  R2cpp_index(idx);
  const int n = idx.size();

  std::vector<double> f_list;
  f_list.reserve((size_t)n * (size_t)(n - 1) / 2);

  for (int i = 0; i < n; ++i) {
    IntegerVector idx_i = idx[i];
    int start_i = idx_i[0];
    int end_i   = idx_i[idx_i.size() - 1];
    arma::uword bi = (arma::uword)idx_i.size();

    SNN.submat(start_i, start_i, end_i, end_i).zeros();

    for (int j = i + 1; j < n; ++j) {
      IntegerVector idx_j = idx[j];
      int start_j = idx_j[0];
      int end_j   = idx_j[idx_j.size() - 1];
      arma::uword bj = (arma::uword)idx_j.size();

      arma::sp_mat sub_view = SNN.submat(start_i, start_j, end_i, end_j);

      std::vector< std::vector< std::pair<arma::uword,double> > > rows(bi);
      rows.shrink_to_fit();

      for (auto it = sub_view.begin(); it != sub_view.end(); ++it) {
        double v = *it;
        if (!std::isfinite(v) || v <= 0) continue;
        rows[it.row()].push_back({it.col(), v});
      }

      std::vector<arma::uword> out_r;
      std::vector<arma::uword> out_c;
      std::vector<double> out_v;
      out_r.reserve(sub_view.n_nonzero);
      out_c.reserve(sub_view.n_nonzero);
      out_v.reserve(sub_view.n_nonzero);

      double f_acc = 0.0;

      for (arma::uword r = 0; r < bi; ++r) {
        auto & vec = rows[r];
        if (vec.empty()) continue;

        if ((int)vec.size() > k) {
          std::nth_element(
            vec.begin(), vec.begin() + k, vec.end(),
            [](const auto &a, const auto &b){ return a.second > b.second; }
          );
          vec.resize(k);
        }

        double s = 0.0;
        for (auto & p : vec) s += p.second;
        if (s <= 0) continue;
        s += 1e-10;


        for (auto & p : vec) {
          double vn = p.second / s;
          if (vn <= 0) continue;
          out_r.push_back(r);
          out_c.push_back(p.first);
          out_v.push_back(vn);
          f_acc += vn;
        }
      }

      if (f_acc <= 0 || out_v.empty()) {
        SNN.submat(start_i, start_j, end_i, end_j).zeros();
        SNN.submat(start_j, start_i, end_j, end_i).zeros();
        f_list.push_back(0.0);
        continue;
      }

      for (double & vv : out_v) vv /= f_acc;
      f_list.push_back(f_acc);

      arma::umat loc(2, out_v.size());
      for (arma::uword t = 0; t < out_v.size(); ++t) {
        loc(0, t) = out_r[t];
        loc(1, t) = out_c[t];
      }
      arma::vec vals(out_v.size());
      for (arma::uword t = 0; t < out_v.size(); ++t) vals[t] = out_v[t];

      arma::sp_mat sub_snn(loc, vals, bi, bj, /*sort_locations=*/true);

      SNN.submat(start_i, start_j, end_i, end_j) = sub_snn;
      SNN.submat(start_j, start_i, end_j, end_i) = sub_snn.t();
    }
  }

  cpp2R_index(idx);

  double fmax = 0.0;
  for (double f : f_list) if (f > fmax) fmax = f;
  if (fmax > 0) SNN *= fmax;

  arma::sp_mat K = SNN;
  arma::vec B_vec = arma::vec(arma::sum(K, 1));

  arma::uvec nz = arma::find(B_vec != 0);
  double mm = 0.0;
  if (nz.n_elem > 0) {
    mm = arma::mean(B_vec.elem(nz)) * (1.0 - lambda);
  }

  return List::create(
    Named("K")     = K,
    Named("B")     = B_vec * lambda,
    Named("value") = mm
  );
}


// [[Rcpp::export]]
arma::vec silhouette_cpp(const arma::vec & labels,
                         const arma::mat & X){

  int n = labels.n_elem;
  arma::vec sil(n, fill::zeros);

  arma::mat dist_matrix = pairwiseEuclideanDistance(X,X);
  dist_matrix.diag().zeros();

  std::map<int, std::vector<int>> clusters;
  for (int i = 0; i < n; ++i) {
    clusters[labels(i)].push_back(i);
  }

  for (int i = 0; i < n; ++i) {
    int current_label = labels(i);
    double a = 0.0;
    double b = std::numeric_limits<double>::max();

    const std::vector<int>& same_cluster = clusters[current_label];
    int count_in_cluster = same_cluster.size() - 1;

    if (count_in_cluster > 0) {
      for (int j : same_cluster) {
        if (i != j) {
          a += dist_matrix(i, j);
        }
      }
      a /= count_in_cluster;
    }

    for (const auto& pair : clusters) {
      int other_label = pair.first;
      if (other_label == current_label) continue;

      const std::vector<int>& other_cluster = pair.second;
      double avg_dist = 0.0;

      for (int j : other_cluster) {
        avg_dist += dist_matrix(i, j);
      }
      avg_dist /= other_cluster.size();

      if (avg_dist < b) {
        b = avg_dist;
      }
    }

    sil(i) = (b - a) / std::max(a, b);
  }

  return sil;

}
