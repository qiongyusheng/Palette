#define ARMA_64BIT_WORD 1
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <unordered_set>
#include <set>
#include <vector>
using namespace Rcpp;
using namespace arma;

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

// [[Rcpp::export]]  //尝试是否可以直接输入sys函数
arma::mat calAv_cpp(const arma::vec & v,
                    const arma::mat & X){
  arma::mat Av = X*v;
  return(Av);
}


// arma::mat GetS(const arma::mat & X,
//                const arma::sp_mat & K,
//                const arma::vec & DVec,
//                double & lambda){
  
//   arma::mat D = arma::zeros<arma::mat>(K.n_cols,K.n_cols);
//   D.diag() = DVec;
//   return(X*(K-lambda*D)*X.t());
// }


// [[Rcpp::export]]
arma::mat GetS(const arma::mat & X,
               const arma::sp_mat & K,
               const arma::vec & DVec,
               double lambda,
               size_t batch_size = 100000) {
  
  // 创建修改后的K矩阵
  arma::sp_mat K_modified = K;
  
  // 优化对角线元素修改
  for (arma::uword i = 0; i < DVec.n_elem; i++) {
    K_modified.diag()[i] -= lambda * DVec(i);
  }
  
  // 自动检测是否需要批处理 - 当K的行数超过batch_size时
  const bool use_batch = (K_modified.n_rows > batch_size);
  
  arma::mat result(X.n_rows, X.n_rows, arma::fill::zeros);
  
  if (use_batch) {
    size_t n_rows = K_modified.n_rows;
    size_t num_batches = (n_rows + batch_size - 1) / batch_size;
    
    for (size_t b = 0; b < num_batches; ++b) {
      size_t start_row = b * batch_size;
      size_t end_row = std::min(start_row + batch_size, n_rows);
      
      // 提取当前批次的K行和对应的X列
      arma::sp_mat K_batch = K_modified.rows(start_row, end_row - 1);
      arma::mat X_batch = X.cols(start_row, end_row - 1);
      
      // 计算当前批次的贡献
      arma::mat tmp = K_batch * X.t();
      result += X_batch * tmp;
    }
  } else {
    // 标准方法
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
  // 创建修改后的K矩阵
  arma::sp_mat K_modified = K -lambda * B;
  
  // 自动检测是否需要批处理 - 当K的行数超过batch_size时
  const bool use_batch = (K_modified.n_rows > batch_size);
  
  arma::mat result(X.n_rows, X.n_rows, arma::fill::zeros);
  
  if (use_batch) {
    size_t n_rows = K_modified.n_rows;
    size_t num_batches = (n_rows + batch_size - 1) / batch_size;
    
    for (size_t b = 0; b < num_batches; ++b) {
      size_t start_row = b * batch_size;
      size_t end_row = std::min(start_row + batch_size, n_rows);
      
      // 提取当前批次的K行和对应的X列
      arma::sp_mat K_batch = K_modified.rows(start_row, end_row - 1);
      arma::mat X_batch = X.cols(start_row, end_row - 1);
      
      // 计算当前批次的贡献
      arma::mat tmp = K_batch * X.t();
      result += X_batch * tmp;
    }
  } else {
    // 标准方法
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
  const size_t n_rows = X.n_rows;  // 低维度（~200）
  const size_t n_cols = X.n_cols;  // 高维度（1e6+）
  arma::mat result(n_rows, n_rows, fill::zeros);
  
  // 优化策略：面向列的分块处理
  const size_t num_batches = (n_cols + batch_size - 1) / batch_size;
  
  // 预计算X转置（行数少时转置代价低）
  const arma::mat Xt = X.t();  // 维度变为（n_cols × n_rows）
  
  // 并行分块处理
#pragma omp parallel for schedule(dynamic)
  for (size_t b = 0; b < num_batches; ++b) {
    // 当前批次范围
    const size_t start = b * batch_size;
    const size_t end = std::min((b+1)*batch_size, n_cols) - 1;
    
    // 提取当前批次的K和B行块
    const arma::sp_mat K_block = K.rows(start, end);
    const arma::sp_mat B_block = B.rows(start, end);
    
    // 分步计算避免中间矩阵
    arma::mat KX = K_block * Xt;  // (batch_size × n_rows)
    arma::mat BX = B_block * Xt;
    arma::mat KBX = KX - BX;
    
    // 获取当前X列块（内存连续访问优化）
    const arma::mat X_block = X.cols(start, end);
    
    // 累加到结果矩阵（线程安全操作）
#pragma omp critical
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
      double dist = norm(A.col(i) - B.col(j), 2); // Calculate Euclidean distance using norm with p=2
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
  
  //  计算得到矩阵和各矩阵的标准化scaler
  arma::sp_mat out1(sim.n_rows, sim.n_cols);
  
  R2cpp_index(idx1);
  R2cpp_index(idx2);
  
  //对行数据计算，找到每组对应的最大值，并随机连接
  int nr = idx1.size();
  for(int i = 0; i < nr; ++i){
    
    arma::rowvec row = groupsim.row(i);
    // 检查该行是否全为0
    if (arma::accu(row != 0) > 0) {
      // 找到该行最大值所在的列索引
      int cross_group = arma::index_max(row);
      arma::uvec r_idx = Rcpp::as<arma::uvec>(idx1[i]);
      arma::uvec c_idx = Rcpp::as<arma::uvec>(idx2[cross_group]);
      
      // 对于不连续的行的每一行,在这些不连续的列上随机取一个列,将该值加1
      for (arma::uword ii = 0; ii < r_idx.n_elem; ++ii) {
        arma::uword row_idx = r_idx(ii);
        arma::uvec col_indices = c_idx;
        arma::uvec col_idx = arma::shuffle(col_indices);
        out1(row_idx, col_idx(0)) += 1;
      }
      
      // 对于不连续的列的每一列,在这些不连续的行上随机取一行,将该值加1
      for (arma::uword jj = 0; jj < c_idx.n_elem; ++jj) {
        arma::uword col_idx = c_idx(jj);
        arma::uvec row_indices = r_idx;
        arma::uvec row_idx = arma::shuffle(row_indices);
        out1(row_idx(0), col_idx) += 1;
      }
    }
  }
  
  //double f1 = arma::norm(out1, "fro");
  int f1 = arma::accu(out1);
  if(f1 != 0){
    out1 /= f1;
  }
  
  arma::sp_mat out2(sim.n_rows, sim.n_cols);
  int nc = idx2.size();
  for(int i = 0; i < nc; ++i){
    
    arma::colvec col = groupsim.col(i);
    // 检查该行是否全为0
    if (arma::accu(col != 0) > 0) {
      // 找到该行最大值所在的列索引
      int cross_group = arma::index_max(col);
      arma::uvec r_idx = Rcpp::as<arma::uvec>(idx1[cross_group]);
      arma::uvec c_idx = Rcpp::as<arma::uvec>(idx2[i]);
      
      // 对于不连续的行的每一行,在这些不连续的列上随机取一个列,将该值加1
      for (arma::uword ii = 0; ii < r_idx.n_elem; ++ii) {
        arma::uword row_idx = r_idx(ii);
        arma::uvec col_indices = c_idx;
        arma::uvec col_idx = arma::shuffle(col_indices);
        out2(row_idx, col_idx(0)) += 1;
      }
      
      // 对于不连续的列的每一列,在这些不连续的行上随机取一行,将该值加1
      for (arma::uword jj = 0; jj < c_idx.n_elem; ++jj) {
        arma::uword col_idx = c_idx(jj);
        arma::uvec row_indices = r_idx;
        arma::uvec row_idx = arma::shuffle(row_indices);
        out2(row_idx(0), col_idx) += 1;
      }
      
    }
    
  }
  
  //double f2 = arma::norm(out2, "fro");
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
  /*
  * X：两个批次数据的各组相似中位数
  * row_mat：X的行表示的那个批次的数据各个组之间的相似中位数
  * col_mat：X的列表示的那个批次的数据各个组之间的相似中位数
  * sim_range：如果两个组相似，那么相似系数不应该低于组自相似值的0.05
  * relative_sim_cutoff：相似系数的皮尔森值
  * 
  */
  
  
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
  //arma::mat X_remain = X%out;
  
  return(X%out);//考虑输出什么
}

// [[Rcpp::export]]
List sKernel_norm(const std::vector<int> & N_list,
                       const int & N,
                       const int & clust,
                       std::vector<std::vector<IntegerVector>> & idx){
 //N_list:每批次细胞个数
 //N：细胞总数
 //clust：总细胞类型个数
 //idx：每批次每种细胞类型索引
  
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
  
  //对每个细胞类型计算kernel
  for(int i = 0; i < clust; ++i){
    arma::sp_mat tmp_out(N,N);
    IntegerVector scaler;//每个细胞类型的缩放因子
    for(int b1 = 0; b1 < (n_batch-1); ++b1){
      
      std::vector<IntegerVector> idx1 = idx[b1];//取b1批次的全部细胞类型的idx
      R2cpp_index(idx1);
      
      IntegerVector tmp = idx1[i];//取b1批次的细胞类型i的idx
      int tmpi = tmp[0];
      //保证b1批次中有该类
      if(tmpi >= 0){
        for(int b2 = b1+1; b2 < n_batch; ++b2){
          
          std::vector<IntegerVector> idx2 = idx[b2];
          R2cpp_index(idx2);
          IntegerVector tmpp = idx2[i];
          int tmpj = tmpp[0];
          //保证b2批次中有该类
          if(tmpj >= 0){
            //取出两批次对应的子矩阵
            arma::sp_mat sub_mat;
            if(b1 == 0){
              sub_mat = tmp_out.submat(0,dim[b2-1],dim[b1]-1,dim[b2]-1);
            }else{
              sub_mat = tmp_out.submat(dim[b1-1],dim[b2-1],dim[b1]-1,dim[b2]-1);
            }
            
            //将对应元素赋值为1
            //double s = 1 - tmp_sim(b1,b2);
             for (size_t rowmat = 0; rowmat < tmp.size(); ++rowmat) {
               int row = tmp[rowmat];
              for (size_t colmat = 0; colmat < tmpp.size(); ++colmat) {
                int col = tmpp[colmat];
                // sub_mat(row, col) = s*1000;
                sub_mat(row, col) = 1;
              }
            }
            
            // int tt = arma::accu(sub_mat);
            int tt = tmp.size() * tmpp.size();
            if(tt != 0){
              sub_mat /= tt;
            }
            scaler.push_back(tt);
            
            //放回原矩阵
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
  
  //int ts = max(scaler_global);
  //out *= ts;
  
  IntegerVector scaler_global2;
  arma::sp_mat out2(N,N);
  //对每个细胞类型计算kernel
  for(int i = 0; i < clust; ++i){
    arma::sp_mat tmp_out(N,N);
    IntegerVector scaler;
    for(int b1 = 0; b1 < n_batch; ++b1){
      
      std::vector<IntegerVector> idx1 = idx[b1];
      R2cpp_index(idx1);
      
      IntegerVector tmp = idx1[i];
      int tmpi = tmp[0];
      //保证b1批次中有该类
      if(tmpi >= 0){
        
        arma::sp_mat sub_mat;
        if(b1 == 0){
          sub_mat = tmp_out.submat(0,0,dim[b1]-1,dim[b1]-1);
        }else{
          sub_mat = tmp_out.submat(dim[b1-1],dim[b1-1],dim[b1]-1,dim[b1]-1);
        }
        //将对应元素赋值为1
        for (size_t rowmat = 0; rowmat < tmp.size(); ++rowmat) {
          int row = tmp[rowmat];
          for (size_t colmat = 0; colmat < tmp.size(); ++colmat) {
            int col = tmp[colmat];
            sub_mat(row, col) = 1;
          }
        }
       
        // int tt = arma::accu(sub_mat);
        int tt = tmp.size() * tmp.size();
        if(tt != 0){
          sub_mat /= tt;
        }
        scaler.push_back(tt);
        
        //放回原矩阵
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
      out(idxD[i], idxD[i]) = BD[idxD[i]];  // 1-based 转为 0-based
    }
  }
    
  out *= ts;
  out2 *= ts;
  return List::create(Named("K") = out, Named("B") = out2);
  
}

arma::mat scaler_assign(const arma::mat & A,
                        const arma::mat & B,
                        const arma::vec & index){
  //A：预测数据
  //B：真实数据
  //index：真实数据的索引
  arma::rowvec la = arma::sqrt(arma::sum(arma::square(A), 0));//计算预测数据的长度
  arma::mat sub_B = B.cols(arma::conv_to<arma::uvec>::from(index));//取真实数据的子集
  arma::rowvec lb = arma::sqrt(arma::sum(arma::square(sub_B), 0));//计算子集的长度
  
  arma::uvec sorted_indices = arma::sort_index(la,"descend");//记录la从大到小排序的索引
  int n = la.size();
  int len = lb.size();
  arma::rowvec sampled(n);
  for(int i = 0; i < n; ++i){
    int idx = arma::randi<arma::uvec>(1, arma::distr_param(0, len - 1))(0);
    
    // Assign the sampled value to the result vector
    sampled(i) = lb(idx);
  }
  arma::rowvec sorted_sample = arma::sort(sampled,"descend");
  arma::rowvec scaler(n,arma::fill::zeros);
  for(int i = 0; i < n; ++i){
    scaler(sorted_indices(i)) = sorted_sample(i);
  }
  
  arma::mat out = A*arma::diagmat(1/la)*arma::diagmat(scaler);
  return(out);
  
}

arma::vec get_index(const arma::mat & A,
                    const arma::mat & B,
                    const int & scaler_k){
  int n = A.n_cols;
  arma::mat dist = pairwiseEuclideanDistance(A,B);
  arma::vec scale_index(n*scaler_k);//创建一个向量记录索引
  int pos = 0;//记录位置
  for(arma::uword i = 0; i < n; ++i){
    arma::rowvec rowi = dist.row(i); 
    arma::uvec new_idx = arma::sort_index(rowi,"ascend");// Partially sort the row
    scale_index.subvec(pos,pos+scaler_k-1) = arma::conv_to<arma::vec>::from(new_idx.subvec(0,scaler_k-1));//new
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
  
  arma::vec scale_index(X.n_cols*k);//创建一个向量记录索引 new
  int pos = 0;//记录位置 new
  
  arma::mat subdsit(n, k);
  for(arma::uword i = 0; i < n; ++i){
    arma::rowvec rowi = dist.row(i); 
    std::partial_sort(rowi.begin(), rowi.begin() + k, rowi.end());// Partially sort the row
    subdsit.row(i) = rowi.subvec(0, k - 1);
    
    arma::rowvec rowii = dist.row(i);
    arma::uvec top_k_indices = arma::sort_index(rowii, "ascend");
    if(do_scaler == true){
      scale_index.subvec(pos,pos+k-1) = arma::conv_to<arma::vec>::from(top_k_indices.subvec(0,k- 1));//new
      pos += k;//new
    }
    //scale_index.subvec(pos,pos+k-1) = arma::conv_to<arma::vec>::from(top_k_indices.subvec(0,k- 1));//new
    //pos += k;//new
  }
  
  arma::mat subdsit1 = exp(-subdsit);
  arma::vec row_sum_w = sum(subdsit1, 1)+1e-10;
  //arma::vec row_sum_w = sum(subdsit1, 0)+1e-10;
  
  //arma::mat weight_mat = subdsit1/arma::diagmat(1 / row_sum_w);
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
      // 使用arma::randn生成标准正态分布的随机数
      arma::colvec samples = arma::randn<arma::colvec>(dim);
      arma::colvec tmp_var = var.col(i);
      arma::colvec rand = samples % tmp_var;
      inferMat.col(i) = tmp + rand;
    }
  }else{
    inferMat += exp;
  }
  
  if(do_scaler){
    // 计算scaler index
    arma::mat input_Z = Z.cols(arma::conv_to<arma::uvec>::from(scale_index));
    arma::vec idx = get_index(input_Z,Z,scaler_k);
    Function unique_r("unique");
    arma::vec sub_idx = Rcpp::as<arma::vec>(unique_r(idx));
    arma::mat Mat = scaler_assign(inferMat,T,sub_idx);
    inferMat = Mat;
  }
  
  //return(M);
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
  
  arma::vec scale_index(X[0].n_cols*k);//创建一个向量记录索引 new
  int pos = 0;//记录位置 new
  
  for(arma::uword i = 0; i < n; ++i){
    arma::rowvec rowi = dist.row(i); 
    std::partial_sort(rowi.begin(), rowi.begin() + k, rowi.end()); // Partially sort the row
    subdsit.row(i) = rowi.subvec(0, k - 1); 
    
    arma::rowvec rowii = dist.row(i);
    arma::uvec top_k_indices = arma::sort_index(rowii, "ascend");
    if(do_scaler == true){
      scale_index.subvec(pos,pos+k-1) = arma::conv_to<arma::vec>::from(top_k_indices.subvec(0,k- 1));//new
      pos += k;//new
    }
    //scale_index.subvec(pos,pos+k-1) = arma::conv_to<arma::vec>::from(top_k_indices.subvec(0,k- 1));//new
    //pos += k;//new
  }
  
  arma::mat subdsit1 = exp(-subdsit);
  arma::vec row_sum_w = sum(subdsit1, 1)+1e-10;
  //arma::mat weight_mat = subdsit1/arma::diagmat(1 / row_sum_w);
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
      // 使用arma::randn生成标准正态分布的随机数
      arma::colvec samples = arma::randn<arma::colvec>(dim);
      arma::colvec tmp_var = var.col(i);
      arma::colvec rand = samples % tmp_var;
      inferMat.col(i) = tmp + rand;
    }
  }else{
    inferMat += exp;
  }
  
  if(do_scaler){
    // 计算scaler index
    arma::mat input_Z = Z.cols(arma::conv_to<arma::uvec>::from(scale_index));
    arma::vec idx = get_index(input_Z,Z,scaler_k);
    Function unique_r("unique");
    arma::vec sub_idx = Rcpp::as<arma::vec>(unique_r(idx));
    arma::mat Mat = scaler_assign(inferMat,T,sub_idx);
    inferMat = Mat;
  }
  
  // return(M);
  return List::create(Named("M") = M,
                      Named("infer") = inferMat);
}

void replace_with_vector(arma::mat& M, const arma::uvec& V) {
  
/*
  if (V.n_elem < M.n_elem) {
    // 向量长度不足,抛出异常
    Rcpp::stop("Vector length is not sufficient to replace all matrix elements");
  }
*/
  
  //arma::uword k = 0;
  for (arma::uword i = 0; i < M.n_rows; i++) {
    for (arma::uword j = 0; j < M.n_cols; j++) {
      int tmp = M(i, j);
      M(i, j) = V(tmp);
      //M(i, j) = V(k++);
      //if (k >= V.n_elem) {
       // // 向量用完,重新从头开始
        //k = 0;
      //}
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
    std::partial_sort(rowi.begin(), rowi.begin() + k, rowi.end()); // Partially sort the row
    subdsit.row(i) = rowi.subvec(0, k - 1); 
  }
  
  arma::mat subdsit1 = exp(-subdsit);
  arma::vec row_sum_w = sum(subdsit1, 1)+1e-10;
  //arma::mat weight_mat = subdsit1/arma::diagmat(1 / row_sum_w);
  arma::mat weight_mat = arma::diagmat(1 / row_sum_w)*subdsit1;

  arma::mat indices(n, k);
  for(arma::uword i = 0; i < n; ++i){
    arma::uvec tmp =  arma::sort_index(dist.row(i));
    tmp = tmp.head(k);
    indices.row(i) = arma::conv_to<arma::rowvec>::from(tmp.t());
    //indices.row(i) = arma::conv_to<arma::irowvec>::from(tmp.t());
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
      // 使用arma::randn生成标准正态分布的随机数
      arma::colvec samples = arma::randn<arma::colvec>(dim);
      arma::colvec tmp_var = var.col(i);
      arma::colvec rand = samples % tmp_var;
      inferMat.col(i) = tmp + rand;
    }
  }else{
    inferMat += exp;
  }

  if(do_scaler){
    // 计算scaler index
    arma::rowvec scale_index = arma::vectorise(indices).t();
    arma::mat input_Z = Z.cols(arma::conv_to<arma::uvec>::from(scale_index));
    arma::vec idx = get_index(input_Z,Z,scaler_k);
    Function unique_r("unique");
    arma::vec sub_idx = Rcpp::as<arma::vec>(unique_r(idx));
    arma::mat Mat = scaler_assign(inferMat,T,sub_idx);
    inferMat = Mat;
  }
  // return(M);
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
  
  const arma::uword n = X[0].n_cols;//待预测细胞数
  int l = Y.size(); //路径数目
  Function cor_r("cor"), cumsum_r("cumsum");
  
  IntegerVector num(l+1);
  num[0] = 0;
  for(int i = 0; i < l; ++i){
    num[i+1] = Y[i][0][0].n_cols;
  }
  IntegerVector acc_num = cumsum_r(num);
  arma::mat dist(n,acc_num[l-1]); //定义一个矩阵，放置所有的相似性
  //计算距离并填充dist
  for(int i = 0; i < l; ++i){
    arma::mat sub_dist(X[i].n_cols,Y[i][0][0].n_cols);
    if(L2 == true){
      sub_dist = 1-Rcpp::as<arma::mat>(cor_r(X[i],Y[i][0][0]));
    }else{
      sub_dist = 1-X[i].t()*Y[i][0][0];
    }
    dist.submat(0,num[i],n-1,num[i+1]-1) = sub_dist;
  }
  
  arma::mat subdsit(n, k);
  arma::sp_mat M(n,dist.n_cols);
  for(arma::uword i = 0; i < n; ++i){
    arma::rowvec rowi = dist.row(i); 
    std::partial_sort(rowi.begin(), rowi.begin() + k, rowi.end()); // Partially sort the row
    subdsit.row(i) = rowi.subvec(0, k - 1);
  }
  
  arma::mat subdsit1 = exp(-subdsit);
  arma::vec row_sum_w = sum(subdsit1, 1)+1e-10;
  //arma::mat weight_mat = subdsit1/arma::diagmat(1 / row_sum_w);
  arma::mat weight_mat = arma::diagmat(1 / row_sum_w)*subdsit1;
  
  arma::mat indices(n, k, arma::fill::ones);
  for(arma::uword i = 0; i < n; ++i){
    arma::uvec tmp =  arma::sort_index(dist.row(i));
    tmp = tmp.head(k);
    indices.row(i) = arma::conv_to<arma::rowvec>::from(tmp.t());
    //indices.row(i) = arma::conv_to<arma::irowvec>::from(tmp.t());
  }
  
  int ll = Y[0].size();//每条路径要走的步数
  std::vector<arma::uvec> IDX(l);
  
  int j = 0;
  while(j < (ll-1)){
    
    for(int i = 0; i < l; ++i){
      arma::mat tmpdist(Y[i][j][1].n_cols,Y[i][j+1][0].n_cols);
      if(L2 == true){
        tmpdist = 1-Y[i][j][1].t()*Y[i][j+1][0];
      }else{
        tmpdist =  pairwiseEuclideanDistance(Y[i][j][1],Y[i][j+1][0]);
      }
      arma::uvec minIndices = arma::index_min(tmpdist, 1);
      if(j != 0){
        arma::uvec tmp =  IDX[i];
        IDX[i] = minIndices.elem(tmp);
      }else{
        IDX[i] = minIndices;
      }
    }
    ++j;
  }
  
  //计算到达模态共有多少个细胞
  IntegerVector na(l+1);
  na[0] = 0;
  for(int i = 0; i < l; ++i){
    na[i+1] = Y[i][ll-1][1].n_cols;
  }
  IntegerVector acc_na = cumsum_r(na);
  //arma::mat W(Y[0][ll-1][1].n_rows,acc_na[l]);
  //调整索引
  for(int i = 0; i < l; ++i){
    IDX[i] += acc_na[i];
    //W.submat(0,acc_na[i],Y[0][ll-1][1].n_rows-1,acc_na[i+1]-1) = Y[i][ll-1][1];
  }
  
  //将不同路径的索引信息拼接起来
  arma::uvec idx(na[l]);
  arma::uword pos = 0;
  for(const arma::uvec& curr_idx : IDX){
    idx.subvec(pos,pos + curr_idx.n_elem - 1) = curr_idx;
    pos += curr_idx.n_elem;
  }
  
  replace_with_vector(indices,idx);
  //arma::mat out(W.n_rows,n);
  
  for(arma::uword i = 0; i < n; ++i){
    arma::uvec sorted_indices = arma::conv_to<arma::uvec>::from(indices.row(i));
    //out.col(i) = arma::mean(W.cols(sorted_indices), 1);
  }
  
  for(arma::uword i = 0; i < n; ++i){
    for(int j = 0; j < k; ++j){
      M(i,indices(i,j)) = weight_mat(i,j);
    }
  }
  
  arma::mat exp = Z * M.t();
  arma::mat inferMat(T.n_rows,n); 
  if(do_var == true){
    arma::mat D = pairwiseEuclideanDistance(exp,T);
    arma::mat var = Get_var(D,exp,T,var_k);
    int dim = T.n_rows;
    for(int i = 0; i < n; ++i){
      arma::colvec tmp = exp.col(i);
      // 使用arma::randn生成标准正态分布的随机数
      arma::colvec samples = arma::randn<arma::colvec>(dim);
      arma::colvec tmp_var = var.col(i);
      arma::colvec rand = samples % tmp_var;
      inferMat.col(i) = tmp + rand;
    }
  }else{
    inferMat += exp;
  }
  
  if(do_scaler == true){
    // 计算scaler index
    arma::rowvec scale_index = arma::vectorise(indices).t();
    arma::mat input_Z = Z.cols(arma::conv_to<arma::uvec>::from(scale_index));
    arma::vec idx_scale = get_index(input_Z,Z,scaler_k);
    Function unique_r("unique");
    arma::vec sub_idx = Rcpp::as<arma::vec>(unique_r(idx_scale));
    arma::mat Mat = scaler_assign(inferMat,T,sub_idx);
    inferMat = Mat;
  }
  
  // return(M);
  return List::create(Named("M") = M,
                      Named("infer") = inferMat);
  
}

arma::sp_mat generateMatchMatrixArma(const CharacterVector& batch1, 
                                     const CharacterVector& batch2) {
  int n1 = batch1.size();
  int n2 = batch2.size();
  if (n1 == 0 || n2 == 0) {
    return arma::sp_mat(n1, n2);
  }
  
  // 构建标签到索引的映射
  std::unordered_map<std::string, std::vector<int>> map1, map2;
  for (int i = 0; i < n1; ++i) {
    map1[as<std::string>(batch1[i])].push_back(i);
  }
  for (int j = 0; j < n2; ++j) {
    map2[as<std::string>(batch2[j])].push_back(j);
  }
  
  // 生成三元组数据
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
  
  // 按列主序排序
  std::vector<size_t> indices(sorted_rows.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
    if (sorted_cols[a] == sorted_cols[b]) {
      return sorted_rows[a] < sorted_rows[b];
    }
    return sorted_cols[a] < sorted_cols[b];
  });
  // 提取排序后的数据
  std::vector<arma::uword> final_rows, final_cols;
  std::vector<double> final_vals;
  for (auto idx : indices) {
    final_rows.push_back(sorted_rows[idx]);
    final_cols.push_back(sorted_cols[idx]);
    final_vals.push_back(sorted_vals[idx]);
  }
  // 转换为 Armadillo 格式
  arma::umat locations(2, final_rows.size());  // 2 行 nnz 列
  locations.row(0) = arma::uvec(final_rows).t(); // 第一行：行索引（转置为行向量）
  locations.row(1) = arma::uvec(final_cols).t(); // 第二行：列索引（转置为行向量）
  arma::vec values(final_vals);
  // 构建稀疏矩阵
  return arma::sp_mat(locations, values, n1, n2, true, true);
}

// [[Rcpp::export]]
arma::sp_mat elementwiseMultiplySp(const arma::sp_mat& A, 
                                   const arma::sp_mat& B) {
  // 确保两个矩阵的维度相同
  if ((A.n_rows != B.n_rows) || (A.n_cols != B.n_cols)) {
    Rcpp::stop("两个矩阵的维度必须相同");
  }
  
  std::vector<arma::uword> rowIndices;
  std::vector<arma::uword> colIndices;
  std::vector<double> values;
  
  // 遍历矩阵A的非零元素
  for (arma::sp_mat::const_iterator it = A.begin(); it != A.end(); ++it) {
    arma::uword i = it.row();
    arma::uword j = it.col();
    double prod = (*it) * B(i, j);
    // 只有乘积非零时才记录
    if (prod != 0.0) {
      rowIndices.push_back(i);
      colIndices.push_back(j);
      values.push_back(prod);
    }
  }
  
  // 将 std::vector 转换为 Armadillo 的向量
  arma::uvec arma_rows(rowIndices);
  arma::uvec arma_cols(colIndices);
  arma::vec arma_vals(values);
  
  // 构造 locations 矩阵，要求为 2 x nnz
  arma::umat locations(2, arma_rows.n_elem);
  locations.row(0) = arma_rows.t();  // 第一行为行索引（转置为行向量）
  locations.row(1) = arma_cols.t();  // 第二行为列索引（转置为行向量）
  
  // 使用 locations 构造稀疏矩阵
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
                   const double & lambda){
  
  R2cpp_index(idx);
  NumericVector f_list;
  int n = idx.size();
  for(int i = 0; i < n; ++i){
    IntegerVector idx_i = idx[i];
    int start_i = idx_i[0];
    int end_i = idx_i[idx_i.size()-1];
    for(int j = 0; j < n; ++j){
      IntegerVector idx_j = idx[j];
      int start_j = idx_j[0];
      int end_j = idx_j[idx_j.size()-1];
      if(i == j){
        SNN.submat(start_i,start_j,end_i,end_j) = arma::sp_mat(idx_i.size(),idx_j.size());
      }else{
        arma::mat tmp = arma::mat(SNN.submat(start_i,start_j,end_i,end_j));
        arma::mat tmpp = tmp;
        
        for(int ii = 0; ii < idx_i.size(); ++ii){
          arma::rowvec row_values = tmp.row(ii);

          if (arma::all(row_values == 0)) {
            continue;
          }

          if(row_values.has_nan()) {
            row_values.replace(arma::datum::nan, 0);
            // Rcpp::stop("detected NaN in SNN graph");
          }

          arma::uvec sort_indices = arma::sort_index(row_values, "descend");
          arma::uvec keep_indices(k); // 存储需要保留的索引
          // 检查前k个元素中有多少个大于0
          
          unsigned int num_nonzero = 0;
          for (unsigned int jj = 0; jj < std::min(sort_indices.n_elem, (arma::uword)k); ++jj) {
            if (row_values(sort_indices(jj)) > 0) {
              keep_indices(num_nonzero) = sort_indices(jj); // 直接赋值
              num_nonzero++;
            }
          }
          
          tmpp.row(ii).zeros();
          for (arma::uword jj = 0; jj < num_nonzero; ++jj) {
            tmpp(ii, keep_indices(jj)) = row_values(keep_indices(jj));
          }
          
        }
          //unsigned int n_keep = (sort_indices.n_elem < k) ? sort_indices.n_elem : k;
          // //arma::uvec keep_indices = sort_indices.head(n_keep);
          //for (arma::uword jj = 0; jj < n_keep; ++jj) {
         //   tmpp(ii, keep_indices(jj)) = row_values(keep_indices(jj));
        //  }
        //}
        
        arma::vec row_sum_w = sum(tmpp, 1)+1e-10;
        arma::mat sub_snn = arma::diagmat(1 / row_sum_w)*tmpp;
        //double f = arma::norm(sub_snn, "fro");
        double f = arma::accu(sub_snn);
        sub_snn /= f;
        f_list.push_back(f);
        SNN.submat(start_i,start_j,end_i,end_j) = arma::sp_mat(sub_snn);
        }
      }
    }
  
  cpp2R_index(idx);
  SNN *= Rcpp::max(f_list);
  arma::sp_mat K = (SNN+SNN.t())*0.5;
  arma::vec B_vec = sum(arma::mat(K), 1);
  arma::vec nonzero_vals = B_vec.elem(arma::find(B_vec != 0));
  double mm = arma::mean(nonzero_vals)*(1-lambda);
  
  return List::create(Named("K") = K,
                      Named("B") = B_vec*lambda,
                      Named("value") = mm);
  //cpp2R_index(idx);
  //arma::sp_mat out = (SNN+SNN.t())*0.5;
  //return(out);
}

// [[Rcpp::export]]
arma::vec silhouette_cpp(const arma::vec & labels,
                         const arma::mat & X){
  // cluster表示每一类细胞的索引
  // X表示数据的低维矩阵
  
  int n = labels.n_elem;
  arma::vec sil(n, fill::zeros);
  
  //计算距离矩阵
  arma::mat dist_matrix = pairwiseEuclideanDistance(X,X);
  dist_matrix.diag().zeros();

  // 创建一个 map，将每个簇的索引存储起来
  std::map<int, std::vector<int>> clusters;
  for (int i = 0; i < n; ++i) {
    clusters[labels(i)].push_back(i);
  }
  
  for (int i = 0; i < n; ++i) {
    int current_label = labels(i);
    double a = 0.0;
    double b = std::numeric_limits<double>::max();
    
    // 计算 a(i)
    const std::vector<int>& same_cluster = clusters[current_label];
    int count_in_cluster = same_cluster.size() - 1; // 排除自己
    
    if (count_in_cluster > 0) {
      for (int j : same_cluster) {
        if (i != j) {
          a += dist_matrix(i, j);
        }
      }
      a /= count_in_cluster;
    }
    
    // 计算 b(i)
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
