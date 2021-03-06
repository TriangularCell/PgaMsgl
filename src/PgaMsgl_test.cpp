// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <Eigen/Dense>

using namespace Rcpp;
using namespace Eigen;
using std::pow; 
using std::abs;
using std::sort;

// [[Rcpp::depends(RcppEigen)]]
//
// Pga for Msgl with L1-L2,1 norm penalty.
//
// [[Rcpp::export]]
List L121_test(SEXP XX, SEXP YY, SEXP B0, SEXP Gm, SEXP mi, SEXP mg, SEXP mc, SEXP Beta, SEXP cutoff)
{
  //  const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
  const MatrixXd X(as<MatrixXd>(XX));
  const MatrixXd Y(as<MatrixXd>(YY));
  const MatrixXd Beta0(as<MatrixXd>(B0));
  const MatrixXd Beta_true(as<MatrixXd>(Beta));
  const MatrixXi G_matr(as<MatrixXi>(Gm));
  const int iter_max = as<int>(mi);
  const int grp_max = as<int>(mg);
  const int coe_max = as<int>(mc);
  const double stop_cutoff = as<double>(cutoff);
  
  //  const int N(X.rows());
  const int p(X.cols());
  const int q(Y.cols());
  const int n_g(G_matr.rows());
  const int n_c=p*q;
  
  MatrixXd beta_old(MatrixXd(p, q).setOnes());
  MatrixXd beta = Beta0;
  double v=pow(X.norm(), -2)*0.5;
  VectorXd tau = VectorXd::Zero(iter_max);
  // VectorXd lambda = VectorXd::Zero(iter_max);
  VectorXd lambda_grpNormalized = VectorXd::Zero(iter_max);
  VectorXd rss = VectorXd::Zero(iter_max); // resisual sum square
  VectorXd rss_relative = VectorXd::Zero(iter_max);
  VectorXd error_relative = VectorXd::Zero(iter_max);
  
  int k=0; // iteration index
  
  while (k < iter_max) 
  {
    checkUserInterrupt();
    beta_old=beta;
    MatrixXd z=beta-2*v*X.adjoint() * (X * beta-Y);
    Map<RowVectorXd> zv(z.data(), z.size());
    VectorXd zsort=zv.cwiseAbs();
    std::sort(zsort.data(),zsort.data()+zsort.size(), std::greater<double>());
    if (coe_max<n_c)
      {
      tau(k)=zsort(coe_max)/v;
      } else {
        tau(k)=zsort(n_c-1)/v;
      }
    
    MatrixXd z_b(MatrixXd(p, q).setZero());
    
    for (int i=0; i<p; i++)
    {
      for (int j=0; j<q; j++)
      {
        z_b(i,j)=(abs(z(i,j))-v*tau(k)>0)?abs(z(i,j))-v*tau(k):0;
      }
    }
    
    VectorXd z_b_g_norm_grpNormalized = VectorXd::Zero(n_g);
    
    for (int g=0; g<n_g; g++)
    {
      int grp_row_start = G_matr(g,0);
      int grp_row_len = int(G_matr(g,1))-int(G_matr(g,0))+1; 
      int grp_col_start = G_matr(g,2); 
      int grp_col_len = int(G_matr(g,3))-int(G_matr(g,2))+1;
      int grp_size = grp_row_len*grp_col_len;
      
      z_b_g_norm_grpNormalized(g) = z_b.block(grp_row_start-1, grp_col_start-1, grp_row_len, grp_col_len).norm()/sqrt(grp_size);
    }
    
    //lambda[k]=sort(z.b.g.norm, decreasing = TRUE)[grp.max+1]/v
    MatrixXd zbsort=z_b_g_norm_grpNormalized;
    Map<RowVectorXd> zbsortv(zbsort.data(), zbsort.size());
    std::sort(zbsortv.data(), zbsortv.data()+zbsortv.size(), std::greater<double>());
    if (grp_max < n_g)
    {
      lambda_grpNormalized(k)=zbsortv(grp_max)/v;
    } else {
      lambda_grpNormalized(k)=zbsortv(n_g-1)/v;
    }
    
    
    for (int g=0; g<n_g; g++)
    {
      int grp_row_start = G_matr(g,0);
      int grp_row_len = int(G_matr(g,1))-int(G_matr(g,0))+1; 
      int grp_col_start = G_matr(g,2); 
      int grp_col_len = int(G_matr(g,3))-int(G_matr(g,2))+1;
      
      MatrixXd z_g=z.block(grp_row_start-1, grp_col_start-1, grp_row_len, grp_col_len);
      MatrixXd z_b_g=z_b.block(grp_row_start-1, grp_col_start-1, grp_row_len, grp_col_len);
      
      if (z_b_g_norm_grpNormalized(g)-v*lambda_grpNormalized(k) > 0)
      {
        MatrixXd z_g_s = z_g.array().sign();
        beta.block(grp_row_start-1, grp_col_start-1, grp_row_len, grp_col_len) << (1-v*lambda_grpNormalized(k)/z_b_g_norm_grpNormalized(g))*z_g_s.cwiseProduct(z_b_g); //a.cwiseProduct(b), dot product of matrix a and matrix b.
      } else {
        beta.block(grp_row_start-1, grp_col_start-1, grp_row_len, grp_col_len) << MatrixXd(grp_row_len, grp_col_len).setZero();
      }
    }
    
    rss(k)=pow((X*beta-Y).norm(), 2);
    rss_relative(k)=pow((X*beta-Y).norm(), 2)/pow(Y.norm(), 2);
    error_relative(k)=(beta-Beta_true).norm()/Beta_true.norm();
    
    k++;
    if ((beta-beta_old).norm()<stop_cutoff){break;};
  }
  return List::create(_["Beta"] = beta,
                      _["Rss"] = rss,
                      _["Tau"] = tau,
                      _["Lambda"] = lambda_grpNormalized,
                      _["iteration.time"] = k,
                      _["Rss_relative"] = rss_relative,
                      _["Error_relative"] = error_relative);
}

// Pga for Msgl with L0-L2,0 norm penalty, .
// Version 1: tau & lambda determined by allowed maximum number of groups and single coefficients at each iteration step.
//
// [[Rcpp::export]]
List L020v1_test(SEXP XX, SEXP YY, SEXP B0, SEXP Gm, SEXP mi, SEXP mg, SEXP mc, SEXP Beta, SEXP cutoff)
{
  const MatrixXd X(as<MatrixXd>(XX));
  const MatrixXd Y(as<MatrixXd>(YY));
  const MatrixXd Beta0(as<MatrixXd>(B0));
  const MatrixXd Beta_true(as<MatrixXd>(Beta));
  const MatrixXi G_matr(as<MatrixXi>(Gm));
  const int iter_max = as<int>(mi);
  const int grp_max = as<int>(mg);
  const int coe_max = as<int>(mc);
  const double stop_cutoff = as<double>(cutoff);
  
  //  const int N(X.rows());
  const int p(X.cols());
  const int q(Y.cols());
  const int n_g(G_matr.rows());
  const int n_c=p*q;
  
  MatrixXd beta_old(MatrixXd(p, q).setOnes());
  MatrixXd beta = Beta0;
  double v=pow(X.norm(), -2)*0.5;
  VectorXd tau = VectorXd::Zero(iter_max);
  VectorXd lambda = VectorXd::Zero(iter_max);
  // VectorXd lambda_grpNormalized = VectorXd::Zero(iter_max);
  VectorXd rss = VectorXd::Zero(iter_max); // resisual sum square
  VectorXd rss_relative = VectorXd::Zero(iter_max); 
  VectorXd error_relative = VectorXd::Zero(iter_max);
  
  int k=0; // iteration index
  
  while (k < iter_max) 
  {
    beta_old=beta;
    MatrixXd z=beta-2*v*X.adjoint() * (X * beta-Y);
    Map<RowVectorXd> zv(z.data(), z.size());
    VectorXd zsort=zv.cwiseAbs();
    std::sort(zsort.data(),zsort.data()+zsort.size(), std::greater<double>());
    if (coe_max<n_c)
    {
      tau(k)=pow(zsort(coe_max), 2)/(2*v);
    } else {
      tau(k)=pow(zsort(n_c-1), 2)/(2*v);
    }
    
    
    MatrixXd z_b(MatrixXd(p, q).setZero());
    
    for (int i=0; i<p; i++)
    {
      for (int j=0; j<q; j++)
      {
        z_b(i,j)=(abs(z(i,j))-sqrt(2*v*tau(k))>0)?z(i,j):0;
      }
    }
    
    //VectorXd z_b_g_norm_grpNormalized = VectorXd::Zero(n_g);
    VectorXd z_b_g_norm = VectorXd::Zero(n_g);
    VectorXd z_b_g_L0 = VectorXd::Zero(n_g);
    VectorXd res = VectorXd::Zero(n_g);
    
    for (int g=0; g<n_g; g++)
    {
      int grp_row_start = G_matr(g,0);
      int grp_row_len = int(G_matr(g,1))-int(G_matr(g,0))+1; 
      int grp_col_start = G_matr(g,2); 
      int grp_col_len = int(G_matr(g,3))-int(G_matr(g,2))+1;
      //int grp_size = grp_row_len*grp_col_len;
      
      //z_b_g_norm_grpNormalized(g) = z_b.block(grp_row_start-1, grp_col_start-1, grp_row_len, grp_col_len).norm()/sqrt(grp_size);
      MatrixXd z_b_g=z_b.block(grp_row_start-1, grp_col_start-1, grp_row_len, grp_col_len);
      z_b_g_norm(g) = z_b_g.norm();
      z_b_g_L0(g)=(z_b_g.array() != 0).count();
      res(g)=pow(z_b_g_norm(g), 2)/(2*v)-(tau(k)*z_b_g_L0(g));
    }
    
    //lambda[k]=sort(z.b.g.norm, decreasing = TRUE)[grp.max+1]/v
    MatrixXd res_sort=res;
    Map<RowVectorXd> res_sortv(res_sort.data(), res_sort.size());
    std::sort(res_sortv.data(), res_sortv.data()+res_sortv.size(), std::greater<double>());
    if (grp_max < n_g)
    {
      lambda(k)=res_sortv(grp_max);
    } else { 
      lambda(k)=res_sortv(n_g-1);
    }
    
    for (int g=0; g<n_g; g++)
    {
      int grp_row_start = G_matr(g,0);
      int grp_row_len = int(G_matr(g,1))-int(G_matr(g,0))+1; 
      int grp_col_start = G_matr(g,2); 
      int grp_col_len = int(G_matr(g,3))-int(G_matr(g,2))+1;
      
      //MatrixXd z_g=z.block(grp_row_start-1, grp_col_start-1, grp_row_len, grp_col_len);
      MatrixXd z_b_g=z_b.block(grp_row_start-1, grp_col_start-1, grp_row_len, grp_col_len);
      
      if (pow(z_b_g_norm(g), 2)-2*v*(lambda(k)+tau(k)*z_b_g_L0(g)) > 1e-15)
      {
        beta.block(grp_row_start-1, grp_col_start-1, grp_row_len, grp_col_len) << z_b_g;
      } else {
        beta.block(grp_row_start-1, grp_col_start-1, grp_row_len, grp_col_len) << MatrixXd(grp_row_len, grp_col_len).setZero();
      }
    }
    
    rss(k)=pow((X*beta-Y).norm(), 2);
    rss_relative(k)=pow((X*beta-Y).norm(), 2)/pow(Y.norm(), 2);
    error_relative(k)=(beta-Beta_true).norm()/Beta_true.norm();
    
    k++;
    if ((beta-beta_old).norm()<stop_cutoff){break;};
  }
  return List::create(_["Beta"] = beta,
                      _["Rss"] = rss,
                      _["Tau"] = tau,
                      _["Lambda"] = lambda,
                      _["iteration.time"] = k,
                      _["Rss_relative"] = rss_relative,
                      _["Error_relative"] = error_relative);
}