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
// Fast Pga for Msgl with L1-L2,1 norm penalty.
//
// [[Rcpp::export]]
List FPga_L121(SEXP XX, SEXP YY, SEXP B0, SEXP Gm, SEXP mi, SEXP mg, SEXP mc)
{
  //  const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
  const MatrixXd X(as<MatrixXd>(XX));
  const MatrixXd Y(as<MatrixXd>(YY));
  const MatrixXd Beta0(as<MatrixXd>(B0));
  const MatrixXi G_matr(as<MatrixXi>(Gm));
  const int iter_max = as<int>(mi);
  const int grp_max = as<int>(mg);
  const int coe_max = as<int>(mc);
  
  //  const int N(X.rows());
  const int p(X.cols());
  const int q(Y.cols());
  const int n_g(G_matr.rows());
  
  MatrixXd beta = Beta0;
  MatrixXd beta_old=beta;
  
  MatrixXd W = Beta0; // Intermaeiate parameter for beta in FISTA
  
  double v=pow(X.norm(), -2)*0.5;
  VectorXd tau = VectorXd::Zero(iter_max);
  VectorXd lambda = VectorXd::Zero(iter_max);
  VectorXd lambda_grpNormalized = VectorXd::Zero(iter_max);
  VectorXd rss = VectorXd::Zero(iter_max); // resisual sum square
  VectorXd rss_relative = VectorXd::Zero(iter_max); 
  
  int k=0; // iteration index
  double t=1; // intermediate parameter
  double t_old=1; 
  
  while (k < iter_max) 
  {
    checkUserInterrupt();
    beta_old=beta;
    t_old=t;
    MatrixXd z=W-2*v*X.adjoint() * (X * W-Y);
    Map<RowVectorXd> zv(z.data(), z.size());
    VectorXd zsort=zv.cwiseAbs();
    std::sort(zsort.data(),zsort.data()+zsort.size(), std::greater<double>());
    tau(k)=zsort(coe_max)/v;
    
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
    lambda_grpNormalized(k)=zbsortv(grp_max)/v;
    
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
    
    k++;
    if ((beta-beta_old).norm()<0.001){break;};
    t=(1+sqrt(1+4*pow(t_old,2)))/2;
    W=beta+(t_old-1)/t*(beta-beta_old);
    
  }
  return List::create(_["Beta"] = beta,
                      _["Rss"] = rss,
                      _["Tau"] = tau,
                      _["Lambda"] = lambda_grpNormalized,
                      _["iteration.time"] = k,
                      _["Rss_relative"] = rss_relative);
}

// Fast Pga for Msgl with L0-L2,0 norm penalty, .
// Version 1: tau & lambda determined by allowed maximum number of groups and single coefficients at each iteration step.
//
// [[Rcpp::export]]
List FPga_L020v1(SEXP XX, SEXP YY, SEXP B0, SEXP Gm, SEXP mi, SEXP mg, SEXP mc)
{
  const MatrixXd X(as<MatrixXd>(XX));
  const MatrixXd Y(as<MatrixXd>(YY));
  const MatrixXd Beta0(as<MatrixXd>(B0));
  const MatrixXi G_matr(as<MatrixXi>(Gm));
  const int iter_max = as<int>(mi);
  const int grp_max = as<int>(mg);
  const int coe_max = as<int>(mc);
  
  //  const int N(X.rows());
  const int p(X.cols());
  const int q(Y.cols());
  const int n_g(G_matr.rows());
  
  MatrixXd beta = Beta0;
  MatrixXd beta_old=beta;
  
  MatrixXd W = Beta0; // Intermaeiate parameter for beta in FISTA
  
  double v=pow(X.norm(), -2)*0.5;
  VectorXd tau = VectorXd::Zero(iter_max);
  VectorXd lambda = VectorXd::Zero(iter_max);
  VectorXd lambda_grpNormalized = VectorXd::Zero(iter_max);
  VectorXd rss = VectorXd::Zero(iter_max); // resisual sum square
  VectorXd rss_relative = VectorXd::Zero(iter_max); 
  
  int k=0; // iteration index
  double t=1; // intermediate parameter
  double t_old=1; 
  
  while (k < iter_max) 
  {
    beta_old=beta;
    t_old=t;
    MatrixXd z=W-2*v*X.adjoint() * (X * W-Y);
    Map<RowVectorXd> zv(z.data(), z.size());
    VectorXd zsort=zv.cwiseAbs();
    std::sort(zsort.data(),zsort.data()+zsort.size(), std::greater<double>());
    tau(k)=pow(zsort(coe_max), 2)/(2*v);
    
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
    lambda(k)=res_sortv(grp_max);
    
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
    
    k++;
    if ((beta-beta_old).norm()<0.001){break;};
    t=(1+sqrt(1+4*pow(t_old,2)))/2;
    W=beta+(t_old-1)/t*(beta-beta_old);   
  }
  return List::create(_["Beta"] = beta,
                      _["Rss"] = rss,
                      _["Tau"] = tau,
                      _["Lambda"] = lambda_grpNormalized,
                      _["iteration.time"] = k,
                      _["Rss_relative"] = rss_relative);
}

// Fast Pga for Msgl with L0-L2,0 norm penalty, .
// Version 2: tau & lambda determined by allowed maximum number of groups and single coefficients at the first iteration step, 
// but decreased with a constant scale.
//
// [[Rcpp::export]]
List FPga_L020v2(SEXP XX, SEXP YY, SEXP B0, SEXP Gm, SEXP mi, SEXP mg, SEXP mc, double minlambda=1e-5, double rlambda=0.98, double mintau=1e-5, double rtau=0.98)
{
  const MatrixXd X(as<MatrixXd>(XX));
  const MatrixXd Y(as<MatrixXd>(YY));
  const MatrixXd Beta0(as<MatrixXd>(B0));
  const MatrixXi G_matr(as<MatrixXi>(Gm));
  const int iter_max = as<int>(mi);
  const int grp_max = as<int>(mg);
  const int coe_max = as<int>(mc);
  
  //  const int N(X.rows());
  const int p(X.cols());
  const int q(Y.cols());
  const int n_g(G_matr.rows());
  
  MatrixXd beta = Beta0;
  MatrixXd beta_old=beta;
  
  MatrixXd W = Beta0; // Intermaeiate parameter for beta in FISTA
  
  double v=pow(X.norm(), -2)*0.5;
  VectorXd tau = VectorXd::Zero(iter_max);
  VectorXd lambda = VectorXd::Zero(iter_max);
  VectorXd lambda_grpNormalized = VectorXd::Zero(iter_max);
  VectorXd rss = VectorXd::Zero(iter_max); // resisual sum square
  VectorXd rss_relative = VectorXd::Zero(iter_max); 
  
  int k=0; // iteration index
  double t=1; // intermediate parameter
  double t_old=1; 
  
  while (k < iter_max) 
  {
    beta_old=beta;
    t_old=t;
    MatrixXd z=W-2*v*X.adjoint() * (X * W-Y);
    if(k==0)
    {
      Map<RowVectorXd> zv(z.data(), z.size());
      VectorXd zsort=zv.cwiseAbs();
      std::sort(zsort.data(),zsort.data()+zsort.size(), std::greater<double>());
      tau(k)=pow(zsort(coe_max), 2)/(2*v);
    } else {
      if(tau(k-1)>mintau)
      {
        tau(k)=tau(k-1)*rtau;
      } else {
        tau(k)=mintau;
      }
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
    
    if(k==0)
    {
      MatrixXd res_sort=res;
      Map<RowVectorXd> res_sortv(res_sort.data(), res_sort.size());
      std::sort(res_sortv.data(), res_sortv.data()+res_sortv.size(), std::greater<double>());
      lambda(k)=res_sortv(grp_max);
    } else {
      if(lambda(k-1)>minlambda)
      {
        lambda(k)=lambda(k-1)*rlambda;
      } else {
        lambda(k)=minlambda;
      }
    }
    
    //lambda[k]=sort(z.b.g.norm, decreasing = TRUE)[grp.max+1]/v
    
    
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
    
    k++;
    if ((beta-beta_old).norm()<0.001){break;};
    t=(1+sqrt(1+4*pow(t_old,2)))/2;
    W=beta+(t_old-1)/t*(beta-beta_old);
  }
  return List::create(_["Beta"] = beta,
                      _["Rss"] = rss,
                      _["Tau"] = tau,
                      _["Lambda"] = lambda_grpNormalized,
                      _["iteration.time"] = k,
                      _["Rss_relative"] = rss_relative);
}
