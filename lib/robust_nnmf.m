function [w,h] = robust_nnmf( x, K )
%% ======================================================================
% robust_nnmf( x, K )
%
% Interface to matlab's nnmf interface, with several options to make the
% method more robust.
%
% Created in 2016 by S. Gliske (sgliske@umich.edu)
% Released under the CC-BY-NC-4.0 License
% http://creativecommons.org/licenses/by-nc/4.0/
%% ======================================================================

  opt = statset('maxiter',10,'display','off');%MaxIter: Maximum number of iterations; create a structure of options
  [w,h] = nnmf(x,K,'rep',100,'opt',opt,'alg','mult');%Perform Initial NNMF with Multiplicative Algorithm
  
  if( isfinite(cond(h)) ) %check if condition number of h is finite (~Inf, i.e., h is not singular matrix); Singular matrices have a condition number of Inf. A singular matrix does not have an inverse; The rank of singular matrix is less than its order (the number of rows or columns). 
      opt = statset('maxiter',1000,'display','off');%Updates the options for NNMF with a higher maximum number of iterations (maxiter, 1000).
      [w,h] = nnmf(x,K,'w0',w,'h0',h,'opt',opt,'alg','als');%w0/h0: Initial value of w/h with the values obtained from the previous NNMF run; Performs NNMF again using the alternating least squares (ALS) algorithm

      if( ~any(w(:)) || ~any(h(:)) ) %checks if all elements in the array w/h are zero: Checks if any elements in w or h are non-zero. If all elements are zero, it evaluates to true.
          opt = statset('maxiter',10,'display','off');
          [w,h] = nnmf(x,K,'rep',1000,'opt',opt,'alg','mult');%Performs NNMF again with the multiplicative update algorithm and 1000 replicates.
      end
  else %If Condition Number is Not Finite  (i.e., h is singular)
      opt = statset('maxiter',10,'display','off');
      [w,h] = nnmf(x,K,'rep',1000,'opt',opt,'alg','mult');%Performs NNMF again with the multiplicative update algorithm and 1000 replicates
  end
end
