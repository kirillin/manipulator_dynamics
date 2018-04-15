// Copyright (C) 2009 - Michael Baudin
// Copyright (C) 2009 - John Burkardt
// This file must be used under the terms of the GNU LGPL license.

function b = nchoosek( n , k )
  //   Returns the binomial number (n,k).
  //
  // Calling Sequence
  //   b = nchoosek ( n , k )
  //
  // Description
  //   Computes the number of k-element subsets of an n-element set.
  //   It is mathematically defined by n!/(k! (n-k)! )
  //   which is equivalent to ( n*(n-1)*...*(n-k+1) ) / ( k*(k-1)*...*1 ).
  //   The implementation, though, uses a robust floating point implementation.
  //
  //   If n<0, k<0 or k>n, then an error is generated.
  //
  // Note about floating point accuracy.
  //   We have factorial(170) ~ 7.e306. 
  //   So n=170 is the greatest integer for which n! can be computed.
  //   If the naive formula b = n!/(k! (n-k)! )
  //   was used directly, the maximum value n for which 
  //   the binomial can be computed is n = 170.
  //   But binomial(171,1) = 1, so there is no reason
  //   to prevent the computation of 1 because an intermediate
  //   result is greater than 1.e308.
  //   This is why the gammaln function, combined with exp, is used instead.
  //
  // Note about rounding for integer inputs.
  //   If n = 4 and k = 1, the gammaln and exp functions 
  //   are accurate only to 1ulp. This leads to the 
  //   result that b = 3.99...99998.
  //   It is the result of the fact that we use the elementary 
  //   functions exp and gammaln.
  //   This is very close to 4, but is not equal to 4.
  //   Assume that you compute c = b (mod 4) and you 
  //   get c = b = 3.99...99998, intead of getting c  = 0.
  //   This is why, when input arguments are integers,
  //   the result is rounded to the nearest integer with the round function.
  //
  // Examples
  // c = nchoosek ( 4 , 1 ) // 4
  // c = nchoosek ( 5 , 0 ) // 1
  // c = nchoosek ( 5 , 1 ) // 5
  // c = nchoosek ( 5 , 2 ) // 10
  // c = nchoosek ( 5 , 3 ) // 10
  // c = nchoosek ( 5 , 4 ) // 5
  // c = nchoosek ( 5 , 5 ) // 1
  // c = nchoosek ( 17 , 18 ) // 0
  // c = nchoosek ( 17 , -1 ) // 0
  // c = nchoosek ( 1.5 , 0.5 ) // 1.5
  // c = nchoosek ( 10000 , 134 ) // 2.050083865024873735e307
  // nchoosek (10,0:10) // [1,10,45,120,210,252,210,120,45,10,1]
  //
  // Authors
  //   Copyright (C) 2009 - Michael Baudin
  //   Copyright (C) 2009 - John Burkardt
  //
  // Bibliography
  //   "Introduction to discrete probabilities with Scilab", Michael Baudin, 2010
  //   http://bugzilla.scilab.org/show_bug.cgi?id=7589, a bug report for the lack of nchoosek function in Scilab.
  //   http://en.wikipedia.org/wiki/Binomial_coefficients
  //   http://wiki.tcl.tk/1755
  
  if ( ( k < 0 ) | ( k > n ) | ( n < 0 ) ) then 
    nstr = strcat(string(n)," ")
    kstr = strcat(string(k)," ")
    errmsg = msprintf ( gettext ( "%s: The parameters n = %s and k = %s are not consistent." ) , ...
      "nchoosek" , nstr , kstr )
      error ( errmsg )
  end
  // If the input are not floating point integers, generates an error.
  if ( or(round(n)<>n) | or(round(k)<>k) ) then
    nstr = strcat(string(n)," ")
    kstr = strcat(string(k)," ")
    errmsg = msprintf ( gettext ( "%s: The parameters n = %s and k = %s are not floating point integers." ) , ...
    "nchoosek" , nstr , kstr )
    error ( errmsg )
  end
  r = gammaln ( n + 1 ) - gammaln (k + 1) - gammaln (n - k + 1)
  b = exp( r )
  // If the input are floating point integers, round the result.
  if ( and(round(n)==n) & and(round(k)==k) ) then
    b = round ( b )
  end
endfunction

function c = nchoosek_evenmorenaive ( n , k )
c = factorial(n)/ factorial(k) / factorial(n-k)
endfunction

function c = nchoosek_naive ( n , k )
c = prod ( n : -1 : n-k+1 )/ prod (1:k)
endfunction


function c = nchooseklog ( n , k )
  if ( ( k < 0 ) | ( k > n ) ) then 
    c = 0
  else
    c = gammaln ( n + 1 ) - gammaln (k + 1) - gammaln (n - k + 1)
  end
endfunction

// An implementation by Samuel Gougeon :
// http://bugzilla.scilab.org/show_bug.cgi?id=7589
// Tries to mimic Matlab's nchoosek.
// See the discussion in the bug report for advantages and 
// limitations of this routine.
// For example, the following script :
//
// n = 10000
// k = n
// combgen(1:n,k);
// nchoosek_re (1:n,k);
//
// makes it crash.


function rep = nchoosek_re (n,k,self)
  //
  // nchoosek(n,k): returns the number C(n,k) of combinations of k >=0
  //      objects among a set of n >=k. n and k are scalar integers
  //                     C(n,k)= n! / (n-k)! / k!
  //      For C(n,k) > 1/%eps , an approximate with relative accuracy
  //       within 2.eps is returned and a warning is displayed.
  //      For C(n,k) > biggest real that Scilab can handle, a warning
  //       is displayed. An approximated result is returned as a string
  //       in exponential notation.
  //
  // nchoose(v,k): returns all combinations of k components of v,
  //      one per line. v may be a vector, a matrix or an hypermatrix 
  //      of any number of dimensions and of any data type.
  //      If v is a matrix or a hypermatrix, linear indices of its 
  //      components are used to build each line of the output
  //     
  // In both cases, all objects are assumed to be distinguishable.
  //
  // EXAMPLES:
  // nchoosek(30,5), tic(); r=nchoosek(1:30,5); toc(), size(r)
  // nchoosek(6,4), r=nchoosek(["a" "b" "c" "d" "e" "f"],4)
  
  // Author (c): Samuel Gougeon, Le Mans - France
  // Date      : 2010-05-24
  // Software  : Scilab
  // License   : CeCILL-B
  
  fname="nchoosek_re"
  if argn(2)==0, head_comments(fname), rep=[], return, end
  
  nn=size(n,"*")  // Number of components of arg#1
  
  if argn(2)==1
    tmp="%s: Wrong number of input arguments: %d or %d expected.\n";
    error(msprintf(gettext(tmp),fname,0,2))
  end

  if ~isdef("self") & nn==1
    if n<1
      tmp="%s: Wrong value for input argument #%d: Must be in the interval %s.\n";
      error(msprintf(gettext(tmp),fname,1,"[1, oo)"))
    end
    if k<0 | k>n
      tmp="%s: Wrong value for input argument #%d: Must be in the interval %s.\n";
      error(msprintf(gettext(tmp),fname,2,"[0, "+string(fix(n))+"]"))
    end
    
    if k==0 | n==k then rep = 1
    else 
      r=sum(log10((n-k+1):n))-sum(log10(1:k))
      if r > -log10(%eps*2)
        tmp="%s: Result is huge! => Approximation within %%eps=%d relative accuracy";
        warning(msprintf(gettext(tmp),fname,%eps))
        rmax=log10(number_properties("huge"))
        if r>rmax
          tmp="%s: Numerical result too big for Scilab (>10^%d) => Output as a string";
          warning(msprintf(gettext(tmp),fname,rmax))
          f=format(), format(16)
          rep=string(10^(r-fix(r)))+"e+"+string(fix(r))
          format(f(2))
        else
          rep=round(10^r)
        end  
      else
        p=ceil(n+(1-k)/2)
        rep=round((prod(p:n)/prod(1:k))*prod((n-k+1):(p-1)))
      end
    end
  // ===============================================================
  else
    // Cn,k combinations of n's components
    if k<0 | k>nn
      tmp="%s: Wrong value for input argument #%d: Must be in the interval %s.\n";
      error(msprintf(gettext(tmp),fname,2,"[0, "+string(fix(nn))+"]"))
    end
    if k==1, if size(n,2)>1, rep=n.', else rep=n, end
    else
      rep=[]  
      if ~isdef("self") | ~self,
        // Working on indices instead of values, to save memory
        N=n // saving original data
        s=2^ceil(log2(ceil(ceil(log2(nn+1))/8)))
        n=1:nn
        if s<=4, n=iconvert(n,10+s), end
      end
        
      for i=1:size(n,"*")-k+1
        tmp=nchoosek_re(n(i+1:$),k-1,%T)
        col=iconvert(ones(size(tmp,1),1),inttype(n))*n(i)
        if rep~=[], rep=[ rep ; col , tmp], 
        else rep=[col , tmp]
        end
        
      end
      
      // Recovering values from combinations of indicesa
      if ~isdef("self") | ~self, rep=matrix(N(rep),size(rep)), end

    end
  end 
endfunction
