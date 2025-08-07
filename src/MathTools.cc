// -*- C++ -*-

#include "MathTools.hh"

#include <cmath>
#include <iostream>
#include <iomanip>

#include "FuncName.hh"

//_____________________________________________________________________________
namespace math
{
//___________________________________________________________________________
G4bool
GaussElim(G4double **a, G4int n, G4double *b, G4int *indx, G4int *ipiv)
{
  G4double big,c,pivinv,sum,c2;
  G4int js,irow=0;

  for(G4int i=0; i<n; ++i) ipiv[i]=0;
  for(G4int i=0; i<n; ++i){
    big=0.0;
    for(G4int j=0; j<n; ++j){
      if(!ipiv[j]){
        if((c=fabs(a[j][i]))>=big){ big=c; irow=j; }
      }
      else if(ipiv[j]>1){
#ifdef DEBUG
        std::cerr << FUNC_NAME << ": Singular Matrix" << std::endl;
#endif
        return false;
      }
    }
    ++(ipiv[irow]); indx[i]=irow;
    if(a[irow][i]==0.0){
#ifdef DEBUG
      std::cerr << FUNC_NAME << ": Singular Matrix" << std::endl;
#endif
      return false;
    }
    pivinv=1.0/a[irow][i];

    for(G4int j=0; j<n; ++j){
      if(ipiv[j]==0){
        c=a[j][i]; a[j][i]=0.0;
        for(G4int k=i+1; k<n; ++k)
          a[j][k]-=a[irow][k]*pivinv*c;
        b[j]-=b[irow]*pivinv*c;
      }
    }
  }

  b[indx[n-1]]/=a[indx[n-1]][n-1];
  for(G4int i=n-2; i>=0; --i){
    sum=b[indx[i]];
    for(G4int j=i+1; j<n; ++j)
      sum-=a[indx[i]][j]*b[indx[j]];
    b[indx[i]]=sum/a[indx[i]][i];
  }
  for(G4int i=0; i<n; ++i){
    if(indx[i]!=i){
      c2=b[i];
      for(G4int j=indx[i]; indx[j]!=j; j=js){
        c=b[j]; b[j]=c2; c2=c; js=indx[j]; indx[j]=j;
      }
    }
  }
  return true;
}

//___________________________________________________________________________
G4bool
GaussJordan(G4double **a, G4int n, G4double *b,
            G4int *indxr, G4int *indxc, G4int *ipiv)
{
  for(G4int j=0; j<n; ++j) ipiv[j]=0;
  for(G4int i=0; i<n; ++i){
    G4double big=0.0;
    G4int irow=-1, icol=-1;
    for(G4int j=0; j<n; ++j)
      if(ipiv[j]!=1)
        for(G4int k=0; k<n; ++k){
          if(ipiv[k]==0){
            if(fabs(a[j][k])>=big){
              big=fabs(a[j][k]);
              irow=j; icol=k;
            }
          }
          else if(ipiv[k]>1){
#ifdef DEBUG
            std::cerr << FUNC_NAME << ": Singular Matrix"
                      << std::endl;
#endif
            return false;
          }
        }
    ++(ipiv[icol]);

    if(irow!=icol){
      for(G4int k=0; k<n; ++k){
        G4double ta=a[irow][k];
        a[irow][k]=a[icol][k];
        a[icol][k]=ta;
      }
      G4double tb=b[irow];
      b[irow]=b[icol];
      b[icol]=tb;
    }

    indxr[i]=irow; indxc[i]=icol;

    if(a[icol][icol]==0.0){
#ifdef DEBUG
      std::cerr << FUNC_NAME << ": Singular Matrix"  << std::endl;
#endif
      return false;
    }
    G4double pivinv=1./a[icol][icol];
    a[icol][icol]=1.;
    for(G4int k=0; k<n; ++k) a[icol][k]*=pivinv;
    b[icol]*=pivinv;
    for(G4int k=0; k<n; ++k){
      if(k!=icol){
        G4double d=a[k][icol];
        a[k][icol]=0.;
        for(G4int l=0; l<n; ++l) a[k][l] -= a[icol][l]*d;
        b[k] -= b[icol]*d;
      }
    }
  }

  for(G4int l=n-1; l>=0; --l){
    if(indxr[l]!=indxc[l]){
      for(G4int k=0; k<n; ++k){
        G4double t=a[k][indxr[l]];
        a[k][indxr[l]]=a[k][indxc[l]];
        a[k][indxc[l]]=t;
      }
    }
  }
  return true;
}

//___________________________________________________________________________
G4bool
InterpolateRatio(G4int n, const G4double *xa, const G4double *ya,
                 G4double *w1, G4double *w2,
                 G4double x, G4double &y, G4double &dy)
{
  G4int i, m, ns=1;
  G4double w, t, hh, h, dd;

  hh=fabs(x-xa[0]);
  for(i=1; i<=n; ++i){
    h=fabs(x-xa[i-1]);
    if(h==0.0) { y=ya[i-1]; dy=0.0; return true; }
    else if(h<hh){ ns=i; hh=h; }
    w1[i-1]=ya[i-1]; w2[i-1]=ya[i-1]*(1.+Epsilon);
  }
  y=ya[ns-1]; ns--;
  for(m=1; m<n; ++m){
    for(i=1; i<=n-m; ++i){
      w=w1[i]-w2[i-1]; h=xa[i+m-1]-x;
      t=(xa[i-1]-x)*w2[i-1]/h; dd=t-w1[i];
      if(dd==0.0){
#ifdef DEBUG
        std::cerr << FUNC_NAME << ": Error" << std::endl;
#endif
        y=Infinity; dy=Infinity; return false;
      }
      dd=w/dd; w2[i-1]=w1[i]*dd; w1[i-1]=t*dd;
    }
    if(2*ns < (n-m)) dy=w1[ns];
    else { dy=w2[ns-1]; ns--; }
    y+=dy;
  }
#if 0
  std::cout << FUNC_NAME << ": x=" << std::setw(10) << x
            << " y=" << std::setw(10) << y << std::endl;
#endif
  return true;
}

//___________________________________________________________________________
G4bool
InterpolatePol(G4int n, const G4double *xa, const G4double *ya,
               G4double *w1, G4double *w2,
               G4double x, G4double &y, G4double &dy)
{
  G4int i, m, ns=1;
  G4double den, dif, dift, ho, hp, w;

  dif=fabs(x-xa[0]);
  for(i=1; i<=n; ++i){
    if((dift=fabs(x-xa[i-1]))<dif){ ns=i; dif=dift; }
    w1[i-1]=w2[i-1]=ya[i-1];
  }
  y=ya[ns-1]; --ns;
  for(m=1; m<n; ++m){
    for(i=1; i<=n-m; ++i){
      ho=xa[i-1]-x; hp=xa[i+m-1]-x;
      w=w1[i]-w2[i-1];
      den=ho-hp;
      if(den==0.0){
#ifdef DEBUG
        std::cerr << FUNC_NAME << ": Error" << std::endl;
#endif
        y=Infinity; dy=Infinity; return false;
      }
      den=w/den;
      w2[i-1]=hp*den; w1[i-1]=ho*den;
    }
    if(2*ns<(n-m)) dy=w1[ns];
    else           { dy=w2[ns-1]; --ns; }
    y+=dy;
  }
  return false;
}

//___________________________________________________________________________
G4bool
SVDksb(G4double **u, const G4double *w, G4double **v,
       G4int m, G4int n, const G4double *b, G4double *x, G4double *wv)
{
  for(G4int j=0; j<n; ++j){
    G4double s=0.0;
    if(w[j]!=0.0){
      for(G4int i=0; i<m; ++i)
        s += u[i][j]*b[i];
      s /= w[j];
    }
    wv[j]=s;
  }
  for(G4int i=0; i<n; ++i){
    G4double s=0.0;
    for(G4int j=0; j<n; ++j)
      s += v[i][j]*wv[j];
    x[i]=s;
  }
  return true;
}

//___________________________________________________________________________
G4bool
SVDcmp(G4double **a, G4int m, G4int n, G4double *w,
       G4double **v, G4double *wv)
{
  G4double g=0.0, scale=0.0, anorm=0.0;
  G4double s, f, h, c;
  G4int nm;

#ifdef DebugPrint
  for(G4int i=0; i<n; ++i) {
    w[i]=wv[i]=0.0;
    for(G4int j=0; j<n; ++j) v[j][i]=0.0;
  }

  {
    std::ios::fmtflags oldFlags = std::cout.flags();
    std::size_t oldPrec = std::cout.precision();
    std::cout.setf(std::ios::scientific);
    std::cout.precision(3);
    std::cout << FUNC_NAME << ": A in SVDcmp 1" <<  std::endl;
    for(G4int ii=0; ii<m; ++ii){
      for(G4int ij=0; ij<n; ++ij){
        std::cout << std::setw(12) << a[ii][ij];
        if(ij!=n-1) std::cout << ",";
      }
      std::cout << std::endl;
    }
    std::cout << FUNC_NAME << ": V in SVDcmp 1" << std::endl;
    for(G4int ii=0; ii<n; ++ii){
      for(G4int ij=0; ij<n; ++ij){
        std::cout << std::setw(12) << v[ii][ij];
        if(ij!=n-1) std::cout << ",";
      }
      std::cout << std::endl;
    }
    std::cout << FUNC_NAME << ": W in SVDcmp 1" << std::endl;
    for(G4int ij=0; ij<n; ++ij){
      std::cout << std::setw(12) << w[ij];
      if(ij!=n-1) std::cout << ",";
    }
    std::cout << std::endl;

    std::cout << FUNC_NAME << ": WV in SVDcmp 1" << std::endl;
    for(G4int ij=0; ij<n; ++ij){
      std::cout << std::setw(12) << wv[ij];
      if(ij!=n-1) std::cout << ",";
    }
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout.flags(oldFlags);
    std::cout.precision(oldPrec);
  }
#endif

  // Householder method
  for(G4int i=0; i<n; ++i){
    wv[i]=scale*g;
    g = scale = 0.0;
    if(i<m){
      for(G4int k=i; k<m; ++k) scale += fabs(a[k][i]);
      if(scale!=0.){
        s = 0;
        for(G4int k=i; k<m; ++k){
          a[k][i] /= scale;
          s += a[k][i]*a[k][i];
        }
        f = a[i][i];
        g = ((f>0.0) ? -sqrt(s) : sqrt(s));
        h = f*g-s;
        a[i][i] = f-g;
        for(G4int j=i+1; j<n; ++j){
          s = 0.0;
          for(G4int k=i; k<m; ++k) s += a[k][i]*a[k][j];
          f = s/h;
          for(G4int k=i; k<m; ++k)  a[k][j] += f*a[k][i];
        }
        for(G4int k=i; k<m; ++k) a[k][i] *= scale;
      }
    }     /* if(i<m) */
    w[i] = scale*g;
    g = s = scale = 0.0;

    if(i<m && i!=n-1){
      for(G4int k=i+1; k<n; ++k) scale += fabs(a[i][k]);
      if(scale!=0.0){
        for(G4int k=i+1; k<n; ++k){
          a[i][k] /= scale;
          s += a[i][k]*a[i][k];
        }
        f = a[i][i+1];
        g = ((f>0.0) ? -sqrt(s) : sqrt(s));
        h = f*g-s;
        a[i][i+1] = f-g;
        for(G4int k=i+1; k<n; ++k) wv[k] = a[i][k]/h;
        for(G4int j=i+1; j<m; ++j){
          s = 0.0;
          for(G4int k=i+1; k<n; ++k) s += a[j][k]*a[i][k];
          for(G4int k=i+1; k<n; ++k) a[j][k] += s*wv[k];
        }
        for(G4int k=i+1; k<n; ++k) a[i][k] *= scale;
      }
    }   /* if(i<m && i!=n-1) */
    G4double tmp=fabs(w[i])+fabs(wv[i]);
    if(tmp>anorm) anorm = tmp;
  }     /* for(G4int i ...) */

#if DebugPrint
  {
    std::ios::fmtflags oldFlags = std::cout.flags();
    std::size_t oldPrec = std::cout.precision();
    std::cout.setf(std::ios::scientific);
    std::cout.precision(3);
    std::cout << FUNC_NAME << ": A in SVDcmp 2" <<  std::endl;
    for(G4int ii=0; ii<m; ++ii){
      for(G4int ij=0; ij<n; ++ij){
        std::cout << std::setw(12) << a[ii][ij];
        if(ij!=n-1) std::cout << ",";
      }
      std::cout << std::endl;
    }
    std::cout << FUNC_NAME << ": V in SVDcmp 2" << std::endl;
    for(G4int ii=0; ii<n; ++ii){
      for(G4int ij=0; ij<n; ++ij){
        std::cout << std::setw(12) << v[ii][ij];
        if(ij!=n-1) std::cout << ",";
      }
      std::cout << std::endl;
    }
    std::cout << FUNC_NAME << ": W in SVDcmp 2" << std::endl;
    for(G4int ij=0; ij<n; ++ij){
      std::cout << std::setw(12) << w[ij];
      if(ij!=n-1) std::cout << ",";
    }
    std::cout << std::endl;

    std::cout << FUNC_NAME << ": WV in SVDcmp 2" << std::endl;
    for(G4int ij=0; ij<n; ++ij){
      std::cout << std::setw(12) << wv[ij];
      if(ij!=n-1) std::cout << ",";
    }
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout.flags(oldFlags);
    std::cout.precision(oldPrec);
  }
#endif

  for(G4int i=n-1; i>=0; --i){
    if(i<n-1){
      if(g!=0.0){
        for(G4int j=i+1; j<n; ++j) v[j][i] = (a[i][j]/a[i][i+1])/g;
        for(G4int j=i+1; j<n; ++j){
          s = 0.0;
          for(G4int k=i+1; k<n; ++k) s += a[i][k]*v[k][j];
          for(G4int k=i+1; k<n; ++k) v[k][j] += s*v[k][i];
        }
      }
      for(G4int j=i+1; j<n; ++j)
        v[i][j] = v[j][i] = 0.0;
    }
    v[i][i]=1.0;  g=wv[i];
  }   /* for(G4int i= ...) */

#ifdef DebugPrint
  {
    std::ios::fmtflags oldFlags = std::cout.flags();
    std::size_t oldPrec = std::cout.precision();
    std::cout.setf(std::ios::scientific);
    std::cout.precision(3);
    std::cout << FUNC_NAME << ": A in SVDcmp 3" <<  std::endl;
    for(G4int ii=0; ii<m; ++ii){
      for(G4int ij=0; ij<n; ++ij){
        std::cout << std::setw(12) << a[ii][ij];
        if(ij!=n-1) std::cout << ",";
      }
      std::cout << std::endl;
    }
    std::cout << FUNC_NAME << ": V in SVDcmp 3" << std::endl;
    for(G4int ii=0; ii<n; ++ii){
      for(G4int ij=0; ij<n; ++ij){
        std::cout << std::setw(12) << v[ii][ij];
        if(ij!=n-1) std::cout << ",";
      }
      std::cout << std::endl;
    }
    std::cout << FUNC_NAME << ": W in SVDcmp 3" << std::endl;
    for(G4int ij=0; ij<n; ++ij){
      std::cout << std::setw(12) << w[ij];
      if(ij!=n-1) std::cout << ",";
    }
    std::cout << std::endl;

    std::cout << FUNC_NAME << ": WV in SVDcmp 3" << std::endl;
    for(G4int ij=0; ij<n; ++ij){
      std::cout << std::setw(12) << wv[ij];
      if(ij!=n-1) std::cout << ",";
    }
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout.flags(oldFlags);
    std::cout.precision(oldPrec);
  }
#endif

  G4int mn = ((m<n) ? m : n);

  for(G4int i=mn-1; i>=0; --i){
    g=w[i];
    for(G4int j=i+1; j<n; ++j) a[i][j]=0.0;
    if(g!=0.0){
      g = 1./g;
      for(G4int j=i+1; j<n; ++j){
        s = 0.0;
        for(G4int k=i+1; k<m; ++k) s += a[k][i]*a[k][j];
        f = (s/a[i][i])*g;
        for(G4int k=i; k<m; ++k) a[k][j] += f*a[k][i];
      }
      for(G4int j=i; j<m; ++j) a[j][i] *= g;
    }
    else
      for(G4int j=i; j<m; ++j) a[j][i] = 0.0;

    a[i][i] += 1.0;
  }   /* for(G4int i= ...) */

#ifdef DebugPrint
  {
    std::ios::fmtflags oldFlags = std::cout.flags();
    std::size_t oldPrec = std::cout.precision();
    std::cout.setf(std::ios::scientific);
    std::cout.precision(3);
    std::cout << FUNC_NAME << ": A in SVDcmp 4" <<  std::endl;
    for(G4int ii=0; ii<m; ++ii){
      for(G4int ij=0; ij<n; ++ij){
        std::cout << std::setw(12) << a[ii][ij];
        if(ij!=n-1) std::cout << ",";
      }
      std::cout << std::endl;
    }
    std::cout << FUNC_NAME << ": V in SVDcmp 4" << std::endl;
    for(G4int ii=0; ii<n; ++ii){
      for(G4int ij=0; ij<n; ++ij){
        std::cout << std::setw(12) << v[ii][ij];
        if(ij!=n-1) std::cout << ",";
      }
      std::cout << std::endl;
    }
    std::cout << FUNC_NAME << ": W in SVDcmp 4" << std::endl;
    for(G4int ij=0; ij<n; ++ij){
      std::cout << std::setw(12) << w[ij];
      if(ij!=n-1) std::cout << ",";
    }
    std::cout << std::endl;

    std::cout << FUNC_NAME << ": WV in SVDcmp 4" << std::endl;
    for(G4int ij=0; ij<n; ++ij){
      std::cout << std::setw(12) << wv[ij];
      if(ij!=n-1) std::cout << ",";
    }
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout.flags(oldFlags);
    std::cout.precision(oldPrec);
  }
#endif

  G4int ll=1;

  for(G4int k=n-1; k>=0; --k){
    for(G4int its=1; its<=30; ++its){
      G4int flag=1; nm=ll;
      for(ll=k; ll>=0; --ll){
        nm=ll-1;
        if(fabs(wv[ll])+anorm == anorm){
          flag=0; break;
        }
        if(fabs(w[nm])+anorm == anorm)
          break;
      }

      if(flag){
        c=0.0; s=1.0;
        for(G4int i=ll; i<=k; ++i){
          f = s*wv[i]; wv[i] *= c;
          if(fabs(f)+anorm == anorm)
            break;
          g=w[i]; h=pythag(f,g); w[i]=h;
          h=1./h; c=g*h; s=-f*h;
          for(G4int j=0; j<m; ++j){
            G4double y=a[j][nm], z=a[j][i];
            a[j][nm]=y*c+z*s; a[j][i]=z*c-y*s;
          }
        }
      }   /* if(flag) */

      G4double z = w[k];
      if(ll==k){
        if(z<0.){
          w[k]=-z;
          for(G4int j=0; j<n; ++j) v[j][k]=-v[j][k];
        }
        break;
      }
#ifdef DEBUG
      if(its==30){
        // 	std::cerr << FUNC_NAME
        // 		  << ": -- no convergence in 30 dvdcmp iterations --"
        // 		  << std::endl;
        return false;
      }
#endif
      nm=k-1;
      G4double x=w[ll], y=w[nm];
      g=wv[nm]; h=wv[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.*h*y);
      g=pythag(f,1.0);
      G4double gtmp = ((f>0.) ? g : -g);
      f=((x-z)*(x+z)+h*((y/(f+gtmp))-h))/x;
      c=s=1.0;
      for(G4int j=ll; j<=nm; ++j){
        g=wv[j+1]; y=w[j+1]; h=s*g; g=c*g;
        z=pythag(f,h);
        wv[j]=z; c=f/z; s=h/z;
        f=x*c+g*s; g=g*c-x*s;
        h=y*s; y=y*c;
        for(G4int jj=0; jj<n; ++jj){
          x=v[jj][j]; z=v[jj][j+1];
          v[jj][j]=x*c+z*s; v[jj][j+1]=z*c-x*s;
        }
        z=pythag(f,h);
        w[j]=z;
        if(z!=0.0){ z=1./z; c=f*z; s=h*z; }
        f=c*g+s*y; x=c*y-s*g;
        for(G4int jj=0; jj<m; ++jj){
          y=a[jj][j]; z=a[jj][j+1];
          a[jj][j]=y*c+z*s; a[jj][j+1]=z*c-y*s;
        }
      }
      wv[ll]=0.0; wv[k]=f; w[k]=x;
    }   /* for(G4int its ...) */
  }     /* for(G4int k= ...) */

#ifdef DebugPrint
  {
    std::ios::fmtflags oldFlags = std::cout.flags();
    std::size_t oldPrec = std::cout.precision();
    std::cout.setf(std::ios::scientific);
    std::cout.precision(3);
    std::cout << FUNC_NAME << ": A in SVDcmp 5" <<  std::endl;
    for(G4int ii=0; ii<m; ++ii){
      for(G4int ij=0; ij<n; ++ij){
        std::cout << std::setw(12) << a[ii][ij];
        if(ij!=n-1) std::cout << ",";
      }
      std::cout << std::endl;
    }
    std::cout << FUNC_NAME << ": V in SVDcmp 5" << std::endl;
    for(G4int ii=0; ii<n; ++ii){
      for(G4int ij=0; ij<n; ++ij){
        std::cout << std::setw(12) << v[ii][ij];
        if(ij!=n-1) std::cout << ",";
      }
      std::cout << std::endl;
    }
    std::cout << FUNC_NAME << ": W in SVDcmp 5" << std::endl;
    for(G4int ij=0; ij<n; ++ij){
      std::cout << std::setw(12) << w[ij];
      if(ij!=n-1) std::cout << ",";
    }
    std::cout << std::endl;

    std::cout << FUNC_NAME << ": WV in SVDcmp 5" << std::endl;
    for(G4int ij=0; ij<n; ++ij){
      std::cout << std::setw(12) << wv[ij];
      if(ij!=n-1) std::cout << ",";
    }
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout.flags(oldFlags);
    std::cout.precision(oldPrec);
  }
#endif

  return true;
}
}
