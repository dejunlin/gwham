#ifndef CG_HPP
#define CG_HPP
/*
 * =====================================================================================
 *
 *       Filename:  cg.hpp
 *
 *    Description:  adapted from mins_ndim.h and mins.h in Numerical Recipes 3rd Edition 
 *
 *        Version:  1.0
 *        Created:  09/21/2014 07:51:40 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include <limits>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <utility>
#include "../typedefs.hpp"
#include "../exception.hpp"
#include "utils.hpp"

using namespace std;

typedef valtype Doub; // default floating type

typedef int Int; // 32 bit integer
typedef unsigned int Uint;

typedef const vector<Doub> VecDoub_I;
typedef vector<Doub> VecDoub;
typedef vector<Doub> VecDoub, VecDoub_O, VecDoub_IO;

typedef bool Bool;

struct Bracketmethod {
	Doub ax,bx,cx,fa,fb,fc;
	template <class T>
	void bracket(const Doub a, const Doub b, T &func)
	{
		const Doub GOLD=1.618034,GLIMIT=100.0,TINY=1.0e-8;
		ax=a; bx=b;
		Doub fu;
		fa=func(ax);
		fb=func(bx);
		if (fb > fa) {
			SWAP(ax,bx);
			SWAP(fb,fa);
		}
		cx=bx+GOLD*(bx-ax);
		fc=func(cx);
		while (fb > fc) {
			Doub r=(bx-ax)*(fb-fc);
			Doub q=(bx-cx)*(fb-fa);
			Doub u=bx-((bx-cx)*q-(bx-ax)*r)/
				(2.0*SIGN(MAX(abs(q-r),TINY),q-r));
			Doub ulim=bx+GLIMIT*(cx-bx);
			if ((bx-u)*(u-cx) > 0.0) {
				fu=func(u);
				if (fu < fc) {
					ax=bx;
					bx=u;
					fa=fb;
					fb=fu;
					return;
				} else if (fu > fb) {
					cx=u;
					fc=fu;
					return;
				}
				u=cx+GOLD*(cx-bx);
				fu=func(u);
			} else if ((cx-u)*(u-ulim) > 0.0) {
				fu=func(u);
				if (fu < fc) {
					shft3(bx,cx,u,u+GOLD*(u-cx));
					shft3(fb,fc,fu,func(u));
				}
			} else if ((u-ulim)*(ulim-cx) >= 0.0) {
				u=ulim;
				fu=func(u);
			} else {
				u=cx+GOLD*(cx-bx);
				fu=func(u);
			}
			shft3(ax,bx,cx,u);
			shft3(fa,fb,fc,fu);
		}
	}
	inline void shft2(Doub &a, Doub &b, const Doub c)
	{
		a=b;
		b=c;
	}
	inline void shft3(Doub &a, Doub &b, Doub &c, const Doub d)
	{
		a=b;
		b=c;
		c=d;
	}
	inline void mov3(Doub &a, Doub &b, Doub &c, const Doub d, const Doub e,
		const Doub f)
	{
		a=d; b=e; c=f;
	}
};

struct Dbrent : Bracketmethod {
	Doub xmin,fmin;
	const Doub tol;
	Dbrent(const Doub toll=3.0e-8) : tol(toll) {}
	template <class T>
	Doub minimize(T &funcd)
	{
		const Int ITMAX=20000;
		const Doub ZEPS=numeric_limits<Doub>::epsilon()*1.0e-3;
		Bool ok1,ok2;
		Doub a,b,d=0.0,d1,d2,du,dv,dw,dx,e=0.0;
		Doub fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
	
		a=(ax < cx ? ax : cx);
		b=(ax > cx ? ax : cx);
		x=w=v=bx;
		fw=fv=fx=funcd(x);
		dw=dv=dx=funcd.df(x);
		for (Int iter=0;iter<ITMAX;iter++) {
			xm=0.5*(a+b);
			tol1=tol*abs(x)+ZEPS;
			tol2=2.0*tol1;
			if (abs(x-xm) <= (tol2-0.5*(b-a))) {
				fmin=fx;
				return xmin=x;
			}
			if (abs(e) > tol1) {
				d1=2.0*(b-a);
				d2=d1;
				if (dw != dx) d1=(w-x)*dx/(dx-dw);
				if (dv != dx) d2=(v-x)*dx/(dx-dv);
				u1=x+d1;
				u2=x+d2;
				ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
				ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
				olde=e;
				e=d;
				if (ok1 || ok2) {
					if (ok1 && ok2)
						d=(abs(d1) < abs(d2) ? d1 : d2);
					else if (ok1)
						d=d1;
					else
						d=d2;
					if (abs(d) <= abs(0.5*olde)) {
						u=x+d;
						if (u-a < tol2 || b-u < tol2)
							d=SIGN(tol1,xm-x);
					} else {
						d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
					}
				} else {
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
				}
			} else {
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
			}
			if (abs(d) >= tol1) {
				u=x+d;
				fu=funcd(u);
			} else {
				u=x+SIGN(tol1,d);
				fu=funcd(u);
				if (fu > fx) {
					fmin=fx;
					return xmin=x;
				}
			}
			du=funcd.df(u);
			if (fu <= fx) {
				if (u >= x) a=x; else b=x;
				mov3(v,fv,dv,w,fw,dw);
				mov3(w,fw,dw,x,fx,dx);
				mov3(x,fx,dx,u,fu,du);
			} else {
				if (u < x) a=u; else b=u;
				if (fu <= fw || w == x) {
					mov3(v,fv,dv,w,fw,dw);
					mov3(w,fw,dw,u,fu,du);
				} else if (fu < fv || v == x || v == w) {
					mov3(v,fv,dv,u,fu,du);
				}
			}
		}
		throw(General_Exception("Too many iterations in routine dbrent"));
	}
};

template <class T>
struct Df1dim {
	const VecDoub &p;
	const VecDoub &xi;
	Int n;
	T &funcd;
	VecDoub xt;
	VecDoub dft;
	Df1dim(VecDoub_I &pp, VecDoub_I &xii, T &funcdd) : p(pp),
		xi(xii), n(pp.size()), funcd(funcdd), xt(n), dft(n) {}
	Doub operator()(const Doub x)
	{
		for (Int j=0;j<n;j++)
			xt[j]=p[j]+x*xi[j];
		return funcd(xt);
	}
	Doub df(const Doub x)
	{
		Doub df1=0.0;
		funcd.df(xt,dft);
		for (Int j=0;j<n;j++)
			df1 += dft[j]*xi[j];
		return df1;
	}
};

template <class T>
struct Dlinemethod {
	VecDoub p;
	VecDoub xi;
	T &func;
	Int n;
	Dlinemethod(T &funcc) : func(funcc) {}
	Doub linmin()
	{
		Doub ax,xx,xmin;
		n=p.size();
		Df1dim<T> df1dim(p,xi,func);
		ax=0.0;
		xx=1.0;
		Dbrent dbrent;
		dbrent.bracket(ax,xx,df1dim);
		xmin=dbrent.minimize(df1dim);
		for (Int j=0;j<n;j++) {
			xi[j] *= xmin;
			p[j] += xi[j];
		}
		return dbrent.fmin;
	}
};

template <class T>
struct Frprmn : Dlinemethod<T> {
	Int iter;
	Doub fret;
	using Dlinemethod<T>::func;
	using Dlinemethod<T>::linmin;
	using Dlinemethod<T>::p;
	using Dlinemethod<T>::xi;
	const Doub ftol;
	Frprmn(T &funcd, const Doub ftoll=3.0e-8) : Dlinemethod<T>(funcd),
		ftol(ftoll) {}
	VecDoub minimize(VecDoub_I &pp)
	{
		const Int ITMAX=2000000000;
		const Doub EPS=1.0e-18;
		const Doub GTOL=1.0e-8;
		Doub gg,dgg;
		Int n=pp.size();
		p=pp;
		VecDoub g(n),h(n);
		xi.resize(n);
		Doub fp=func(p);
		func.df(p,xi);
		for (Int j=0;j<n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j];
		}
		for (Int its=0;its<ITMAX;its++) {
		        if(!(its % 100)) {
			  cout << "#At iteration: " << its << " states are:\n";
			  for(uint i = 0; i < p.size(); ++i) cout << "#\t" << i << "\t" << p[i] << endl;
			  cout.flush();
			}
			iter=its;
			fret=linmin();
			if (2.0*abs(fret-fp) <= ftol*(abs(fret)+abs(fp)+EPS)) {
			  cout << "#Convergence met at iteration " << its << " because the difference in F(p) between 2 consecutive iteration is smaller than " << ftol*(abs(fret)+abs(fp)+EPS)/2 << endl;
				return p;
			}
			fp=fret;
			func.df(p,xi);
			Doub test=0.0;
			Doub den=MAX(fp, Doub{1.0});
			for (Int j=0;j<n;j++) {
				Doub temp=abs(xi[j])*MAX(abs(p[j]),Doub{1.0})/den;
				if (temp > test) test=temp;
			}
			if (test < GTOL) { 
			  cout << "#Convergence met at iteration " << its << " because gradient is smaller than " << GTOL << endl;
			  return p;
			}
			dgg=gg=0.0;
			for (Int j=0;j<n;j++) {
				gg += g[j]*g[j];
//			  dgg += xi[j]*xi[j];
				dgg += (xi[j]+g[j])*xi[j];
			}
			if (gg == 0.0) {
			  cout << "#Convergence met at iteration " << its << " because the search direction vector is exactly zero" << endl;
			  return p;
			}
			Doub gam=dgg/gg;
			for (Int j=0;j<n;j++) {
				g[j] = -xi[j];
				xi[j]=h[j]=g[j]+gam*h[j];
			}
		}
		throw(General_Exception("Too many iterations in frprmn"));
	}
};

#endif
