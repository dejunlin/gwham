#ifndef BFGS_HPP
#define BFGS_HPP
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
#include "matrix.hpp"
#include "utils.hpp"

using namespace std;

typedef valtype Doub; // default floating type
typedef NRmatrix<Doub> MatDoub, MatDoub_O, MatDoub_IO;

typedef int Int; // 32 bit integer
typedef unsigned int Uint;

typedef const vector<Doub> VecDoub_I;
typedef vector<Doub> VecDoub;
typedef vector<Doub> VecDoub, VecDoub_O, VecDoub_IO;

typedef bool Bool;

template <class T>
void lnsrch(VecDoub_I &xold, const Doub fold, VecDoub_I &g, VecDoub_IO &p,
VecDoub_O &x, Doub &f, const Doub stpmax, Bool &check, T &func) {
	const Doub ALF=1.0e-4, TOLX=numeric_limits<Doub>::epsilon();
	Doub a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
	Doub rhs1,rhs2,slope=0.0,sum=0.0,temp,test,tmplam;
	Int i,n=xold.size();
	check=false;
	for (i=0;i<n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=0;i<n;i++)
			p[i] *= stpmax/sum;
	for (i=0;i<n;i++)
		slope += g[i]*p[i];
	if (slope >= 0.0) throw("Roundoff problem in lnsrch.");
	test=0.0;
	for (i=0;i<n;i++) {
		temp=abs(p[i])/MAX(abs(xold[i]),Doub{1.0});
		if (temp > test) test=temp;
	}
	alamin=TOLX/test;
	alam=1.0;
	for (;;) {
		for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
		f=func(x);
		if (alam < alamin) {
			for (i=0;i<n;i++) x[i]=xold[i];
			check=true;
			return;
		} else if (f <= fold+ALF*alam*slope) return;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(f-fold-slope));
			else {
				rhs1=f-fold-alam*slope;
				rhs2=f2-fold-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc < 0.0) tmplam=0.5*alam;
					else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
					else tmplam=-slope/(b+sqrt(disc));
				}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = f;
		alam=MAX(tmplam,0.1*alam);
	}
}

template <class T>
void dfpmin(VecDoub_IO &p, const Doub gtol, Int &iter, Doub &fret, T &funcd)
{
	const Int ITMAX=200;
	//const Doub EPS=numeric_limits<Doub>::epsilon();
	const Doub EPS=gtol;
	const Doub TOLX=4*EPS,STPMX=100.0;
	Bool check;
	Doub den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
	Int n=p.size();
	VecDoub dg(n),g(n),hdg(n),pnew(n),xi(n);
	MatDoub hessin(n,n);
	fp=funcd(p);
	funcd.df(p,g);
	for (Int i=0;i<n;i++) {
		for (Int j=0;j<n;j++) hessin[i][j]=0.0;
		hessin[i][i]=1.0;
		xi[i] = -g[i];
		sum += p[i]*p[i];
	}
	stpmax=STPMX*MAX(sqrt(sum),Doub(n));
	for (Int its=0;its<ITMAX;its++) {
		iter=its;
		lnsrch(p,fp,g,xi,pnew,fret,stpmax,check,funcd);
		fp=fret;
		for (Int i=0;i<n;i++) {
			xi[i]=pnew[i]-p[i];
			p[i]=pnew[i];
		}
		test=0.0;
		for (Int i=0;i<n;i++) {
			temp=abs(xi[i])/MAX(abs(p[i]),Doub{1.0});
			if (temp > test) test=temp;
		}
		if (test < TOLX)
			return;
		for (Int i=0;i<n;i++) dg[i]=g[i];
		funcd.df(p,g);
		test=0.0;
		den=MAX(fret,Doub{1.0});
		for (Int i=0;i<n;i++) {
			temp=abs(g[i])*MAX(abs(p[i]),Doub{1.0})/den;
			if (temp > test) test=temp;
		}
		if (test < gtol)
			return;
		for (Int i=0;i<n;i++)
			dg[i]=g[i]-dg[i];
		for (Int i=0;i<n;i++) {
			hdg[i]=0.0;
			for (Int j=0;j<n;j++) hdg[i] += hessin[i][j]*dg[j];
		}
		fac=fae=sumdg=sumxi=0.0;
		for (Int i=0;i<n;i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}
		if (fac > sqrt(EPS*sumdg*sumxi)) {
			fac=1.0/fac;
			fad=1.0/fae;
			for (Int i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
			for (Int i=0;i<n;i++) {
				for (Int j=i;j<n;j++) {
					hessin[i][j] += fac*xi[i]*xi[j]
						-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
					hessin[j][i]=hessin[i][j];
				}
			}
		}
		for (Int i=0;i<n;i++) {
			xi[i]=0.0;
			for (Int j=0;j<n;j++) xi[i] -= hessin[i][j]*g[j];
		}
	}
	throw("too many iterations in dfpmin");
}

#endif
