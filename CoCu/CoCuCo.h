#ifndef COCUCO_H
#define COCUCO_H

#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>

using namespace std;
using namespace Eigen;
typedef Matrix2d m2d;
typedef vector<Matrix2d, aligned_allocator<Matrix2d>> vm2d;

double gmean(double x, double y){
	double gmean;
	if (x == y)
		gmean = x;
	else if (x*y > 0){
		if (x < 0)
			gmean=-sqrt(x*y);
		if (x > 0)
		gmean = sqrt(x*y);
      	}
	else
       		gmean=(x+y)/2.;
 	return gmean;
}

Matrix<complex<double>, 9, 9> eint1(vector<double> vec, double x, double y, double z)
{
//     THIS ROUTINE IS SPIN DEPENDENT :-
//     IT IS WRITTEN FOR s,p and d BANDS ONLY
      double sss, sps, pps, ppp, ss, ps, pp, ds, dp, dd;
      sss = vec[0]; sps = vec[1]; pps = vec[2]; ppp = vec[3];
      ss = vec[4]; ps = vec[5]; pp = vec[6]; ds = vec[7];
      dp = vec[8]; dd = vec[9];
      Matrix<complex<double>, 9, 9> b;
      double xx=x*x;
      double xy=x*y;
      double yy=y*y;
      double yz=y*z;
      double zz=z*z;
      double zx=z*x;
      double xxyy=xx*yy;
      double yyzz=yy*zz;
      double zzxx=zz*xx;
      double aux=pps-ppp;
      double r3=sqrt(3.);
      double aux1=r3*ss;
      double f8=3.*zz-1.;
      double f1=xx+yy;
      double f2=xx-yy;
      double f3=zz-.5*f1;
      double g1=1.5*f2*ds;
      double g2=r3*f3*ds;
      b(0,0)=sss;
      b(0,1)=x*sps;
      b(0,2)=y*sps;
      b(0,3)=z*sps;
      b(1,0)=-b(0,1);
      b(1,1)=xx*pps+(1.-xx)*ppp;
      b(1,2)=xy*aux;
      b(1,3)=zx*aux;
      b(2,0)=-b(0,2);
      b(2,1)=b(1,2);
      b(2,2)=yy*pps+(1.-yy)*ppp;
      b(2,3)=yz*aux;
      b(3,0)=-b(0,3);
      b(3,1)=b(1,3);
      b(3,2)=b(2,3);
      b(3,3)=zz*pps+(1.-zz)*ppp;
      b(0,4)=xy*aux1;
      b(0,5)=yz*aux1;
      b(0,6)=zx*aux1;
      b(0,7)=.5*f2*aux1;
      b(0,8)=.5*f8*ss;
      b(4,0)=b(0,4);
      b(5,0)=b(0,5);
      b(6,0)=b(0,6);
      b(7,0)=b(0,7);
      b(8,0)=b(0,8);
      double f4=.5*r3*f2*ps;
      double f5=.5*f8*ps;
      double aux2=r3*xx*ps+(1.-2.*xx)*pp;
      b(1,4)=aux2*y;
      b(1,5)=(r3*ps-2.*pp)*xy*z;
      b(1,6)=aux2*z;
      b(1,7)=(f4+(1.-f2)*pp)*x;
      b(1,8)=(f5-r3*zz*pp)*x;
      double aux3=(r3*yy*ps+(1.-2.*yy)*pp);
      b(2,4)=aux3*x;
      b(2,5)=aux3*z;
      b(2,6)=b(1,5);
      b(2,7)=(f4-(1.+f2)*pp)*y;
      b(2,8)=(f5-r3*zz*pp)*y;
      double aux4=r3*zz*ps+(1.-2.*zz)*pp;
      b(3,4)=b(1,5);
      b(3,5)=aux4*y;
      b(3,6)=aux4*x;
      b(3,7)=(f4-f2*pp)*z;
      b(3,8)=(f5+r3*f1*pp)*z;
      b(4,1)=-b(1,4);
      b(5,1)=-b(1,5);
      b(6,1)=-b(1,6);
      b(7,1)=-b(1,7);
      b(8,1)=-b(1,8);
      b(4,2)=-b(2,4);
      b(5,2)=-b(2,5);
      b(6,2)=-b(2,6);
      b(7,2)=-b(2,7);
      b(8,2)=-b(2,8);
      b(4,3)=-b(3,4);
      b(5,3)=-b(3,5);
      b(6,3)=-b(3,6);
      b(7,3)=-b(3,7);
      b(8,3)=-b(3,8);
      b(4,4)=3.*xxyy*ds+(f1-4.*xxyy)*dp+(zz+xxyy)*dd;
      b(4,5)=(3.*yy*ds+(1.-4.*yy)*dp+(yy-1.)*dd)*zx;
      b(4,6)=(3.*xx*ds+(1.-4.*xx)*dp+(xx-1.)*dd)*yz;
      b(4,7)=(g1-2.*f2*dp+.5*f2*dd)*xy;
      b(4,8)=(g2-2.0*r3*zz*dp+.5*r3*(1.+zz)*dd)*xy;
      b(5,4)=b(4,5);
      b(5,5)=3.*yyzz*ds+(yy+zz-4.0*yyzz)*dp+(xx+yyzz)*dd;
      b(5,6)=(3.0*zz*ds+(1.0-4.0*zz)*dp+(zz-1.0)*dd)*xy;
      b(5,7)=(g1-(1.0+2.0*f2)*dp+(1.0+.5*f2)*dd)*yz;
      b(5,8)=(g2+r3*(f1-zz)*dp-.5*r3*f1*dd)*yz;
      b(6,4)=b(4,6);
      b(6,5)=b(5,6);
      b(6,6)=3.*zzxx*ds+(zz+xx-4.*zzxx)*dp+(yy+zzxx)*dd;
      b(6,7)=(g1+(1.-2.*f2)*dp-(1.-.5*f2)*dd)*zx;
      b(6,8)=(g2+r3*(f1-zz)*dp-.5*r3*f1*dd)*zx;
      b(7,4)=b(4,7);
      b(7,5)=b(5,7);
      b(7,6)=b(6,7);
      b(7,7)=.75*f2*f2*ds+(f1-f2*f2)*dp+(zz+.25*f2*f2)*dd;
      b(7,8)=.5*f2*g2-r3*zz*f2*dp+.25*r3*(1.+zz)*f2*dd;
      b(8,4)=b(4,8);
      b(8,5)=b(5,8);
      b(8,6)=b(6,8);
      b(8,7)=b(7,8);
      b(8,8)=f3*f3*ds+3.*zz*f1*dp+.75*f1*f1*dd;
      return b;
}

Matrix<complex<double>, 9, 9> U(int numat, int isp){

      const double cshift = .575530 - .715751;
      const double delta = 1.113608939931278e-1;
      const double vex = delta*0.5;
      double s0, p0, d0t, d0e;
      //Co
      if (numat == 1){
        s0 =  1.12946 + cshift; // on-site
        p0 =  1.75262 + cshift;
        d0t =  0.5*(0.60547 + 0.60445) + cshift - (2*isp-1)*vex;
        d0e =  0.5*(0.60547 + 0.60445) + cshift - (2*isp-1)*vex;
      }

      if (numat == 2){
	//Cu
        s0 =  0.79466;
        p0 =  1.35351;
        d0t =  0.5*(0.37307 + 0.37180);
        d0e =  0.5*(0.37307 + 0.37180);
      }
      Matrix<complex<double>, 9, 9> result;
      result = Matrix<complex<double>, 9, 9>::Zero();
      result(0,0) = s0;
      result(1,1) = p0;
      result(2,2) = p0;
      result(3,3) = p0;
      result(4,4) = d0t;
      result(5,5) = d0t;
      result(6,6) = d0t;
      result(7,7) = d0e;
      result(8,8) = d0e;
      return result;
}


vector<double> param(int numat, int numnn){
//     THIS ROUTINE IS ATOM DEPENDENT :-
//     -----------------------------------------------------------------
//     The first index in the tight binding parameter arrays refers to
//         1: Bulk Co
//         2: Bulk Cu

//         Interatomic hoppings at end
//     -----------------------------------------------------------------
	vector<double> result;
	result.reserve(10);
	double sss, sps, pps, ppp, sds, pds, pdp, dds, ddp, ddd;

	if (numat == 1){
//     Co :
//     first n.n.

		if (numnn == 1){
			sss = -0.09043;   //  same atom hopping
			sps =  0.13649;
			pps =  0.23748;
			ppp = -0.00142;
			sds = -0.03806;
			pds = -0.04069;
			pdp =  0.02797;
			dds = -0.04213;
			ddp =  0.02976;
			ddd = -0.00684;
		}

//     second n.n.

		if (numnn == 2){
			sss = -0.00337;
			sps =  0.00135;
			pps =  0.02849;
			ppp =  0.01099;
			sds = -0.01119;
			pds = -0.01061;
			pdp =  0.01134;
			dds = -0.00759;
			ddp =  0.00495;
			ddd = -0.00016;
		}
	}

//     -----------------------------------------------------------------
	if (numat == 2){
//     Cu:

//     first n.n.

		if (numnn == 1){
			sss = -0.07518;
			sps =  0.11571;
			pps =  0.19669;
			ppp =  0.01940;
			sds = -0.03107;
			pds = -0.03289;
			pdp =  0.01753;
			dds = -0.02566;
			ddp =  0.01800;
			ddd = -0.00408;
		}

//     second n.n.

		if (numnn == 2){
			sss = -0.00092;
			sps =  0.01221;
			pps =  0.05389;
			ppp =  0.00846;
			sds = -0.00852;
			pds = -0.00536;
			pdp =  0.00321;
			dds = -0.00451;
			ddp =  0.00241;
			ddd = -0.00029;
		}
	}

	result.emplace_back(sss);
	result.emplace_back(sps);
	result.emplace_back(pps);
	result.emplace_back(ppp);
	result.emplace_back(sds);
	result.emplace_back(pds);
	result.emplace_back(pdp);
	result.emplace_back(dds);
	result.emplace_back(ddp);
	result.emplace_back(ddd);

	return result;
}
#endif
