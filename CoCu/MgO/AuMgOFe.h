#ifndef COCUCO_H
#define COCUCO_H

#include <cmath>
#include <Eigen/Dense>
#include <Eigen/StdVector>

using namespace std;
using namespace Eigen;

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

double gmean(double x, double y, double atxdist, double atydist, double atxydist, double order){
	double potx, poty;
	potx = x*pow(atxdist, order);
	poty = y*pow(atydist, order);
	double potxy = gmean(potx, poty);
	double gmean = potxy/(1.*pow(atxydist, order));
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

string species(int numat){
	string atom;
	if (numat == 1)
		atom = "Au";
	else if (numat == 2)
		atom = "Fe";
	else if (numat == 3)
		atom = "Mg";
	else if (numat == 4)
		atom = "O";
	else
		cout<<"Selection not recognised!"<<endl;
	return atom;
}

Matrix<complex<double>, 9, 9> U(int numat, int spin){

      double ryd=13.6058;
      double s0, p0, d0t, d0e;
      /* double ef=0.57553;//TODO isn't this the Papa Fermi level of Cu..? */
      double ef=0.7466;//fermi level of Fe as calculated by me
      double Au_ef = 0.5074;// calculated by me
//    dee=2.1         //   approximately reproduces Yuassas results
      double dee=3.5;         //   classic
      double xmgo_shift=-dee/ryd;
      //TODO Papa Fermi level of Au = 0.5380
      if (numat == 1){
      //Au
	      s0 = 0.56220 - Au_ef;
	      p0 = 1.27897 - Au_ef;
	      d0t = 0.26097 - Au_ef;
	      d0e = 0.25309 - Au_ef;
      }
      else if (numat == 2){//TODO need to figure out relative Fermi shifts
      	if (spin == 0){
      //BULK Fe up:
      	      s0 =  1.13516 - ef;
       	      p0 =  1.81739 - ef;
     	      d0t = 0.64840 - ef;
      	      d0e = 0.62960 - ef;
	}
	else if (spin == 1){
      //BULK Fe down:
	      s0 =  1.14481 - ef;
      	      p0 =  1.80769 - ef;
      	      d0t = 0.78456 - ef;
      	      d0e = 0.75661 - ef;
      	}
      }
      //Mg :
      else if (numat == 3){	
	      s0 = 9.88/ryd + xmgo_shift;
	      p0 = 100. + xmgo_shift;
	      d0t = 100. + xmgo_shift;
	      d0e = 100. + xmgo_shift;
      }
      //O :
      else if (numat == 4){
	      s0 = 100. + xmgo_shift;
	      p0 = -2.03/ryd + xmgo_shift;
	      d0t = 100. + xmgo_shift;
	      d0e = 100. + xmgo_shift;
      }
      else 
	      cout<<"Warning, selection not recognised in function 'U'"<<endl;
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


vector<double> param(int numat, int numnn, int spin){
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
        double ryd=13.6058;

	if (numat == 1){
//     Au :
//     first n.n.

		if (numnn == 1){
			sss = -0.06680;
			pps =  0.17866;
			ppp = -0.01645;
			dds = -0.04971;
			ddp =  0.02624;
			ddd = -0.00457;
			sps =  0.09721;
			sds = -0.04722;
			pds = -0.06399;
			pdp =  0.01896;
		}

//     second n.n.

		if (numnn == 2){
			sss =  0.00277;
			pps =  0.03707;
			ppp = -0.01025;
			dds = -0.00305;
			ddp =  0.00240;
			ddd = -0.00057;
			sps =  0.00261;
			sds = -0.00784;
			pds = -0.00762;
			pdp =  0.00470;
		}
	}

//     -----------------------------------------------------------------
	if (numat == 2){
		if (spin == 0){
//     Fe UP:

//     first n.n.

        		if (numnn == 1){
        		      sss = -0.12950;
        		      pps =  0.25741;
        		      ppp =  0.02422;
        		      dds = -0.04541;
        		      ddp =  0.02714;
        		      ddd = -0.00260;
        		      sps =  0.17363;
        		      sds = -0.06115;
        		      pds = -0.08485;
        		      pdp =  0.01778;
        		}
        
        //     second n.n.
        
        		if (numnn == 2){
        		      sss = -0.02915;
        		      pps =  0.16827;
        		      ppp =  0.04112;
        		      dds = -0.02713;
        		      ddp =  0.00589;
        		      ddd =  0.00060;
        		      sps =  0.06571;
        		      sds = -0.03560;
        		      pds = -0.05473;
        		      pdp = -0.00280;
        		}
        
        //     third n.n.
        
        		if (numnn == 3){
        		      sss =  0.01595;
        		      pps = -0.04985;
        		      ppp =  0.01796;
        		      dds =  0.00112;
        		      ddp =  0.00034;
        		      ddd = -0.00056;
        		      sps = -0.02477;
        		      sds = -0.00073;
        		      pds = -0.00082;
        		      pdp = -0.00241;
        		}
		}
		if (spin == 1){
        //     Fe DOWN:
        
        //     first n.n.
        
        		if (numnn == 1){
        		      sss = -0.13243;
        		      pps =  0.25911;
        		      ppp =  0.02653;
        		      dds = -0.05266;
        		      ddp =  0.03276;
        		      ddd = -0.00286;
        		      sps =  0.17278;
        		      sds = -0.07145;
        		      pds = -0.09702;
        		      pdp =  0.02129;
        		}
        
        //     second n.n.
        
        		if (numnn == 2){
        		      sss = -0.03003;
        		      pps =  0.18256;
        		      ppp =  0.03703;
        		      dds = -0.03396;
        		      ddp =  0.00581;
        		      ddd =  0.00114;
        		      sps =  0.07159;
        		      sds = -0.04075;
        		      pds = -0.06522;
        		      pdp = -0.00467;
        		}
        
        //     third n.n.
        
        		if (numnn == 3){
        		      sss =  0.01589;
        		      pps = -0.04253;
        		      ppp =  0.01538;
        		      dds =  0.00233;
        		      ddp =  0.00013;
        		      ddd = -0.00060;
        		      sps = -0.02306;
        		      sds =  0.00016;
        		      pds =  0.00222;
        		      pdp = -0.00351;
        		}
        	}
	}


//     -----------------------------------------------------------------
	if (numat == 3){
//     Mg:

//     first n.n.

		if (numnn == 1){
		        sss = 0.;
		        sps = 1.1/ryd;
		        pps = 0.;
		        ppp = 0.;
		        sds = 0.;
		        pds = 0.;
		        pdp = 0.;
		        dds = 0.;
		        ddp = 0.;
		        ddd = 0.;
		}

//     second n.n.

		if (numnn == 2){
		        sss =-0.18/ryd;
		        sps = 0.;
		        pps = 0.;
		        ppp = 0.;
		        sds = 0.;
		        pds = 0.;
		        pdp = 0.;
		        dds = 0.;
		        ddp = 0.;
		        ddd = 0.;
		}

//     third n.n.

		if (numnn == 3){
		        sss = 0.;
		        sps = 0.89/ryd;
		        pps = 0.;
		        ppp = 0.;
		        sds = 0.;
		        pds = 0.;
		        pdp = 0.;
		        dds = 0.;
		        ddp = 0.;
		        ddd = 0.;
		}
	}

//     -----------------------------------------------------------------
	if (numat == 4){
//     O:

//     first n.n.

		if (numnn == 1){
		        sss = 0.;
		        sps = 1.1/ryd;
		        pps = 0.;
		        ppp = 0.;
		        sds = 0.;
		        pds = 0.;
		        pdp = 0.;
		        dds = 0.;
		        ddp = 0.;
		        ddd = 0.;
		}

//     second n.n.

		if (numnn == 2){
		        sss = 0.;
		        sps = 0.;
		        pps = 0.65/ryd;
		        ppp =-0.07/ryd;
		        sds = 0.;
		        pds = 0.;
		        pdp = 0.;
		        dds = 0.;
		        ddp = 0.;
		        ddd = 0.;
		}

//     third n.n.

		if (numnn == 3){
		        sss = 0.;
		        sps = 0.89/ryd;
		        pps = 0.;
		        ppp = 0.;
		        sds = 0.;
		        pds = 0.;
		        pdp = 0.;
		        dds = 0.;
		        ddp = 0.;
		        ddd = 0.;
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
