#include "UncoupledConstrainedMixture.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FEUncoupledConstrainedMixture, FEUncoupledMaterial)
	// #################### Elastin ####################
	ADD_PARAMETER(m_E_phi    , FE_RANGE_CLOSED(0.0, 1.0), "E_phi"); // Elastin mass fraction
	ADD_PARAMETER(m_E_c    , FE_RANGE_GREATER_OR_EQUAL(0.0), "E_c"); // Elastin shear modulus
	ADD_PARAMETER(m_E_G[0]    , FE_RANGE_GREATER_OR_EQUAL(1.0), "E_G_z"); // Z pre-stretch
	ADD_PARAMETER(m_E_G[1]    , FE_RANGE_GREATER_OR_EQUAL(1.0), "E_G_theta");// theta pre-stretch. radial pre-stetch is calculated by incompressibility of E_G
	// #################### Elastin ####################

	// #################### Collagen ####################
	ADD_PARAMETER(m_C_phi    , FE_RANGE_CLOSED(0.0, 1.0), "C_phi"); // Collagen mass fraction. TODO: add individual collagen fibers mass fractions.
	ADD_PARAMETER(m_C_k1   , FE_RANGE_GREATER_OR_EQUAL(0.0), "C_k1");
	ADD_PARAMETER(m_C_k2   , FE_RANGE_GREATER_OR_EQUAL(0.0), "C_k2");
	ADD_PARAMETER(m_C_kappa, FE_RANGE_CLOSED(0.0, 1.0/3.0), "C_kappa");
	ADD_PARAMETER(m_C_gdeg , "C_gamma");
	ADD_PARAMETER(m_C_G    , FE_RANGE_GREATER_OR_EQUAL(1.0), "C_G"); // pre-stretch collagen (in fiber direction)
	//ADD_PARAMETER(m_C_G[0]    , FE_RANGE_GREATER_OR_EQUAL(1.0), "C_G_z"); // Z pre-stretch
	//ADD_PARAMETER(m_C_G[1]    , FE_RANGE_GREATER_OR_EQUAL(1.0), "C_G_theta");// theta pre-stretch. radial pre-stetch is calculated by incompressibility of E_G
	// #################### Collagen ####################

	// #################### Pre-stretch deposition time ####################
	ADD_PARAMETER(m_G_Time    , FE_RANGE_GREATER_OR_EQUAL(0.0), "G_Time");
	// #################### Pre-stretch deposition time ####################

END_FECORE_CLASS();



//-----------------------------------------------------------------------------
//! Calculates the deviatoric stress
mat3ds FEUncoupledConstrainedMixture::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	return (NeoHookeDevStress(mp)+CollagenDevStress(mp));
	//return NeoHookeDevStress(mp);
}

//-----------------------------------------------------------------------------
//! Calculates the deviatoric tangent
tens4ds FEUncoupledConstrainedMixture::DevTangent(FEMaterialPoint& mp)
{

	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
  return (NeoHookeDevTangent(mp)+CollagenDevTangent(mp));
	//return NeoHookeDevTangent(mp);
}

mat3ds FEUncoupledConstrainedMixture::NeoHookeDevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	double time = GetFEModel()->GetTime().currentTime;
  double E_phi = m_E_phi(mp);
	double E_c = m_E_c(mp);
	double E_G_z = m_E_G[0](mp); //Copy the material properties.
	double E_G_t = m_E_G[1](mp);
	double E_G_r = 1.0/(E_G_z*E_G_t);
	double G_Time = m_G_Time;
	//Form the stretch deposition tensor E_G with a diagonal matrix.
	//| E_G_r				0						0 		|
	//|   0			E_G_theta				0 		|
	//|		0					0						E_G_z |
	mat3dd I(1);
	mat3dd E_G = (time <= G_Time) ? (mat3dd(E_G_r,E_G_t,E_G_z)-I)*time/G_Time+I : mat3dd(E_G_r,E_G_t,E_G_z);

	mat3d Q = GetLocalCS(mp); //Get material axes.
	mat3d E_G_rot = (Q*E_G*Q.transpose()); //Rotate the pre-stretch tensor respect to the material axes.

	// deformation gradient tensor
	mat3d F = pt.m_F;

	// determinant of deformation gradient
	double J = pt.m_J;

  // Elastin left Cauchy-Green strain tensor
	mat3ds b_E = (F*E_G_rot*E_G_rot.transpose()*F.transpose()).sym();

	// Evaluate the distortional deformation gradient
	// double Jm13 = pow(J, -1. / 3.);
	// mat3d F = pt.m_F*Jm13;

	// calculate deviatoric left Cauchy-Green tensor: b = F*Ft
	// mat3ds b = pt.LeftCauchyGreen();

	// Evaluate the ground matrix stress
	//mat3ds s = E_G*c*I*E_G*b;
	//mat3ds s = (E_G_rot*E_c*I*E_G_rot*b).sym();

	return E_phi*E_c*pow(J,-5.0/3.0)*(b_E-mat3dd(b_E.tr()/3.));

}


tens4ds FEUncoupledConstrainedMixture::NeoHookeDevTangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		double time = GetFEModel()->GetTime().currentTime;
    double E_phi = m_E_phi(mp);
		double E_c = m_E_c(mp);
		double E_G_z = m_E_G[0](mp);
		double E_G_t = m_E_G[1](mp);
		double E_G_r = 1.0/(E_G_z*E_G_t);
		double G_Time = m_G_Time;

		mat3dd I(1); // Identity
		mat3dd E_G = (time <= G_Time) ? (mat3dd(E_G_r,E_G_t,E_G_z)-I)*time/G_Time+I : mat3dd(E_G_r,E_G_t,E_G_z);
		mat3d Q = GetLocalCS(mp);
		mat3d E_G_rot = (Q*E_G*Q.transpose());

		// deformation gradient tensor
		mat3d F = pt.m_F;

		// determinant of deformation gradient
		double J = pt.m_J;

	  // Elastin left Cauchy-Green strain tensor
		mat3ds b_E = (F*E_G_rot*E_G_rot.transpose()*F.transpose()).sym();

		// trace of b_E
		double IbE=b_E.tr();

		//mat3dd I(1); // Identity

		tens4ds IxI=dyad1s(I);
		tens4ds I4=dyad4s(I);
		tens4ds bExI=dyad1s(b_E,I); // = bE x I + I x bE

    // Evaluate the ground matrix stress
    return (I4*IbE -bExI +IxI*(IbE/3.0))*(2.0*E_c*E_phi*pow(J,-5.0/3.0)/3.0);

		//feLog("S\n");
		//feLog("%f %f %f\n",s.xx(),s.xy(),s.xz());
		//feLog("%f %f %f\n",s.xy(),s.yy(),s.yz());
		//feLog("%f %f %f\n",s.xz(),s.yz(),s.zz());

		//mat3ds sbar = s.dev();

    // Evaluate the elasticity tensor
    //tens4ds IxI = dyad1s(I);
    //tens4ds I4  = dyad4s(I);
    //tens4ds ce;
		//ce.zero();
		//ce = ((I4 - IxI/3.)*s.tr()-dyad1s(sbar,I))*(2./3.);
		//ce = ce.pp(E_G_rot);
		// TODO: check this tangent tensor.
		// With big E_G values, the tangent tensor test fails.
    //return ce / J;
}



mat3ds FEUncoupledConstrainedMixture::CollagenDevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	double time = GetFEModel()->GetTime().currentTime;
  double C_phi = m_C_phi(mp);
	double k1 = m_C_k1(mp); //Copy the material properties.
	double k2 = m_C_k2(mp);
	double kappa = m_C_kappa(mp);
	double g = m_C_gdeg(mp)*PI/180;

  double C_G = m_C_G(mp);
	//double C_G_z = m_C_G[0](mp); //Copy the pre-stretch parameters.
	//double C_G_t = m_C_G[1](mp);
	//double C_G_r = 1.0/(C_G_z*C_G_t);
	double G_Time = m_G_Time;
	// deposition stretch matrix respect to the fiber direction
  mat3dd I(1);
	feLog("C_G: %f\n",m_C_G(mp));
	mat3ds C_G_fib = (time <= G_Time) ? (mat3ds(1/sqrt(C_G),C_G,1/sqrt(C_G),0,0,0)-I)*time/G_Time+I : mat3ds(1/sqrt(C_G),C_G,1/sqrt(C_G),0,0,0);
	feLog("C_G_fib\n");
	feLog("%f %f %f\n",C_G_fib.xx(),C_G_fib.yy(),C_G_fib.zz());
	mat3d T0;
	T0.zero();
	T0[0][0]=1.0;   T0[0][1]=0.0;     T0[0][2]=0;
	T0[1][0]=0.0;   T0[1][1]=sin(g);  T0[1][2]=-cos(g);
	T0[2][0]=0.0;   T0[2][1]=cos(g);  T0[2][2]=sin(g);
	mat3d T1;
	T1.zero();
	T1[0][0]=1.0; T1[0][1]=0.0;     T1[0][2]=0.0;
	T1[1][0]=0.0; T1[1][1]=-sin(g); T1[1][2]=-cos(g);
	T1[2][0]=0.0; T1[2][1]=cos(g);  T1[2][2]=-sin(g);
	mat3d C_G_fib0 = (T0*C_G_fib*T0.transpose());
	mat3d C_G_fib1 = (T1*C_G_fib*T1.transpose());
	feLog("C_G_fib0\n");
	feLog("%f %f %f\n",C_G_fib0[0][0],C_G_fib0[0][1],C_G_fib0[0][2]);
	feLog("%f %f %f\n",C_G_fib0[1][0],C_G_fib0[1][1],C_G_fib0[1][2]);
	feLog("%f %f %f\n",C_G_fib0[2][0],C_G_fib0[2][1],C_G_fib0[2][2]);
	feLog("C_G_fib1\n");
	feLog("%f %f %f\n",C_G_fib1[0][0],C_G_fib1[0][1],C_G_fib1[0][2]);
	feLog("%f %f %f\n",C_G_fib1[1][0],C_G_fib1[1][1],C_G_fib1[1][2]);
	feLog("%f %f %f\n",C_G_fib1[2][0],C_G_fib1[2][1],C_G_fib1[2][2]);
	//Form the stretch deposition tensor C_G with a diagonal matrix.
	//| C_G_r				0						0 		|
	//|   0			C_G_theta				0 		|
	//|		0					0						C_G_z |
	//mat3dd I(1);
	//mat3dd C_G_fib = (time <= G_Time) ? (mat3dd(C_G_r,C_G_t,C_G_z)-I)*time/G_Time+I : mat3dd(C_G_r,C_G_t,C_G_z);

	mat3d Q = GetLocalCS(mp); //Get material axes.
	mat3d C_G_rot0 = (Q*C_G_fib0*Q.transpose()); //Rotate the pre-stretch tensor respect to the material axes.
	mat3d C_G_rot1 = (Q*C_G_fib1*Q.transpose()); //Rotate the pre-stretch tensor respect to the material axes.

	feLog("C_G_rot0\n");
	feLog("%f %f %f\n",C_G_rot0[0][0],C_G_rot0[0][1],C_G_rot0[0][2]);
	feLog("%f %f %f\n",C_G_rot0[1][0],C_G_rot0[1][1],C_G_rot0[1][2]);
	feLog("%f %f %f\n",C_G_rot0[2][0],C_G_rot0[2][1],C_G_rot0[2][2]);
	feLog("C_G_rot1\n");
	feLog("%f %f %f\n",C_G_rot1[0][0],C_G_rot1[0][1],C_G_rot1[0][2]);
	feLog("%f %f %f\n",C_G_rot1[1][0],C_G_rot1[1][1],C_G_rot1[1][2]);
	feLog("%f %f %f\n",C_G_rot1[2][0],C_G_rot1[2][1],C_G_rot1[2][2]);

	// determinant of deformation gradient
	double J = pt.m_J;

	// Evaluate the distortional deformation gradient
  double Jm13 = pow(J, -1. / 3.);
  mat3d F = pt.m_F*Jm13;

	// Collagen deformation gradient
	mat3d F_C0 = F*C_G_rot0;
	mat3d F_C1 = F*C_G_rot1;

  // Collagen left Cauchy-Green strain tensor
  mat3ds b_C0 = (F_C0*F_C0.transpose()).sym();
	mat3ds b_C1 = (F_C1*F_C1.transpose()).sym();

  // Collagen right Cauchy-Green strain tensor
  mat3ds C_C0 = (F_C0.transpose()*F_C0).sym();
	mat3ds C_C1 = (F_C1.transpose()*F_C1).sym();

	feLog("C_C0\n");
	feLog("%f %f %f\n",C_C0.xx(),C_C0.xy(),C_C0.xz());
	feLog("%f %f %f\n",C_C0.xy(),C_C0.yy(),C_C0.yz());
	feLog("%f %f %f\n",C_C0.xz(),C_C0.yz(),C_C0.zz());

	// trace of b_C
	double I1bC0=b_C0.tr();
	double I1bC1=b_C1.tr();

	// calculate deviatoric left Cauchy-Green tensor: b = F*Ft
//	mat3ds b = pt.DevLeftCauchyGreen();
//	mat3ds C = pt.DevRightCauchyGreen();
//	double I1 = C.tr();

	// Copy the local element basis directions to n
	// According to the r/theta/z, the diagonal fibers are located in the plane formed by e2/e3
	vec3d n[2];
	n[0].x = Q[0][2]; n[0].y = Q[1][2]; n[0].z = Q[2][2]; //Z direction
	n[1].x = Q[0][1]; n[1].y = Q[1][1]; n[1].z = Q[2][1]; // Theta direction

	// Evaluate the structural direction in the current configuration
	double cg = cos(g); double sg = sin(g);
	vec3d ar[2],a[2];
	ar[0] = n[0]*cg + n[1]*sg; a[0] = F_C0*ar[0];
	ar[1] = n[0]*cg - n[1]*sg; a[1] = F_C1*ar[1];

	double I4bC = ar[0]*(C_C0*ar[0]);
	double E0 = ((kappa*(I1bC0) + (1-3*kappa)*(I4bC)))-1.0;
	mat3ds s;
	s.zero();
	mat3ds h0;
	if (E0 >= 0) {
			h0 = kappa*b_C0 + (1-3*kappa)*dyad(a[0]);
			s += h0*(2.*k1*E0*exp(k2*E0*E0));
	}
	double I6bC = ar[1]*(C_C1*ar[1]);
	double E1 = ((kappa*(I1bC1) + (1-3*kappa)*(I6bC)))-1.0;
	mat3ds h1;
	if (E1 >= 0) {
			h1 = kappa*b_C1 + (1-3*kappa)*dyad(a[1]);
			s += h1*(2.*k1*E1*exp(k2*E1*E1));
	}
	feLog("S\n");
	feLog("%f %f %f\n",s.xx(),s.xy(),s.xz());
	feLog("%f %f %f\n",s.xy(),s.yy(),s.yz());
	feLog("%f %f %f\n",s.xz(),s.yz(),s.zz());
	//change return C_phi*s.dev()/J;
	return C_phi*s.dev()/J;
}



tens4ds FEUncoupledConstrainedMixture::CollagenDevTangent(FEMaterialPoint& mp)
{
 FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
 double time = GetFEModel()->GetTime().currentTime;
 double C_phi = m_C_phi(mp);
 double k1 = m_C_k1(mp); //Copy the material properties.
 double k2 = m_C_k2(mp);
 double kappa = m_C_kappa(mp);
 double g = m_C_gdeg(mp)*PI/180;

 double C_G = m_C_G(mp);
 //double C_G_z = m_C_G[0](mp); //Copy the pre-stretch parameters.
 //double C_G_t = m_C_G[1](mp);
 //double C_G_r = 1.0/(C_G_z*C_G_t);
 double G_Time = m_G_Time;

 mat3dd I(1);
 mat3dd C_G_fib = (time <= G_Time) ? (mat3dd(1/sqrt(C_G),C_G,1/sqrt(C_G))-I)*time/G_Time+I : mat3dd(1/sqrt(C_G),C_G,1/sqrt(C_G));
 mat3d T0;
 T0.zero();
 T0[0][0]=1.0;   T0[0][1]=0.0;     T0[0][2]=0.0;
 T0[1][0]=0.0;   T0[1][1]=sin(g);  T0[1][2]=-cos(g);
 T0[2][0]=0.0;   T0[2][1]=cos(g);  T0[2][2]=sin(g);
 mat3d T1;
 T1.zero();
 T1[0][0]=1.0; T1[0][1]=0.0;     T1[0][2]=0.0;
 T1[1][0]=0.0; T1[1][1]=-sin(g); T1[1][2]=-cos(g);
 T1[2][0]=0.0; T1[2][1]=cos(g);  T1[2][2]=-sin(g);

 mat3ds C_G_fib0 = (T0*C_G_fib*T0.transpose()).sym();
 mat3ds C_G_fib1 = (T1*C_G_fib*T1.transpose()).sym();

 mat3d Q = GetLocalCS(mp); //Get material axes.
 mat3ds C_G_rot0 = (Q*C_G_fib0*Q.transpose()).sym(); //Rotate the pre-stretch tensor respect to the material axes.
 mat3ds C_G_rot1 = (Q*C_G_fib1*Q.transpose()).sym(); //Rotate the pre-stretch tensor respect to the material axes.


 // determinant of deformation gradient
 double J = pt.m_J;

 // Evaluate the distortional deformation gradient
 double Jm13 = pow(J, -1. / 3.);
 mat3d F = pt.m_F*Jm13;

 // Collagen deformation gradient
 mat3d F_C0 = F*C_G_rot0;
 mat3d F_C1 = F*C_G_rot1;

 // Collagen left Cauchy-Green strain tensor
 mat3ds b_C0 = (F_C0*F_C0.transpose()).sym();
 mat3ds b_C1 = (F_C1*F_C1.transpose()).sym();

 // Collagen right Cauchy-Green strain tensor
 mat3ds C_C0 = (F_C0.transpose()*F_C0).sym();
 mat3ds C_C1 = (F_C1.transpose()*F_C1).sym();

 // trace of b_C
 double I1bC0=b_C0.tr();
 double I1bC1=b_C1.tr();

 // Copy the local element basis directions to n
 // According to the r/theta/z, the diagonal fibers are located in the plane formed by e2/e3
 vec3d n[2];
 n[0].x = Q[0][2]; n[0].y = Q[1][2]; n[0].z = Q[2][2]; //Z direction
 n[1].x = Q[0][1]; n[1].y = Q[1][1]; n[1].z = Q[2][1]; // Theta direction

 // Evaluate the structural direction in the current configuration
 double cg = cos(g); double sg = sin(g);
 vec3d ar[2],a[2];
 ar[0] = n[0]*cg + n[1]*sg; a[0] = F_C0*ar[0];
 ar[1] = n[0]*cg - n[1]*sg; a[1] = F_C1*ar[1];

 double I4bC = ar[0]*(C_C0*ar[0]);
 double E0 = ((kappa*(I1bC0) + (1-3*kappa)*(I4bC)))-1.0;
 mat3ds s;
 s.zero();
 mat3ds h0;
 if (E0 >= 0) {
 		h0 = kappa*b_C0 + (1-3*kappa)*dyad(a[0]);
 		s += h0*(2.*k1*E0*exp(k2*E0*E0));
 }
 double I6bC = ar[1]*(C_C1*ar[1]);
 double E1 = ((kappa*(I1bC1) + (1-3*kappa)*(I6bC)))-1.0;
 mat3ds h1;
 if (E1 >= 0) {
 		h1 = kappa*b_C1 + (1-3*kappa)*dyad(a[1]);
 		s += h1*(2.*k1*E1*exp(k2*E1*E1));
 }

 mat3ds sbar = s.dev();

 // Evaluate the elasticity tensor
 tens4ds IxI = dyad1s(I);
 tens4ds I4  = dyad4s(I);
 tens4ds ce  = ((I4 - IxI/3.)*s.tr()-dyad1s(sbar,I))*(2.0/3.0);
     if (E0 >= 0) ce += dyad1s(h0.dev())*(4.*k1*(1 + 2 * k2*E0*E0)*exp(k2*E0*E0));
 		 if (E1 >= 0) ce += dyad1s(h1.dev())*(4.*k1*(1 + 2 * k2*E1*E1)*exp(k2*E1*E1));

 return C_phi*ce / J;
}

//		mat3dd I(1);
//		mat3dd E_G = (time <= G_Time) ? (mat3dd(C_G_r,C_G_t,C_G_z)-I)*time/G_Time+I : mat3dd(C_G_r,C_G_t,C_G_z);
//		mat3d Q = GetLocalCS(mp);
//		mat3d C_G_rot = (Q*E_G*Q.transpose());

    // determinant of deformation gradient
//    double J = pt.m_J;
    // Evaluate the distortional deformation gradient
//	  double Jm13 = pow(J, -1. / 3.);
//    mat3d F = pt.m_F*Jm13;

    // calculate deviatoric left Cauchy-Green tensor: b = F*Ft
//		mat3ds b = pt.DevLeftCauchyGreen();;
//    mat3ds C = pt.DevRightCauchyGreen();
//    double I1 = C.tr();

		// Copy the local element basis directions to n
//		vec3d n[2];
//		n[0].x = Q[0][2]; n[0].y = Q[1][2]; n[0].z = Q[2][2]; //Z direction
//		n[1].x = Q[0][1]; n[1].y = Q[1][1]; n[1].z = Q[2][1]; // Theta direction

		// Evaluate the structural direction in the current configuration
//		double cg = cos(g); double sg = sin(g);
//		vec3d ar[2],a[2];
//		ar[0] = n[0]*cg + n[1]*sg; a[0] = F*ar[0];
//		ar[1] = n[0]*cg - n[1]*sg; a[1] = F*ar[1];
//		double I40 = ar[0]*(C*ar[0]);

//		double E0 = kappa*(I1-3) + (1-3*kappa)*(I40-1);
//		mat3ds h0;
//		if (E0 >= 0) {
//				h0 = kappa*b + (1-3*kappa)*dyad(a[0]);
//		}

//		double I41 = ar[1]*(C*ar[1]);

//		double E1 = kappa*(I1-3) + (1-3*kappa)*(I41-1);
//		mat3ds h1;
//		if (E1 >= 0) {
//				h1 = kappa*b + (1-3*kappa)*dyad(a[1]);
//		}
    // Evaluate the elasticity tensor
//		tens4ds ce;
//		ce.zero(); // Zero eleasticity tensor to avoid null values.
//		if (E0 >= 0) ce += dyad1s(h0.dev())*(4.*k1*(1 + 2 * k2*E0*E0)*exp(k2*E0*E0));
//		if (E1 >= 0) ce += dyad1s(h1.dev())*(4.*k1*(1 + 2 * k2*E1*E1)*exp(k2*E1*E1));
		//ce = ce.pp(C_G_rot);
//    return ce / J;
//}
///////////////////////////////////////////////////
// TODO: Add the parameter validation, specially over the mass fractions.

// bool FEUncoupledConstrainedMixture::Init()
// {
// 	if (FEUncoupledConstrainedMixture::Init() == false) return false;
// 	return true;
// }
//
// bool FEUncoupledConstrainedMixture::Validate()
// {
// 	if (FEUncoupledConstrainedMixture::Validate() == false) return false;
//
// 	//TODO: How to validate parameters?
// 	// if( (m_E_phi+m_C_phi) > 1.0){
// 	// 	return false;
// 	// 	double msum = m_E_phi+m_C_phi;
// 	// 	feLogError("Error in mass fraction.phiSum = %f.", msum);
// 	// }
//
// 	return true;
// }
