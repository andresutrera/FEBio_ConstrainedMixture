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
//	ADD_PARAMETER(m_C_phi    , FE_RANGE_CLOSED(0.0, 1.0), "C_phi"); // Collagen mass fraction. TODO: add individual collagen fibers mass fractions.
//	ADD_PARAMETER(m_C_k1   , FE_RANGE_GREATER_OR_EQUAL(0.0), "C_k1");
//	ADD_PARAMETER(m_C_k2   , FE_RANGE_GREATER_OR_EQUAL(0.0), "C_k2");
//	ADD_PARAMETER(m_C_kappa, FE_RANGE_CLOSED(0.0, 1.0/3.0), "C_kappa");
//	ADD_PARAMETER(m_C_gdeg , "C_gamma");
//	ADD_PARAMETER(m_C_G[0]    , FE_RANGE_GREATER_OR_EQUAL(1.0), "C_G_z"); // Z pre-stretch
//	ADD_PARAMETER(m_C_G[1]    , FE_RANGE_GREATER_OR_EQUAL(1.0), "C_G_theta");// theta pre-stretch. radial pre-stetch is calculated by incompressibility of E_G
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
	//return (NeoHookeDevStress(mp)+CollagenDevStress(mp));
	return NeoHookeDevStress(mp);
}

//-----------------------------------------------------------------------------
//! Calculates the deviatoric tangent
tens4ds FEUncoupledConstrainedMixture::DevTangent(FEMaterialPoint& mp)
{

	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
  //return (NeoHookeDevTangent(mp)+CollagenDevTangent(mp));
	return NeoHookeDevTangent(mp);
}

mat3ds FEUncoupledConstrainedMixture::NeoHookeDevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	double time = GetFEModel()->GetTime().currentTime;
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

	return (E_c*pow(J,-5.0/3.0)*(b_E-mat3dd(b_E.tr()/3.)));

}


tens4ds FEUncoupledConstrainedMixture::NeoHookeDevTangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		double time = GetFEModel()->GetTime().currentTime;
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
    return (I4*IbE -bExI +IxI*(IbE/3.0))*(2.0*E_c*pow(J,-5.0/3.0)/3.0);

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



//mat3ds FEUncoupledConstrainedMixture::CollagenDevStress(FEMaterialPoint& mp)
//{
//	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
//	double time = GetFEModel()->GetTime().currentTime;

//	double k1 = m_C_k1(mp);
//	double k2 = m_C_k2(mp);
//	double kappa = m_C_kappa(mp);
//	double g = m_C_gdeg(mp)*PI/180;

//	double C_G_z = m_C_G[0](mp); //Copy the material properties.
//	double C_G_t = m_C_G[1](mp);
//	double C_G_r = 1.0/(C_G_z*C_G_t);
//	double G_Time = m_G_Time;
	//Form the stretch deposition tensor E_G with a diagonal matrix.
	//| C_G_r				0						0 		|
	//|   0			C_G_theta				0 		|
	//|		0					0						C_G_z |
//	mat3dd I(1);
//	mat3dd C_G = (time <= G_Time) ? (mat3dd(C_G_r,C_G_t,C_G_z)-I)*time/G_Time+I : mat3dd(C_G_r,C_G_t,C_G_z);

//	mat3d Q = GetLocalCS(mp); //Get material axes.
//	mat3d C_G_rot = (Q*C_G*Q.transpose()); //Rotate the pre-stretch tensor respect to the material axes.

	// determinant of deformation gradient
//	double J = pt.m_J;
	// Evaluate the distortional deformation gradient
//	double Jm13 = pow(J, -1. / 3.);
//	mat3d F = pt.m_F*Jm13;

	// calculate deviatoric left Cauchy-Green tensor: b = F*Ft
//	mat3ds b = pt.DevLeftCauchyGreen();
//	mat3ds C = pt.DevRightCauchyGreen();
//	double I1 = C.tr();

	// Copy the local element basis directions to n
	// According to the r/theta/z, the diagonal fibers are located in the plane formed by e2/e3
//	vec3d n[2];
//	n[0].x = Q[0][2]; n[0].y = Q[1][2]; n[0].z = Q[2][2]; //Z direction
//	n[1].x = Q[0][1]; n[1].y = Q[1][1]; n[1].z = Q[2][1]; // Theta direction

	// Evaluate the structural direction in the current configuration
//	double cg = cos(g); double sg = sin(g);
//	vec3d ar[2],a[2];
//	ar[0] = n[0]*cg + n[1]*sg; a[0] = F*ar[0];
//	ar[1] = n[0]*cg - n[1]*sg; a[1] = F*ar[1];

//	double I40 = ar[0]*(C*ar[0]);
//	double E0 = kappa*(I1-3) + (1-3*kappa)*(I40-1);
//	mat3ds s;
//	s.zero();
//	if (E0 >= 0) {
//			mat3ds h0 = kappa*b + (1-3*kappa)*dyad(a[0]);
//			s += h0*(2.*k1*E0*exp(k2*E0*E0));
//	}
//	double I41 = ar[1]*(C*ar[1]);
//	double E1 = kappa*(I1-3) + (1-3*kappa)*(I41-1);
//	if (E1 >= 0) {
//			mat3ds h1 = kappa*b + (1-3*kappa)*dyad(a[1]);
//			s += h1*(2.*k1*E1*exp(k2*E1*E1));
//	}

//	return s.dev()/J;
//}



//tens4ds FEUncoupledConstrainedMixture::CollagenDevTangent(FEMaterialPoint& mp)
//{
//    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
//		double time = GetFEModel()->GetTime().currentTime;

//		double k1 = m_C_k1(mp);
//		double k2 = m_C_k2(mp);
//		double kappa = m_C_kappa(mp);
//		double g = m_C_gdeg(mp)*PI/180;

//		double C_G_z = m_C_G[0](mp); //Copy the material properties.
//		double C_G_t = m_C_G[1](mp);
//		double C_G_r = 1.0/(C_G_z*C_G_t);
//		double G_Time = m_G_Time;

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
