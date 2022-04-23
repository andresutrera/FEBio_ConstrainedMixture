#include "UncoupledConstrainedMixture.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FEUncoupledConstrainedMixture, FEUncoupledMaterial)
	// #################### F ####################
	ADD_PARAMETER(m_F_r    , FE_RANGE_GREATER_OR_EQUAL(0.0), "F_r");
	ADD_PARAMETER(m_F_theta    , FE_RANGE_GREATER_OR_EQUAL(0.0), "F_theta");
	ADD_PARAMETER(m_F_z    , FE_RANGE_GREATER_OR_EQUAL(0.0), "F_z");
	// #################### F ####################

	// #################### Elastin ####################
	ADD_PARAMETER(m_E_phi    , FE_RANGE_CLOSED(0.0, 1.0), "E_phi"); // Elastin mass fraction
	ADD_PARAMETER(m_E_c    , FE_RANGE_GREATER_OR_EQUAL(0.0), "E_c"); // Elastin shear modulus
	ADD_PARAMETER(m_E_G[0]    , FE_RANGE_GREATER_OR_EQUAL(0.0), "E_G_z"); // Z pre-stretch
	ADD_PARAMETER(m_E_G[1]    , FE_RANGE_GREATER_OR_EQUAL(0.0), "E_G_theta");// theta pre-stretch. radial pre-stetch is calculated by incompressibility of E_G
	// #################### Elastin ####################

	// #################### Collagen ####################
	ADD_PARAMETER(m_C_phi    , FE_RANGE_CLOSED(0.0, 1.0), "C_phi"); // Collagen mass fraction. TODO: add individual collagen fibers mass fractions.
	ADD_PARAMETER(m_C_k1   , FE_RANGE_GREATER_OR_EQUAL(0.0), "C_k1");
	ADD_PARAMETER(m_C_k2   , FE_RANGE_GREATER_OR_EQUAL(0.0), "C_k2");
	ADD_PARAMETER(m_C_kappa, FE_RANGE_CLOSED(0.0, 1.0/3.0), "C_kappa");
	ADD_PARAMETER(m_C_gdeg , "C_gamma");
	ADD_PARAMETER(m_C_G    , FE_RANGE_GREATER_OR_EQUAL(0.0), "C_G"); // pre-stretch collagen (in fiber direction)
	// #################### Collagen ####################

	// #################### Smooth Muscle ####################
	ADD_PARAMETER(m_SM_phi    , FE_RANGE_CLOSED(0.0, 1.0), "SM_phi"); // Collagen mass fraction. TODO: add individual collagen fibers mass fractions.
	ADD_PARAMETER(m_SM_k1   , FE_RANGE_GREATER_OR_EQUAL(0.0), "SM_k1");
	ADD_PARAMETER(m_SM_k2   , FE_RANGE_GREATER_OR_EQUAL(0.0), "SM_k2");
	ADD_PARAMETER(m_SM_kappa, FE_RANGE_CLOSED(0.0, 1.0/3.0), "SM_kappa");
	ADD_PARAMETER(m_SM_gdeg , "SM_gamma");
	ADD_PARAMETER(m_SM_G    , FE_RANGE_GREATER_OR_EQUAL(0.0), "SM_G"); // pre-stretch collagen (in fiber direction)
	// #################### Smooth Muscle ####################

	// #################### Pre-stretch deposition time ####################
	ADD_PARAMETER(m_G_Time    , FE_RANGE_GREATER_OR_EQUAL(0.0), "G_Time");
	// #################### Pre-stretch deposition time ####################

END_FECORE_CLASS();



//-----------------------------------------------------------------------------
//! Calculates the global deviatoric stress
mat3ds FEUncoupledConstrainedMixture::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	return (NeoHookeDevStress(mp)+CollagenDevStress(mp)+SmoothMuscleDevStress(mp));
}

//-----------------------------------------------------------------------------
//! Calculates the global deviatoric tangent
tens4ds FEUncoupledConstrainedMixture::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
  return (NeoHookeDevTangent(mp)+CollagenDevTangent(mp)+SmoothMuscleDevTangent(mp));
}

//-----------------------------------------------------------------------------
//! Calculates the elastin deviatoric stress (Neo-Hookean material)
mat3ds FEUncoupledConstrainedMixture::NeoHookeDevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	double time = GetFEModel()->GetTime().currentTime;
  double E_phi = m_E_phi(mp); //Copy the mass fraction.
	double E_c = m_E_c(mp); //Copy the material properties.
	double E_G_z = m_E_G[0](mp); //Copy the deposition stretch values.
	double E_G_t = m_E_G[1](mp);
	double E_G_r = 1.0/(E_G_z*E_G_t); //Defined by incompressibility condition
	double G_Time = m_G_Time; //Define the overall time where an "incremental" deposition stretch is applied

	double F_r = m_F_r(mp);
	double F_theta = m_F_theta(mp);
	double F_z = m_F_z(mp);

	//Form the stretch deposition tensor E_G with a diagonal matrix.
	//| E_G_r				0						0 		|
	//|   0			E_G_theta				0 		|
	//|		0					0						E_G_z |
	mat3dd I(1); // Identity

	//Current deposition stretch tensor
	mat3dd E_G = (time <= G_Time) ? (mat3dd(E_G_r,E_G_t,E_G_z)-I)*time/G_Time+I : mat3dd(E_G_r,E_G_t,E_G_z);
	mat3dd F_ro = (time <= G_Time) ? (mat3dd(F_r,F_theta,F_z)-I)*time/G_Time+I : mat3dd(F_r,F_theta,F_z);

	E_G = E_G*F_ro;

	mat3d Q = GetLocalCS(mp); //Get material axes.
	mat3d E_G_rot = (Q*E_G*Q.transpose()); //Rotate the pre-stretch tensor respect to the material axes.

	// deformation gradient tensor
	mat3d F = pt.m_F;

	// determinant of deformation gradient
	double J = pt.m_J;

  // Elastin left Cauchy-Green strain tensor (F_E*F_E^T)
	mat3ds b_E = (F*E_G_rot*E_G_rot.transpose()*F.transpose()).sym();

	return E_phi*E_c*pow(J,-5.0/3.0)*(b_E-mat3dd(b_E.tr()/3.));

}

//-----------------------------------------------------------------------------
//! Calculates the elastin deviatoric tangent (Neo-Hookean material)
tens4ds FEUncoupledConstrainedMixture::NeoHookeDevTangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		double time = GetFEModel()->GetTime().currentTime;
    double E_phi = m_E_phi(mp); //Copy the mass fraction.
		double E_c = m_E_c(mp); //Copy the material properties.
		double E_G_z = m_E_G[0](mp); //Copy the deposition stretch values.
		double E_G_t = m_E_G[1](mp);
		double E_G_r = 1.0/(E_G_z*E_G_t); //Defined by incompressibility condition
		double G_Time = m_G_Time; //Define the overall time where an "incremental" deposition stretch is applied

		double F_r = m_F_r(mp);
		double F_theta = m_F_theta(mp);
		double F_z = m_F_z(mp);

		//Form the stretch deposition tensor E_G with a diagonal matrix.
		//| E_G_r				0						0 		|
		//|   0			E_G_theta				0 		|
		//|		0					0						E_G_z |
		mat3dd I(1); // Identity

		//Current deposition stretch tensor
		mat3dd E_G = (time <= G_Time) ? (mat3dd(E_G_r,E_G_t,E_G_z)-I)*time/G_Time+I : mat3dd(E_G_r,E_G_t,E_G_z);
		mat3dd F_ro = (time <= G_Time) ? (mat3dd(F_r,F_theta,F_z)-I)*time/G_Time+I : mat3dd(F_r,F_theta,F_z);

		E_G = E_G*F_ro;
		mat3d Q = GetLocalCS(mp); //Get material axes
		mat3d E_G_rot = (Q*E_G*Q.transpose()); //Rotate the pre-stretch tensor respect to the material axes.

		// deformation gradient tensor
		mat3d F = pt.m_F;

		// determinant of deformation gradient
		double J = pt.m_J;

	  // Elastin left Cauchy-Green strain tensor (F_E*F_E^T)
		mat3ds b_E = (F*E_G_rot*E_G_rot.transpose()*F.transpose()).sym();

		// trace of b_E
		double IbE=b_E.tr();

		tens4ds IxI=dyad1s(I); // I x I
		tens4ds I4=dyad4s(I); // =-dC^-1/dC
		tens4ds bExI=dyad1s(b_E,I); // = bE x I + I x bE

    // Evaluate the ground matrix stress
    return (I4*IbE -bExI +IxI*(IbE/3.0))*(2.0*E_c*E_phi*pow(J,-5.0/3.0)/3.0);
}

//-----------------------------------------------------------------------------
//! Calculates the collagen deviatoric stress (Holzapfel-Gasser-Ogden material)
mat3ds FEUncoupledConstrainedMixture::CollagenDevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	double time = GetFEModel()->GetTime().currentTime;
  double C_phi = m_C_phi(mp); //Copy the mass fraction.
	double k1 = m_C_k1(mp); //Copy the material properties.
	double k2 = m_C_k2(mp);
	double kappa = m_C_kappa(mp);
	double g = m_C_gdeg(mp)*PI/180;

	double F_r = m_F_r(mp);
	double F_theta = m_F_theta(mp);
	double F_z = m_F_z(mp);

	double C_G = m_C_G(mp); //Copy the deposition stretch value.
	double G_Time = m_G_Time; //Define the overall time where an "incremental" deposition stretch is applied

  mat3dd I(1); // Identity

	//Current deposition stretch tensor oriented respect to an auxiliar system of coordinates.
	//One of the axes of this system of coordinates is coincident with the fiber direction
  mat3ds C_G_fib = (time <= G_Time) ? (mat3ds(1/sqrt(C_G),C_G,1/sqrt(C_G),0,0,0)-I)*time/G_Time+I : mat3ds(1/sqrt(C_G),C_G,1/sqrt(C_G),0,0,0);
  double cg = cos(g); double sg = sin(g);

	//Rotation matrix from the auxiliar system of coordinates to cylindrical coordinate system (r,theta,zeta)
  //Family 0
	mat3d T0;
	T0.zero();
	T0[0][0]=1.0;   T0[0][1]=0.0; T0[0][2]=0.0;
	T0[1][0]=0.0;   T0[1][1]=sg;  T0[1][2]=-cg;
	T0[2][0]=0.0;   T0[2][1]=cg;  T0[2][2]=sg;

	//Family 1
	mat3d T1;
	T1.zero();
	T1[0][0]=1.0; T1[0][1]=0.0; T1[0][2]=0.0;
	T1[1][0]=0.0; T1[1][1]=-sg; T1[1][2]=-cg;
	T1[2][0]=0.0; T1[2][1]=cg;  T1[2][2]=-sg;

	//Current deposition stretch tensor oriented respect to cylindrical coordinates (r,theta,zeta).
	mat3d C_G_fib0 = (T0*C_G_fib*T0.transpose());
	mat3d C_G_fib1 = (T1*C_G_fib*T1.transpose());

	mat3dd F_ro = (time <= G_Time) ? (mat3dd(F_r,F_theta,F_z)-I)*time/G_Time+I : mat3dd(F_r,F_theta,F_z);
	C_G_fib0 = C_G_fib0*F_ro;
	C_G_fib1 = C_G_fib1*F_ro;

	mat3d Q = GetLocalCS(mp); //Get material axes.

  //Current deposition stretch tensor oriented respect to cartesian coordinates (x,y,z).
	mat3d C_G_rot0 = (Q*C_G_fib0*Q.transpose());
	mat3d C_G_rot1 = (Q*C_G_fib1*Q.transpose());

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
	vec3d ar[2],a[2];
	ar[0] = n[0]*cg + n[1]*sg; a[0] = F_C0*ar[0];
	ar[1] = n[0]*cg - n[1]*sg; a[1] = F_C1*ar[1];

	double I4bC = ar[0]*(C_C0*ar[0]);
	if(I4bC < 1.0) {
		I4bC = 1.0;
	}
	double E0 = ((kappa*(I1bC0) + (1-3*kappa)*(I4bC)))-1.0;
	mat3ds s;
	s.zero();
	mat3ds h0;
	h0 = kappa*b_C0 + (1-3*kappa)*dyad(a[0]);
	s += h0*(2.*k1*E0*exp(k2*E0*E0));

	double I6bC = ar[1]*(C_C1*ar[1]);
	if(I6bC < 1.0) {
		I6bC = 1.0;
	}
	double E1 = ((kappa*(I1bC1) + (1-3*kappa)*(I6bC)))-1.0;
	mat3ds h1;
	h1 = kappa*b_C1 + (1-3*kappa)*dyad(a[1]);
	s += h1*(2.*k1*E1*exp(k2*E1*E1));

	return C_phi*s.dev()/J;
}

//-----------------------------------------------------------------------------
//! Calculates the collagen deviatoric tangent (Holzapfel-Gasser-Ogden material)
tens4ds FEUncoupledConstrainedMixture::CollagenDevTangent(FEMaterialPoint& mp)
{
 FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
 double time = GetFEModel()->GetTime().currentTime;
 double C_phi = m_C_phi(mp); //Copy the mass fraction.
 double k1 = m_C_k1(mp); //Copy the material properties.
 double k2 = m_C_k2(mp);
 double kappa = m_C_kappa(mp);
 double g = m_C_gdeg(mp)*PI/180;

 double F_r = m_F_r(mp);
 double F_theta = m_F_theta(mp);
 double F_z = m_F_z(mp);

 double C_G = m_C_G(mp); //Copy the deposition stretch value.
 double G_Time = m_G_Time; //Define the overall time where an "incremental" deposition stretch is applied

 mat3dd I(1); // Identity

 //Current deposition stretch tensor oriented respect to an auxiliar system of coordinates.
 //One of the axes of this system of coordinates is coincident with the fiber direction
 mat3dd C_G_fib = (time <= G_Time) ? (mat3dd(1/sqrt(C_G),C_G,1/sqrt(C_G))-I)*time/G_Time+I : mat3dd(1/sqrt(C_G),C_G,1/sqrt(C_G));
 double cg = cos(g); double sg = sin(g);

 //Rotation matrix from the auxiliar system of coordinates to cylindrical coordinate system (r,theta,zeta)
 //Family 0
 mat3d T0;
 T0.zero();
 T0[0][0]=1.0;   T0[0][1]=0.0;     T0[0][2]=0.0;
 T0[1][0]=0.0;   T0[1][1]=sg;  T0[1][2]=-cg;
 T0[2][0]=0.0;   T0[2][1]=cg;  T0[2][2]=sg;

 //Family 1
 mat3d T1;
 T1.zero();
 T1[0][0]=1.0; T1[0][1]=0.0;     T1[0][2]=0.0;
 T1[1][0]=0.0; T1[1][1]=-sg; T1[1][2]=-cg;
 T1[2][0]=0.0; T1[2][1]=cg;  T1[2][2]=-sg;

 //Current deposition stretch tensor oriented respect to cylindrical coordinates (r,theta,zeta).
 mat3ds C_G_fib0 = (T0*C_G_fib*T0.transpose()).sym();
 mat3ds C_G_fib1 = (T1*C_G_fib*T1.transpose()).sym();

 mat3dd F_ro = (time <= G_Time) ? (mat3dd(F_r,F_theta,F_z)-I)*time/G_Time+I : mat3dd(F_r,F_theta,F_z);
 C_G_fib0 = C_G_fib0*F_ro;
 C_G_fib1 = C_G_fib1*F_ro;

 mat3d Q = GetLocalCS(mp); //Get material axes.

 //Current deposition stretch tensor oriented respect to cartesian coordinates (x,y,z).
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
 vec3d ar[2],a[2];
 ar[0] = n[0]*cg + n[1]*sg; a[0] = F_C0*ar[0];
 ar[1] = n[0]*cg - n[1]*sg; a[1] = F_C1*ar[1];

 double I4bC = ar[0]*(C_C0*ar[0]);
 if(I4bC < 1.0) {
	 I4bC = 1.0;
 }
 double E0 = ((kappa*(I1bC0) + (1-3*kappa)*(I4bC)))-1.0;
 mat3ds s;
 s.zero();
 mat3ds h0;
 h0 = kappa*b_C0 + (1-3*kappa)*dyad(a[0]);
 s += h0*(2.*k1*E0*exp(k2*E0*E0));

 double I6bC = ar[1]*(C_C1*ar[1]);
 if(I6bC < 1.0) {
	 I6bC = 1.0;
 }
 double E1 = ((kappa*(I1bC1) + (1-3*kappa)*(I6bC)))-1.0;
 mat3ds h1;
 h1 = kappa*b_C1 + (1-3*kappa)*dyad(a[1]);
 s += h1*(2.*k1*E1*exp(k2*E1*E1));


 mat3ds sbar = s.dev();

 // Evaluate the elasticity tensor
 tens4ds IxI = dyad1s(I);
 tens4ds I4  = dyad4s(I);
 tens4ds ce  = ((I4 - IxI/3.)*s.tr()-dyad1s(sbar,I))*(2.0/3.0);
     if (E0 >= 0) ce += dyad1s(h0.dev())*(4.*k1*(1 + 2 * k2*E0*E0)*exp(k2*E0*E0));
 		 if (E1 >= 0) ce += dyad1s(h1.dev())*(4.*k1*(1 + 2 * k2*E1*E1)*exp(k2*E1*E1));

 return C_phi*ce / J;
}

//######################## SMOOTH MUSCLE ###############################
//! Calculates the SM deviatoric stress (Holzapfel-Gasser-Ogden material)
mat3ds FEUncoupledConstrainedMixture::SmoothMuscleDevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	double time = GetFEModel()->GetTime().currentTime;
  double SM_phi = m_SM_phi(mp); //Copy the mass fraction.
	double k1 = m_SM_k1(mp); //Copy the material properties.
	double k2 = m_SM_k2(mp);
	double kappa = m_SM_kappa(mp);
	double g = m_SM_gdeg(mp)*PI/180;

	double F_r = m_F_r(mp);
  double F_theta = m_F_theta(mp);
  double F_z = m_F_z(mp);

	double SM_G = m_SM_G(mp); //Copy the deposition stretch value.
	double G_Time = m_G_Time; //Define the overall time where an "incremental" deposition stretch is applied

  mat3dd I(1); // Identity

	//Current deposition stretch tensor oriented respect to an auxiliar system of coordinates.
	//One of the axes of this system of coordinates is coincident with the fiber direction
  mat3ds SM_G_fib = (time <= G_Time) ? (mat3ds(1/sqrt(SM_G),SM_G,1/sqrt(SM_G),0,0,0)-I)*time/G_Time+I : mat3ds(1/sqrt(SM_G),SM_G,1/sqrt(SM_G),0,0,0);
  double cg = cos(g); double sg = sin(g);

	//Rotation matrix from the auxiliar system of coordinates to cylindrical coordinate system (r,theta,zeta)
  //Family 0
	mat3d T0;
	T0.zero();
	T0[0][0]=1.0;   T0[0][1]=0.0; T0[0][2]=0.0;
	T0[1][0]=0.0;   T0[1][1]=sg;  T0[1][2]=-cg;
	T0[2][0]=0.0;   T0[2][1]=cg;  T0[2][2]=sg;

	//Family 1
	mat3d T1;
	T1.zero();
	T1[0][0]=1.0; T1[0][1]=0.0; T1[0][2]=0.0;
	T1[1][0]=0.0; T1[1][1]=-sg; T1[1][2]=-cg;
	T1[2][0]=0.0; T1[2][1]=cg;  T1[2][2]=-sg;

	//Current deposition stretch tensor oriented respect to cylindrical coordinates (r,theta,zeta).
	mat3d SM_G_fib0 = (T0*SM_G_fib*T0.transpose());
	mat3d SM_G_fib1 = (T1*SM_G_fib*T1.transpose());

	mat3dd F_ro = (time <= G_Time) ? (mat3dd(F_r,F_theta,F_z)-I)*time/G_Time+I : mat3dd(F_r,F_theta,F_z);
  SM_G_fib0 = SM_G_fib0*F_ro;
  SM_G_fib1 = SM_G_fib1*F_ro;

	mat3d Q = GetLocalCS(mp); //Get material axes.

  //Current deposition stretch tensor oriented respect to cartesian coordinates (x,y,z).
	mat3d SM_G_rot0 = (Q*SM_G_fib0*Q.transpose());
	mat3d SM_G_rot1 = (Q*SM_G_fib1*Q.transpose());

	// determinant of deformation gradient
	double J = pt.m_J;

	// Evaluate the distortional deformation gradient
  double Jm13 = pow(J, -1. / 3.);
  mat3d F = pt.m_F*Jm13;

	// Collagen deformation gradient
	mat3d F_C0 = F*SM_G_rot0;
	mat3d F_C1 = F*SM_G_rot1;

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
	vec3d ar[2],a[2];
	ar[0] = n[0]*cg + n[1]*sg; a[0] = F_C0*ar[0];
	ar[1] = n[0]*cg - n[1]*sg; a[1] = F_C1*ar[1];

	double I4bC = ar[0]*(C_C0*ar[0]);
	if(I4bC < 1.0) {
		I4bC = 1.0;
	}
	double E0 = ((kappa*(I1bC0) + (1-3*kappa)*(I4bC)))-1.0;
	mat3ds s;
	s.zero();
	mat3ds h0;
	h0 = kappa*b_C0 + (1-3*kappa)*dyad(a[0]);
	s += h0*(2.*k1*E0*exp(k2*E0*E0));

	double I6bC = ar[1]*(C_C1*ar[1]);
	if(I6bC < 1.0) {
		I6bC = 1.0;
	}
	double E1 = ((kappa*(I1bC1) + (1-3*kappa)*(I6bC)))-1.0;
	mat3ds h1;
	h1 = kappa*b_C1 + (1-3*kappa)*dyad(a[1]);
	s += h1*(2.*k1*E1*exp(k2*E1*E1));

	return SM_phi*s.dev()/J;
}

//-----------------------------------------------------------------------------
//! Calculates the SM deviatoric tangent (Holzapfel-Gasser-Ogden material)
tens4ds FEUncoupledConstrainedMixture::SmoothMuscleDevTangent(FEMaterialPoint& mp)
{
 FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
 double time = GetFEModel()->GetTime().currentTime;
 double SM_phi = m_SM_phi(mp); //Copy the mass fraction.
 double k1 = m_SM_k1(mp); //Copy the material properties.
 double k2 = m_SM_k2(mp);
 double kappa = m_SM_kappa(mp);
 double g = m_SM_gdeg(mp)*PI/180;

 double F_r = m_F_r(mp);
 double F_theta = m_F_theta(mp);
 double F_z = m_F_z(mp);

 double SM_G = m_SM_G(mp); //Copy the deposition stretch value.
 double G_Time = m_G_Time; //Define the overall time where an "incremental" deposition stretch is applied

 mat3dd I(1); // Identity

 //Current deposition stretch tensor oriented respect to an auxiliar system of coordinates.
 //One of the axes of this system of coordinates is coincident with the fiber direction
 mat3dd SM_G_fib = (time <= G_Time) ? (mat3dd(1/sqrt(SM_G),SM_G,1/sqrt(SM_G))-I)*time/G_Time+I : mat3dd(1/sqrt(SM_G),SM_G,1/sqrt(SM_G));
 double cg = cos(g); double sg = sin(g);

 //Rotation matrix from the auxiliar system of coordinates to cylindrical coordinate system (r,theta,zeta)
 //Family 0
 mat3d T0;
 T0.zero();
 T0[0][0]=1.0;   T0[0][1]=0.0; T0[0][2]=0.0;
 T0[1][0]=0.0;   T0[1][1]=sg;  T0[1][2]=-cg;
 T0[2][0]=0.0;   T0[2][1]=cg;  T0[2][2]=sg;

 //Family 1
 mat3d T1;
 T1.zero();
 T1[0][0]=1.0; T1[0][1]=0.0; T1[0][2]=0.0;
 T1[1][0]=0.0; T1[1][1]=-sg; T1[1][2]=-cg;
 T1[2][0]=0.0; T1[2][1]=cg;  T1[2][2]=-sg;

 //Current deposition stretch tensor oriented respect to cylindrical coordinates (r,theta,zeta).
 mat3ds SM_G_fib0 = (T0*SM_G_fib*T0.transpose()).sym();
 mat3ds SM_G_fib1 = (T1*SM_G_fib*T1.transpose()).sym();

 mat3dd F_ro = (time <= G_Time) ? (mat3dd(F_r,F_theta,F_z)-I)*time/G_Time+I : mat3dd(F_r,F_theta,F_z);
 SM_G_fib0 = SM_G_fib0*F_ro;
 SM_G_fib1 = SM_G_fib1*F_ro;

 mat3d Q = GetLocalCS(mp); //Get material axes.

 //Current deposition stretch tensor oriented respect to cartesian coordinates (x,y,z).
 mat3ds SM_G_rot0 = (Q*SM_G_fib0*Q.transpose()).sym(); //Rotate the pre-stretch tensor respect to the material axes.
 mat3ds SM_G_rot1 = (Q*SM_G_fib1*Q.transpose()).sym(); //Rotate the pre-stretch tensor respect to the material axes.


 // determinant of deformation gradient
 double J = pt.m_J;

 // Evaluate the distortional deformation gradient
 double Jm13 = pow(J, -1. / 3.);
 mat3d F = pt.m_F*Jm13;

 // Collagen deformation gradient
 mat3d F_C0 = F*SM_G_rot0;
 mat3d F_C1 = F*SM_G_rot1;

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
 vec3d ar[2],a[2];
 ar[0] = n[0]*cg + n[1]*sg; a[0] = F_C0*ar[0];
 ar[1] = n[0]*cg - n[1]*sg; a[1] = F_C1*ar[1];

 double I4bC = ar[0]*(C_C0*ar[0]);
 if(I4bC < 1.0) {
	 I4bC = 1.0;
 }
 double E0 = ((kappa*(I1bC0) + (1-3*kappa)*(I4bC)))-1.0;
 mat3ds s;
 s.zero();
 mat3ds h0;
 h0 = kappa*b_C0 + (1-3*kappa)*dyad(a[0]);
 s += h0*(2.*k1*E0*exp(k2*E0*E0));

 double I6bC = ar[1]*(C_C1*ar[1]);
 if(I6bC < 1.0) {
	 I6bC = 1.0;
 }
 double E1 = ((kappa*(I1bC1) + (1-3*kappa)*(I6bC)))-1.0;
 mat3ds h1;
 h1 = kappa*b_C1 + (1-3*kappa)*dyad(a[1]);
 s += h1*(2.*k1*E1*exp(k2*E1*E1));
 

 mat3ds sbar = s.dev();

 // Evaluate the elasticity tensor
 tens4ds IxI = dyad1s(I);
 tens4ds I4  = dyad4s(I);
 tens4ds ce  = ((I4 - IxI/3.)*s.tr()-dyad1s(sbar,I))*(2.0/3.0);
     if (E0 >= 0) ce += dyad1s(h0.dev())*(4.*k1*(1 + 2 * k2*E0*E0)*exp(k2*E0*E0));
 		 if (E1 >= 0) ce += dyad1s(h1.dev())*(4.*k1*(1 + 2 * k2*E1*E1)*exp(k2*E1*E1));

 return SM_phi*ce / J;
}


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
