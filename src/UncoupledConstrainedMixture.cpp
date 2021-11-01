#include "UncoupledConstrainedMixture.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FEUncoupledConstrainedMixture, FEUncoupledMaterial)
	ADD_PARAMETER(m_c    , FE_RANGE_GREATER_OR_EQUAL(0.0), "c");
	ADD_PARAMETER(m_E_G[0]    , FE_RANGE_GREATER_OR_EQUAL(1.0), "E_G_z"); // Z pre-stretch
	ADD_PARAMETER(m_E_G[1]    , FE_RANGE_GREATER_OR_EQUAL(1.0), "E_G_theta");// theta pre-stretch. radial pre-stetch is calculated by incompressibility of E_G
	ADD_PARAMETER(m_G_Time    , FE_RANGE_GREATER_OR_EQUAL(0.0), "G_Time");

END_FECORE_CLASS();


//-----------------------------------------------------------------------------
//! Calculates the deviatoric stress
mat3ds FEUncoupledConstrainedMixture::DevStress(FEMaterialPoint& mp)
{

	return NeoHookeDevStress(mp);
}

//-----------------------------------------------------------------------------
//! Calculates the deviatoric tangent
tens4ds FEUncoupledConstrainedMixture::DevTangent(FEMaterialPoint& mp)
{
	  return NeoHookeDevTangent(mp);
}

mat3ds FEUncoupledConstrainedMixture::NeoHookeDevStress(FEMaterialPoint& mp)
{
	mat3dd I(1);
	double time = GetFEModel()->GetTime().currentTime;
	double c = m_c(mp);
	double E_G_z = m_E_G[0](mp); //Copy the material properties.
	double E_G_t = m_E_G[1](mp);
	double E_G_r = 1.0/(E_G_z*E_G_t);
	double G_Time = m_G_Time;
	//Form the stretch deposition tensor E_G with a diagonal matrix.
	//| E_G_r				0						0 		|
	//|   0			E_G_theta				0 		|
	//|		0					0						E_G_z |
	mat3dd E_G = (time <= G_Time) ? (mat3dd(E_G_r,E_G_t,E_G_z)-I)*time/G_Time+I : mat3dd(E_G_r,E_G_t,E_G_z);

	mat3d Q = GetLocalCS(mp); //Get material axes.
	mat3d E_G_rot = (Q*E_G*Q.transpose()); //Rotate the pre-stretch tensor respect to the material axes.

	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// determinant of deformation gradient
	double J = pt.m_J;


	// Evaluate the distortional deformation gradient
	double Jm13 = pow(J, -1. / 3.);
	mat3d F = pt.m_F*Jm13;

	// calculate deviatoric left Cauchy-Green tensor: b = F*Ft
	mat3ds b = pt.DevLeftCauchyGreen();

	// Evaluate the ground matrix stress
	//mat3ds s = E_G*c*I*E_G*b;
	mat3ds s = (E_G_rot*c*I*E_G_rot*b).sym();


	return s.dev() / J;
}


tens4ds FEUncoupledConstrainedMixture::NeoHookeDevTangent(FEMaterialPoint& mp)
{
	  mat3dd I(1);
		double time = GetFEModel()->GetTime().currentTime;
    double c = m_c(mp);
		double E_G_z = m_E_G[0](mp);
		double E_G_t = m_E_G[1](mp);
		double E_G_r = 1.0/(E_G_z*E_G_t);
		double G_Time = m_G_Time;

		mat3dd E_G = (time <= G_Time) ? (mat3dd(E_G_r,E_G_t,E_G_z)-I)*time/G_Time+I : mat3dd(E_G_r,E_G_t,E_G_z);
		mat3d Q = GetLocalCS(mp);
		mat3d E_G_rot = (Q*E_G*Q.transpose());

    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

    // determinant of deformation gradient
    double J = pt.m_J;

    // Evaluate the distortional deformation gradient
	  double Jm13 = pow(J, -1. / 3.);
    mat3d F = pt.m_F*Jm13;

    // calculate deviatoric left Cauchy-Green tensor: b = F*Ft
    mat3ds b = pt.DevLeftCauchyGreen();

    // Evaluate the ground matrix stress
    mat3ds s = (E_G_rot*c*I*E_G_rot*b).sym();

		mat3ds sbar = s.dev();

    // Evaluate the elasticity tensor
    tens4ds IxI = dyad1s(I);
    tens4ds I4  = dyad4s(I);
    tens4ds ce = ((I4 - IxI/3.)*s.tr()-dyad1s(sbar,I))*(2./3.);
		ce = ce.pp(E_G_rot);
    return ce / J;
}
