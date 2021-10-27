#include "UncoupledConstrainedMixture.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
//-----------------------------------------------------------------------------
// These macros define the material parameter list for the FENeoHookeanPI material.
//
// The BEGIN_FECORE_CLASS macro takes the material class and its base class as
// parameters.
//
// The ADD_PARAMETER macro defines the actual parameter. It takes three parameters:
// - the variable as defined in the class.
// - a range speficier which defines the valid range of the parameter
// - a string that defines the name of the parameter as it will appear in the input file.
//
// The END_PARAMETER_LIST macro just defines the end of the parameter list.
// define the material parameters
BEGIN_FECORE_CLASS(FEUncoupledConstrainedMixture, FEUncoupledMaterial)
	ADD_PARAMETER(m_c    , FE_RANGE_GREATER_OR_EQUAL(0.0), "c");
	ADD_PARAMETER(m_E_G[0]    , FE_RANGE_GREATER_OR_EQUAL(1.0), "E_G_x");
	ADD_PARAMETER(m_E_G[1]    , FE_RANGE_GREATER_OR_EQUAL(1.0), "E_G_y");
	ADD_PARAMETER(m_E_G_Time    , FE_RANGE_GREATER_OR_EQUAL(0.0), "E_G_Time");

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
	double E_G_x = m_E_G[0](mp);
	double E_G_y = m_E_G[1](mp);
	double E_G_z = 1.0/(E_G_x*E_G_y);
	double E_G_Time = m_E_G_Time;
	mat3dd E_G = (time <= E_G_Time) ? (mat3dd(E_G_x,E_G_y,E_G_z)-I)*time/E_G_Time+I : mat3dd(E_G_x,E_G_y,E_G_z);

	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// determinant of deformation gradient
	double J = pt.m_J;


	// Evaluate the distortional deformation gradient
	double Jm13 = pow(J, -1. / 3.);
	mat3d F = pt.m_F*Jm13;

	// calculate deviatoric left Cauchy-Green tensor: b = F*Ft
	mat3ds b = pt.DevLeftCauchyGreen();
	mat3ds C = pt.DevRightCauchyGreen();
	double I1 = C.tr();

	// Evaluate the ground matrix stress
	mat3ds s = E_G*c*I*E_G*b;


	return s.dev() / J;
}


tens4ds FEUncoupledConstrainedMixture::NeoHookeDevTangent(FEMaterialPoint& mp)
{
	  mat3dd I(1);
		double time = GetFEModel()->GetTime().currentTime;
    double c = m_c(mp);
		double E_G_x = m_E_G[0](mp);
		double E_G_y = m_E_G[1](mp);
		double E_G_z = 1.0/(E_G_x*E_G_y);
		double E_G_Time = m_E_G_Time;

		mat3dd E_G = (time <= E_G_Time) ? (mat3dd(E_G_x,E_G_y,E_G_z)-I)*time/E_G_Time+I : mat3dd(E_G_x,E_G_y,E_G_z);

    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

    // determinant of deformation gradient
    double J = pt.m_J;

    // Evaluate the distortional deformation gradient
	  double Jm13 = pow(J, -1. / 3.);
    mat3d F = pt.m_F*Jm13;

    // calculate deviatoric left Cauchy-Green tensor: b = F*Ft
    mat3ds b = pt.DevLeftCauchyGreen();
    mat3ds C = pt.DevRightCauchyGreen();
    double I1 = C.tr();


    // Evaluate the ground matrix stress
    mat3ds s = E_G*c*I*E_G*b;


	mat3ds sbar = s.dev();

    // Evaluate the elasticity tensor

    tens4ds IxI = dyad1s(I);
    tens4ds I4  = dyad4s(I);
    tens4ds ce = ((I4 - IxI/3.)*s.tr()-dyad1s(sbar,I))*(2./3.);
		ce = ce.pp(E_G);
    return ce / J;
}
