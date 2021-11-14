#pragma once

#include <FEBioMech/FEUncoupledMaterial.h>
#include <FECore/FEModelParam.h>

class FEUncoupledConstrainedMixture : public FEUncoupledMaterial
{
public:
	FEUncoupledConstrainedMixture(FEModel* pfem) : FEUncoupledMaterial(pfem) { m_npmodel = 3; }

	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;
	//
	// bool Init();
	// bool Validate();


	DECLARE_FECORE_CLASS();

private:
	// #################### Elastin ####################
	FEParamDouble  m_E_phi; //Elastin mass fraction
	FEParamDouble  m_E_c;         // neo-Hookean c coefficient
	FEParamDouble  m_E_G[2];		//Pre-stretch tensor
	// #################### Elastin ####################

	// #################### Collagen ####################
	//FEParamDouble  m_C_phi; //Collagen mass fraction
	//FEParamDouble  m_C_k1, m_C_k2, m_C_kappa, m_C_gdeg;	//Fiber constants
	//FEParamDouble  m_C_G[2]; //Pre-stretch tensor
	// #################### Collagen ####################

	// #################### Time ####################
	double	m_G_Time;	// Elastin desposition time
	// #################### Collagen ####################


	// ########################################
	mat3ds NeoHookeDevStress(FEMaterialPoint& mp);
	tens4ds NeoHookeDevTangent(FEMaterialPoint& mp);
	//mat3ds CollagenDevStress(FEMaterialPoint& mp);
	//tens4ds CollagenDevTangent(FEMaterialPoint& mp);


};
