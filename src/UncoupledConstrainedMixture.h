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

	DECLARE_FECORE_CLASS();

private:
	FEParamDouble  m_c;         // neo-Hookean c coefficient
	FEParamDouble  m_E_G[2];
	double	m_G_Time;	// Elastin desposition time
	mat3ds NeoHookeDevStress(FEMaterialPoint& mp);
	tens4ds NeoHookeDevTangent(FEMaterialPoint& mp);



};
