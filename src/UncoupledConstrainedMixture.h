#pragma once
//=============================================================================
// This plugin example illustrates how to create a new material.
// It requires FEBio 3.0 (or up)
//
// Author : Steve Maas
// Copyright (c) 2015 - 2020
// All rights reserved
//
//=============================================================================

//-----------------------------------------------------------------------------
// We need to include this file since our new material class will inherit from
// FEElasticMaterial, which is defined in this include file.
#include <FEBioMech/FEUncoupledMaterial.h>
#include <FECore/FEModelParam.h>
//-----------------------------------------------------------------------------
// This material class implements a neo-Hookean constitutive model.
// Since it is a (compressible, coupled) hyper-elatic material, it needs to inherit
// from FEElasticMaterial.
class FEUncoupledConstrainedMixture : public FEUncoupledMaterial
{
public:
	FEUncoupledConstrainedMixture(FEModel* pfem) : FEUncoupledMaterial(pfem) { m_npmodel = 3; }

	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;


	// declare parameter list
	DECLARE_FECORE_CLASS();

private:
	// The neo-Hookean material defines two material parameters.
	// They are defined here, but the class also needs to let the framework
	// know that these parameters exist. This is necessary to define the parameters
	// in the FEBio input file. This is done by declaring the DECLARE_FECORE_CLASS() below.
	FEParamDouble  m_c;         // neo-Hookean c coefficient
	FEParamDouble  m_E_G[2];
	double	m_E_G_Time;	// Elastin desposition time
	mat3ds NeoHookeDevStress(FEMaterialPoint& mp);
	tens4ds NeoHookeDevTangent(FEMaterialPoint& mp);



};
