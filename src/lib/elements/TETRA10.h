// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "elements/element.h"
#include "elements/isoparametric.h"
#include "FEStorage.h"
#include "solidmech.h"

namespace nla3d {
//4-node 2D QUAD  element for steady thermal analysis
class ElementTETRA10 : public ElementIsoParamTETRA10 {
  public:
  ElementTETRA10 () {
    intOrder = 2;
    type = ElementType::TETRA10;
  }
  void pre();
  void buildK();
  void update();
  void makeB (uint16 nPoint, math::Mat<6,30> &B);
  void makeC (math::MatSym<6> &C);
  // Elastic module
  double E = 0.0;
  // Poissons coef.
  double my = 0.0;

  // stresses in the element (calculated after the solving of the global equation system in
  // update() function.
  //stress[M_XX], stress[M_YY], stress[M_ZZ], stress[M_XY], stress[M_YZ], stress[M_XZ]
  math::Vec<6> stress; // Cauchy stresses

  //strains[M_XX], strains[M_YY], strains[M_ZZ], strains[M_XY], strains[M_YZ], strains[M_XZ]
  math::Vec<6> strains;

  double vol;

  //postproc procedures
  bool getScalar(double* scalar, scalarQuery code, uint16 gp, const double scale);
  bool getTensor(math::MatSym<3>* tensor, tensorQuery code, uint16 gp, const double scale);

};
} // nla3d namespace
