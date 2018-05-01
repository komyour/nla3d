#include "elements/TRIANGLE_THERMO.h"

namespace nla3d {

ElementTRIANGLE_THERMO::ElementTRIANGLE_THERMO () {
  type = ElementType::TRIANGLE_THERMO;
}

void ElementTRIANGLE_THERMO::pre () {
  for (uint16 i = 0; i < getNNodes(); i++) {
    storage->addNodeDof(getNodeNumber(i), {Dof::TEMP});
  }
}

void ElementTRIANGLE_THERMO::buildK() {
  math::Mat<3,3> matS(1. , storage->getNode(getNodeNumber(0)).pos[0] , storage->getNode(getNodeNumber(0)).pos[1], 
                      1. , storage->getNode(getNodeNumber(1)).pos[0] , storage->getNode(getNodeNumber(1)).pos[1], 
                      1. , storage->getNode(getNodeNumber(2)).pos[0] , storage->getNode(getNodeNumber(2)).pos[1]) ;
  square = matS.det() / 2.;

  math::MatSym<3> matKe;
  matKe.zero();

  math::Mat<2,3> matB;
  matB.zero();

  math::MatSym<2> matC;
  matC.zero();

  makeC(matC);
  makeB(matB);
  math::matBTDBprod(matB, matC, square, matKe);
  assembleK(matKe, {Dof::TEMP});
}

void ElementTRIANGLE_THERMO::update () {
  math::Mat<2,3> matB;
  matB.zero();

  math::MatSym<2> matC;
  matC.zero();

  makeC(matC);
  makeB(matB);

  math::Vec<3> U;
  for (uint16 i = 0; i < getNNodes(); i++) {
    U[i] = storage->getNodeDofSolution(getNodeNumber(i), Dof::TEMP);
  }

  flux.zero();
  math::matBVprod(matB, U, k, flux);
}

void ElementTRIANGLE_THERMO::makeB(math::Mat<2,3> &B)
{
  double *B_L = B.ptr();
  double b[3], c[3];

  b[0] = storage->getNode(getNodeNumber(1)).pos[1] - storage->getNode(getNodeNumber(2)).pos[1];
  b[1] = storage->getNode(getNodeNumber(2)).pos[1] - storage->getNode(getNodeNumber(0)).pos[1];
  b[2] = storage->getNode(getNodeNumber(0)).pos[1] - storage->getNode(getNodeNumber(1)).pos[1];
  c[0] = storage->getNode(getNodeNumber(2)).pos[0] - storage->getNode(getNodeNumber(1)).pos[0];
  c[1] = storage->getNode(getNodeNumber(0)).pos[0] - storage->getNode(getNodeNumber(2)).pos[0];
  c[2] = storage->getNode(getNodeNumber(1)).pos[0] - storage->getNode(getNodeNumber(0)).pos[0];


  const double A = 1./2./square;
  B_L[0] = b[0]*A;
  B_L[1] = b[1]*A;
  B_L[2] = b[2]*A;
  B_L[3] = c[0]*A;
  B_L[4] = c[1]*A;
  B_L[5] = c[2]*A;
}

void ElementTRIANGLE_THERMO::makeC (math::MatSym<2> &C) {
  C.comp(0,0) = -k;
  C.comp(1,1) = -k;
}

bool ElementTRIANGLE_THERMO::getScalar(double* scalar, scalarQuery query, uint16 gp, const double scale) {
  if (query == scalarQuery::VOL){
     *scalar += square*scale;
    return true;
  }
  return false;
}

bool  ElementTRIANGLE_THERMO::getVector(math::Vec<2>* vector, vectorQuery query, uint16 gp, const double scale) {
  switch (query) {
    case vectorQuery::FLUX:
      *vector += flux*scale;
      return true;
    case vectorQuery::GRADT:
      *vector += flux*(scale/k);
      return true;
  }  
  return false;
}

} //namespace nla3d
