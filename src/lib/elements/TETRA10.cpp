// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "elements/TETRA10.h"

namespace nla3d {
using namespace math;


void ElementTETRA10::pre() {
  if (det.size()==0) {
    makeJacob();
  }

  // register element equations
  for (uint16 i = 0; i < getNNodes(); i++) {
    storage->addNodeDof(getNodeNumber(i), {Dof::UX, Dof::UY, Dof::UZ});
  }
}

void ElementTETRA10::buildK() {
  MatSym<30> Ke;  // element stiff. matrix
  Ke.zero();

  Vec<10> Fe;
  Fe.zero();

  // matC is 3d elastic  matrix
  math::MatSym<6> matC;
  matC.zero();
  Mat<6,30> matB;
  makeC(matC);
  // build Ke
  double dWt; //Gaussian quadrature weight
  for (uint16 np=0; np < nOfIntPoints(); np++) {
    dWt = intWeight(np);
    makeB(np, matB);
    matBTDBprod(matB, matC, dWt / det[np], Ke);
  }// loop over integration points
  assembleK(Ke, {Dof::UX, Dof::UY, Dof::UZ});
}

void ElementTETRA10::makeB (uint16 np, Mat<6,30> &B)
{
  double *B_NL = B.ptr();
  for (uint16 i=0; i < 10; i++) {
    B_NL[0*30+(i*3+0)] += NiXj[np][i][0] / det[np];
    B_NL[1*30+(i*3+1)] += NiXj[np][i][1] / det[np];
    B_NL[2*30+(i*3+2)] += NiXj[np][i][2] / det[np];
    B_NL[3*30+(i*3+0)] += NiXj[np][i][1] / det[np];
    B_NL[3*30+(i*3+1)] += NiXj[np][i][0] / det[np];
    B_NL[4*30+(i*3+1)] += NiXj[np][i][1] / det[np];
    B_NL[4*30+(i*3+2)] += NiXj[np][i][2] / det[np];
    B_NL[5*30+(i*3+0)] += NiXj[np][i][2] / det[np];
    B_NL[5*30+(i*3+2)] += NiXj[np][i][0] / det[np];
  }
}

// after solution it's handy to calculate stresses, strains and other stuff in elements.
void ElementTETRA10::update () {
  // matB is strain matrix
  math::Mat<6,30> matB;
  matB.zero();

  // matC is 3d elastic  matrix
  math::MatSym<6> matC;
  matC.zero();

  // fill here matC
  makeC(matC);
  for (uint16 np=0; np < nOfIntPoints(); np++) {  
    // fill here matB
    makeB(np, matB);
    // get nodal solutions from storage
    math::Vec<30> U;
    for (uint16 i = 0; i < getNNodes(); i++) {
      U[i*3 + 0] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UX);
      U[i*3 + 1] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UY);
      U[i*3 + 2] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UZ);
    }
    
    //restore strains
    strains.zero();
    math::matBVprod(matB, U, -1.0, strains);
    
    stress.zero();
    math::matBVprod(matC, strains, 1.0, stress);
  }
}

void ElementTETRA10::makeC (math::MatSym<6> &C) {
  const double A = E/((1.+my)*(1.-2.*my));

  C.comp(0,0) = (1.-my)*A;
  C.comp(0,1) = my*A;
  C.comp(0,2) = my*A;
  C.comp(1,1) = (1.-my)*A;
  C.comp(1,2) = my*A;
  C.comp(2,2) = (1.-my)*A;

  C.comp(3,3) = (1./2.-my)*A;
  C.comp(4,4) = (1./2.-my)*A;
  C.comp(5,5) = (1./2.-my)*A;
}

bool ElementTETRA10::getScalar(double* scalar, scalarQuery query, uint16 gp, const double scale) {
  return false;
}

bool  ElementTETRA10::getTensor(math::MatSym<3>* tensor, tensorQuery query, uint16 gp, const double scale) {
  if (query == tensorQuery::C){
      tensor->comp(0,0) += strains[0];
      tensor->comp(1,1) += strains[1];
      tensor->comp(2,2) += strains[2];
      tensor->comp(0,1) += strains[3];
      tensor->comp(1,2) += strains[4];
      tensor->comp(0,2) += strains[5];
      return true;
  }
  if (query == tensorQuery::E){
    tensor->comp(0,0) += stress[0];
    tensor->comp(1,1) += stress[1];
    tensor->comp(2,2) += stress[2];
    tensor->comp(0,1) += stress[3];
    tensor->comp(1,2) += stress[4];
    tensor->comp(0,2) += stress[5];
    return true;
  }
  
  return false;
}
} //namespace nla3d

