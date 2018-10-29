// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "elements/ElementFactory.h"
#include "elements/PLANE41.h"
#include "elements/SOLID81.h"
#include "elements/TRUSS3.h"
#include "elements/TETRA0.h"
#include "elements/TETRA1.h"
#include "elements/QUADTH.h"
#include "elements/TRIANGLE4.h"
#include "elements/TRIANGLE_THERMO.h"
#include "elements/INTER0.h"
#include "elements/INTER3.h"

namespace nla3d {


ElementType ElementFactory::elName2elType (std::string elName) {
  for (uint16 i = 0; i < (uint16) ElementType::UNDEFINED; i++) {
    if (elName.compare(elTypeLabels[i]) == 0) {
      return (ElementType) i;
    }
  }
  return ElementType::UNDEFINED;
}


void ElementFactory::createElements (ElementType elId, const uint32 n, std::vector<Element*>& ptr) {
  if (elId == ElementType::UNDEFINED) {
    LOG(FATAL) << "Element type is undefined";
  }
  if (n == 0) {
    return;
  }
  // reserve memory in the vector for newly created elements
  ptr.reserve(ptr.size() + n);

  switch (elId) {
      case ElementType::PLANE41:
        for (uint32 i = 0; i < n; i++) {
          ptr.push_back(new ElementPLANE41());
        }
        break;
      case ElementType::SOLID81:
        for (uint32 i = 0; i < n; i++) {
          ptr.push_back(new ElementSOLID81());
        }
        break;
      case ElementType::TRUSS3:
        for (uint32 i = 0; i < n; i++) {
          ptr.push_back(new ElementTRUSS3());
        }
        break;
      case ElementType::TRIANGLE4:
        for (uint32 i = 0; i < n; i++) {
          ptr.push_back(new ElementTRIANGLE4());
        }
        break;
      case ElementType::TETRA0:
        for (uint32 i = 0; i < n; i++) {
          ptr.push_back(new ElementTETRA0());
        }
        break;
      case ElementType::TETRA1:
        for (uint32 i = 0; i < n; i++) {
          ptr.push_back(new ElementTETRA1());
        }
        break;
      case ElementType::QUADTH:
        for (uint32 i = 0; i < n; i++) {
          ptr.push_back(new ElementQUADTH());
        }
        break;
      case ElementType::TRIANGLE_THERMO:
        for (uint32 i = 0; i < n; i++) {
          ptr.push_back(new ElementTRIANGLE_THERMO());
        }
        break;
      case ElementType::SurfaceLINETH:
        for (uint32 i = 0; i < n; i++) {
          ptr.push_back(new SurfaceLINETH());
        }
        break;
      case ElementType::INTER0:
        for (uint32 i = 0; i < n; i++) {
          ptr.push_back(new ElementINTER0());
        }
        break;
      case ElementType::INTER3:
        for (uint32 i = 0; i < n; i++) {
          ptr.push_back(new ElementINTER3());
        }
        break;
      default:
        LOG(ERROR) << "Don't have an element with id " << (uint16) elId;
    }
}

} // namespace nla3d
