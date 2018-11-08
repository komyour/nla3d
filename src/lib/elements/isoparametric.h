// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"

// Formulations for isoparametric FE for different shapes. Isoparametic formulation means usage of the
// same shape functions for geometry interpolation and for field (displacement, for instance)
//   interpolation.
// These derived classes specify next things:
// 1. Specify form functions
// 2. Specify integration rule (Gauss quadrature)
// 3. Specify field derivatives in local and global coordinates

namespace nla3d {

// structures for convenient keeping of quadrature constants
struct QuadPt1D {
  double r;
  double w;
};


struct QuadPt2D {
  double r;
  double s;
  double w;
};


struct QuadPt3D {
  double r;
  double s;
  double t;
  double w;
};

struct QuadVol2D {
  double l1;
  double l2;
  double l3;
  double w;
};

struct QuadVol3D {
  double l1;
  double l2;
  double l3;
  double l4;
  double w;
};


// gauss quadrature for 1D line from (-1) to (1)
// 1st order
static const QuadPt1D _quad_line_o1[] = {
  {+0.0000000000000000, +2.0000000000000000}
};

// 2nd order
static const QuadPt1D _quad_line_o2[] = {
  {-0.5773502691896258, +1.0000000000000000},
  {+0.5773502691896258, +1.0000000000000000}
};

// 3rd order
static const QuadPt1D _quad_line_o3[] = {
  {-0.7745966692414834, +0.5555555555555556},
  {+0.0000000000000000, +0.8888888888888888},
  {+0.7745966692414834, +0.5555555555555556}
};

// array of number of quadrature points in integration scheme
static const uint16 _np_line[] = {
  sizeof(_quad_line_o1) / sizeof(QuadPt1D),
  sizeof(_quad_line_o2) / sizeof(QuadPt1D),
  sizeof(_quad_line_o3) / sizeof(QuadPt1D)
};

static const QuadPt1D* _table_line[] = {
  _quad_line_o1,
  _quad_line_o2,
  _quad_line_o3
};


// gauss quadrature for 2D rect from (-1,-1) to (1,1)
// 1st order
static const QuadPt2D _quad_quad_o1[] = {
  {+0.0000000000000000, +0.0000000000000000, +4.0000000000000000}
};

// 2nd order
static const QuadPt2D _quad_quad_o2[] = {
  {-0.5773502691896258, -0.5773502691896258, +1.0000000000000000},
  {+0.5773502691896258, -0.5773502691896258, +1.0000000000000000},
  {-0.5773502691896258, +0.5773502691896258, +1.0000000000000000},
  {+0.5773502691896258, +0.5773502691896258, +1.0000000000000000}
};

// 3rd order
static const QuadPt2D _quad_quad_o3[] = {
  {-0.7745966692414834, -0.7745966692414834, +0.3086419753086420},
  {+0.0000000000000000, -0.7745966692414834, +0.4938271604938271},
  {+0.7745966692414834, -0.7745966692414834, +0.3086419753086420},

  {-0.7745966692414834, +0.0000000000000000, +0.4938271604938271},
  {+0.0000000000000000, +0.0000000000000000, +0.7901234567901234},
  {+0.7745966692414834, +0.0000000000000000, +0.4938271604938271},

  {-0.7745966692414834, +0.7745966692414834, +0.3086419753086420},
  {+0.0000000000000000, +0.7745966692414834, +0.4938271604938271},
  {+0.7745966692414834, +0.7745966692414834, +0.3086419753086420}
};


// array of number of quadrature points in integration scheme
static const uint16 _np_quad[] = {
  sizeof(_quad_quad_o1) / sizeof(QuadPt2D),
  sizeof(_quad_quad_o2) / sizeof(QuadPt2D),
  sizeof(_quad_quad_o3) / sizeof(QuadPt2D)
};


static const QuadPt2D* _table_quad[] = {
  _quad_quad_o1,
  _quad_quad_o2,
  _quad_quad_o3
};


// gauss quadrature for 3D brick from (-1,-1,-1) to (1,1,1)
// 1st order
static const QuadPt3D _quad_hexahedron_o1[] = {
  {+0.0000000000000000, +0.0000000000000000, +0.0000000000000000, +8.0000000000000000}
};


// 2nd order
static const QuadPt3D _quad_hexahedron_o2[] = {
  {-0.5773502691896258, -0.5773502691896258, -0.5773502691896258, +1.0000000000000000},
  {+0.5773502691896258, -0.5773502691896258, -0.5773502691896258, +1.0000000000000000},
  {-0.5773502691896258, +0.5773502691896258, -0.5773502691896258, +1.0000000000000000},
  {+0.5773502691896258, +0.5773502691896258, -0.5773502691896258, +1.0000000000000000},

  {-0.5773502691896258, -0.5773502691896258, +0.5773502691896258, +1.0000000000000000},
  {+0.5773502691896258, -0.5773502691896258, +0.5773502691896258, +1.0000000000000000},
  {-0.5773502691896258, +0.5773502691896258, +0.5773502691896258, +1.0000000000000000},
  {+0.5773502691896258, +0.5773502691896258, +0.5773502691896258, +1.0000000000000000}
};


// 3rd order
static const QuadPt3D _quad_hexahedron_o3[] = {
  {-0.7745966692414834, -0.7745966692414834, -0.7745966692414834, +0.1714677640603567},
  {+0.0000000000000000, -0.7745966692414834, -0.7745966692414834, +0.2743484224965707},
  {+0.7745966692414834, -0.7745966692414834, -0.7745966692414834, +0.1714677640603567},

  {-0.7745966692414834, +0.0000000000000000, -0.7745966692414834, +0.2743484224965707},
  {+0.0000000000000000, +0.0000000000000000, -0.7745966692414834, +0.4389574759945130},
  {+0.7745966692414834, +0.0000000000000000, -0.7745966692414834, +0.2743484224965707},

  {-0.7745966692414834, +0.7745966692414834, -0.7745966692414834, +0.1714677640603567},
  {+0.0000000000000000, +0.7745966692414834, -0.7745966692414834, +0.2743484224965707},
  {+0.7745966692414834, +0.7745966692414834, -0.7745966692414834, +0.1714677640603567},


  {-0.7745966692414834, -0.7745966692414834, +0.0000000000000000, +0.2743484224965707},
  {+0.0000000000000000, -0.7745966692414834, +0.0000000000000000, +0.4389574759945130},
  {+0.7745966692414834, -0.7745966692414834, +0.0000000000000000, +0.2743484224965707},

  {-0.7745966692414834, +0.0000000000000000, +0.0000000000000000, +0.4389574759945130},
  {+0.0000000000000000, +0.0000000000000000, +0.0000000000000000, +0.7023319615912207},
  {+0.7745966692414834, +0.0000000000000000, +0.0000000000000000, +0.4389574759945130},

  {-0.7745966692414834, +0.7745966692414834, +0.0000000000000000, +0.2743484224965707},
  {+0.0000000000000000, +0.7745966692414834, +0.0000000000000000, +0.4389574759945130},
  {+0.7745966692414834, +0.7745966692414834, +0.0000000000000000, +0.2743484224965707},


  {-0.7745966692414834, -0.7745966692414834, +0.7745966692414834, +0.1714677640603567},
  {+0.0000000000000000, -0.7745966692414834, +0.7745966692414834, +0.2743484224965707},
  {+0.7745966692414834, -0.7745966692414834, +0.7745966692414834, +0.1714677640603567},

  {-0.7745966692414834, +0.0000000000000000, +0.7745966692414834, +0.2743484224965707},
  {+0.0000000000000000, +0.0000000000000000, +0.7745966692414834, +0.4389574759945130},
  {+0.7745966692414834, +0.0000000000000000, +0.7745966692414834, +0.2743484224965707},

  {-0.7745966692414834, +0.7745966692414834, +0.7745966692414834, +0.1714677640603567},
  {+0.0000000000000000, +0.7745966692414834, +0.7745966692414834, +0.2743484224965707},
  {+0.7745966692414834, +0.7745966692414834, +0.7745966692414834, +0.1714677640603567}
};

static const uint16 _np_hexahedron[] = {
  sizeof(_quad_hexahedron_o1) / sizeof(QuadPt3D),
  sizeof(_quad_hexahedron_o2) / sizeof(QuadPt3D),
  sizeof(_quad_hexahedron_o3) / sizeof(QuadPt3D)
};


static const QuadPt3D* _table_hexahedron[] = {
  _quad_hexahedron_o1,
  _quad_hexahedron_o2,
  _quad_hexahedron_o3
};


// gauss quadrature for tetra from (0,0,0,0) to (1,1,1,1)
// 1st order
static const QuadVol3D _quad_tetra_o1[] = {
  {0.25, 0.25, 0.25, 0.25, 1.00}
};
// 2nd order
static const QuadVol3D _quad_tetra_o2[] = {
  {0.5854102, 0.1381966, 0.1381966, 0.1381966, 0.25},
  {0.1381966, 0.5854102, 0.1381966, 0.1381966, 0.25},
  {0.1381966, 0.1381966, 0.5854102, 0.1381966, 0.25},
  {0.1381966, 0.1381966, 0.1381966, 0.5854102, 0.25}
};

// 3d order
static const QuadVol3D _quad_tetra_o3[] = {
  {0.25, 0.25, 0.25, 0.25, -0.8},
  {0.33, 0.66, 0.66, 0.66, 0.45},
  {0.66, 0.33, 0.66, 0.66, 0.45},
  {0.66, 0.66, 0.33, 0.66, 0.45},
  {0.66, 0.66, 0.66, 0.33, 0.45}
};

static const uint16 _np_tetra[] = {
  sizeof(_quad_tetra_o1) / sizeof(QuadVol3D),
  sizeof(_quad_tetra_o2) / sizeof(QuadVol3D),
  sizeof(_quad_tetra_o3) / sizeof(QuadVol3D)
};


static const QuadVol3D* _table_tetra[] = {
  _quad_tetra_o1,
  _quad_tetra_o2,
  _quad_tetra_o3
};

class ElementIsoParamLINE : public ElementLINE {
  public:
    double det; //Jacobian

    void makeJacob();

    double intWeight(uint16 np);
    double volume();
    void np2rst(uint16 np, double *xi); //by number of gauss point find local coordinates
    uint16 nOfIntPoints();

    // get form function values in local point (r, s)
    math::Vec<2> formFunc(double r);
    // get form function values in integration point np
    math::Vec<2> formFunc(uint16 np);
    math::Mat<2, 1> formFuncDeriv(double r);

  protected:
    uint16 i_int = 0; // index of integration scheme
};

class ElementIsoParamQUAD : public ElementQUAD {
  public:
    std::vector<math::Mat<4, 2> > NiXj; //derivates form function / local coordinates
    std::vector<double> det;  //Jacobian

    double sideDet[4]; //Jacobian for side integration

    // function to calculate all staff for isoparametric FE
    void makeJacob(); 
    double intWeight(uint16 np);
    double volume();
    void np2rst(uint16 np, double *xi); //by number of gauss point find local coordinates
    uint16 nOfIntPoints();

    // get form function values in local point (r, s)
    math::Vec<4> formFunc(double r, double s);
    // get form function values in integration point np
    math::Vec<4> formFunc(uint16 np);
    math::Mat<4, 2> formFuncDeriv(double r, double s);

  protected:
    uint16 i_int = 0; // index of integration scheme
};

class ElementIsoParamTETRA10 : public ElementQUADTETRA {
  public:
    std::vector<math::Mat<10, 4> > NiXj; //derivates form function / local coordinates
    std::vector<double> det;  //Jacobian

    // function to calculate all staff for isoparametric FE
    void makeJacob(); 
    double intWeight(uint16 np);
    double volume();
    void np2rst(uint16 np, double *xi); //by number of gauss point find local coordinates
    uint16 nOfIntPoints();
    // get form function values in local point (l1, l2, l3, l4)
    math::Vec<10> formFunc(double l1, double l2, double l3, double l4);
    // get form function values in integration point np
    math::Vec<10> formFunc(uint16 np);
    // get form function derivatives (vs. r,s,t) in local point (l1,l2,l3,l4)
    math::Mat<10, 4> formFuncDeriv(double l1, double l2, double l3, double l4);

  protected:
    uint16 i_int = 0; // index of integration scheme
};

class ElementIsoParamHEXAHEDRON : public ElementHEXAHEDRON {
  public:
    std::vector<math::Mat<8, 3> > NiXj; //derivates form function / local coordinates
    std::vector<double> det;  //Jacobian

    // function to calculate all staff for isoparametric FE
    void makeJacob(); 
    double intWeight(uint16 np);
    double volume();
    void np2rst(uint16 np, double *xi); //by number of gauss point find local coordinates
    uint16 nOfIntPoints();
    // get form function values in local point (r, s, t)
    math::Vec<8> formFunc(double r, double s, double t);
    // get form function values in integration point np
    math::Vec<8> formFunc(uint16 np);
    // get form function derivatives (vs. r,s,t) in local point (r,s,t)
    math::Mat<8, 3> formFuncDeriv(double r, double s, double t);

  protected:
    uint16 i_int = 0; // index of integration scheme
};


inline double ElementIsoParamLINE::intWeight(uint16 np) {
  assert(np < _np_line[i_int]);
  assert(det > 0.0);

  return _table_line[i_int][np].w * det;
}


inline void ElementIsoParamLINE::np2rst (uint16 np, double *v) {
  assert(np < _np_line[i_int]);
  v[0] = _table_line[i_int][np].r;
}


inline uint16 ElementIsoParamLINE::nOfIntPoints() {
  return _np_line[i_int];
}


inline double ElementIsoParamQUAD::intWeight(uint16 np) {
  assert(np < _np_quad[i_int]);
  assert(det.size() != 0);

  return _table_quad[i_int][np].w * det[np];
}


inline void ElementIsoParamQUAD::np2rst (uint16 np, double *v) {
  assert(np < _np_quad[i_int]);
  v[0] = _table_quad[i_int][np].r;
  v[1] = _table_quad[i_int][np].s;
}


inline uint16 ElementIsoParamQUAD::nOfIntPoints() {
  return _np_quad[i_int];
}

inline double ElementIsoParamHEXAHEDRON::intWeight(uint16 np) {
  assert(np < _np_hexahedron[i_int]);
  assert(det.size() != 0);

  return _table_hexahedron[i_int][np].w * det[np];
}

inline void ElementIsoParamHEXAHEDRON::np2rst (uint16 np, double *v) {
  assert(np < _np_hexahedron[i_int]);
  v[0] = _table_hexahedron[i_int][np].r;
  v[1] = _table_hexahedron[i_int][np].s;
  v[2] = _table_hexahedron[i_int][np].t;
}


inline uint16 ElementIsoParamHEXAHEDRON::nOfIntPoints() {
  return _np_hexahedron[i_int];
}

inline double ElementIsoParamTETRA10::intWeight(uint16 np) {
  assert(np < _np_tetra[i_int]);
  assert(det.size() != 0);

  return _table_tetra[i_int][np].w * det[np];
}


inline void ElementIsoParamTETRA10::np2rst (uint16 np, double *v) {
  assert(np < _np_tetra[i_int]);
  v[0] = _table_tetra[i_int][np].l1;
  v[1] = _table_tetra[i_int][np].l2;
  v[2] = _table_tetra[i_int][np].l3;
  v[4] = _table_tetra[i_int][np].l4;
}


inline uint16 ElementIsoParamTETRA10::nOfIntPoints() {
  return _np_tetra[i_int];
}

} // namespace nla3d
