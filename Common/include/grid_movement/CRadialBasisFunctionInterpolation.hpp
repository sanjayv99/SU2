/*!
 * \file CRadialBasisFunctionInterpolation.hpp
 * \brief Headers of the CRadialBasisFunctionInterpolation class.
 * \author F. van Steen
 * \version 8.0.1 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
#include "CVolumetricMovement.hpp"
#include "CRadialBasisFunctionNode.hpp"
#include "../../include/toolboxes/CSymmetricMatrix.hpp"
#include "../../include/toolboxes/geometry_toolbox.hpp"

#include <unordered_set>
#include "../../include/adt/CADTPointsOnlyClass.hpp" // TODO -   is also in .cpp file

/*!
 * \class CRadialBasisFunctionInterpolation
 * \brief Class for moving the volumetric numerical grid using Radial Basis Function interpolation.
 * \author F. van Steen
 */

class CRadialBasisFunctionInterpolation : public CVolumetricMovement {
protected:
  
  vector<su2double> CtrlNodeDeformation;  /*!< \brief Control Node Deformation.*/ 

  vector<su2double> InterpCoeff;          /*!< \brief Control node interpolation coefficients.*/

  unsigned long nCtrlNodesGlobal{0};      /*!< \brief Total number of control nodes.*/
  unsigned long nCtrlNodesLocal{0};      /*!< \brief Total number of local control nodes.*/
  su2activematrix CtrlCoords;             /*!< \brief Coordinates of the control nodes.*/

  su2double MaxErrorGlobal{0.0};          /*!< \brief Maximum error data reduction algorithm.*/
  vector<su2double> sensitivity_update;

  vector<string> CtrlTypeVec;                                     /*!< \brief This vector contains the control nodes at that moment */  
  vector<CRadialBasisFunctionNode*> nodes;                        // vector containing all boundary nodes

  unordered_map<string, vector<unsigned long>> node_type_indices; // map containing the indices for the different type of nodes
  unordered_set<unsigned long> control_node_indices;              // in case of DR this contains the control node indices
  unordered_map<string, unordered_set<unsigned long>> ctrl_nodes_type;  // map containing the control nodes of the different types in case of DR
  
  vector<unsigned short> PeriodicAxis{0,0,0};
  vector<su2double> PeriodicLength{0,0,0};
  su2double PeriodicAngle{0.0};
  unsigned short RotationalAxis;

  bool IsCylindrical = false; 

  vector<CRadialBasisFunctionNode*> per_nodes;

  
public:

  /*!
  * \brief Constructor of the class.
  */
  CRadialBasisFunctionInterpolation(CGeometry* geometry, CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CRadialBasisFunctionInterpolation(void) override;

  /*!
   * \brief Grid deformation using the spring analogy method.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] UpdateGeo - Update geometry.
   * \param[in] Derivative - Compute the derivative (disabled by default). Does not actually deform the grid if enabled.
   */
  void SetVolume_Deformation(CGeometry* geometry, CConfig* config, bool UpdateGeo, bool Derivative,
                                                bool ForwardProjectionDerivative);
  
  /*!
  * \brief Selecting unique set of boundary nodes based on marker information.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  */
  void SetBoundNodes(CGeometry* geometry, CConfig* config);

  /*!
  * \brief Selecting internal nodes for the volumetric deformation.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem. 
  * \param[in] internalNode - Internal nodes.
  */
  void SetInternalNodes(CGeometry* geometry, CConfig* config, vector<unsigned long>& internalNodes);

  /*!
  * \brief Assigning the control nodes.
  * \param[in] config -Definition of the particular problem.
  * */

  void SetCtrlNodes(CConfig* config);

  /*!
  * \brief Solving the RBF system to obtain the interpolation coefficients.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] type - Type of radial basis function.
  * \param[in] radius - Support radius of the radial basis function.
  */

  void SolveRBF_System(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, bool Derivative, vector<unsigned long>& internalNodes, bool ForwardProjectionDerivative);

  /*!
  * \brief Obtaining the interpolation coefficients of the control nodes.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] type - Type of radial basis function.
  * \param[in] radius - Support radius of the radial basis function.
  */

  void GetInterpCoeffs(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, bool Derivative, vector<unsigned long>& internalNodes, bool ForwardProjectionDerivative);


  /*!
  * \brief Gathering of all control node coordinates.
  * \param[in] geometry - Geometrical definition of the problem.
  */
  void SetCtrlNodeCoords(CGeometry* geometry, CConfig* config);

  /*!
  * \brief Build the deformation vector with surface displacements of the control nodes.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  */
  void SetBoundaryDisplacements(CGeometry* geometry, CConfig* config);

  void SetCtrlNodeDerivatives(CGeometry* geometry, CConfig* config, bool ForwardProjectionDerivative);
 
  void SetInternalNodeDerivatives(CGeometry* geometry, CConfig* config, vector<unsigned long>& internalNodes, bool ForwardProjectionDerivative);
 
  /*!
  * \brief Computation of the interpolation matrix and inverting in.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] type - Type of radial basis function.
  * \param[in] radius - Support radius of the radial basis function.
  * \param[in] invInterpMat - Inverse of the interpolation matrix.
  */
  void GetInterpMat_sequential(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius, CSymmetricMatrix& invInterpMat);
  void GetInterpMat_parallel(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius, CSymmetricMatrix& interpMat);
  void GetInvInterpMat(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, su2passivematrix& invInterpMat);
  /*!
  * \brief Computation of interpolation coefficients
  * \param[in] invInterpMat - Inverse of interpolation matrix
  */
  void ComputeInterpCoeffs(su2passivematrix& invInterpMat);

  /*!
  * \brief Finding initial data reduction control nodes based on maximum deformation.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] maxErrorNodeLocal - Local maximum error node.
  * \param[in] maxErrorLocal - Local maximum error.
  */
  void GetInitMaxErrorNode(CGeometry* geometry, CConfig* config, bool Derivative, unsigned long& maxErrorNodeLocal, su2double& maxErrorLocal);

  /*! 
  * \brief Addition of control node to the reduced set.
  * \param[in] maxErrorNode - Node with maximum error to be added.
  */
  void AddControlNode(unsigned long maxErrorNode);

  /*! 
  * \brief Compute global number of control nodes.
  */
  void Get_nCtrlNodesGlobal(CConfig* config);

  /*! 
  * \brief Compute interpolation error.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] type - Type of radial basis function.
  * \param[in] radius - Support radius of the radial basis function.
  * \param[in] maxErrorNodeLocal - Local maximum error node.
  * \param[in] maxErrorLocal - Local maximum error.
  */
  void GetInterpError(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, bool Derivative, unsigned long& maxErrorNodeLocal, su2double& maxErrorLocal);

  /*! 
  * \brief Compute error of single node.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] type - Type of radial basis function.
  * \param[in] radius - Support radius of the radial basis function.
  * \param[in] iNode - Local node in consideration.
  * \param[in] localError - Local error.
  */
  void GetNodalError(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, CRadialBasisFunctionNode* iNode, bool Derivative, su2double* localError);

  /*!
  * \brief Updating the grid coordinates.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] type - Type of radial basis function.
  * \param[in] radius - Support radius of the radial basis function.
  * \param[in] internalNodes - Internal nodes.
  */
  void UpdateGridCoord(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, const vector<unsigned long>& internalNodes);

  /*!
  * \brief Updating the internal node coordinates.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] type - Type of radial basis function.
  * \param[in] radius - Support radius of the radial basis function.
  * \param[in] internalNodes - Internal nodes.
  */
  void UpdateInternalCoords(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius, const vector<unsigned long>& internalNodes);

  /*!
  * \brief Updating the internal node coordinates.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] type - Type of radial basis function.
  * \param[in] radius - Support radius of the radial basis function.
  */
  void UpdateBoundCoords(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius);

  /*! 
  * \brief Apply correction to the nonzero error boundary nodes.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] type - Type of radial basis function.
  * \param[in] internalNodes - Internal nodes.
  */
  void SetCorrection(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const vector<unsigned long>& internalNodes);

  /*!
  * \brief Custom comparison function, for sorting the CRadialBasisFunctionNode objects based on their index.
  * \param[in] a - First considered Radial Basis Function Node.
  * \param[in] b - Second considered Radial Basis Function Node.
  * \return True if index of a is smaller than index of b.
  */
  inline static bool HasSmallerIndex(CRadialBasisFunctionNode* a, CRadialBasisFunctionNode* b){
    return a->GetIndex() < b->GetIndex();
  }

  /*!
  * \brief Custom equality function, for obtaining a unique set of CRadialBasisFunctionNode objects.
  * \param[in] a - First considered Radial Basis Function Node.
  * \param[in] b - Second considered Radial Basis Function Node.
  * \return True if index of a and b are equal.
  */
  inline static bool HasEqualIndex(CRadialBasisFunctionNode* a, CRadialBasisFunctionNode* b){
    return a->GetIndex() == b->GetIndex();
  }

  void UpdateGridCoord_Derivatives(CGeometry* geometry, CConfig* config, bool ForwardProjectionDerivative);

  template <typename T>
  vector<unsigned long> ConvertToVector(const T& container) {
    return vector<unsigned long>(container.begin(), container.end());
  }

  template <typename T>
  inline vector<unsigned long> GetIndices(unordered_map<string, T>& ctrl_nodes, const string& type) {
    return ConvertToVector(ctrl_nodes[type]);
  }
  
  void ComputeSensitivity(CGeometry* geometry,  const RADIAL_BASIS& type, const su2double radius, su2passivematrix &invInterpMat, vector<unsigned long>& internalNodes);
  
  void SetPeriodicVars(CConfig* config);

  inline su2double GetDistance(const su2double* a, const su2double*b) const {
    su2double d(0);   
    // dist = sqrt(pow(m.coords_polar_cylindrical(node1,0),2) + pow(m.coords_polar_cylindrical(node2,0),2) -2*m.coords_polar_cylindrical(node1,0)*m.coords_polar_cylindrical(node2,0)*cos(m.periodic_length/M_PI*sin( (m.coords_polar_cylindrical(node2,1)-m.coords_polar_cylindrical(node1,1))*M_PI/m.periodic_length)) + pow(m.coords_polar_cylindrical(node2,2) - m.coords_polar_cylindrical(node1,2),2) );
    if (IsCylindrical) {
      d = a[0] * a[0] + b[0]*b[0] - 2*a[0]*b[0]*cos(fabs(PeriodicAngle)/PI_NUMBER * sin((a[1]-b[1])*PI_NUMBER/fabs(PeriodicAngle)));
      if (nDim == 3){
        d += (a[2]-b[2]) * (a[2]-b[2]);
      }
    } else {
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        // su2double diff = boolPeriodic[iDim] ? (per_length[iDim]/PI_NUMBER * sin( (a[iDim] - b[iDim]) * PI_NUMBER / per_length[iDim])) : (a[iDim] - b[iDim]);
        su2double diff;
        if (PeriodicAxis[iDim]) {
            diff = PeriodicLength[iDim] / PI_NUMBER * sin((a[iDim] - b[iDim]) * PI_NUMBER / PeriodicLength[iDim]);
        } else {
            diff = a[iDim] - b[iDim];
        }
        d += diff * diff;
      }  
    }

   

    return sqrt(d);
  }

  void Cart_to_Cyl(CGeometry* geometry, CConfig* config);
  void delta_Cart_to_cyl(const su2double* init_coord_cart, const su2double* var_coord_cart, su2double* var_coord_cyl);
  void delta_cyl_to_Cart(const su2double* init_coord_cart, su2double* var_coord);

  inline void cart_to_cyl(const su2double* coord, su2double* cyl_coord) const {
    cyl_coord[0] = GeometryToolbox::Norm(2, coord);
    cyl_coord[1] = atan2(coord[1], coord[0]);

    if (nDim == 3){
      cyl_coord[2] = coord[2];
    }
  }

  inline void cyl_to_cart(const su2double* coord, su2double* cart_coord) const {
    cart_coord[0] = coord[0]*cos(coord[1]);
    cart_coord[1] = coord[0]*sin(coord[1]);

    if (nDim == 3){
      cart_coord[2] = coord[2];
    }
  }
  
  void GetNodalDeformation(CGeometry* geometry, CRadialBasisFunctionNode* iNode, su2double* varCoord);
  su2double ComputeDistance(const su2double* ctrlCoords, const su2double* targetCoords);

  unique_ptr<CADTPointsOnlyClass> CreateADT(CGeometry* geometry, const string& type);
  void ProjectBoundNodes(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, const string& nodetype, CADTPointsOnlyClass* BoundADT);
  void ApplyRBF(const su2double* coord, const RADIAL_BASIS& type, const su2double radius, su2double* new_coord);
  void ApplyProjection(CGeometry* geometry, CConfig* config, unsigned short iMarker, unsigned long pointID, su2double* coord, bool Edge3D, su2double* new_coord);
};