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
#include "../../include/adt/CADTPointsOnlyClass.hpp" 

/*!
 * \class CRadialBasisFunctionInterpolation
 * \brief Class for moving the volumetric numerical grid using Radial Basis Function interpolation.
 * \author F. van Steen
 */
class CRadialBasisFunctionInterpolation : public CVolumetricMovement {
protected:
  bool DataReduction;
  vector<su2double> CtrlNodeDeformation;  /*!< \brief Control Node Deformation.*/  // TODO -  change name such that it can also contain sensitivity

  vector<su2double> InterpCoeff;          /*!< \brief Control node interpolation coefficients.*/

  unsigned long nCtrlNodesGlobal{0};      /*!< \brief Total number of control nodes.*/
  unsigned long nCtrlNodesLocal{0};      /*!< \brief Total number of local control nodes.*/
  su2activematrix CtrlCoords;             /*!< \brief Coordinates of the control nodes.*/
  vector<su2double> CtrlNormals;
  

  su2double MaxErrorGlobal{0.0};          /*!< \brief Maximum error data reduction algorithm.*/
  vector<su2double> sensitivityUpdate;

  vector<CRadialBasisFunctionNode::NODETYPE> CtrlTypeVec;                                     /*!< \brief This vector contains the control nodes at that moment */  
  vector<CRadialBasisFunctionNode*> nodes;                        // vector containing all boundary nodes

  // unordered_map<unsigned short, vector<unsigned long>> InflationLayerEdgeNodes;

  unordered_map<unsigned short, vector<unsigned long>> layerNodes;
  int nIter;

  struct NodeTypeHash {
    std::size_t operator()(CRadialBasisFunctionNode::NODETYPE t) const noexcept {
      using U = typename std::underlying_type<CRadialBasisFunctionNode::NODETYPE>::type;
      return std::hash<U>{}(static_cast<U>(t));
    }
  };

  unordered_map<CRadialBasisFunctionNode::NODETYPE, unordered_map<unsigned short, unordered_set<unsigned long>>, NodeTypeHash> node_indices; // map containing the indices for the different type of nodes
  unordered_map<CRadialBasisFunctionNode::NODETYPE, unordered_map<unsigned short, unordered_set<unsigned long>>, NodeTypeHash> per_node_indices;
  unordered_map<CRadialBasisFunctionNode::NODETYPE, unordered_map<unsigned short, unordered_set<unsigned long>>, NodeTypeHash> ctrl_node_indices;

  vector<CRadialBasisFunctionNode*> per_nodes; // periodic nodes

  vector<su2double> PeriodicLength{0,0,0};
  su2double PeriodicAngle{0.0};
  int RotationalAxis = -1;

  bool IsCylindrical = false; 
  su2double dataReductionTolerance{0.0};
  
  static constexpr passivedouble NORMAL_THRESHOLD = M_PI/180.0; // // TODO -  add explanation 

  enum class RHS_Data {DISPLACEMENT, SENSITIVITY};

  vector<unsigned long> primaryMarker;
  vector<unsigned long> secondaryMarker;

  unordered_map<unsigned short, unsigned short> sec2prim;
  
  bool sharpEdge = false;


  bool PreserveIL = false;


  
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
  * \param[in] Derivative - Compute the derivative (disabled by default). Does not actually deform the grid if enabled.
  * \return True if surface correction has to be performed in case of data reduction.
  */
  const bool SetBoundNodes(CGeometry* geometry, CConfig* config, const bool Derivative); //done

  /*!
  * \brief Determines whether given node is the primary periodic node. 
  * \param[in] geometry - Geometrical definition of the problem. 
  * \param[in] config - Definition of the particular problem. 
  * \param[in] nodeIndex - Index of node. 
  * \return True if primary periodic node.
  */
  const bool isPrimaryPeriodicNode(CGeometry* geometry, CConfig* config, const unsigned long nodeIndex) const; //done

  /*!
  * \brief Selecting internal nodes for the volumetric deformation.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem. 
  * \param[in] internalNodes - Internal nodes.
  */
  void SetInternalNodes(CGeometry* geometry, CConfig* config, vector<unsigned long>& internalNodes); //DONE

  /*!
  * \brief Solving the RBF system to obtain the interpolation coefficients or sensitivity update.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] type - Type of radial basis function.
  * \param[in] radius - Support radius of the radial basis function.
  * \param[in] Derivative - Compute the derivative (disabled by default). Does not actually deform the grid if enabled.
  * \param[in] internalNodes - Internal nodes of the problem.
  * \param[in] ForwardProjectionDerivative - Forward computation of the derivatives.
  * \param[in] Screen_Output - determines if text is written to screen.
  */
  void SolveRBF_System(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, const bool Derivative, const vector<unsigned long>& internalNodes, const bool ForwardProjectionDerivative, const bool Screen_Output); // DONE

  void SolveRBF_System_IL(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const bool Derivative, const bool ForwardProjectionDerivative, const bool Screen_Output);

  /*!
  * \brief Obtaining the interpolation coefficients of the control nodes.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] type - Type of radial basis function.
  * \param[in] radius - Support radius of the radial basis function.
  */
  void ProjectSlidingNodes(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, const bool dataReduction, const CRadialBasisFunctionNode::NODETYPE nodeType);
  
  
  /*!
  * \brief Finds initial max error node and establishes data reduction tolerance.
  * \param[in] geometry - Geometrical definition of the problem. 
  * \param[in] config - Definition of the particular problem. 
  * \param[in] Derivative - Compute the derivative (disabled by default). Does not actually deform the grid if enabled.
  * \param[in] maxErrorNodeLocal - Local node index with maximum error. 
  * \param[in] maxErrorLocal - Maximum local error. 
  */
  void InitializeDataReduction(CGeometry* geometry, CConfig* config, bool Derivative, unsigned long& maxErrorNodeLocal, su2double& maxErrorLocal); //DONE
  void InitializeDataReduction(CGeometry* geometry, const bool Derivative, const unordered_set<unsigned long>& nodeSet, unsigned long& maxErrorNodeLocal, su2double& maxErrorLocal); 

  /*!
  * \brief Computes the interpolation coefficents of RBF interpolation system. 
  * \param[in] geometry - Geometrical definition of the problem. 
  * \param[in] config - Definition of the particular problem. 
  * \param[in] type - Type of radial basis function. 
  * \param[in] radius - Support radius of the radial basis function. 
  * \param[in] internalNode - Internal nodes to which the surface deformation is interpolated. 
  * \param[in] ForwardProjectionDerivative - Forward computation of the derivatives.
  * \param[in] rhs - Determines whether right hand side vector of RBF system contains displacement or sensitivity. 
  * \param[in] invInterpMatrix - Inverse of the interpolation matrix.
  */
  void GetInterpolationCoefficients(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius,  bool ForwardProjectionDerivative, RHS_Data rhs, su2passivematrix& invInterpMatrix);
   
  /*!
  * \brief Gathering of the control node coordinates.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  */
  void SetCtrlNodeCoords(CGeometry* geometry, CConfig* config); //DONE

  /*!
  * \brief Build the deformation vector with surface displacements of the control nodes for the rhs of the RBF system.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  */
  void SetBoundaryDisplacements(CGeometry* geometry, CConfig* config); 
  
  /*!
  * \brief Build the sensitivity vector for the rhs of the RBF system. 
  * \param[in] geometry - Geometrical definition of the problem. 
  * \param[in] config - Definition of the particular problem. 
  * \param[in] ForwardProjectionDerivative - Forward computation of the derivatives.
  */
  void SetCtrlNodeDerivatives(CGeometry* geometry, CConfig* config, bool ForwardProjectionDerivative); //DONE
  
  /*!
  * \brief Fill RHS vector of RBF system with internal node derivative values.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] internalNodes - Nodes considered as internal nodes. 
  * \param[in] ForwardProjectionDerivative - Forward computation of the derivatives. 
  */
  void SetInternalNodeDerivatives(CGeometry* geometry, CConfig* config, const vector<unsigned long>& internalNodes, bool ForwardProjectionDerivative); //DONE
 
  /*!
  * \brief Computation of the interpolation matrix sequentially.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] type - Type of radial basis function.
  * \param[in] radius - Support radius of the radial basis function.
  * \param[in] interpMat - Interpolation matrix.
  */
  void GetInterpolationMatrixSequential(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius, CSymmetricMatrix& interpMat) const ; 

  /*!
  * \brief Computation of the interpolation matrix in parallel.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] type - Type of radial basis function.
  * \param[in] radius - Support radius of the radial basis function.
  * \param[in] interpMat - Interpolation matrix.
  */
  void GetInterpolationMatrixParallel(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius, CSymmetricMatrix& interpMat) const; //DONE

  /*!
  * \breif Computation of the interpolation matrix and inverting it.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] type - Type of radial basis function.
  * \param[in] radius - Support radius of the radial basis function.
  * \param[in] invInterpMat - Inverted interpolation matrix. 
  */
  void GetInverseInterpolationMatrix(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius, su2passivematrix& invInterpMat) const;
  
  /*!
  * \brief Computation of interpolation coefficients
  * \param[in] invInterpMat - Inverse of interpolation matrix
  */
  void ComputeInterpolationCoefficients(const su2passivematrix& invInterpMat); //DONE

  /*! 
  * \brief Addition of a control node to the reduced set of control ndoes.
  * \param[in] config - Definition of particular problem. 
  * \param[in] maxErrorNode - Index of node with maximum error to be added.
  */
  void AddControlNode(CConfig* config, unsigned long maxErrorNode); //DONE

  /*! 
  * \brief Compute global number of control nodes.
  * \param[in] config - Definition of the particular problem.
  */
  void Get_nCtrlNodesGlobal(CConfig* config); //DONE

  /*! 
  * \brief Compute interpolation error.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] type - Type of radial basis function.
  * \param[in] radius - Support radius of the radial basis function.
  * \param[in] maxErrorNodeLocal - Local maximum error node.
  * \param[in] maxErrorLocal - Local maximum error.
  */
  void GetInterpError(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, bool Derivative, const vector<CRadialBasisFunctionNode::NODETYPE>& ctrltypes, unsigned long& maxErrorNodeLocal, su2double& maxErrorLocal);

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
  void UpdateGridCoord(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, const vector<unsigned long>& internalNodes, bool surfaceCorrection);

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
  void UpdateBoundCoords_IL(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius);

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

  /*!
  * \brief Compute sensitivity update
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] type - Type of radial basis function.
  * \param[in] radius - Support radius of the radial basis function.
  * \param[in] invInterpMat - Inverted interpolation matrix. 
  * \param[in] internalNodes - Nodes considered as internal nodes.
  */  
  void ComputeSensitivity(CGeometry* geometry,  const RADIAL_BASIS& type, const su2double radius, const su2passivematrix &invInterpMat, const vector<unsigned long>& internalNodes); //DONE

  /*!
  * \brief Sets required periodic parameters for the RBF interpolation.
  * \param[in] config - Definition of particular problem.
  */
  void SetPeriodicVars(CConfig* config); //done

  inline su2double GetDistance(const su2double* a, const su2double*b) const {
    su2double d(0);   

    if (IsCylindrical) {
      d = a[0] * a[0] + b[0]*b[0] - 2*a[0]*b[0]*cos(PeriodicAngle/PI_NUMBER * sin((a[1]-b[1])*PI_NUMBER/PeriodicAngle));
      if (nDim == 3){
        d += (a[2]-b[2]) * (a[2]-b[2]);
      }
    } else {
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        su2double diff = a[iDim] - b[iDim];

        if (SU2_TYPE::GetValue(PeriodicLength[iDim]) != 0.0) {
            diff = PeriodicLength[iDim] / PI_NUMBER * sin(diff * PI_NUMBER / PeriodicLength[iDim]);
        } 

        d += diff * diff;
      }  
    }
    return sqrt(d);
  }

  /*!
  * \brief Transforming all boundary nodes coordinates to cylindrical coordinates.
  * \param[in] geometry - Geometrical definition of the problem. 
  */
  void TransformBoundaryNodesToCylindricalCoords(CGeometry* geometry); //DONE

  void CartDispToCyl(const su2double* init_coord_cart, const su2double* var_coord_cart, su2double* var_coord_cyl) const;
  void CylDispToCart(const su2double* init_coord_cart, su2double* var_coord) const;

  inline void CartToCyl(const su2double* coord, su2double* cyl_coord) const {
    int idx1 = (RotationalAxis + 1) % 3;
    int idx2 = (RotationalAxis + 2) % 3;

    cyl_coord[0] = sqrt(coord[idx1] * coord[idx1] + coord[idx2] * coord[idx2]);
    cyl_coord[1] = atan2(coord[idx2], coord[idx1]);

    if (nDim == 3){
      cyl_coord[2] = coord[RotationalAxis];
    }
  }

  inline void CylToCart(const su2double* coord, su2double* cart_coord) const {
    int idx1 = (RotationalAxis + 1) % 3;
    int idx2 = (RotationalAxis + 2) % 3;

    cart_coord[idx1] = coord[0]*cos(coord[1]);
    cart_coord[idx2] = coord[0]*sin(coord[1]);

    if (nDim == 3){
      cart_coord[RotationalAxis] = coord[2];
    }
  }

  /*!
  * \brief Obtains deformation of single node.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] iNode - Node for which the deformation is obtained.
  * \param[in] varCoord - Resulting variation in coordinates.
  */  
  void GetNodalDeformation(CGeometry* geometry, CConfig* config, const CRadialBasisFunctionNode* const iNode, su2double* varCoord) const; 
  
  su2double ComputeDistance(const su2double* ctrlCoords, const su2double* targetCoords);

  void SetCorrectionSurface(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type);



  void GetPeriodicNodeErrors(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, bool Derivative);



  /*!
  * \brief Generating AD tree with specific node type for single marker.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] targetNodes - Nodes of specific type and global marker.
  * \param[in] markerLocal - Index of local marker.
  * \param[in] isPeriodic - Determines whether ADT will include periodic images.
  * \return Pointer to ADT object. 
  */
  unique_ptr<CADTPointsOnlyClass> CreateADT(CGeometry* geometry, const unordered_set<unsigned long>& targetNodes, const short markerLocal, bool isPeriodic, bool sequentialIDs) const; //DONE

  /*!
  * \brief Applies Radial Basis Function interpolation to update node coordinates.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] type - Type of radial basis function.
  * \param[in] radius - Support radius of the radial basis function.
  * \param[in] node - Node whose coordinates are updated.
  */
  void ApplyRBF(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius, CRadialBasisFunctionNode* const node);

  /*!
  * \brief Finds and sets nearest node and rank for given node and AD tree.
  * \param[in] BoundADT - AD tree used for nearest neighbor query.
  * \param[in] node - Node for nearest neighbor query.
  */
  void GetNearestNode(CADTPointsOnlyClass* const BoundADT, CRadialBasisFunctionNode* const node);

  inline void GetNearestNode(CADTPointsOnlyClass* const BoundADT, const su2double* coord, unsigned long& near_point_id, int& near_rank_id, su2double& dist) const {
    BoundADT->DetermineNearestNode(coord, dist, near_point_id, near_rank_id);
  }

  /*! 
  * \brief Sets nearest local normals for a RBF node
  * \param[in] geometry - Geometrical definition of the problem. 
  * \param[in] config - Definition of the particular problem
  * \param[in] marker - Local index of the considered marker. 
  * \param[in] node - RBF control node.
  */
  void SetLocalNormal(CGeometry* geometry, CConfig* config, CRadialBasisFunctionNode* const node);

  /*!
  * \brief Applies the projection for a RBF node
  * \param[in] geometry - Geometrical definition of the problem. 
  * \param[in] marker - Local index of the considered marker.
  * \param[in] node - RBF control node.
  * \param[in] projection - Resulting projection vector.
  */
  void ApplyProjection(CGeometry* geometry, const short marker, const CRadialBasisFunctionNode* const node, su2double* const projection) const;

  /*!
  * \brief Adds normal component of the distance vector to the projection vector.
  * \param[in] distanceVector - Displacement vector of the considered points.
  * \param[in] normal - Normal vector of nearest boundary node.
  * \param[in] projection - Projection vector to which components are added.
  */
  void AddProjectionComponent(const su2double* const distanceVector, const su2double* const normal, su2double* const projection) const; 

  /*!
  * \brief Exchanges nearest node data (normals and coordinates) between MPI ranks.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] targetNodes - Set of target nodes of specific type and marker.
  * \param[in] marker - Local index of the considered marker.
  * \return Non-local nearest node data.
  */
  
  vector<su2double> ExchangeNearestNodeData(CGeometry* geometry, CConfig* config, const unsigned short marker, const unordered_set<unsigned long>& targetSet,const vector<CRadialBasisFunctionNode*>& targetNodes, const CRadialBasisFunctionNode::NODETYPE nodetype) const;

  /*!
  * \brief Assigns nearest node data collected from other MPI ranks
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] responseRecvBuffer - non-local nearest node data. 
  * \param[in] targetNodes - Set of target nodes of specific type and marker whose nearest node data is set or updated
  */

  void SetNearestNodeData(CGeometry* geometry, CConfig* config, const vector<su2double>& responseRecvBuffer, const unordered_set<unsigned long>& targetSet, const vector<CRadialBasisFunctionNode*>& targetNodes);

  /*!
  * \brief Updates the coordinate variation for a RBF node.
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem.
  * \param[in] node - RBF node whose coordinate variation is updated. 
  * \param[in] projection - Projection vector of the node. 
  */  

  void UpdateVarCoord(CGeometry* geometry, CConfig* config, const CRadialBasisFunctionNode* const node, const su2double* projection) const;

  /*!
  * \brief Sets small components of normalized vector to zero based on threshold. 
  * \param[in] vec - Vector to be pruned on entry, pruned vector on exit.
  */

  inline void PruneVector(su2double* const vec) const {

    const su2double normalMagnitude = GeometryToolbox::Norm(nDim, vec);

    for (auto iDim = 0u; iDim < nDim; iDim++) {
      if (fabs(vec[iDim]) / normalMagnitude < NORMAL_THRESHOLD) {
        vec[iDim] = 0.0;// TODO -  add comment on why normal_threshold is as it is. 
      }    
    }
  }

  /*!
  * \brief  Obtains the normal(s) of a given vertex
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] config - Definition of the particular problem. 
  * \param[in] marker - Local marker index.
  * \param[in] nodetype - Type of RBF node.
  * \param[in] nearestVertex - Nearest boundary vertex, whose normals are obtained.
  * \param[in] normal - Resulting normal vector of the vertex. 
  * \param[in] secondaryNormal - Resulting secondary normal vector of the vertex (only for 3D and edge nodetypes). */

  void GetVertexNormals(CGeometry* geometry, CConfig* config, const unsigned short markerLocal, const CRadialBasisFunctionNode::NODETYPE nodetype, CVertex* const nearestVertex, su2double* normal, su2double* secondaryNormal) const;
       
  /*!
  * \brief Returns nearest vertex after ADT query. 
  * \param[in] geometry - Geometrical definition of the problem.
  * \param[in] pointID - ID of the nearest node in ADT query.
  * \param[in] marker - Local marker index. 
  */
  CVertex* GetNearestVertex(CGeometry* geometry, const unsigned long pointID, const unsigned short marker) const;

  /*!
  * \brief Obtains coordinates of the nearest vertex 
  * \param[in] nearestVertex - Nearest vertex.
  * \param[in] nearestCoord - Nearest coordinate. 
  */
  void GetNearestCoord( CVertex* const nearestVertex, su2double* const nearestCoord ) const;

  void SetInternalNodesDerivative(CGeometry* geometry, CConfig* config, vector<unsigned long>& internalNodes);

  void CheckSharpEdge(CGeometry* geometry, unsigned long iNode);
  unsigned long retrieveIndex(unsigned long index_in);
  void SetCorrection_IL(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const vector<unsigned long>& internalNodes, unsigned short marker);
  bool Opps(CGeometry* geometry, unordered_set<unsigned long> &nodeIndices);
  su2double GetRbfWeight(const su2double* normal1, const su2double* normal2) const ;
};


