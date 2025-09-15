/*!
 * \file CRadialBasisFunctionInterpolation.cpp
 * \brief Subroutines for moving mesh volume elements using Radial Basis Function interpolation.
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

#include "../../include/grid_movement/CRadialBasisFunctionInterpolation.hpp"
#include "../../include/interface_interpolation/CRadialBasisFunction.hpp"

using NODETYPE = CRadialBasisFunctionNode::NODETYPE;
constexpr passivedouble CRadialBasisFunctionInterpolation::NORMAL_THRESHOLD; 

CRadialBasisFunctionInterpolation::CRadialBasisFunctionInterpolation(CGeometry* geometry, CConfig* config) : CVolumetricMovement(geometry) {}

CRadialBasisFunctionInterpolation::~CRadialBasisFunctionInterpolation() = default;


void CRadialBasisFunctionInterpolation::SetVolume_Deformation(CGeometry* geometry, CConfig* config, bool UpdateGeo, bool Derivative,
                                                bool ForwardProjectionDerivative){


  // TODO -  Debug MPI 
  // {
  //   if (rank == 0){
  //     int i = 0;
  //     while(0 == i){
  //       sleep(1);
  //     }
  //   }
  // }
  // if (rank == MASTER_NODE) cout << "this is the development code!" << endl;
  

  /*--- Retrieve type of RBF and its support radius ---*/ 
  const auto kindRBF = config->GetKindRadialBasisFunction();
  const su2double radius = config->GetRadialBasisFunctionParameter();
  
  su2double MinVolume, MaxVolume;

  /*--- Retrieving number of deformation steps and screen output from config ---*/

  const auto Nonlinear_Iter = config->GetGridDef_Nonlinear_Iter();
  auto Screen_Output = config->GetDeform_Output();
  
  /*--- Disable the screen output if we're running SU2_CFD ---*/
  if (config->GetKind_SU2() == SU2_COMPONENT::SU2_CFD && !Derivative) Screen_Output = false;
  if (config->GetSmoothGradient()) Screen_Output = true;

  /*--- Setting periodic variables if neccessary ---*/                                               
  if (config->GetnMarker_Periodic() != 0) SetPeriodicVars(config);

  /*--- Determining the boundary and internal nodes. Setting the control nodes. ---*/ 
  const bool surfaceCorrection = SetBoundNodes(geometry, config, Derivative);
  
  vector<unsigned long> internalNodes; 
  SetInternalNodes(geometry, config, internalNodes); 
  
  /*--- Looping over the number of deformation iterations ---*/
  for (auto iNonlinear_Iter = 0ul; iNonlinear_Iter < Nonlinear_Iter; iNonlinear_Iter++) {
    
    /*--- Compute min volume in the entire mesh. ---*/
    ComputeDeforming_Element_Volume(geometry, MinVolume, MaxVolume, Screen_Output);
    if (rank == MASTER_NODE && Screen_Output)
      cout << "Min. volume: " << MinVolume << ", max. volume: " << MaxVolume << "." << endl;
    

    /*--- Solving the RBF system, resulting in the interpolation coefficients ---*/
    SolveRBF_System(geometry, config, kindRBF, radius, Derivative, internalNodes, ForwardProjectionDerivative, Screen_Output);
   
    /*--- Updating the coordinates of the grid ---*/
    if (!Derivative) {
      UpdateGridCoord(geometry, config, kindRBF, radius, internalNodes, surfaceCorrection);
    } else {
      UpdateGridCoord_Derivatives(geometry, config, ForwardProjectionDerivative);
    }

    if(UpdateGeo){
      UpdateDualGrid(geometry, config);
    }
    
    if (!Derivative){
      /*--- Check for failed deformation (negative volumes). ---*/

      ComputeDeforming_Element_Volume(geometry, MinVolume, MaxVolume, Screen_Output);

      /*--- Calculate amount of nonconvex elements ---*/

      ComputenNonconvexElements(geometry, Screen_Output);
    }
   
    if (rank == MASTER_NODE && Screen_Output) {
      cout << "Non-linear iter.: " << iNonlinear_Iter + 1 << "/" << Nonlinear_Iter << ". ";
      if (nDim == 2)
        cout << "Min. area: " << MinVolume <<  "." << endl;
      else
        cout << "Min. volume: " << MinVolume <<  "." << endl;
    }  
  }  


  // deallocate memory
  for (auto node : nodes) {
      delete node;
  }
  nodes.clear();

  for (auto node : per_nodes) {
      delete node;
  }
  per_nodes.clear();
}

// BUG -   there is likely an error in how the projection is done when there is more than a single deformation step
void CRadialBasisFunctionInterpolation::SolveRBF_System(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, const bool Derivative, const vector<unsigned long>& internalNodes, const bool ForwardProjectionDerivative, const bool Screen_Output){

  const auto dataReduction = config->GetRBF_DataReduction();
  unsigned long maxErrorNodeLocal;
  unsigned short greedyIterations = 0;
  su2double maxErrorLocal{0};
  const RHS_Data rhs = (Derivative && dataReduction) ?  RHS_Data::SENSITIVITY : RHS_Data::DISPLACEMENT;
  vector <NODETYPE> slideCtrlTypes = dataReduction ? vector<NODETYPE>{NODETYPE::EDGE} : vector<NODETYPE>{NODETYPE::EDGE, NODETYPE::SURFACE};
  bool RBF_StageSurface = dataReduction ? false : true;
  su2passivematrix invInterpMat;

  if (IsCylindrical) TransformBoundaryNodesToCylindricalCoords(geometry);

  if (dataReduction && nCtrlNodesLocal == 0) { 
    InitializeDataReduction(geometry, config, Derivative, maxErrorNodeLocal, maxErrorLocal);
  } 

  /*--- While the maximum error is above the tolerance, data reduction algorithm is continued. ---*/
  while(((MaxErrorGlobal > dataReductionTolerance)  || greedyIterations == 0)){ 
    CtrlTypeVec.resize(1);
    CtrlTypeVec[0] = NODETYPE::DISPLACED;
      
    /*--- In case of a nonzero local error, control nodes are added ---*/
    if(maxErrorLocal> 0){
      AddControlNode(config, maxErrorNodeLocal);
    }
    
    GetInterpolationCoefficients(geometry, config, type, radius, internalNodes, ForwardProjectionDerivative, rhs, invInterpMat);

    for (const auto &iType : slideCtrlTypes) {
      ProjectSlidingNodes(geometry, config, type, radius, dataReduction, iType);    

      CtrlTypeVec.push_back(iType);
      GetInterpolationCoefficients(geometry, config, type, radius, internalNodes, ForwardProjectionDerivative, rhs, invInterpMat); 
    }

    GetInterpError(geometry, config, type, radius, Derivative, slideCtrlTypes, maxErrorNodeLocal, maxErrorLocal);     

    
    if (!RBF_StageSurface && MaxErrorGlobal < dataReductionTolerance && nDim == 3) {      
      MaxErrorGlobal = dataReductionTolerance + std::numeric_limits<su2double>::epsilon();
      slideCtrlTypes.push_back(NODETYPE::SURFACE);
      RBF_StageSurface = true;
    }
    
    greedyIterations++;

    if (rank == MASTER_NODE && Screen_Output && greedyIterations % 1 == 0 ) {
      cout << "Greedy iteration: " << greedyIterations
           << ". Max error: " << MaxErrorGlobal 
           << ". Global nr. of ctrl nodes: "  << nCtrlNodesGlobal << endl;
    } 
  } 
  GetPeriodicNodeErrors(geometry, config, type, radius, Derivative);
// TODO -  remove
  // cout << "rank: " << rank << " local nCtrl: " << nCtrlNodesLocal << " global nCtrl: " << nCtrlNodesGlobal <<  " total points: " << geometry->GetnPoint() << " global n points: " << geometry->GetGlobal_nPoint() << " internal nodes: " << internalNodes.size() << endl; 
  
  /*--- Once the data reduction tolerance is reached the error for the periodic nodes has to be determined. ---*/
  if (Derivative){
    SetInternalNodeDerivatives(geometry, config, internalNodes, ForwardProjectionDerivative);
    ComputeSensitivity(geometry, type, radius, invInterpMat, internalNodes); 
  }
}

void CRadialBasisFunctionInterpolation::ProjectSlidingNodes(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, const bool dataReduction, const CRadialBasisFunctionNode::NODETYPE nodeType) {
  
  const auto& markerToNodeSet = dataReduction ? ctrl_node_indices[nodeType] : node_indices[nodeType];
  
  for (const auto& iMarker : markerToNodeSet){
    /*--- Obtaining the local marker ---*/
    const auto markerGlobal = iMarker.first;
    const auto markerLocal = config->GetMarker_Local(markerGlobal);

    /*--- Create global ADT for points of specific marker and nodetype.
            For periodic domains also periodic images are included, 
            since periodic distance is not supported. ---*/
    const auto& targetNodesADT = node_indices[nodeType][markerGlobal];
    const bool isPeriodic = config->GetnMarker_Periodic() != 0;
    auto BoundADT = CreateADT(geometry, targetNodesADT, markerLocal, isPeriodic);

    /*--- Sliding nodes ---*/
    const auto& targetNodes = iMarker.second;

    /*--- Apply RBF deformation and obtain nearest boundary node after deformation ---*/ 
    for (const auto iNode : targetNodes) {
      auto* const node = nodes[iNode];
      ApplyRBF(geometry, type, radius, node); // TODO -  start here
      GetNearestNode(BoundADT.get(), node);
    }

    /*--- Communicate nearest node data among ranks and store it ---*/
    auto response_recv_buffer = ExchangeNearestNodeData(geometry, config, markerLocal, targetNodes, nodes,  nodeType);
    SetNearestNodeData(geometry, config, response_recv_buffer, targetNodes, nodes);
    
    /*--- Apply projection and set resulting position as coordinate variation ---*/      
    for (const auto iNode : targetNodes) {
      const auto* const node = nodes[iNode];
      su2double projection[3] = {0.0, 0.0, 0.0};
      ApplyProjection(geometry, markerLocal, node, projection);
      UpdateVarCoord(geometry, config, node, projection);
    }
  }
}

void CRadialBasisFunctionInterpolation::InitializeDataReduction(CGeometry* geometry, CConfig* config, const bool Derivative, unsigned long& maxErrorNodeLocal, su2double& maxErrorLocal ) {
  /*--- The initial maximum error node is found based on maximum displacement 
          or sensitivity of the displaced boundary nodes. 
          The data reduction tolerance is then obtained as a fraction 
          (specified in config file of the maximum error. ---*/
  
  maxErrorLocal = 0.0;
  su2double normSquaredDeformation = 0;
  passivedouble sens_i = 0.0;

  const auto& markers = node_indices[NODETYPE::DISPLACED]; 
  for (const auto& iMarker : markers) {
    for (const auto& iNode : iMarker.second) {
      normSquaredDeformation = 0.0;

      /*--- Compute to squared norm of the deformation ---*/
      if (Derivative) {
        const auto iPoint = nodes[iNode]->GetIndex();
        for (auto iDim = 0u; iDim < nDim; iDim++) {
          sens_i = SU2_TYPE::GetValue(geometry->GetSensitivity(iPoint, iDim));
          normSquaredDeformation += sens_i * sens_i;
        }
      } else {
        normSquaredDeformation = GeometryToolbox::SquaredNorm(nDim, geometry->vertex[nodes[iNode]->GetMarker()][nodes[iNode]->GetVertex()]->GetVarCoord());    
      }
      
      /*--- In case squared norm deformation is larger than the error, update the error ---*/
      if(normSquaredDeformation > maxErrorLocal) {
        maxErrorLocal = normSquaredDeformation;
        maxErrorNodeLocal = iNode;
      }
    }
  }

  /*--- Account for the possibility of applying the deformation in multiple steps ---*/
  maxErrorLocal = sqrt(maxErrorLocal) / ((su2double)config->GetGridDef_Nonlinear_Iter());

  SU2_MPI::Allreduce(&maxErrorLocal, &MaxErrorGlobal, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());  
    
  /*--- Error tolerance for the data reduction tolerance ---*/
  dataReductionTolerance = config->GetRBF_DataRedTolerance() * MaxErrorGlobal; 
  if (rank == MASTER_NODE) cout << "DATA REDUCTION TOLERANCE: " << dataReductionTolerance << "\n" << endl; // TODO -  placed elsewhere in the SU2 output
}

void CRadialBasisFunctionInterpolation::GetInterpolationCoefficients(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, const vector<unsigned long>& internalNodes, bool ForwardProjectionDerivative, RHS_Data rhs, su2passivematrix& invInterpMat){
  
  /*---Update global number of control nodes ---*/
  Get_nCtrlNodesGlobal(config);

  /*--- Gather and distribute control node coordinates ---*/
  SetCtrlNodeCoords(geometry, config);

  /*--- Set right hand side vector ---*/
  if (rhs == RHS_Data::DISPLACEMENT) { 
    SetBoundaryDisplacements(geometry, config);
  }else{
    SetCtrlNodeDerivatives(geometry, config, ForwardProjectionDerivative);
  }

  /*--- Computate the (inverse) interpolation matrix. ---*/
  GetInverseInterpolationMatrix(geometry, type, radius, invInterpMat);

  /*--- Compute the interpolation coefficients. ---*/
  ComputeInterpolationCoefficients(invInterpMat);
}

void CRadialBasisFunctionInterpolation::ComputeSensitivity(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius, const su2passivematrix &invInterpMat, const vector<unsigned long>& internalNodes){
  /*--- The sensitivity update for the control nodes is given by: 
          delta_sens_c = inv(Phi_cc) * Phi_i,c^T * sens_i = inv(Phi_cc) * Phi_ci * sens_i
          Where Phi are the interpolation matrices and c and i represent the control and internal nodes.
          This computation is split into two parts, the calculation of Phi_ci * sens_i resulting in a intermediate product vector,
          and inv(Phi_cc) * intermediateProduct.   ---*/

  /*--- Compute Phi_ci * sens_i, result is stored in intermediateProduct ---*/
  const auto vectorSize = nCtrlNodesGlobal * nDim;
  vector<su2double> intermediateProductLocal(vectorSize, 0.0);

  for (auto iNode = 0ul; iNode < nCtrlNodesGlobal; iNode++){
    const auto iOffset = iNode * nDim;
    for (auto iDim =0u; iDim < nDim; iDim++){      
      for (auto jNode =0ul; jNode < internalNodes.size(); jNode++){
        const su2double dist = GetDistance(CtrlCoords[iOffset], geometry->nodes->GetCoord(internalNodes[jNode]));
        const su2double rbf_eval = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));
        intermediateProductLocal[iOffset + iDim] += rbf_eval * CtrlNodeDeformation[jNode*nDim+iDim];
      }
    }
  }

  /*--- Sum contributions of all ranks ---*/
  vector<su2double> intermediateProduct(vectorSize);
  SU2_MPI::Allreduce(intermediateProductLocal.data(), intermediateProduct.data(), vectorSize, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());


  /*--- Compute the sensitivity update: inv(Phi_cc) * IntermediateProduct ---*/
  vector<su2double> sensUpdateLocal(vectorSize, 0.0);  
  for (auto iNode = rank; iNode < nCtrlNodesGlobal; iNode += size) {
    const auto iOffset = iNode * nDim;
    for (auto iDim =0u; iDim < nDim; iDim++) {
      for (auto jNode = 0ul; jNode < nCtrlNodesGlobal; jNode++) {
        sensUpdateLocal[iOffset + iDim] += invInterpMat(iNode, jNode) * intermediateProduct[jNode*nDim+iDim];
      }
    }
  }

  /*--- Sum contributions of all ranks ---*/
  sensitivityUpdate.resize(vectorSize);
  SU2_MPI::Allreduce(sensUpdateLocal.data(), sensitivityUpdate.data(), vectorSize, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
}



const bool CRadialBasisFunctionInterpolation::SetBoundNodes(CGeometry* geometry, CConfig* config, const bool Derivative){
  bool surfaceCorrection = false;
  
  auto add_node = [&](vector<CRadialBasisFunctionNode*>& list, unsigned long iNode, unsigned int iMarker, unsigned long iVertex, NODETYPE type) {
        list.push_back(new CRadialBasisFunctionNode(iNode, iMarker, iVertex));
        list.back()->SetNodeType(type);
    };
  
  if (Derivative) {
    for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
      if (!config->GetMarker_All_SendRecv(iMarker) && !config->GetMarker_All_PerBound(iMarker)) {

        for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) {

          const auto iNode = geometry->vertex[iMarker][iVertex]->GetNode(); 

          if (!geometry->nodes->GetDomain(iNode)) continue;

          if(geometry->nodes->GetPeriodicBoundary(iNode)){
            auto& target_list = isPrimaryPeriodicNode(geometry, config, iNode) ? nodes : per_nodes;
            add_node(target_list, iNode, iMarker, iVertex, NODETYPE::DISPLACED);
          }else{
            add_node(nodes, iNode, iMarker, iVertex, NODETYPE::DISPLACED);
          }
          
        }
      }
    }
  } else {


      for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {

        /*--- If the marker is not a send/receive or periodic marker ---*/
        if (config->GetMarker_All_SendRecv(iMarker) || config->GetMarker_All_PerBound(iMarker)) 
          continue;

        const bool isSliding = config->GetMarker_All_Deform_Mesh_Slide(iMarker);
      
        for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          
          /*--- Get index of node, if not part of domain then skip ---*/
          const auto iNode = geometry->vertex[iMarker][iVertex]->GetNode(); 
          if (!geometry->nodes->GetDomain(iNode)) continue;
          
          if (isSliding) {
            /*--- Count how many boundaries include this specific node ---*/
            unsigned short nVertex = 0;          
            for (auto jMarker = 0u; jMarker < config->GetnMarker_All(); jMarker++){
              if (geometry->nodes->GetVertex(iNode, jMarker) != -1 && !config->GetMarker_All_SendRecv(jMarker)){
                nVertex++;
              }
            }

            if (nVertex == 1) {
              add_node(nodes, iNode, iMarker, iVertex, nDim == 2 ? NODETYPE::EDGE : NODETYPE::SURFACE);
            } else if (nVertex == nDim - 1) {

               /*--- If part of 2 markers (only for 3D). These are nodes on the interface of two markers. 
                    For a periodic marker these are considered as sliding surface nodes, else as edge nodes. ---*/

              if (geometry->nodes->GetPeriodicBoundary(iNode)) { 
                auto& target_list = isPrimaryPeriodicNode(geometry, config, iNode) ? nodes : per_nodes;
                add_node(target_list, iNode, iMarker, iVertex, NODETYPE::SURFACE);
              } else {

                 /*--- Check if node is part of a moving boundary (has prescribed displacement). 
                      If yes, then it's skipped and that node is included through the moving marker.
                      But it also means that there is a hole in the domain boundary. 
                      This needs to be accounted for when doing the correction for data reduction. ---*/             
                
                bool isMoving = false;
                for (auto jMarker = 0u; jMarker < config->GetnMarker_All(); jMarker++){
                  if ( config->GetMarker_All_DV(jMarker) && geometry->nodes->GetVertex(iNode, jMarker) != -1) {
                    surfaceCorrection = true;
                    isMoving = true;
                    break;
                  }
                }
                
                /*--- If not part of moving boundary then it's simply an edge node ---*/
                if (!isMoving) {
                  add_node(nodes, iNode, iMarker, iVertex, NODETYPE::EDGE);
                }
              }
            } else {
                /*--- These are nodes that are the vertices of the domain. 
                If periodic then they are classified as sliding edge nodes, else as displaced nodes. */

                if(geometry->nodes->GetPeriodicBoundary(iNode)){                  
                  auto& target_list = isPrimaryPeriodicNode(geometry, config, iNode) ? nodes : per_nodes;
                  add_node(target_list, iNode, iMarker, iVertex, NODETYPE::EDGE);
                } else {
                  add_node(nodes, iNode, iMarker, iVertex, NODETYPE::DISPLACED);
                }
            }
          } else {
              if(geometry->nodes->GetPeriodicBoundary(iNode)){
                auto& target_list = isPrimaryPeriodicNode(geometry, config, iNode) ? nodes : per_nodes;
                add_node(target_list, iNode, iMarker, iVertex, NODETYPE::DISPLACED);
              } else {
                  add_node(nodes, iNode, iMarker, iVertex, NODETYPE::DISPLACED);
              }
          }
        }
      }
  }
      

  // ensuring that all global markers are present in the mappings. 
  for (auto iMarker = 0u; iMarker < config->GetnMarker_CfgFile(); iMarker++) {
    if (config->GetMarker_CfgFile_Deform_Mesh_Slide(config->GetMarker_CfgFile_TagBound(iMarker))) {

      node_indices[NODETYPE::EDGE][iMarker];
      node_indices[NODETYPE::SURFACE][iMarker];

      per_node_indices[NODETYPE::EDGE][iMarker];
      per_node_indices[NODETYPE::SURFACE][iMarker];

      if (config->GetRBF_DataReduction()) {
        ctrl_node_indices[NODETYPE::EDGE][iMarker];
        ctrl_node_indices[NODETYPE::SURFACE][iMarker];
      }
    }
  }

  // sort and make unique
  stable_sort(nodes.begin(), nodes.end(), HasSmallerIndex);
  nodes.resize(distance(nodes.begin(), unique(nodes.begin(), nodes.end(), HasEqualIndex)));

  stable_sort(per_nodes.begin(), per_nodes.end(), HasSmallerIndex);
  per_nodes.resize(distance(per_nodes.begin(), unique(per_nodes.begin(), per_nodes.end(), HasEqualIndex)));

  // filling the mappings
  for(unsigned long x = 0ul; x < nodes.size(); x++){

    auto marker = nodes[x]->GetMarker(); // this is the local marker index
    auto tag = config->GetMarker_All_TagBound(marker);
    auto global_marker = config->GetMarker_CfgFile_TagBound(tag);

    auto type = nodes[x]->GetNodeType();
    node_indices[type][global_marker].insert(x);
  }

  for (unsigned long x = 0ul; x < per_nodes.size(); x++) {

    auto marker = per_nodes[x]->GetMarker();
    auto tag = config->GetMarker_All_TagBound(marker);
    auto global_marker = config->GetMarker_CfgFile_TagBound(tag);
    
    auto type = per_nodes[x]->GetNodeType();
    per_node_indices[type][global_marker].insert(x);
  }
  // TODO -  debug output 
  
  ofstream out("edge"+to_string(rank)+".txt");
  auto idx = node_indices[NODETYPE::EDGE];
  for (auto& pair : idx) {
      const unsigned short& marker = pair.first;      // marker (key)
      const auto& indices = pair.second;  // indices (value)

      for (auto i : indices) {
          out << geometry->nodes->GetGlobalIndex(nodes[i]->GetIndex())  /*<< "\t" << nodes[i]->GetMarker() << "\t" << marker*/ << endl;
      }
  }
  out.close();

  // ofstream out2("surface"+to_string(rank)+".txt");
  // auto idx2 = node_indices[NODETYPE::SURFACE];
  // for (auto& pair : idx2) {
  //     const unsigned short& marker = pair.first;      // marker (key)
  //     const std::vector<unsigned long>& indices = pair.second;  // indices (value)

  //     for (auto i : indices) {
  //         out2 << geometry->nodes->GetGlobalIndex(nodes[i]->GetIndex()) << endl;
  //     }
  // }
  // out2.close();


  ofstream out3("disp"+to_string(rank)+".txt");
  auto idx3 = node_indices[NODETYPE::DISPLACED];
  for (auto& pair : idx3) {
      const unsigned short& marker = pair.first;      // marker (key)
      const auto& indices = pair.second;  // indices (value)

      for (auto i : indices) {
          out3 << geometry->nodes->GetGlobalIndex(nodes[i]->GetIndex()) << endl;
      }
  }
  out3.close();

  // ofstream out4("nodes"+to_string(rank)+".txt");
  // for (auto i : nodes){
  //   out4 << geometry->nodes->GetGlobalIndex(i->GetIndex()) << /*"\t" <<  i->getNodetype()  <<*/  endl;
  // }
  // out4.close();

  ofstream out5("per_nodes"+to_string(rank)+".txt");
  for (auto i : per_nodes){
    out5 << geometry->nodes->GetGlobalIndex(i->GetIndex()) << /*"\t" <<  i->getNodetype()  <<*/  endl;
  }
  out5.close();
  return surfaceCorrection;
}

const bool CRadialBasisFunctionInterpolation::isPrimaryPeriodicNode(CGeometry* geometry, CConfig* config, const unsigned long nodeIndex) const {
  bool isPrimaryNode = false;
  
  /*--- Finding the corresponding periodic marker index ---*/
  unsigned short perMarker;              
  for (auto jMarker = 0u; jMarker < config->GetnMarker_All(); jMarker++){
    if (geometry->nodes->GetVertex(nodeIndex, jMarker) != -1 && config->GetMarker_All_PerBound(jMarker)) {

      /*--- If node is on periodic marker with the smaller periodic index it is considered as primary periodic node ---*/
      const auto idx_PeriodicMarker = config->GetMarker_Periodic(config->GetMarker_All_TagBound(jMarker));
      if(idx_PeriodicMarker < config->GetMarker_Periodic_Donor2(idx_PeriodicMarker)){
        return isPrimaryNode = true;
      }
    }
  }
  return isPrimaryNode;
}

void CRadialBasisFunctionInterpolation::GetInverseInterpolationMatrix(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius, su2passivematrix& invInterpMat) const {
  CSymmetricMatrix interpMat;

  /*---  Build the interpolation matrix in parallel or sequentially ---*/
  #ifdef HAVE_MPI
    GetInterpolationMatrixParallel(geometry, type, radius, interpMat);    
  #else
    GetInterpolationMatrixSequential(geometry, type, radius, interpMat);    
  #endif
  
  /*--- Check if the kernel is symmetric positive definite ---*/
  const bool kernelIsSPD = (type == RADIAL_BASIS::WENDLAND_C2) || (type == RADIAL_BASIS::GAUSSIAN) ||
                          (type == RADIAL_BASIS::INV_MULTI_QUADRIC);

  /*--- Inverting matrix and transferring data to output variable ---*/
  interpMat.Invert(kernelIsSPD); 
  invInterpMat = interpMat.StealData();
}

void CRadialBasisFunctionInterpolation::GetInterpolationMatrixParallel(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius, CSymmetricMatrix& interpMat) const{
  vector<su2double> rbfValsGlobal, rbfValsLocal;
  unsigned long startRow, endRow;

  /*--- Since the matrix is symmetric only lower half needs to be computed.
          Number of elements of the lower half matrix: n(n+1)/2
          Must be a double for the ceil operation. ---*/
  const passivedouble nLowerTriangle = (nCtrlNodesGlobal*(nCtrlNodesGlobal+1))/2;

  /*--- Average number of RBF evaluations per process ---*/
  const su2double nProcess = ceil(nLowerTriangle/size);
  
  /*--- To balance the number of RBF evaluation over the processes, a start and end row of the global matrix are assigned for each process.
          The start row is be found with: row(row+1)/2 = rank*nProcess -> row = (-1 + sqrt(1+8*rank*nProcess))/2
          The end row is found with:      row(row+1)/2 = (rank+1)*nProcess -> row = (-1 + sqrt(1+8*(rank+1)*nProcess))/2
          Checks are applied to ensure that the row values are not out of bounds. ---*/
  startRow = ceil((-1 + sqrt(1+8*nProcess*(rank))) / 2);
  endRow = ceil((-1 + sqrt(1+8*nProcess*(rank+1))) / 2);
  if (startRow > nCtrlNodesGlobal) startRow = nCtrlNodesGlobal;
  if (endRow > nCtrlNodesGlobal)  endRow = nCtrlNodesGlobal;  
  
  /*--- Number of RBF evaluations to be carried out ---*/
  const int nElems = (endRow*(endRow+1) - startRow*(startRow+1))/2;

  /*--- Obtaining local RBF evaluations ---*/
  rbfValsLocal.resize(nElems);
  auto iValue = 0ul;
  for (auto iRow = startRow; iRow < endRow; iRow++) {
    const auto offset = iRow * nDim;
    for (auto iCol = 0ul; iCol <= iRow; iCol++) {
      su2double dist = GetDistance(CtrlCoords[offset], CtrlCoords[iCol * nDim]);      
      rbfValsLocal[iValue++] = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));
    }
  }

  /*--- Distributing the local RBF evaluations ---*/
  vector<int> recvCnts(size), disp(size);
  rbfValsGlobal.resize(nLowerTriangle);
  
  SU2_MPI::Allgather(&nElems, 1, MPI_INT, recvCnts.data(), 1, MPI_INT, SU2_MPI::GetComm());

  disp[0] = 0;
  for (int i = 1; i < size; i++) {
      disp[i] = disp[i-1] + recvCnts[i-1];
  }

  SU2_MPI::Allgatherv(rbfValsLocal.data(), nElems, MPI_DOUBLE, rbfValsGlobal.data(), recvCnts.data(), disp.data(), MPI_DOUBLE, SU2_MPI::GetComm());


  /*--- Filling lower part of the interpolation matrix ---*/
  interpMat.Initialize(nCtrlNodesGlobal);
  iValue = 0;
  for (auto iRow = 0ul; iRow < nCtrlNodesGlobal; iRow++) {
    for (auto iCol = 0ul; iCol <= iRow; iCol++) {    
      interpMat.Set(iRow, iCol, SU2_TYPE::GetValue(rbfValsGlobal[iValue++]));      
    }
  }
}

void CRadialBasisFunctionInterpolation::GetInterpolationMatrixSequential(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius, CSymmetricMatrix& interpMat) const {
  interpMat.Initialize(nCtrlNodesGlobal);
    
  /*--- Construction of the interpolation matrix. 
      Since this matrix is symmetric only upper halve has to be considered.
      Evaluates the distance between control nodes and fills matrix with RBF evaluations ---*/  
  
  for (auto iNode = 0ul; iNode < nCtrlNodesGlobal; iNode++) {
    const auto iOffset = iNode * nDim;    
    for (auto jNode = iNode; jNode < nCtrlNodesGlobal; jNode++) {
      const auto jOffset = jNode * nDim;
      const su2double dist = GetDistance(CtrlCoords[iOffset], CtrlCoords[jOffset]);
      interpMat(iNode, jNode) = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));
    }
  }  
}


void CRadialBasisFunctionInterpolation::SetBoundaryDisplacements(CGeometry* geometry, CConfig* config){
  const SU2_COMPONENT Kind_SU2 = config->GetKind_SU2();
  const bool dataRed = config->GetRBF_DataReduction();
  unsigned long ctrlNodeIndex = 0;

  CtrlNodeDeformation.resize(nCtrlNodesLocal*nDim, 0.0);   
  
  for (const auto& iType : CtrlTypeVec) {
    const auto& markerToNodeSet = dataRed ? ctrl_node_indices[iType] : node_indices[iType];
    
    for (const auto& iMarker : markerToNodeSet){
      const auto& markerLocal = config->GetMarker_Local(iMarker.first); 
      
      /*--- If the marker is moving then its nodal displacement is added to the deformation vector,
              else it is set to zero ---*/
      const bool isMoving = ((config->GetMarker_All_Moving(markerLocal) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_CFD)) ||
        ((config->GetMarker_All_DV(markerLocal) == YES || config->GetMarker_All_Deform_Mesh_Slide(markerLocal) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_DEF)) ||
        ((config->GetDirectDiff() == D_DESIGN) && (Kind_SU2 == SU2_COMPONENT::SU2_CFD) &&
        (config->GetMarker_All_DV(markerLocal) == YES)) /*NOTE: This feature has not been tested for RBF interpolation*/ ||
        ((config->GetMarker_All_DV(markerLocal) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_DOT));

      for (const auto iNode : iMarker.second) {
        auto baseIndex = ctrlNodeIndex * nDim;
        if (isMoving) {
          su2double varCoord[3] = {0.0, 0.0, 0.0};
          GetNodalDeformation(geometry, config, nodes[iNode], varCoord);
          
          for (auto iDim = 0u; iDim < nDim; iDim++){
            CtrlNodeDeformation[baseIndex + iDim] = SU2_TYPE::GetValue(varCoord[iDim]);
          }
        } else {
            for (auto iDim = 0u; iDim < nDim; iDim++) {
              CtrlNodeDeformation[baseIndex + iDim] = 0.0;
            }
        }
        ctrlNodeIndex++;
      }      
    }
  }
}

void CRadialBasisFunctionInterpolation::GetNodalDeformation(CGeometry* geometry, CConfig* config, const CRadialBasisFunctionNode* const iNode, su2double* varCoord) const {
  /*--- If requested (no by default) impose the surface deflections in
    increments and solve the grid deformation with
    successive small deformations. ---*/
  const su2double VarIncrement = 1.0 / ((su2double)config->GetGridDef_Nonlinear_Iter());

  const su2double* varCoordCartOrig = geometry->vertex[iNode->GetMarker()][iNode->GetVertex()]->GetVarCoord();
  su2double varCoordCart[3] = {0.0, 0.0, 0.0};
  for (auto iDim = 0u; iDim < nDim; iDim++) {
    varCoordCart[iDim] = varCoordCartOrig[iDim]*VarIncrement;
  }

  if (IsCylindrical) {
    const su2double* coord = geometry->nodes->GetCoord(iNode->GetIndex());
    CartDispToCyl(coord, varCoordCart, varCoord);
  } else{
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      varCoord[iDim] = varCoordCart[iDim];
    }
  }
}

void CRadialBasisFunctionInterpolation::SetInternalNodeDerivatives(CGeometry* geometry, CConfig* config, const vector<unsigned long>& internalNodes, bool ForwardProjectionDerivative) {
  const SU2_COMPONENT Kind_SU2 = config->GetKind_SU2();
  
  /*--- Early exit for unsupported features ---*/
  if ((Kind_SU2 == SU2_COMPONENT::SU2_DOT) && !ForwardProjectionDerivative) {
    CtrlNodeDeformation.resize(internalNodes.size() * nDim, 0.0); 
  } else {
    SU2_MPI::Error("Missing feature in RBF interpolation", CURRENT_FUNCTION);
  }

  /*--- Set sensitivity values in the RHS vector of the RBF system ---*/
  for (auto iNode = 0ul; iNode < internalNodes.size(); iNode++) {
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      CtrlNodeDeformation[iNode * nDim + iDim] = SU2_TYPE::GetValue(geometry->GetSensitivity(internalNodes[iNode], iDim));
    }
  }
}

void CRadialBasisFunctionInterpolation::SetCtrlNodeDerivatives(CGeometry* geometry, CConfig* config, const bool ForwardProjectionDerivative){
  const SU2_COMPONENT Kind_SU2 = config->GetKind_SU2();
  unsigned long ctrlNodeIndex = 0;
  const bool dataRed = config->GetRBF_DataReduction();

  /*--- Early exit for unsupported features ---*/
  if (!(Kind_SU2 == SU2_COMPONENT::SU2_DOT && !ForwardProjectionDerivative)) {
    SU2_MPI::Error("Missing feature in RBF interpolation", CURRENT_FUNCTION);
  }

  CtrlNodeDeformation.resize(nCtrlNodesLocal * nDim, 0.0);
  
  for (auto iType : CtrlTypeVec) {
    const auto& markerToNodeSet = dataRed ? ctrl_node_indices[iType] : node_indices[iType];

    for (const auto& iMarker : markerToNodeSet) {
      for (const auto iNode : iMarker.second) {
        const auto nodeIndex = nodes[iNode]->GetIndex();
        for (auto iDim = 0u; iDim < nDim; iDim++) {
          CtrlNodeDeformation[ctrlNodeIndex * nDim + iDim] = SU2_TYPE::GetValue(geometry->GetSensitivity(nodeIndex, iDim));
        }
        ctrlNodeIndex++;
      }
    }
  }
  
}

void CRadialBasisFunctionInterpolation::SetInternalNodes(CGeometry* geometry, CConfig* config, vector<unsigned long>& internalNodes) const { 

  /*--- helper function to check if node is in the nodes vector ---*/
  auto is_in_nodes = [&](unsigned long iNode) {
        return find_if(nodes.begin(), nodes.end(),
            [&](CRadialBasisFunctionNode* n) { return n->GetIndex() == iNode; }) != nodes.end();
    };

  /*--- Add all non-boundary nodes and periodic nodes not already in nodes vector ---*/
  for (auto iNode = 0ul; iNode < geometry->GetnPoint(); iNode++) {    
    if (!geometry->nodes->GetBoundary(iNode)) {
      internalNodes.push_back(iNode);
    } else if (geometry->nodes->GetPeriodicBoundary(iNode) && !is_in_nodes(iNode)) {      
      internalNodes.push_back(iNode);
    }
  }  
  
  /*--- In case of a parallel computation, the nodes on the send/receive markers are included as internal nodes
          if they are not already a boundary node with known deformation ---*/
  #ifdef HAVE_MPI
    for (auto iMarker = 0u; iMarker < geometry->GetnMarker(); iMarker++) { 
      if (config->GetMarker_All_SendRecv(iMarker)) { 
        for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) { 
          const auto iNode = geometry->vertex[iMarker][iVertex]->GetNode();

          /*--- if not among the boundary nodes ---*/
          if (!is_in_nodes(iNode)) {
            internalNodes.push_back(iNode);
          }             
        }
      }
    }
  #endif

  /*--- sorting of the local indices and obtain unique set ---*/
  sort(internalNodes.begin(), internalNodes.end());
  internalNodes.resize(std::distance(internalNodes.begin(), unique(internalNodes.begin(), internalNodes.end())));

  // ofstream out4("int_nodes"+to_string(rank)+".txt");
  // for (auto i : internalNodes){
  //   out4 << geometry->nodes->GetGlobalIndex(i) << /*"\t" <<  i->getNodetype()  <<*/  endl;
  // }
  // out4.close();
}

void CRadialBasisFunctionInterpolation::ComputeInterpolationCoefficients(const su2passivematrix& invInterpMat) {
  /*--- The interpolation coefficients are found from the resulting matrix vector multiplication: Coeff_control = Phi_inv * Def_control.
          Which for a single coefficient is a summation from j=1 to nCtrlNodesGlobal: Coeff_i = sum( InterpMat_i,j * d_j ).
          The deformation vector is scattered accros the ranks and therefore each rank computes its contribution,
          and finally these contributions are summed. ---*/

  vector<su2double> localInterpCoeffSum(nDim * nCtrlNodesGlobal, 0.0);

  /*--- Distribute local control node sizes and compute rank specific starting index ---*/
  vector<unsigned long> localControlNodeSizes(size);
  SU2_MPI::Allgather(&nCtrlNodesLocal, 1, MPI_UNSIGNED_LONG, localControlNodeSizes.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm()); 
  
  unsigned long rankStartIndex = 0;
  for (auto iProc = 0; iProc < rank; iProc++){
    rankStartIndex += localControlNodeSizes[iProc];
  }

  /*--- Obtain local contributions ---*/
  for (auto iNode = 0ul; iNode < nCtrlNodesLocal; iNode++) {
    const auto iNodeOffset = iNode * nDim;
    const auto column = rankStartIndex+iNode;
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      const su2double deformation = CtrlNodeDeformation[iNodeOffset+iDim];
      for (auto jNode = 0ul; jNode < nCtrlNodesGlobal; jNode++ ) {
        localInterpCoeffSum[jNode*nDim + iDim] += invInterpMat(jNode, column) * deformation;
      }
    }
  }
  
  /*--- Sum local contributions to obtain interpolation coefficients ---*/
  InterpCoeff.resize(nDim * nCtrlNodesGlobal);
  #ifdef HAVE_MPI
    SU2_MPI::Allreduce(localInterpCoeffSum.data(), InterpCoeff.data(), localInterpCoeffSum.size(), MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  #else
    InterpCoeff = move(localInterpCoeffSum);
  #endif
}

void CRadialBasisFunctionInterpolation::UpdateGridCoord(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, const vector<unsigned long>& internalNodes, bool surfaceCorrection){

  if(rank == MASTER_NODE){
    cout << "\nupdating the grid coordinates" << endl;
  }
  
  /*--- Update of internal node coordinates ---*/
  UpdateInternalCoords(geometry, type, radius, internalNodes);
  
  /*--- Update of boundary node coordinates ---*/
  UpdateBoundCoords(geometry, config, type, radius);   

  /*--- In case of data reduction, perform the correction for nonzero error nodes ---*/
  if(config->GetRBF_DataReduction() && nodes.size() -  nCtrlNodesLocal > 0 ){
    if (surfaceCorrection) {
      SetCorrectionSurface(geometry, config, type);
    }
    SetCorrection(geometry, config, type, internalNodes); 
  }
}

void CRadialBasisFunctionInterpolation::UpdateGridCoord_Derivatives(CGeometry* geometry, CConfig* config, bool ForwardProjectionDerivative){
  
  
  SU2_COMPONENT Kind_SU2 = config->GetKind_SU2();

  if ((Kind_SU2 == SU2_COMPONENT::SU2_DOT) && !ForwardProjectionDerivative) {
  
    unsigned long nCtrlNode = nCtrlNodesLocal;
    unsigned long nCtrlNodes[size];
    SU2_MPI::Allgather(&nCtrlNode, 1, MPI_UNSIGNED_LONG, nCtrlNodes, 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm()); 

    unsigned long start_idx = 0;
    for (auto iProc = 0; iProc < rank; iProc++){
      start_idx += nCtrlNodes[iProc]; 
    }
    
    unsigned long idx = start_idx;
    for (auto type : CtrlTypeVec) {
      // auto ctrl_idx = config->GetRBF_DataReduction() ? GetIndices(ctrl_nodes_type, type) : GetIndices(node_indices, type);
      auto test = config->GetRBF_DataReduction()
      ? ctrl_node_indices[type]
      : node_indices[type];
    
      for (auto mark : test){
        for (auto idx_i : mark.second) {
          auto iMarker = nodes[idx_i]->GetMarker();

          if (config->GetSolid_Wall(iMarker) || (config->GetMarker_All_DV(iMarker) == YES)) {
            auto iPoint = nodes[idx_i]->GetIndex();
            for (auto iDim = 0u; iDim < nDim; iDim++) {
              // summation of current sensitivity and the computed update
              su2double sens_new =  geometry->GetSensitivity(iPoint, iDim) + sensitivityUpdate[idx * nDim + iDim];
              geometry->SetSensitivity(iPoint, iDim, sens_new);
            }
          }
          idx++;
        }
      }
    }
  } else {
    SU2_MPI::Error("Missing feature in RBF interpolation", CURRENT_FUNCTION);
  }

}


void CRadialBasisFunctionInterpolation::UpdateInternalCoords(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius, const vector<unsigned long>& internalNodes){
   /*--- Vector for storing the coordinate variation ---*/
  su2double var_coord[nDim]{0.0};
  
  /*--- Loop over the internal nodes ---*/
  for(auto iNode = 0ul; iNode < internalNodes.size(); iNode++){
    auto targetCoords = geometry->nodes->GetCoord(internalNodes[iNode]);
    /*--- Loop for contribution of each control node ---*/
    for(auto jNode = 0ul; jNode < nCtrlNodesGlobal; jNode++){

      /*--- Determine distance between considered internal and control node ---*/
      su2double dist = ComputeDistance(CtrlCoords[jNode*nDim], targetCoords);

      /*--- Evaluate RBF based on distance ---*/
      su2double rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));
    
      /*--- Add contribution to total coordinate variation ---*/
      for(auto iDim = 0u; iDim < nDim; iDim++){
        var_coord[iDim] += rbf*InterpCoeff[jNode * nDim + iDim];
      }
    }
    
    if (IsCylindrical) CylDispToCart(targetCoords, var_coord);
    
    /*--- Apply the coordinate variation and resetting the var_coord vector to zero ---*/
    for(auto iDim = 0u; iDim < nDim; iDim++){
      geometry->nodes->AddCoord(internalNodes[iNode], iDim, var_coord[iDim]);
      var_coord[iDim] = 0;
    } 
    
  }  
}

su2double CRadialBasisFunctionInterpolation::ComputeDistance(const su2double* ctrlCoords, const su2double* targetCoords) {
  su2double dist;
  if (IsCylindrical) {
    su2double cyl_coords[nDim];
    CartToCyl(targetCoords, cyl_coords);  
    dist = GetDistance(ctrlCoords, cyl_coords);
  } else {
    dist = GetDistance(ctrlCoords, targetCoords);
  }
  return dist;
}

void CRadialBasisFunctionInterpolation::UpdateBoundCoords(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius){
  
  /*--- Vector for storing the coordinate variation ---*/
  su2double var_coord[nDim]{0.0};
  
  /*--- In case of data reduction, the non-control boundary nodes are treated as if they where internal nodes ---*/
  if(config->GetRBF_DataReduction()){
    
    for (auto nodetype : CtrlTypeVec) {
      /*--- Looping over the non selected boundary nodes ---*/
      
      auto markers = node_indices[nodetype];
      for (auto mark : markers){

        for (auto iNode : mark.second) {

          if (!nodes[iNode]->GetControl()){
          
            auto targetCoords = geometry->nodes->GetCoord(nodes[iNode]->GetIndex());

            /*--- Finding contribution of each control node ---*/
            for( auto jNode = 0ul; jNode <  nCtrlNodesGlobal; jNode++){
              
              /*--- Distance of non-selected boundary node to control node ---*/
              su2double dist = ComputeDistance(CtrlCoords[jNode*nDim], targetCoords);

              /*--- Evaluation of the radial basis function based on the distance ---*/
              su2double rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));

              /*--- Computing and add the resulting coordinate variation ---*/
              for(auto iDim = 0u; iDim < nDim; iDim++){
                var_coord[iDim] += rbf*InterpCoeff[jNode * nDim + iDim];
              }
            }

          
            if (IsCylindrical) {
              CylDispToCart(geometry->nodes->GetCoord(nodes[iNode]->GetIndex()), var_coord);
            }
            
            /*--- Applying the coordinate variation and resetting the var_coord vector*/
            for(auto iDim = 0u; iDim < nDim; iDim++){
              geometry->nodes->AddCoord(nodes[iNode]->GetIndex(), iDim, var_coord[iDim]);
              var_coord[iDim] = 0;
            }
          
          }
        }
      }
    }
  }

  const su2double VarIncrement = 1.0 / ((su2double)config->GetGridDef_Nonlinear_Iter());

  unsigned long idx = 0;
  for (auto type : CtrlTypeVec) {
    auto markers = config->GetRBF_DataReduction()
    ? ctrl_node_indices[type]
    : node_indices[type];

    for (const auto& iMarker : markers) {
      auto ctrl_idx = iMarker.second;
      for (auto idx_i : ctrl_idx) {
        auto varCoord = geometry->vertex[nodes[idx_i]->GetMarker()][nodes[idx_i]->GetVertex()]->GetVarCoord();

        for (auto iDim = 0u; iDim < nDim; iDim++) {
          geometry->nodes->AddCoord(nodes[idx_i]->GetIndex(), iDim, varCoord[iDim] * VarIncrement);   
        }
      }
    }
  }
}


void CRadialBasisFunctionInterpolation::SetCtrlNodeCoords(CGeometry* geometry, CConfig* config){
  unsigned long index = 0;
  unsigned long localCoordsSize = nCtrlNodesLocal * nDim;
  vector<su2double> localCoords(localCoordsSize, 0.0);  
  const bool dataRed = config->GetRBF_DataReduction();

  /*--- Gathering local control node coordinates, cylindrical if necessary. ---*/  
  for (const auto& iType : CtrlTypeVec) {
    const auto& markerToNodeSet = dataRed ? ctrl_node_indices[iType] : node_indices[iType];
    for (const auto& iMarker  : markerToNodeSet) {
      for (const auto iNode : iMarker.second) {    
        const su2double* coord = IsCylindrical ? nodes[iNode]->GetCylCoord() : geometry->nodes->GetCoord(nodes[iNode]->GetIndex());
        for (auto iDim = 0u; iDim < nDim; iDim++) {
          localCoords[index++] = coord[iDim];
        }
      }
    }
  }

  /*--- Distributing global control node coordinates among all processes ---*/
  CtrlCoords.resize(nCtrlNodesGlobal*nDim);

  vector<int> localCoordsSizes(size);
  SU2_MPI::Allgather(&localCoordsSize, 1, MPI_INT, localCoordsSizes.data(), 1, MPI_INT, SU2_MPI::GetComm()); 

  vector<int> disps(size, 0);
  for(auto iProc = 1; iProc < SU2_MPI::GetSize(); iProc++){
    disps[iProc] = disps[iProc-1]+localCoordsSizes[iProc-1];
  }  
  
  SU2_MPI::Allgatherv(localCoords.data(), localCoordsSize, MPI_DOUBLE, CtrlCoords.data(), localCoordsSizes.data(), disps.data(), MPI_DOUBLE, SU2_MPI::GetComm()); 
};


void CRadialBasisFunctionInterpolation::GetInterpError(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, bool Derivative, const vector<NODETYPE>& slideCtrlTypes,  unsigned long& maxErrorNodeLocal, su2double& maxErrorLocal){
  if (!config->GetRBF_DataReduction()) {
    maxErrorLocal = 0;
    return;
  }
  
  /*--- Array containing the local error ---*/
  su2double localError[nDim];

  /*--- Magnitude of the local maximum error ---*/
  maxErrorLocal = 0.0;
  // cout << "errors: " << endl;
  /*--- Loop over non-selected boundary nodes ---*/
  NODETYPE errortype; 
  auto nodes_disp = node_indices[NODETYPE::DISPLACED];
  for (auto mark : nodes_disp) {
    // if (rank == MASTER_NODE) cout << rank << " " << mark.first << " " << mark.second.size() << endl;
    for (auto idx : mark.second) {
      if(!nodes[idx]->GetControl()){

        /*--- Compute nodal error ---*/
        GetNodalError(geometry, config, type, radius, nodes[idx], Derivative, localError);
        // cout << localError[0] << " " << localError[1] << endl;
        /*--- Setting error ---*/
        nodes[idx]->SetError(localError, nDim);
        
        /*--- Compute error magnitude and update local maximum error if necessary ---*/
        su2double errorMagnitude = GeometryToolbox::Norm(nDim, localError);
        if(errorMagnitude > maxErrorLocal){
          maxErrorLocal = errorMagnitude;
          maxErrorNodeLocal = idx;
          errortype = NODETYPE::DISPLACED;
        }
      }
    }
  }

  
  for (auto i_type : slideCtrlTypes) {
    
    auto nodes_edge = node_indices[i_type];

    for (const auto& mark : nodes_edge) {
      
      auto markerGlobal = mark.first;
      const auto& markerLocal = config->GetMarker_Local(markerGlobal);

      const auto targetNodes = mark.second;

      const auto& targetNodesADT = node_indices[i_type][markerGlobal];
      const bool isPeriodic = config->GetnMarker_Periodic() != 0;
      auto BoundADT = CreateADT(geometry, targetNodesADT, markerLocal, isPeriodic);
      
      for (auto x : targetNodes) {
        const auto iNode = nodes[x];
        ApplyRBF(geometry, type, radius, iNode);
      }

      // function that stores the nearest node id and rank id of the rbf nodes; 
      for (auto x : targetNodes) {
        const auto iNode = nodes[x];
        GetNearestNode(BoundADT.get(), iNode);
      }
      
      auto response_recv_buffer = ExchangeNearestNodeData(geometry, config, markerLocal, targetNodes, nodes, i_type);
      SetNearestNodeData(geometry, config, response_recv_buffer, targetNodes, nodes);
      

      //now the normals and the nearest node coords are available -> do projection        
      for (auto x : targetNodes) {
        auto* iNode = nodes[x];
        su2double projection[3] = {0.0, 0.0, 0.0};
        ApplyProjection(geometry, markerLocal, iNode, projection);

        if(IsCylindrical){
          su2double temp[nDim];
          auto new_coord = iNode->GetNewCoord();
          CylToCart(new_coord, temp);
          CylDispToCart(temp, projection);
        }

        iNode->SetError(projection, nDim);
      }

      const auto& indices = mark.second;  
      for (const auto idx : indices) {
        if(!nodes[idx]->GetControl()){
          auto localError = nodes[idx]->GetError();

          su2double errorMagnitude = GeometryToolbox::Norm(nDim, localError);
          if(errorMagnitude > maxErrorLocal){
            maxErrorLocal = errorMagnitude;
            maxErrorNodeLocal = idx;
            errortype = nodes[idx]->GetNodeType();
          }
        }
      }
    }
  }

  SU2_MPI::Allreduce(&maxErrorLocal, &MaxErrorGlobal, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
}

void CRadialBasisFunctionInterpolation::GetNodalError(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, CRadialBasisFunctionNode* iNode, bool Derivative, su2double* localError){ 
  
  /*--- If requested (no by default) impose the surface deflections in increments ---*/
  const su2double VarIncrement = 1.0 / ((su2double)config->GetGridDef_Nonlinear_Iter());
  
  /*--- If node is part of a moving boundary then the error is defined as the difference
           between the found and prescribed displacements. Thus, here the displacement is substracted from the error ---*/

  // cout << "node: " << iNode->GetIndex() << endl;
  if (Derivative) {
    for (auto iDim =0u; iDim < nDim; iDim++) {
      localError[iDim] = -geometry->GetSensitivity(iNode->GetIndex(), iDim) * VarIncrement;
    }
    // cout << localError[0] << " " << localError[1] << endl;
  }else{  
    auto disp_true = geometry->vertex[iNode->GetMarker()][iNode->GetVertex()]->GetVarCoord();
    for(auto iDim = 0u; iDim < nDim; iDim++){
      localError[iDim] = -disp_true[iDim] * VarIncrement;
    }
  }


  su2double disp_interp[nDim] = {0.0};
  /*--- Resulting displacement from the RBF interpolation is added to the error ---*/ 

  auto targetCoord = geometry->nodes->GetCoord(iNode->GetIndex());
  // cout << "target coord: " << targetCoord[0] << " " << targetCoord[1] << endl;
  /*--- Finding contribution of each control node ---*/
  for(auto jNode = 0ul; jNode < nCtrlNodesGlobal; jNode++){

    /*--- Distance between non-selected boundary node and control node ---*/
    su2double dist = ComputeDistance(CtrlCoords[jNode*nDim], targetCoord);

    /*--- Evaluation of Radial Basis Function ---*/
    su2double rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));

    /*--- Add contribution to error ---*/
    for(auto iDim = 0u; iDim < nDim; iDim++){
      // localError[iDim] += rbf*InterpCoeff[jNode*nDim + iDim];
      disp_interp[iDim] += rbf*InterpCoeff[jNode*nDim + iDim];
    }
  }
  // cout << "disp_interp: " << disp_interp[0] << " " << disp_interp[1] << endl;

  if (IsCylindrical) {
    CylDispToCart(geometry->nodes->GetCoord(iNode->GetIndex()),  disp_interp);
  }

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    localError[iDim] += disp_interp[iDim];
  }
  // cout << "final error: " << localError[0] << " " << localError[1] << endl;
}

void CRadialBasisFunctionInterpolation::SetCorrection(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const vector<unsigned long>& internalNodes){
  /*--- The non-selected control nodes still have a nonzero error once the maximum error falls below the data reduction tolerance. 
          This error is applied as correction and interpolated into the volumetric mesh for internal nodes that fall within the correction radius.
          To evaluate whether an internal node falls within the correction radius an AD tree is constructed of the boundary nodes,
          making it possible to determine the distance to the nearest boundary node. 
          The selected control nodes are included in the AD tree to avoid problems with hight aspect ratio cells near the considered geometry. ---*/
  
  
  /*--- Correction Radius, equal to maximum error times a prescribed constant ---*/
  const su2double CorrectionRadius = config->GetRBF_DataRedCorrectionFactor()*MaxErrorGlobal;


  /*--- Construction of the AD tree consisting of the selected boundary nodes ---*/

  /*--- Number of selected boundary nodes ---*/
  const auto nVertexNonPer = nodes.size();
  const auto nVertexPer   = per_nodes.size();
  auto nVertexBound = nVertexNonPer + nVertexPer;

  /*--- Vector storing the coordinates of the boundary nodes ---*/
  vector<su2double> Coord_bound(nDim*nVertexBound);

  /*--- Vector storing the IDs of the boundary nodes ---*/
  vector<unsigned long> PointIDs(nVertexBound);   

  /*--- Storing boundary node information ---*/

  for(auto iVertex = 0ul; iVertex < nVertexNonPer; iVertex++){
      PointIDs[iVertex] = iVertex;
      auto iNode = nodes[iVertex]->GetIndex();
      auto coord = geometry->nodes->GetCoord(iNode);
      for(auto iDim = 0u; iDim < nDim; iDim++){
        Coord_bound[iVertex*nDim + iDim] = coord[iDim]; 
      }
  }


  for (auto iVertex = 0ul; iVertex < nVertexPer; iVertex++){
    PointIDs[iVertex + nVertexNonPer] = nVertexNonPer + iVertex;
    auto iNode = per_nodes[iVertex]->GetIndex();
    auto coord = geometry->nodes->GetCoord(iNode);
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      Coord_bound[(iVertex + nVertexNonPer) * nDim + iDim] = coord[iDim];
    }
  }
  
  
  /*--- Construction of AD tree ---*/
  CADTPointsOnlyClass BoundADT(nDim, nVertexBound, Coord_bound.data(), PointIDs.data(), false);

  /*--- ID of nearest boundary node ---*/
  unsigned long pointID;
  /*--- Distance to nearest boundary node ---*/
  su2double dist;
  /*--- rank of nearest boundary node ---*/
  int rankID;

  /*--- Interpolation of the correction to the internal nodes that fall within the correction radius ---*/
  for(auto iNode = 0ul; iNode < internalNodes.size(); iNode++){
    
    /*--- Find nearest node ---*/
    BoundADT.DetermineNearestNode(geometry->nodes->GetCoord(internalNodes[iNode]), dist, pointID, rankID);  

    /*--- evaluate RBF ---*/
    su2double rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, CorrectionRadius, dist));
    
    /*--- Get error of nearest node ---*/
    auto err = (pointID < nodes.size()) ? nodes[pointID]->GetError() : per_nodes[pointID-nodes.size()]->GetError();
    
    /*--- Apply correction to the internal node ---*/
    for(auto iDim = 0u; iDim < nDim; iDim++){      
      geometry->nodes->AddCoord(internalNodes[iNode], iDim, -rbf*err[iDim]); 
    }
  }


  /*--- Applying the correction to the non-selected boundary nodes ---*/
  for(auto iNode : nodes){
    if (!iNode->GetControl()) {
      auto err = iNode->GetError();
      for(auto iDim = 0u; iDim < nDim; iDim++){
        geometry->nodes->AddCoord(iNode->GetIndex(), iDim, -err[iDim]);
      }
    }
  }
}

void CRadialBasisFunctionInterpolation::SetCorrectionSurface(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type){
  // required for updated normal
  geometry->SetBoundControlVolume(config, UPDATE);

  /*--- Correction Radius, equal to maximum error times a prescribed constant ---*/
  const su2double CorrectionRadius = config->GetRBF_DataRedCorrectionFactor()*MaxErrorGlobal;

  

  // consideration for the sliding surface nodes that would fall fdwithin the correction radius from a node with a prescribed deformation


  // making ADT of all the moving nodes
  unsigned long nDisplaced = 0;
  auto disp_map = node_indices[NODETYPE::DISPLACED];
  for (auto& pair : disp_map){
     nDisplaced += pair.second.size(); // conservative size;
  }

  vector<su2double> coord_moving(nDisplaced*nDim);
  vector<unsigned long> ID_moving(nDisplaced);

  unsigned long cnt = 0;
  for(auto iVertex = 0ul; iVertex < nodes.size(); iVertex++){
    if (nodes[iVertex]->GetNodeType() == NODETYPE::DISPLACED /*&& !nodes[iVertex]->GetControl()*/) {

      auto coord = geometry->nodes->GetCoord(nodes[iVertex]->GetIndex());
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        coord_moving[cnt*nDim+iDim] = coord[iDim];
      }
      ID_moving[cnt] = iVertex;
      cnt++;
    } 
  }

  coord_moving.resize(nDim*cnt);
  ID_moving.resize(cnt);
  
  CADTPointsOnlyClass MovingADT(nDim, cnt, coord_moving.data(), ID_moving.data(), false);
  
  /*--- ID of nearest boundary node ---*/
  unsigned long pointID;
  /*--- Distance to nearest boundary node ---*/
  su2double dist;
  /*--- rank of nearest boundary node ---*/
  int rankID;



  auto surf_map = node_indices[NODETYPE::SURFACE];
  for (auto& pair : surf_map) {
    auto node_vec = pair.second;
    for (auto iNode : node_vec) {

      // for each surface node find the closest displaced node
      MovingADT.DetermineNearestNode(geometry->nodes->GetCoord(nodes[iNode]->GetIndex()), dist, pointID, rankID);  
      

      if (nodes[pointID]->GetControl()){
        continue;
      }

      /*--- evaluate RBF ---*/
      su2double rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, CorrectionRadius, dist));

      if (rbf > 0) {
        
        // obtain error of the displaced node
        auto err = nodes[pointID]->GetError();
       
        // coord and normal of the surface node
        auto coord = geometry->nodes->GetCoord(nodes[iNode]->GetIndex());
        auto normal = geometry->vertex[nodes[iNode]->GetMarker()][nodes[iNode]->GetVertex()]->GetNormal();
        
        // magnitude of the surface normal
        auto mag = GeometryToolbox::Norm(nDim, normal);
        for (auto iDim = 0u; iDim < nDim; iDim++){
          if (abs(normal[iDim])/mag < 0.017){
            normal[iDim] = 0;
          }
        }

        mag = GeometryToolbox::Norm(nDim, normal);
        for (auto iDim = 0u; iDim < nDim; iDim++){  
            normal[iDim] = normal[iDim]/mag;
        }
        
        // dot product of error and cylindrical normal
        auto dp = GeometryToolbox::DotProduct(nDim, err, normal);
        
        // add this error to the surface node as it will be applied later on
        for(auto iDim = 0u; iDim < nDim; iDim++){
          err[iDim] -= dp*normal[iDim];
          nodes[iNode]->AddError(err[iDim]*rbf, iDim);
        }  
      }
    }
  }

}

void CRadialBasisFunctionInterpolation::AddControlNode(CConfig* config, unsigned long maxErrorNode){ 
  auto* node = nodes[maxErrorNode];
  node->setControl();

  vector<su2double> zero(nDim, 0.0);
  node->SetError(zero.data(), nDim);

  auto marker = node->GetMarker(); // this is the local marker index
  const auto& tag = config->GetMarker_All_TagBound(marker);
  auto global_marker = config->GetMarker_CfgFile_TagBound(tag);
  ctrl_node_indices[node->GetNodeType()][global_marker].insert(maxErrorNode);
}


void CRadialBasisFunctionInterpolation::Get_nCtrlNodesGlobal(CConfig* config){
  nCtrlNodesLocal = 0;  

  /*--- Summation of control nodes size of all markers ---*/
  for (const auto& iType : CtrlTypeVec) {
    const auto& markerToNodeSet = config->GetRBF_DataReduction() ? ctrl_node_indices[iType] : node_indices[iType];
    for (const auto& iMarker : markerToNodeSet) {
      nCtrlNodesLocal += iMarker.second.size();
    } 
  } 

  SU2_MPI::Allreduce(&nCtrlNodesLocal, &nCtrlNodesGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
}


void CRadialBasisFunctionInterpolation::SetPeriodicVars(CConfig* config){
  
  /*--- Counter number of periodic angles encountered ---*/
  unsigned short rotationalAngleCnt = 0;
  
  /*--- Loop over periodic markers and find rotational periodic parameters ---*/
  for (auto iMarker = 0u; iMarker < config->GetnMarker_Periodic(); iMarker++) {

    /*--- Obtain periodic angle and periodic axis of rotation ---*/
    const su2double* periodicAngles = config->GetPeriodic_RotAngles(iMarker);
    for (auto iDim = 0u; iDim < 3; iDim ++) {

      if (periodicAngles[iDim] != 0.0) {

        /*--- For a 2D domain the rotational axis is perpendicular to the 2D plane (z-axis) ---*/
        if (nDim == 2 && iDim != Z_DIR) {
          SU2_MPI::Error("Rotational periodicity can only be applied around the third coordinate axis for a 2D domain.", CURRENT_FUNCTION);
        } 
        
        /*--- Store Rotational axis, periodic angle. 
                Present errors in case of multiple rotational axes or periodic angles. ---*/
        if (rotationalAngleCnt == 0) {
          RotationalAxis = iDim; 
          PeriodicAngle = fabs(periodicAngles[iDim]);
          IsCylindrical = true;
        } else if (RotationalAxis != iDim) {
          SU2_MPI::Error("Only a single rotationally periodic angle can be provided with MARKER_PERIODIC.", CURRENT_FUNCTION);
        } else if (fabs(periodicAngles[iDim]) != PeriodicAngle) {
          SU2_MPI::Error("Two different values of periodic angles detected in MARKER_PERIODIC.", CURRENT_FUNCTION);
        }
        rotationalAngleCnt++;
      }
    }
  }

  /*--- Loop over periodic markers to find translational periodic parameters ---*/
  for (auto iMarker = 0u; iMarker < config->GetnMarker_Periodic(); iMarker++) {

    /*--- Store periodic lengths ---*/
    const su2double* periodicTranslation = config->GetPeriodic_Translation(iMarker);
    for (auto iDim = 0u; iDim < 3; iDim++) {

        if (periodicTranslation[iDim] != 0.0) {

          /*--- For 2D, periodicity cannot be rotational and translational.
                  For 3D domain with rotational periodicity, the domain can only be translationally periodic along the rotational periodic axis. ---*/
          if (IsCylindrical) {
            if (nDim == 2) {
              SU2_MPI::Error("A 2D domain cannot have both rotational and translational periodicity simultaneously.", CURRENT_FUNCTION);
            } else if (RotationalAxis != iDim) {
              SU2_MPI::Error("Periodic translation can only be applied along the rotational axis.", CURRENT_FUNCTION);
            }
          }

          /*--- Save periodic axes and lengths ---*/
          PeriodicLength[iDim] = fabs(periodicTranslation[iDim]);
        }
    }
  }
}

void CRadialBasisFunctionInterpolation::TransformBoundaryNodesToCylindricalCoords(CGeometry* geometry){
  
  auto transform = [&](const vector<CRadialBasisFunctionNode*>& nodes) {
    for (auto* iNode : nodes) {
      su2double cyl_coord[nDim];
      CartToCyl(geometry->nodes->GetCoord(iNode->GetIndex()), cyl_coord);
      iNode->SetCylCoord(cyl_coord, nDim);
    }
  };

  transform(nodes);
  transform(per_nodes);
}

void CRadialBasisFunctionInterpolation::CartDispToCyl(const su2double* init_coord_cart, const su2double* var_coord_cart, su2double* var_coord_cyl) const {
  // new Carthesian coord
  su2double new_coord[3];
  for (auto iDim = 0u; iDim < nDim; iDim++){
    new_coord[iDim] = init_coord_cart[iDim] + var_coord_cart[iDim];
  }

  su2double init_coord_cyl[3];
  CartToCyl(init_coord_cart, init_coord_cyl);
  CartToCyl(new_coord, var_coord_cyl);

  for (auto iDim = 0u; iDim < nDim; iDim++){
    var_coord_cyl[iDim] -= init_coord_cyl[iDim]; 
  }
}

void CRadialBasisFunctionInterpolation::CylDispToCart(const su2double* init_coord_cart, su2double* var_coord) const {
  // transform initial coordinates to cylindrical
  su2double init_coord_cyl[nDim];
  CartToCyl(init_coord_cart, init_coord_cyl);

  // find new position in cylindrical coordinates
  su2double new_coord_cyl[nDim];
  for (auto iDim = 0u; iDim < nDim; iDim++){
    new_coord_cyl[iDim] = init_coord_cyl[iDim] + var_coord[iDim];
  }

  // transform new position into cartesian coordinates
  CylToCart(new_coord_cyl, var_coord);

  // subtract initial coordinate position to find the deltaf
  for (auto iDim = 0u; iDim < nDim; iDim++){
    var_coord[iDim] -= init_coord_cart[iDim]; 
  }
}


void CRadialBasisFunctionInterpolation::GetPeriodicNodeErrors(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, bool Derivative) {

  if (!config->GetRBF_DataReduction() || Derivative) return;

  const auto& markers = per_node_indices[NODETYPE::DISPLACED];
  for (const auto& mark : markers) {
      for (auto idx : mark.second) {
          su2double localError[nDim];
          GetNodalError(geometry, config, type, radius, per_nodes[idx], Derivative, localError);
          per_nodes[idx]->SetError(localError, nDim);
      }
  }

  for (auto i = 1u; i < CtrlTypeVec.size(); i++) {
    /*// TODO -  const*/ auto& iType = CtrlTypeVec[i];
    const auto& markers = per_node_indices[iType];
    for (const auto& mark : markers) {

      auto markerGlobal = mark.first;
      const auto& markerLocal = config->GetMarker_Local(markerGlobal);

      const auto targetNodes = mark.second;

      const auto& targetNodesADT = node_indices[iType][markerGlobal];
      auto BoundADT = CreateADT(geometry, targetNodesADT, markerLocal, config->GetnMarker_Periodic() != 0);
              


      for (auto x : targetNodes) {
        const auto iNode = per_nodes[x];
        ApplyRBF(geometry, type, radius, iNode);
      }

      // function that stores the nearest node id and rank id of the rbf nodes; 
      for (auto x : targetNodes) {
        const auto iNode = per_nodes[x];
        GetNearestNode(BoundADT.get(), iNode);
      }
      
      auto response_recv_buffer = ExchangeNearestNodeData(geometry, config, markerLocal, targetNodes, per_nodes, iType);
      SetNearestNodeData(geometry, config, response_recv_buffer, targetNodes, per_nodes);
      

      //now the normals and the nearest node coords are available -> do projection        
      for (auto x : targetNodes) {
        const auto iNode = per_nodes[x];
        su2double projection[3] = {0.0, 0.0, 0.0};
        ApplyProjection(geometry, markerLocal, iNode, projection);

        if(IsCylindrical){
          su2double temp[nDim];
          auto new_coord = iNode->GetNewCoord();
          CylToCart(new_coord, temp);
          CylDispToCart(temp, projection);
        }

        iNode->SetError(projection, nDim);
      }
    }
  }
}


unique_ptr<CADTPointsOnlyClass> CRadialBasisFunctionInterpolation::CreateADT(CGeometry* geometry, const unordered_set<unsigned long>& targetNodes, const short markerLocal, bool isPeriodic) const {
  const unsigned long size = targetNodes.size();

  /*--- Early exit in case the global marker is not part of the local domain or there are no target nodes. ---*/
  if (markerLocal == -1 || size == 0) {
    return unique_ptr<CADTPointsOnlyClass>( new CADTPointsOnlyClass(nDim, 0ul, nullptr, nullptr, true));
  }
  
  /*--- In case of periodicity, 3 images are included: original + two periodic images. ---*/
  const int nImages = isPeriodic ? 3 : 1; 

  vector<su2double> CoordBound(size * nDim * nImages);
  vector<unsigned long> PointIDs(size * nImages);
  const auto nVertex = geometry->GetnVertex(markerLocal); 
  
  /*--- Obtaining the periodic translations ---*/
  const su2double angular_periodic_shift[3] = {0.0, +PeriodicAngle, -PeriodicAngle};  
  su2double translation_periodic_shift[3][3];
  for (auto iDim = 0u; iDim < nDim; iDim++) {
    translation_periodic_shift[0][iDim] = 0.0;
    translation_periodic_shift[1][iDim] = +PeriodicLength[iDim];
    translation_periodic_shift[2][iDim] = -PeriodicLength[iDim]; 
  }

  /*--- Add coordinates and point IDs to respective vectors. In case of periodicity multiple images are added. 
          They are indexed as: 0 - original image, 1 - positive periodic shift, 2 - negative periodic shift. 
          The vertex and image are coded into the PointID as: vertex + image * nVertex. 
          They can later be retrieved as remainder and quotient: vertex = PointID % nVertex, image = PointID / nVertex. ---*/
  unsigned long nodeCounter = 0;
  for (const auto iNode : targetNodes) {
    const auto* node = nodes[iNode];
    for (auto jImage = 0u; jImage < nImages; jImage++) {
      PointIDs[ nodeCounter + size * jImage] = node->GetVertex() + jImage * nVertex;
    }
    
    const su2double* coord = IsCylindrical ? node->GetCylCoord() : geometry->nodes->GetCoord(node->GetIndex());
    
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      for (auto jImage = 0; jImage < nImages; jImage++) {        
        // NOTE for angular periodicity the rotational axis is always the theta/second coordinate
        CoordBound[(nodeCounter + jImage * size)*nDim + iDim] = coord[iDim] + translation_periodic_shift[jImage][iDim] + ((iDim == 1) ? angular_periodic_shift[jImage] : 0.0); 
      }
    }
    nodeCounter++;
  }

  /*--- Construction of AD tree ---*/
  return unique_ptr<CADTPointsOnlyClass>( new CADTPointsOnlyClass(nDim, size*nImages, CoordBound.data(), PointIDs.data(), true));
}

void CRadialBasisFunctionInterpolation::ApplyRBF(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius, CRadialBasisFunctionNode* const node) {
  // TODO -  add comments
  const su2double* coord = IsCylindrical ? node->GetCylCoord() : geometry->nodes->GetCoord(node->GetIndex()); 

  // copy coord to new_coord
  node->SetNewCoord(coord, nDim);

  // find contribution of each control node 
  for(auto jNode = 0ul; jNode < nCtrlNodesGlobal; jNode++){
    
    const unsigned long offset = jNode*nDim;

    /*--- Determine distance between considered internal and control node ---*/
    const su2double dist = GetDistance(CtrlCoords[offset], coord);
    
    /*--- Evaluate RBF based on distance ---*/
    const auto rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));

    /*--- Add contribution to new coordinates -- -*/
    const su2double* coeff = &InterpCoeff[offset];
    for(auto iDim = 0u; iDim < nDim; iDim++){
      node->AddNewCoord(rbf * coeff[iDim], iDim);
    } 
  }  
}

void CRadialBasisFunctionInterpolation::GetNearestNode(CADTPointsOnlyClass* const BoundADT, CRadialBasisFunctionNode* const node) {
  // TODO -   add comments
  const su2double* coord = node->GetNewCoord();
  
  su2double dist;
  unsigned long near_point_id;
  int near_rank_id; 
  
  BoundADT->DetermineNearestNode(coord, dist, near_point_id, near_rank_id);
  node->SetNearestNode(near_point_id);
  node->SetNearestRank(near_rank_id);
}


CVertex* CRadialBasisFunctionInterpolation::GetNearestVertex(CGeometry* geometry, const unsigned long pointID, const unsigned short marker) const {

  const auto nVertex = geometry->GetnVertex(marker);
  const auto nearestVertexIndex = pointID % nVertex;

  return geometry->vertex[marker][nearestVertexIndex];
}

void CRadialBasisFunctionInterpolation::GetNearestCoord( CVertex* const nearestVertex, su2double* const nearestCoord ) const {
  const su2double* cartCoord = nearestVertex->GetCoord();
  if (IsCylindrical) {
      CartToCyl(cartCoord, nearestCoord);
  } else {
      copy_n(cartCoord, nDim, nearestCoord);
  }
}

void CRadialBasisFunctionInterpolation::SetLocalNormal(CGeometry* geometry, CConfig* config, CRadialBasisFunctionNode* const node) {
  // TODO -   add comments
    const auto marker = node->GetMarker();
    const auto pointID = node->GetNearestNode();
    // corresponding vertex
    auto nearestVertex = GetNearestVertex(geometry, pointID, marker);
         
    // piece of code to obtain the nearestCoord
    su2double nearestCoord[3] = {0.0, 0.0, 0.0};
    GetNearestCoord(nearestVertex, nearestCoord);

    su2double normal[3] = {0.0,0.0,0.0};
    su2double secondaryNormal[3] = {0.0,0.0,0.0};

    const auto nodetype = node->GetNodeType();
    GetVertexNormals(geometry, config, marker, nodetype, nearestVertex, normal, secondaryNormal);

    if (nDim == 3 && nodetype == NODETYPE::EDGE){
      node->SetNearestNormal2(secondaryNormal, nDim);
    }

    node->SetNearestNormal(normal, nDim);
    node->SetNearestCoord(nearestCoord, nDim);  
}

void CRadialBasisFunctionInterpolation::ApplyProjection(CGeometry* geometry, const short marker, const CRadialBasisFunctionNode* const node, su2double* const projection) const {

  const auto nVertex = geometry->GetnVertex(marker);
  const auto pointID = node->GetNearestNode();
  const auto nearestVertexIndex = pointID%nVertex; 
  const auto nearestImg = pointID / nVertex;

  const su2double* nearestCoord = node->GetNearestCoord();
  const su2double* normal = node->GetNearestNormal();
  const su2double* newCoord = node->GetNewCoord();

  
  su2double periodicCoord[3] = {0.0, 0.0, 0.0};
  
  const su2double shift = (nearestImg == 1) ? + 1.0 : (nearestImg == 2) ? -1.0 : 0.0;

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    periodicCoord[iDim] = nearestCoord[iDim] + shift * PeriodicLength[iDim];
  }
  
  if (nearestVertexIndex >= nVertex) { // TODO -   can I do the check when storing these values, instead of where I am using them?
    SU2_MPI::Error("Error: nearest node index is out of marker bounds", CURRENT_FUNCTION);
  }

  su2double distanceVector[3];
  GeometryToolbox::Distance(nDim, newCoord, periodicCoord, distanceVector);

  AddProjectionComponent(distanceVector, normal, projection);

  if(nDim == 3 && node->GetNodeType() == NODETYPE::EDGE) {
    const su2double* secNormal = node->GetNearestNormal2();

    AddProjectionComponent(distanceVector, secNormal, projection);
  }  
}

void CRadialBasisFunctionInterpolation::AddProjectionComponent(const su2double* const distanceVector, const su2double* const normal, su2double* const projection) const {
  // get dot product
  const su2double normalComponent = GeometryToolbox::DotProduct(nDim, normal, distanceVector);

  const su2double normalMagnitudeSq = GeometryToolbox::SquaredNorm(nDim, normal);
  const su2double normalizedNormalComponent = normalComponent / normalMagnitudeSq;
  // get projection
  for(auto iDim = 0u; iDim < nDim; iDim++){ 
    projection[iDim] += normal[iDim] * normalizedNormalComponent;
  }
}

vector<su2double> CRadialBasisFunctionInterpolation::ExchangeNearestNodeData(CGeometry* geometry, CConfig* config, const unsigned short marker, const unordered_set<unsigned long>& targetSet, const vector<CRadialBasisFunctionNode*>& targetNodes, const NODETYPE nodetype) const {

  // make vectors based on the rank that is nearest in order to send/receive blocks of information from one rank to other.

  // Maps each MPI rank to a list of nodes whose nearest neighbor resides on that rank
  unordered_map<int, vector<unsigned long>> queryMap;

  // loop through the local control nodes to fill the query map
  // for (const auto* iNode : targetNodes) {
  for (const auto index : targetSet) {
    const auto* iNode = targetNodes[index];

    // if nearest rank is own rank then set the local normals
    // else add to query map
    if (iNode->GetNearestRank() != rank) {
      queryMap[iNode->GetNearestRank()].push_back(iNode->GetNearestNode());
    }
  }

  // finding the number of send and receive counts
  vector<int> sendCounts(size, 0), recvCounts(size, 0);
  for (const auto& rankEntry : queryMap) {
    const auto& destRank = rankEntry.first;
    const auto& queryNodes  = rankEntry.second;
    sendCounts[destRank] = queryNodes.size();
  }
  // sending to each rank the nodes for which data is required to be send back
  SU2_MPI::Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, SU2_MPI::GetComm());

  // setting up data structures for sending the query nodes to their residing rank
  vector<unsigned long> sendBuffer, recvBuffer;
  vector<int> sendDisps(size), recvDisps(size);

  auto totalRecv = accumulate(recvCounts.begin(), recvCounts.end(), 0);
  auto totalSend = accumulate(sendCounts.begin(), sendCounts.end(), 0);

  recvBuffer.resize(totalRecv);
  sendBuffer.reserve(totalSend);

  // iterate over the ranks
  for ( auto iRank = 0, disp = 0; iRank < size; iRank++) {

    // set starting index of the rank
    sendDisps[iRank] = disp;

    // if iRank exists in map
    if (queryMap.count(iRank)) {
      // add querynodes at end of sendBuffer vector
      const auto& queryNodes = queryMap[iRank];
      sendBuffer.insert(sendBuffer.end(), queryNodes.begin(), queryNodes.end());

      // update starting index
      disp += queryNodes.size();
    }
  }

  // Get starting indices of the receiving data
  recvDisps[0] = 0;
  for (auto iRank = 1; iRank < size; iRank++) {
    recvDisps[iRank] = recvDisps[iRank - 1] + recvCounts[iRank - 1];
  }

  // sending the query data to all ranks
  SU2_MPI::Alltoallv(sendBuffer.data(), sendCounts.data(), sendDisps.data(), MPI_UNSIGNED_LONG,
                recvBuffer.data(), recvCounts.data(), recvDisps.data(), MPI_UNSIGNED_LONG,
                SU2_MPI::GetComm());


  // sending information back
  vector<su2double> responseSendBuffer;


  // iterating over the received query nodes
  int dataSize = (nDim == 3 && nodetype == NODETYPE::EDGE) ? 3 : 2;
  responseSendBuffer.resize(recvBuffer.size()*dataSize*nDim);

  auto cnt = 0;
  for (const auto x : recvBuffer) {
      // corresponding vertex
      const auto nearestVertex = GetNearestVertex(geometry, x, marker);
            
      // piece of code to obtain the nearestCoord
      su2double nearestCoord[3] = {0.0, 0.0, 0.0}, normal[3] = {0.0, 0.0, 0.0}, secondaryNormal[3] = {0.0, 0.0, 0.0};
  
      GetNearestCoord(nearestVertex, nearestCoord);

      GetVertexNormals(geometry, config, marker, nodetype, nearestVertex, normal, secondaryNormal);
      
    for (int i = 0; i < nDim; ++i)
        responseSendBuffer[cnt++] = normal[i];

    for (int i = 0; i < nDim; ++i)
        responseSendBuffer[cnt++] = nearestCoord[i];
      
    if (nDim == 3 && nodetype == NODETYPE::EDGE){
      for (int i = 0; i < nDim; i++)
      responseSendBuffer[cnt++] = secondaryNormal[i];
    }   
  }

  
  // obtaining the response send and receive counts and starting indices
  vector<int> responseSendCounts(size), responseSendDisps(size);
  vector<int> responseRecvCounts(size), responseRecvDisps(size);

  

  for (int i = 0; i < size; ++i) {
      responseSendCounts[i] = recvCounts[i] * dataSize * nDim;
      responseRecvCounts[i] = sendCounts[i] * dataSize * nDim;
  }

  responseSendDisps[0] = 0;
  responseRecvDisps[0] = 0;
  for (int i = 1; i < size; ++i) {
      responseSendDisps[i] = responseSendDisps[i - 1] + responseSendCounts[i - 1];
      responseRecvDisps[i] = responseRecvDisps[i - 1] + responseRecvCounts[i - 1];
  }

  // Allocate receive buffer
  std::vector<su2double> response_recv_buffer(
      accumulate(responseRecvCounts.begin(), responseRecvCounts.end(), 0)
  );

  SU2_MPI::Alltoallv(responseSendBuffer.data(), responseSendCounts.data(), responseSendDisps.data(), MPI_DOUBLE,
              response_recv_buffer.data(), responseRecvCounts.data(), responseRecvDisps.data(), MPI_DOUBLE,
              SU2_MPI::GetComm());

  return response_recv_buffer;
}

void CRadialBasisFunctionInterpolation::GetVertexNormals(CGeometry* geometry, CConfig* config, const unsigned short marker, const NODETYPE type, CVertex* const nearestVertex, su2double* normal, su2double* secondaryNormal ) const {

  const su2double* cartNormal = nearestVertex->GetNormal();
  const su2double* cartCoord = nearestVertex->GetCoord();

  if (IsCylindrical) {
      CartDispToCyl(cartCoord, cartNormal, normal);
  } else {
      copy_n(cartNormal, nDim, normal);
  }

  PruneVector(normal);
  
  
  if (nDim == 3 && type == NODETYPE::EDGE){

    const auto nearestNodeIndex = nearestVertex->GetNode();
    const auto nMarkers = config->GetnMarker_All();      
    unsigned short neighborBoundaryMarker = nMarkers; 
    for( auto iMarker = 0u; iMarker < nMarkers; iMarker++){
      if(geometry->nodes->GetVertex(nearestNodeIndex, iMarker) != -1 && iMarker != marker && !config->GetMarker_All_PerBound(iMarker)){
        neighborBoundaryMarker = iMarker;
        break;
      } 
    }

    if (neighborBoundaryMarker != nMarkers) {
      const auto neighborBoundaryVertex = geometry->nodes->GetVertex(nearestNodeIndex, neighborBoundaryMarker);

      const su2double* const secondaryCartNormal = geometry->vertex[neighborBoundaryMarker][neighborBoundaryVertex]->GetNormal();
      
      if (IsCylindrical){
        CartDispToCyl(cartCoord, secondaryCartNormal, secondaryNormal);
      } else{
        copy_n(secondaryCartNormal, nDim, secondaryNormal);
      }

      PruneVector(secondaryNormal);
    }
  }
}

void CRadialBasisFunctionInterpolation::SetNearestNodeData(CGeometry* geometry, CConfig* config, const vector<su2double>& responseRecvBuffer, const unordered_set<unsigned long>& targetSet, const vector<CRadialBasisFunctionNode*>& targetNodes) {

  if (targetSet.empty()) return;

  auto offset = 0;
  const auto dataSize = (nDim == 3 && targetNodes[*targetSet.begin()]->GetNodeType() == NODETYPE::EDGE) ? 3 : 2;

  su2double temp[3] = {0.0, 0.0, 0.0};

  for (const auto index : targetSet) {
     auto* iNode = targetNodes[index];

  // for (auto* iNode : targetNodes) {
  // // for (auto index : targetSet) {
  // //   auto* iNode = targetNodes[index];
    if (iNode->GetNearestRank() == rank) {
      SetLocalNormal(geometry, config, iNode); 
    } else {
        
        for (int iDim = 0; iDim < nDim; iDim++) {
          temp[iDim] = responseRecvBuffer[offset + iDim];
        }
        iNode->SetNearestNormal(temp, nDim);

        for (int iDim = 0; iDim < nDim; iDim++) {           
          temp[iDim] = responseRecvBuffer[offset + nDim + iDim];
        }
        iNode->SetNearestCoord(temp, nDim);

        if (dataSize == 3) {
          for (int iDim = 0; iDim < nDim; iDim++) {
            temp[iDim] = responseRecvBuffer[offset + 2 * nDim + iDim];
          }
          iNode->SetNearestNormal2(temp, nDim);
        }

      offset += dataSize * nDim;
    }
  }
}

void CRadialBasisFunctionInterpolation::UpdateVarCoord(CGeometry* geometry, CConfig* config, const CRadialBasisFunctionNode* const node, const su2double* const projection) const {
  su2double coord[nDim];
  const su2double* cart_coord = geometry->nodes->GetCoord(node->GetIndex());
  const short marker = node->GetMarker();

  if (IsCylindrical) {
      CartToCyl(cart_coord, coord);
  } else {
      std::copy_n(cart_coord, nDim, coord);
  }

  const su2double* new_coord = node->GetNewCoord();
  su2double coordVar[3] = {0.0, 0.0, 0.0};

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    coordVar[iDim] = new_coord[iDim] - projection[iDim] - coord[iDim]; // contains coordinate variation
  }
  
  if(IsCylindrical){
    CylDispToCart(cart_coord, coordVar);
  }

  const auto iNonLinear_iter = (config->GetGridDef_Nonlinear_Iter());
  for (auto iDim = 0u; iDim < nDim; iDim++){
    coordVar[iDim] = iNonLinear_iter * coordVar[iDim];
  }
  
  geometry->vertex[marker][node->GetVertex()]->SetVarCoord(coordVar);
}