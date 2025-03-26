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
#include "../../include/toolboxes/geometry_toolbox.hpp"
#include "../../include/adt/CADTPointsOnlyClass.hpp"


CRadialBasisFunctionInterpolation::CRadialBasisFunctionInterpolation(CGeometry* geometry, CConfig* config) : CVolumetricMovement(geometry) {}

CRadialBasisFunctionInterpolation::~CRadialBasisFunctionInterpolation() = default;

void CRadialBasisFunctionInterpolation::SetVolume_Deformation(CGeometry* geometry, CConfig* config, bool UpdateGeo, bool Derivative,
                                                bool ForwardProjectionDerivative){



  // {
  //   if (rank==MASTER_NODE){
  //     int i = 0;
  //     while(0==i){
  //       sleep(5);
  //     }
  //   }
  // }

  
  // SU2_MPI::Barrier(SU2_MPI::GetComm());
  // SU2_MPI::Abort(SU2_MPI::GetComm(), 0);
  

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


  /*--- Determining the boundary and internal nodes. Setting the control nodes. ---*/ 
  SetBoundNodes(geometry, config);
  
  vector<unsigned long> internalNodes; 
  SetInternalNodes(geometry, config, internalNodes); 
  
  SetCtrlNodes(config);

  /*--- Looping over the number of deformation iterations ---*/
  for (auto iNonlinear_Iter = 0ul; iNonlinear_Iter < Nonlinear_Iter; iNonlinear_Iter++) {
    
    /*--- Compute min volume in the entire mesh. ---*/

    ComputeDeforming_Element_Volume(geometry, MinVolume, MaxVolume, Screen_Output);
    if (rank == MASTER_NODE && Screen_Output)
      cout << "Min. volume: " << MinVolume << ", max. volume: " << MaxVolume << "." << endl;
    

    /*--- Solving the RBF system, resulting in the interpolation coefficients ---*/
    SolveRBF_System(geometry, config, kindRBF, radius, Derivative, internalNodes, ForwardProjectionDerivative);
   
    /*--- Updating the coordinates of the grid ---*/
    if (!Derivative) {
      UpdateGridCoord(geometry, config, kindRBF, radius, internalNodes);
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
}

void CRadialBasisFunctionInterpolation::SolveRBF_System(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, bool Derivative, vector<unsigned long>& internalNodes, bool ForwardProjectionDerivative){
  
  /*--- In case of data reduction an iterative greedy algorithm is applied 
          to perform the interpolation with a reduced set of control nodes.
          Otherwise with a full set of control nodes. ---*/
  
  if(config->GetRBF_DataReduction()){
    
    /*--- Local maximum error node and corresponding maximum error  ---*/
    unsigned long maxErrorNodeLocal;
    su2double maxErrorLocal{0};

    /*--- Obtaining the initial maximum error nodes, which are found based on the maximum applied deformation. */
    if(ControlNodes->empty()){
      GetInitMaxErrorNode(geometry, config, Derivative, maxErrorNodeLocal, maxErrorLocal); 
      SU2_MPI::Allreduce(&maxErrorLocal, &MaxErrorGlobal, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
    }

    /*--- Number of greedy iterations. ---*/
    unsigned short greedyIter = 0;
    
    /*--- Error tolerance for the data reduction tolerance ---*/
    const su2double dataReductionTolerance = config->GetRBF_DataRedTolerance() * MaxErrorGlobal; 

    su2passivematrix invInterpMat;
    /*--- While the maximum error is above the tolerance, data reduction algorithm is continued. ---*/
    while(MaxErrorGlobal > dataReductionTolerance || greedyIter == 0){ 
      
      /*--- In case of a nonzero local error, control nodes are added ---*/
      if(maxErrorLocal> 0){
        AddControlNode(maxErrorNodeLocal);
      }

      /*--- Obtaining the global number of control nodes. ---*/
      Get_nCtrlNodesGlobal();
      
      /*--- Obtaining the control nodes coordinates and distributing over all processes. ---*/
      SetCtrlNodeCoords(geometry);

      /*--- Obtaining the deformation of the control nodes. ---*/
      if (!Derivative){
        SetBoundaryDisplacements(geometry, config);
      }else{
        SetCtrlNodeDerivatives(geometry, config, ForwardProjectionDerivative);
      } 

      /*--- Computation of the (inverse) interpolation matrix. ---*/
      
      GetInvInterpMat(geometry, type, radius, invInterpMat);
      
      /*--- Obtaining the interpolation coefficients. ---*/
      ComputeInterpCoeffs(invInterpMat);

      /*--- Determining the interpolation error, of the non-control boundary nodes. ---*/
      GetInterpError(geometry, config, type, radius, Derivative, maxErrorNodeLocal, maxErrorLocal); 
      SU2_MPI::Allreduce(&maxErrorLocal, &MaxErrorGlobal, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());

      if(rank == MASTER_NODE) cout << "Greedy iteration: " << greedyIter << ". Max error: " << MaxErrorGlobal << ". Global nr. of ctrl nodes: "  << nCtrlNodesGlobal << "\n" << endl;
      
      //TODO  debug output
      if (Derivative){
        ofstream res_out("greedy_out.txt");
        for (auto x =0ul; x < ControlNodes->size(); x++){
          res_out << (*ControlNodes)[x]->GetIndex() << "\t" << SU2_TYPE::GetValue(CtrlNodeDeformation[x*nDim]) << "\t" << SU2_TYPE::GetValue(CtrlNodeDeformation[x*nDim+1]) << endl;
        }

        for (auto x : BoundNodes){
          res_out << x->GetIndex() << "\t" << SU2_TYPE::GetValue(geometry->GetSensitivity(x->GetIndex(), 0))  +  x->GetError()[0] << "\t" <<  SU2_TYPE::GetValue(geometry->GetSensitivity(x->GetIndex(), 1)) + x->GetError()[1] << endl;
        }
        res_out.close();
      }
      
      greedyIter++;

    } 

    if (Derivative){
      SetInternalNodeDerivatives(geometry, config, internalNodes, ForwardProjectionDerivative);
      ComputeSensitivity(geometry, config, type, radius, invInterpMat, internalNodes); //TODO Investigate possibility of doing correction for the non-selected nodes. 
    }
  }else{
    /*--- Obtaining the interpolation coefficients. ---*/
    GetInterpCoeffs(geometry, config, type, radius, Derivative, internalNodes, ForwardProjectionDerivative);
  }
}

void CRadialBasisFunctionInterpolation::GetInterpCoeffs(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, bool Derivative, vector<unsigned long>& internalNodes, bool ForwardProjectionDerivative){
  
  /*--- Obtaining the control nodes coordinates and distributing over all processes. ---*/
  SetCtrlNodeCoords(geometry);

  /*--- Obtaining the deformation of the control nodes. ---*/
  if (!Derivative){
    SetBoundaryDisplacements(geometry, config);
  }else{
    SetInternalNodeDerivatives(geometry, config, internalNodes, ForwardProjectionDerivative);
  } 

  /*--- Computation of the (inverse) interpolation matrix. ---*/
  su2passivematrix invInterpMat;
  GetInvInterpMat(geometry, type, radius, invInterpMat);

  if(!Derivative){
    /*--- Obtaining the interpolation coefficients. ---*/
    ComputeInterpCoeffs(invInterpMat);
  }else{
    ComputeSensitivity(geometry, config, type, radius, invInterpMat, internalNodes);
  }
}

void CRadialBasisFunctionInterpolation::ComputeSensitivity(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius,su2passivematrix &invInterpMat, vector<unsigned long>& internalNodes){

  //from here
  vector<su2double> inter_res(nCtrlNodesGlobal * nDim, 0.0);

  //loop over all global control nodes
  for (auto iNode = 0ul; iNode < nCtrlNodesGlobal; iNode++){
    for (auto iDim =0u; iDim < nDim; iDim++){
      
      // loop over the local internal nodes
      for (auto jNode =0ul; jNode < internalNodes.size(); jNode++){

        auto dist = GeometryToolbox::Distance(nDim, CtrlCoords[iNode*nDim], geometry->nodes->GetCoord(internalNodes[jNode]));
        auto rbf_eval = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));
        
        // CtrlNodeDeformation vector contains the internal node sensitivities. 
        // Phi_c,i * d_i
        inter_res[iNode*nDim+iDim] += rbf_eval * CtrlNodeDeformation[jNode*nDim+iDim];
      }
    }
  }

  // Summing contributions of all processes
  vector<su2double> inter_res_sum(nCtrlNodesGlobal * nDim);
  SU2_MPI::Allreduce(inter_res.data(), inter_res_sum.data(), nCtrlNodesGlobal*nDim, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

  // arrays storing global and local sensitivity deltas
  sensitivity_update.resize(nCtrlNodesGlobal * nDim);
  vector<su2double> sens_update_local(nCtrlNodesGlobal*nDim, 0.0);
  
  // inv(Phi_cc) *inter_res
  for (auto iNode = rank; iNode < nCtrlNodesGlobal; iNode += size){
    for (auto iDim =0u; iDim < nDim; iDim++){
      for (auto jNode = 0ul; jNode < nCtrlNodesGlobal; jNode++){
        sens_update_local[iNode*nDim+iDim] += invInterpMat(iNode, jNode) * inter_res_sum[jNode*nDim+iDim];
      }
    }
  }

  SU2_MPI::Allreduce(sens_update_local.data(), sensitivity_update.data(), nCtrlNodesGlobal*nDim, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
}

void CRadialBasisFunctionInterpolation::SetBoundNodes(CGeometry* geometry, CConfig* config){
  
  /*--- Storing of the local node, marker and vertex information of the boundary nodes ---*/

  /*--- Looping over the markers ---*/
  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {

    /*--- Checking if not internal or send/receive marker ---*/
    if (!config->GetMarker_All_Deform_Mesh_Internal(iMarker) && !config->GetMarker_All_SendRecv(iMarker)) {

      /*--- Looping over the vertices of marker ---*/
      for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        /*--- Node in consideration ---*/
        auto iNode = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Check whether node is part of the subdomain and not shared with a receiving marker (for parallel computation)*/
        if (geometry->nodes->GetDomain(iNode)) {
          BoundNodes.push_back(new CRadialBasisFunctionNode(iNode, iMarker, iVertex));        
        }        
      }
    }
  }

  /*--- Sorting of the boundary nodes based on their index ---*/
  sort(BoundNodes.begin(), BoundNodes.end(), HasSmallerIndex);

  /*--- Obtaining unique set ---*/
  BoundNodes.resize(std::distance(BoundNodes.begin(), unique(BoundNodes.begin(), BoundNodes.end(), HasEqualIndex)));
}

void CRadialBasisFunctionInterpolation::SetCtrlNodes(CConfig* config){
  
  /*--- Assigning the control nodes based on whether data reduction is applied or not. ---*/
  if(config->GetRBF_DataReduction()){

    /*--- Control nodes are an empty set ---*/
    ControlNodes = &ReducedControlNodes;
  }else{

    /*--- Control nodes are the boundary nodes ---*/
    ControlNodes = &BoundNodes;
  }

  /*--- Obtaining the total number of control nodes. ---*/
  Get_nCtrlNodesGlobal();

};

void CRadialBasisFunctionInterpolation::GetInvInterpMat(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius, su2passivematrix& invInterpMat) {
  
  CSymmetricMatrix interpMat;

  #ifdef HAVE_MPI
    GetInterpMat_parallel(geometry, type, radius, interpMat);    
  #else
    GetInterpMat_sequential(geometry, type, radius, interpMat);    
  #endif
  
  const bool kernelIsSPD = (type == RADIAL_BASIS::WENDLAND_C2) || (type == RADIAL_BASIS::GAUSSIAN) ||
                          (type == RADIAL_BASIS::INV_MULTI_QUADRIC);

  interpMat.Invert(kernelIsSPD); 

  invInterpMat = interpMat.StealData();

  if (rank == MASTER_NODE){
    ofstream out;
    out.open("rbf_mat_inv.txt");
    if (out.is_open()) {
      for (auto row_i = 0ul; row_i < nCtrlNodesGlobal; row_i++){
        for (auto col_i = 0ul; col_i < nCtrlNodesGlobal; col_i++){
          out << invInterpMat(row_i, col_i) << "\t";
        }  
        out << endl;
      }
      out.close();
    } 
  }
}

void CRadialBasisFunctionInterpolation::GetInterpMat_parallel(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius, CSymmetricMatrix& interpMat){
  
  /*--- Initialization of the interpolation matrix ---*/
  interpMat.Initialize(nCtrlNodesGlobal);

  // Total number of elems in the lower triangular matrix
  unsigned long N_lowerTriangle = (nCtrlNodesGlobal*(nCtrlNodesGlobal+1))/2;

  // Average number of elements per process
  unsigned long N_perProcess = ceil(N_lowerTriangle/size);

  // For balancing the number of elements per process, the start and end rows are determined.
  // The number of elements in a lower triangle matrix of size n x n is given by: T_n = n(n+1)/2
  // Starting and ending row is determined by solving this equations for the number of elements (rank*N_perProcess) up to that row:
  // row(row+1)/2 = rank*N_perProcess. Using quadratic formula results in: row = (-1 + sqrt(1+8*N_perProcess*rank))/2.
  // Ceil is used to obtain integer numbers.

  unsigned long start_row = ceil((-1 + sqrt(1+8*N_perProcess*(rank))) / 2);
  unsigned long end_row = ceil((-1 + sqrt(1+8*N_perProcess*(rank+1))) / 2);
  if (end_row > nCtrlNodesGlobal)  end_row = nCtrlNodesGlobal;  
  
  // Number of elements to be evaluated
  int nr_elems = (end_row*(end_row+1) - start_row*(start_row+1))/2;
  
  // Finding RBF evaluations
  vector<su2double> rbf_vals(nr_elems);
  unsigned long cnt = 0;

  for (auto row_i = start_row; row_i < end_row; row_i++){
    for (auto col_i = 0ul; col_i <= row_i; col_i++){
      auto dist = GeometryToolbox::Distance(nDim, CtrlCoords[row_i*nDim], CtrlCoords[col_i*nDim]);
      rbf_vals[cnt++] = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));
    }
  }
  // Gathering the number of elements on all nodes
  int recv_cnts[size];
  SU2_MPI::Allgather(&nr_elems, 1, MPI_INT, recv_cnts, 1, MPI_INT, SU2_MPI::GetComm());

  // Determining the starting index/displacement for each process
  int disp[size];
  disp[0] = 0;
  for (int i = 1; i < size; i++) {
      disp[i] = disp[i-1] + recv_cnts[i-1];
  }

  // global rbf evaluations are stored in rbf_vals_all
  vector<su2double> rbf_vals_all(N_lowerTriangle);

  // Allgather the RBF evaluations
  SU2_MPI::Allgatherv(rbf_vals.data(), nr_elems, MPI_DOUBLE, rbf_vals_all.data(), recv_cnts, disp, MPI_DOUBLE, SU2_MPI::GetComm());

  // fill the lower triangular part of the interpolation matrix  
  cnt = 0;
  for (auto row_i = 0ul; row_i < nCtrlNodesGlobal; row_i++) {
    for (auto col_i = 0ul; col_i <= row_i; col_i++) {
      interpMat.Set(row_i, col_i, SU2_TYPE::GetValue(rbf_vals_all[cnt++]));
    }
  }
}

void CRadialBasisFunctionInterpolation::GetInterpMat_sequential(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius, CSymmetricMatrix& interpMat){

  /*--- Initialization of the interpolation matrix ---*/
  interpMat.Initialize(nCtrlNodesGlobal);    
    
  /*--- Construction of the interpolation matrix. 
      Since this matrix is symmetric only upper halve has to be considered ---*/

  /*--- Looping over the target nodes ---*/
  for( auto iNode = 0ul; iNode < nCtrlNodesGlobal; iNode++ ){

    /*--- Looping over the control nodes ---*/
    for ( auto jNode = iNode; jNode < nCtrlNodesGlobal; jNode++){
      
      /*--- Distance between nodes ---*/
      auto dist = GeometryToolbox::Distance(nDim, CtrlCoords[iNode*nDim], CtrlCoords[jNode*nDim]);   

      /*--- Evaluation of RBF ---*/
      interpMat(iNode, jNode) = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));
    }
  }
}


void CRadialBasisFunctionInterpolation::SetBoundaryDisplacements(CGeometry* geometry, CConfig* config){

  
  /* --- Initialization of the deformation vector ---*/
  CtrlNodeDeformation.resize(ControlNodes->size()*nDim, 0.0); 

  /*--- If requested (no by default) impose the surface deflections in
    increments and solve the grid deformation with
    successive small deformations. ---*/
  const su2double VarIncrement = 1.0 / ((su2double)config->GetGridDef_Nonlinear_Iter());

  SU2_COMPONENT Kind_SU2 = config->GetKind_SU2();

  /*--- Loop over the control nodes ---*/
  for (auto iNode = 0ul; iNode < ControlNodes->size(); iNode++) {
    
     /*--- Setting nonzero displacement of the moving markers, else setting zero displacement for static markers---*/
    auto iMarker = (*ControlNodes)[iNode]->GetMarker();
    if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_CFD)) ||
        ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_DEF)) ||
        ((config->GetDirectDiff() == D_DESIGN) && (Kind_SU2 == SU2_COMPONENT::SU2_CFD) &&
         (config->GetMarker_All_DV(iMarker) == YES)) /*NOTE: This feature has not been tested for RBF interpolation*/ ||
        ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_DOT))) {
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        CtrlNodeDeformation[iNode*nDim + iDim] = SU2_TYPE::GetValue(geometry->vertex[(*ControlNodes)[iNode]->GetMarker()][(*ControlNodes)[iNode]->GetVertex()]->GetVarCoord()[iDim] * VarIncrement);
      }
    }
  }
}

void CRadialBasisFunctionInterpolation::SetInternalNodeDerivatives(CGeometry* geometry, CConfig* config, vector<unsigned long>& internalNodes, bool ForwardProjectionDerivative) {
  
  SU2_COMPONENT Kind_SU2 = config->GetKind_SU2();
  
  if ((Kind_SU2 == SU2_COMPONENT::SU2_DOT) && !ForwardProjectionDerivative) {
    CtrlNodeDeformation.resize(internalNodes.size() * nDim, 0.0); 

    for (auto iNode = 0ul; iNode < internalNodes.size(); iNode++){
      for (auto iDim =0u; iDim < nDim; iDim++){
        CtrlNodeDeformation[iNode * nDim + iDim] = SU2_TYPE::GetValue(geometry->GetSensitivity(internalNodes[iNode], iDim));
      }
    }
  } else {
    SU2_MPI::Error("Missing feature in RBF interpolation", CURRENT_FUNCTION);
  }
}

void CRadialBasisFunctionInterpolation::SetCtrlNodeDerivatives(CGeometry* geometry, CConfig* config, bool ForwardProjectionDerivative){
  
  SU2_COMPONENT Kind_SU2 = config->GetKind_SU2();

  if ((Kind_SU2 == SU2_COMPONENT::SU2_DOT) && !ForwardProjectionDerivative) {
    CtrlNodeDeformation.resize(ControlNodes->size() * nDim, 0.0);

    for (auto iNode = 0ul; iNode < ControlNodes->size(); iNode++){
      if ((config->GetMarker_All_DV((*ControlNodes)[iNode]->GetMarker()) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_DOT)) {  
        for (auto iDim = 0u; iDim < nDim; iDim++){
          CtrlNodeDeformation[iNode * nDim + iDim] = SU2_TYPE::GetValue(geometry->GetSensitivity((*ControlNodes)[iNode]->GetIndex(), iDim));
        }
      }
    }
  } else {
    SU2_MPI::Error("Missing feature in RBF interpolation", CURRENT_FUNCTION);
  }
}

void CRadialBasisFunctionInterpolation::SetInternalNodes(CGeometry* geometry, CConfig* config, vector<unsigned long>& internalNodes){ 

  /*--- Looping over all nodes and check if part of domain and not on boundary ---*/
  for (auto iNode = 0ul; iNode < geometry->GetnPoint(); iNode++) {    
    if (!geometry->nodes->GetBoundary(iNode)) {
      internalNodes.push_back(iNode);
    }   
  }  

  /*--- Adding nodes on markers considered as internal nodes ---*/
  for (auto iMarker = 0u; iMarker < geometry->GetnMarker(); iMarker++){

    /*--- Check if marker is considered as internal nodes ---*/
    if(config->GetMarker_All_Deform_Mesh_Internal(iMarker)){
      
      /*--- Loop over marker vertices ---*/
      for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) { 

        /*--- Local node index ---*/
        auto iNode = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- if not among the boundary nodes ---*/
        if (find_if (BoundNodes.begin(), BoundNodes.end(), [&](CRadialBasisFunctionNode* i){return i->GetIndex() == iNode;}) == BoundNodes.end()) {
          internalNodes.push_back(iNode);
        }            
      }
    }
  }


  /*--- In case of a parallel computation, the nodes on the send/receive markers are included as internal nodes
          if they are not already a boundary node with known deformation ---*/

  #ifdef HAVE_MPI
    /*--- Looping over the markers ---*/
    for (auto iMarker = 0u; iMarker < geometry->GetnMarker(); iMarker++) { 

      /*--- If send or receive marker ---*/
      if (config->GetMarker_All_SendRecv(iMarker)) { 

        /*--- Loop over marker vertices ---*/
        for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) { 
          
          /*--- Local node index ---*/
          auto iNode = geometry->vertex[iMarker][iVertex]->GetNode();

          /*--- if not among the boundary nodes ---*/
          if (find_if (BoundNodes.begin(), BoundNodes.end(), [&](CRadialBasisFunctionNode* i){return i->GetIndex() == iNode;}) == BoundNodes.end()) {
            internalNodes.push_back(iNode);
          }             
        }
      }
    }

    /*--- sorting of the local indices ---*/
    sort(internalNodes.begin(), internalNodes.end());

    /*--- Obtaining unique set of internal nodes ---*/
    internalNodes.resize(std::distance(internalNodes.begin(), unique(internalNodes.begin(), internalNodes.end())));
  #endif

}

void CRadialBasisFunctionInterpolation::ComputeInterpCoeffs(su2passivematrix& invInterpMat) {

  /*--- resizing the interpolation coefficient vector ---*/
  InterpCoeff.resize(nDim * nCtrlNodesGlobal);

  /*--- Each process computes a portion of the result ---*/
  vector<su2double> interpCoeffLocal(nDim * nCtrlNodesGlobal, 0.0);


  unsigned long nCtrlNode = ControlNodes->size();
  unsigned long nCtrlNodes[size];
  SU2_MPI::Allgather(&nCtrlNode, 1, MPI_UNSIGNED_LONG, nCtrlNodes, 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm()); 
  
  unsigned long start_idx = 0;
  for (auto iProc = 0; iProc < rank; iProc++){
    start_idx += nCtrlNodes[iProc];
  }

  for (auto iNode = 0; iNode < ControlNodes->size(); iNode++){
    for (auto iDim = 0u; iDim < nDim; iDim++){
      for (auto jNode = 0; jNode < nCtrlNodesGlobal; jNode++ ){
        // interpCoeffLocal has sums of global vector, CtrlNodeDeformation vector should be the local one
        interpCoeffLocal[jNode*nDim + iDim] += invInterpMat(jNode, iNode + start_idx) * CtrlNodeDeformation[iNode*nDim+iDim];
      }
    }
  }
  
  /*--- Gather the results from all processes on the master node ---*/
  #ifdef HAVE_MPI
    SU2_MPI::Allreduce(interpCoeffLocal.data(), InterpCoeff.data(), interpCoeffLocal.size(), MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  #else
    InterpCoeff = move(interpCoeffLocal);
  #endif
}

void CRadialBasisFunctionInterpolation::UpdateGridCoord(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, const vector<unsigned long>& internalNodes){
  
  if(rank == MASTER_NODE){
    cout << "updating the grid coordinates" << endl;
  }
  
  /*--- Update of internal node coordinates ---*/
  UpdateInternalCoords(geometry, type, radius, internalNodes);

  /*--- Update of boundary node coordinates ---*/
  UpdateBoundCoords(geometry, config, type, radius);   

  /*--- In case of data reduction, perform the correction for nonzero error nodes ---*/
  if(config->GetRBF_DataReduction() && BoundNodes.size() > 0){
    SetCorrection(geometry, config, type, internalNodes); 
  }
}

void CRadialBasisFunctionInterpolation::UpdateGridCoord_Derivatives(CGeometry* geometry, CConfig* config, bool ForwardProjectionDerivative){
  
  
  SU2_COMPONENT Kind_SU2 = config->GetKind_SU2();

  if ((Kind_SU2 == SU2_COMPONENT::SU2_DOT) && !ForwardProjectionDerivative) {
  
    unsigned long nCtrlNode = ControlNodes->size();
    unsigned long nCtrlNodes[size];
    SU2_MPI::Allgather(&nCtrlNode, 1, MPI_UNSIGNED_LONG, nCtrlNodes, 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm()); 

    unsigned long start_idx = 0;
    for (auto iProc = 0; iProc < rank; iProc++){
      start_idx += nCtrlNodes[iProc]*nDim;
    }
    
    for (auto iNode = 0ul; iNode < ControlNodes->size(); iNode++) {
      auto iMarker = (*ControlNodes)[iNode]->GetMarker();
      if (config->GetSolid_Wall(iMarker) || (config->GetMarker_All_DV(iMarker) == YES)) {
              
        auto iPoint = (*ControlNodes)[iNode]->GetIndex();
        for (auto iDim = 0u; iDim < nDim; iDim++) {      
          
          // summation of current sensitivity and the computed update
          su2double sens_new =  geometry->GetSensitivity(iPoint, iDim) + sensitivity_update[start_idx + iNode*nDim+iDim];
          geometry->SetSensitivity(iPoint, iDim, sens_new);
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
    
    /*--- Loop for contribution of each control node ---*/
    for(auto jNode = 0ul; jNode < nCtrlNodesGlobal; jNode++){

      /*--- Determine distance between considered internal and control node ---*/
      auto dist = GeometryToolbox::Distance(nDim, CtrlCoords[jNode*nDim], geometry->nodes->GetCoord(internalNodes[iNode]));

      /*--- Evaluate RBF based on distance ---*/
      auto rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));
      
      /*--- Add contribution to total coordinate variation ---*/
      for(auto iDim = 0u; iDim < nDim; iDim++){
        var_coord[iDim] += rbf*InterpCoeff[jNode * nDim + iDim];
      }
    }

    /*--- Apply the coordinate variation and resetting the var_coord vector to zero ---*/
    for(auto iDim = 0u; iDim < nDim; iDim++){
      geometry->nodes->AddCoord(internalNodes[iNode], iDim, var_coord[iDim]);
      var_coord[iDim] = 0;
    } 
  }  
}

void CRadialBasisFunctionInterpolation::UpdateBoundCoords(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius){
  
  /*--- Vector for storing the coordinate variation ---*/
  su2double var_coord[nDim]{0.0};
  
  /*--- In case of data reduction, the non-control boundary nodes are treated as if they where internal nodes ---*/
  if(config->GetRBF_DataReduction()){

    /*--- Looping over the non selected boundary nodes ---*/
    for(auto iNode = 0ul; iNode < BoundNodes.size(); iNode++){
      
      /*--- Finding contribution of each control node ---*/
      for( auto jNode = 0ul; jNode <  nCtrlNodesGlobal; jNode++){
        
        /*--- Distance of non-selected boundary node to control node ---*/
        auto dist = GeometryToolbox::Distance(nDim, CtrlCoords[jNode*nDim], geometry->nodes->GetCoord(BoundNodes[iNode]->GetIndex()));
        
        /*--- Evaluation of the radial basis function based on the distance ---*/
        auto rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));

        /*--- Computing and add the resulting coordinate variation ---*/
        for(auto iDim = 0u; iDim < nDim; iDim++){
          var_coord[iDim] += rbf*InterpCoeff[jNode * nDim + iDim];
        }
      }

      /*--- Applying the coordinate variation and resetting the var_coord vector*/
      for(auto iDim = 0u; iDim < nDim; iDim++){
        geometry->nodes->AddCoord(BoundNodes[iNode]->GetIndex(), iDim, var_coord[iDim]);
        var_coord[iDim] = 0;
      }

    }
  }

  /*--- Applying the surface deformation, which are stored in the deformation vector ---*/
  for(auto jNode = 0ul; jNode < ControlNodes->size(); jNode++){ 
    if(config->GetMarker_All_Moving((*ControlNodes)[jNode]->GetMarker()) || config->GetMarker_All_DV((*ControlNodes)[jNode]->GetMarker())){
      for(auto iDim = 0u; iDim < nDim; iDim++){
          geometry->nodes->AddCoord((*ControlNodes)[jNode]->GetIndex(), iDim, CtrlNodeDeformation[jNode*nDim + iDim]); 
      }
    }
  }
}

void CRadialBasisFunctionInterpolation::GetInitMaxErrorNode(CGeometry* geometry, CConfig* config, bool Derivative,  unsigned long& maxErrorNodeLocal, su2double& maxErrorLocal){

  /*--- Set max error to zero ---*/
  maxErrorLocal = 0.0;

  /*--- Loop over the nodes ---*/  
  for (auto iNode = 0ul; iNode < BoundNodes.size(); iNode++){

    su2double normSquaredDeformation = 0;
    /*--- Compute to squared norm of the deformation ---*/
    if (Derivative){
      auto iPoint = BoundNodes[iNode]->GetIndex();
      for (auto iDim = 0u; iDim < nDim; iDim++){
        passivedouble sens_i = SU2_TYPE::GetValue(geometry->GetSensitivity(iPoint, iDim));
        normSquaredDeformation += sens_i * sens_i;
      }
    } else {
      normSquaredDeformation = GeometryToolbox::SquaredNorm(nDim, geometry->vertex[BoundNodes[iNode]->GetMarker()][BoundNodes[iNode]->GetVertex()]->GetVarCoord());    
    }
    
    /*--- In case squared norm deformation is larger than the error, update the error ---*/
    if(normSquaredDeformation > maxErrorLocal){
      maxErrorLocal = normSquaredDeformation;
      maxErrorNodeLocal = iNode;
    }
  }

  /*--- Account for the possibility of applying the deformation in multiple steps ---*/
  maxErrorLocal = sqrt(maxErrorLocal) / ((su2double)config->GetGridDef_Nonlinear_Iter());
}


void CRadialBasisFunctionInterpolation::SetCtrlNodeCoords(CGeometry* geometry){
  /*--- The coordinates of all control nodes are made available on all processes ---*/
  
  /*--- resizing the matrix containing the global control node coordinates ---*/
  CtrlCoords.resize(nCtrlNodesGlobal*nDim);
  
  /*--- Array containing the local control node coordinates ---*/ 
  su2double localCoords[nDim*ControlNodes->size()];
  
  /*--- Storing local control node coordinates ---*/
  for(auto iNode = 0ul; iNode < ControlNodes->size(); iNode++){
    auto coord = geometry->nodes->GetCoord((*ControlNodes)[iNode]->GetIndex());  
    for ( auto iDim = 0u ; iDim < nDim; iDim++ ){
      localCoords[ iNode * nDim + iDim ] = coord[iDim];
    }
  }

  /*--- Gathering local control node coordinate sizes on all processes. ---*/
  int LocalCoordsSizes[size];
  int localCoordsSize = nDim*ControlNodes->size();
  SU2_MPI::Allgather(&localCoordsSize, 1, MPI_INT, LocalCoordsSizes, 1, MPI_INT, SU2_MPI::GetComm()); 

  /*--- Array containing the starting indices for the allgatherv operation */
  int disps[SU2_MPI::GetSize()] = {0};    

  for(auto iProc = 1; iProc < SU2_MPI::GetSize(); iProc++){
    disps[iProc] = disps[iProc-1]+LocalCoordsSizes[iProc-1];
  }
  
  /*--- Distributing global control node coordinates among all processes ---*/
  SU2_MPI::Allgatherv(&localCoords, localCoordsSize, MPI_DOUBLE, CtrlCoords.data(), LocalCoordsSizes, disps, MPI_DOUBLE, SU2_MPI::GetComm()); 
};


void CRadialBasisFunctionInterpolation::GetInterpError(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, bool Derivative, unsigned long& maxErrorNodeLocal, su2double& maxErrorLocal){
  /*--- Array containing the local error ---*/
  su2double localError[nDim];

  /*--- Magnitude of the local maximum error ---*/
  maxErrorLocal = 0.0;

  /*--- Loop over non-selected boundary nodes ---*/
  for(auto iNode = 0ul; iNode < BoundNodes.size(); iNode++){

    /*--- Compute nodal error ---*/
    GetNodalError(geometry, config, type, radius, iNode, Derivative, localError);

    /*--- Setting error ---*/
    BoundNodes[iNode]->SetError(localError, nDim);
    
    /*--- Compute error magnitude and update local maximum error if necessary ---*/
    su2double errorMagnitude = GeometryToolbox::Norm(nDim, localError);
    if(errorMagnitude > maxErrorLocal){
      maxErrorLocal = errorMagnitude;
      maxErrorNodeLocal = iNode;
    }
  }  
}

void CRadialBasisFunctionInterpolation::GetNodalError(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, unsigned long iNode, bool Derivative, su2double* localError){ 
  
  /*--- If requested (no by default) impose the surface deflections in increments ---*/
  const su2double VarIncrement = 1.0 / ((su2double)config->GetGridDef_Nonlinear_Iter());
  
  /*--- If node is part of a moving boundary then the error is defined as the difference
           between the found and prescribed displacements. Thus, here the displacement is substracted from the error ---*/

  if (Derivative) {
    for (auto iDim =0u; iDim < nDim; iDim++) {
      localError[iDim] = -geometry->GetSensitivity(BoundNodes[iNode]->GetIndex(), iDim) * VarIncrement;
    }
  }else{
  
    if(config->GetMarker_All_Moving(BoundNodes[iNode]->GetMarker())){
      auto displacement = geometry->vertex[BoundNodes[iNode]->GetMarker()][BoundNodes[iNode]->GetVertex()]->GetVarCoord();

      for(auto iDim = 0u; iDim < nDim; iDim++){
        localError[iDim] = -displacement[iDim] * VarIncrement;
      }
    }else{
      for(auto iDim = 0u; iDim < nDim; iDim++){
        localError[iDim] = 0; 
      }
    }
  }


  /*--- Resulting displacement from the RBF interpolation is added to the error ---*/

  /*--- Finding contribution of each control node ---*/
  for(auto jNode = 0ul; jNode < nCtrlNodesGlobal; jNode++){

    /*--- Distance between non-selected boundary node and control node ---*/
    auto dist = GeometryToolbox::Distance(nDim, CtrlCoords[jNode *nDim], geometry->nodes->GetCoord(BoundNodes[iNode]->GetIndex()));

    /*--- Evaluation of Radial Basis Function ---*/
    auto rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));

    /*--- Add contribution to error ---*/
    for(auto iDim = 0u; iDim < nDim; iDim++){
      localError[iDim] += rbf*InterpCoeff[jNode*nDim + iDim];
    }
  }
}

void CRadialBasisFunctionInterpolation::SetCorrection(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const vector<unsigned long>& internalNodes){

  /*--- The non-selected control nodes still have a nonzero error once the maximum error falls below the data reduction tolerance. 
          This error is applied as correction and interpolated into the volumetric mesh for internal nodes that fall within the correction radius.
          To evaluate whether an internal node falls within the correction radius an AD tree is constructed of the boundary nodes,
          making it possible to determine the distance to the nearest boundary node. ---*/

  /*--- Construction of the AD tree consisting of the non-selected boundary nodes ---*/

  /*--- Number of non-selected boundary nodes ---*/
  const unsigned long nVertexBound = BoundNodes.size();
  
  /*--- Vector storing the coordinates of the boundary nodes ---*/
  vector<su2double> Coord_bound(nDim*nVertexBound);

  /*--- Vector storing the IDs of the boundary nodes ---*/
  vector<unsigned long> PointIDs(nVertexBound);

  /*--- Correction Radius, equal to maximum error times a prescribed constant ---*/
  const su2double CorrectionRadius = config->GetRBF_DataRedCorrectionFactor()*MaxErrorGlobal; 

  /*--- Storing boundary node information ---*/
  unsigned long i = 0;
  unsigned long j = 0;
  for(auto iVertex = 0ul; iVertex < nVertexBound; iVertex++){
    auto iNode = BoundNodes[iVertex]->GetIndex();
    PointIDs[i++] = iVertex;
    for(auto iDim = 0u; iDim < nDim; iDim++){
      Coord_bound[j++] = geometry->nodes->GetCoord(iNode, iDim);
    }
  }

  /*--- Construction of AD tree ---*/
  CADTPointsOnlyClass BoundADT(nDim, nVertexBound, Coord_bound.data(), PointIDs.data(), true);

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

    /*--- Get error of nearest node ---*/
    auto err = BoundNodes[pointID]->GetError();

    /*--- evaluate RBF ---*/
    auto rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, CorrectionRadius, dist));

    /*--- Apply correction to the internal node ---*/
    for(auto iDim = 0u; iDim < nDim; iDim++){
      geometry->nodes->AddCoord(internalNodes[iNode], iDim, -rbf*err[iDim]);
    }
  }

  /*--- Applying the correction to the non-selected boundary nodes ---*/
  for(auto iNode = 0ul; iNode < BoundNodes.size(); iNode++){
    auto err =  BoundNodes[iNode]->GetError();
    for(auto iDim = 0u; iDim < nDim; iDim++){
      geometry->nodes->AddCoord(BoundNodes[iNode]->GetIndex(), iDim, -err[iDim]);
    }
  }
}


void CRadialBasisFunctionInterpolation::AddControlNode(unsigned long maxErrorNode){
  /*--- Addition of node to the reduced set of control nodes ---*/
  ReducedControlNodes.push_back(move(BoundNodes[maxErrorNode]));

  /*--- Removal of node among the non-selected boundary nodes ---*/
  BoundNodes.erase(BoundNodes.begin()+maxErrorNode);
}


void CRadialBasisFunctionInterpolation::Get_nCtrlNodesGlobal(){
  /*--- Determining the global number of control nodes ---*/

  /*--- Local number of control nodes ---*/
  auto local_nControlNodes = ControlNodes->size();
  
  /*--- Summation of local number of control nodes ---*/
  SU2_MPI::Allreduce(&local_nControlNodes, &nCtrlNodesGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
}