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


  // TODO -  Debug MPI 
  {
    if (rank==0){
      int i = 0;
      while(0==i){
        sleep(1);
      }
    }
  }
  

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


  if (config->GetnMarker_Periodic() != 0) SetPeriodicVars(config);

  /*--- Determining the boundary and internal nodes. Setting the control nodes. ---*/ 
  bool surfaceCorrection = false;
  SetBoundNodes(geometry, config, Derivative, surfaceCorrection);
  
  vector<unsigned long> internalNodes; 
  SetInternalNodes(geometry, config, internalNodes); 
  

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

void CRadialBasisFunctionInterpolation::SolveRBF_System(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, bool Derivative, vector<unsigned long>& internalNodes, bool ForwardProjectionDerivative){
  
  if (IsCylindrical) TransformBoundNodesToCylindrical(geometry, config);
  SetCtrlNodes(config);

  /*--- In case of data reduction an iterative greedy algorithm is applied 
          to perform the interpolation with a reduced set of control nodes.
          Otherwise with a full set of control nodes. ---*/
  
  if(config->GetRBF_DataReduction()){
    /*--- Local maximum error node and corresponding maximum error  ---*/
    unsigned long maxErrorNodeLocal;
    su2double maxErrorLocal{0};

    /*--- Obtaining the initial maximum error nodes, which are found based on the maximum applied deformation. */
    if(nCtrlNodesGlobal == 0){  
      GetInitMaxErrorNode(geometry, config, Derivative, maxErrorNodeLocal, maxErrorLocal); 
      SU2_MPI::Allreduce(&maxErrorLocal, &MaxErrorGlobal, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
    }

    /*--- Number of greedy iterations. ---*/
    unsigned short greedyIter = 0;
    
    /*--- Error tolerance for the data reduction tolerance ---*/
    if (dataReductionTolerance == 0.0) {
      dataReductionTolerance = config->GetRBF_DataRedTolerance() * MaxErrorGlobal; 
      cout << "DATA REDUCTION TOLERANCE: " << dataReductionTolerance << endl;
    }
    // const su2double 

    su2passivematrix invInterpMat;

    vector <string> ctrltypes = {"edge"};

    /*--- While the maximum error is above the tolerance, data reduction algorithm is continued. ---*/
    while(((MaxErrorGlobal > dataReductionTolerance)  || greedyIter == 0)){ 
      CtrlTypeVec.resize(1);

      /*--- In case of a nonzero local error, control nodes are added ---*/
      if(maxErrorLocal> 0){
        AddControlNode(maxErrorNodeLocal);
      }

      /*--- Obtaining the global number of control nodes. ---*/
      Get_nCtrlNodesGlobal(config);
      
      /*--- Obtaining the control nodes coordinates and distributing over all processes. ---*/
      SetCtrlNodeCoords(geometry, config);

      /*--- Obtaining the deformation of the control nodes. ---*/
      if (!Derivative){
        SetBoundaryDisplacements(geometry, config);
      }else{
        SetCtrlNodeDerivatives(geometry, config, ForwardProjectionDerivative);
      } 

      /*--- Computation of the (inverse) interpolation matrix. ---*/      
      GetInvInterpMat(geometry, config, type, radius, invInterpMat);
      
      /*--- Obtaining the interpolation coefficients. ---*/
      ComputeInterpCoeffs(invInterpMat);

      
      for (const auto &itype : ctrltypes) {
        
        auto marker = ctrl_node_indices[itype];
        for (auto mark : marker){

          if (std::find(CtrlTypeVec.begin(), CtrlTypeVec.end(), itype) == CtrlTypeVec.end()) {
            CtrlTypeVec.push_back(itype);
          }
                    
          ProjectBoundNodes(geometry, config, type, radius, itype, false, nodes, mark);

          }
        
        /*--- Obtaining the global number of control nodes. ---*/
        Get_nCtrlNodesGlobal(config);
        
        /*--- Obtaining the control nodes coordinates and distributing over all processes. ---*/
        SetCtrlNodeCoords(geometry, config);
  
        /*--- Obtaining the deformation of the control nodes. ---*/
        if (!Derivative){
          SetBoundaryDisplacements(geometry, config);
        }else{
          SetCtrlNodeDerivatives(geometry, config, ForwardProjectionDerivative);
        } 
  
        /*--- Computation of the (inverse) interpolation matrix. ---*/
        
        GetInvInterpMat(geometry, config, type, radius, invInterpMat);
        
        /*--- Obtaining the interpolation coefficients. ---*/
        ComputeInterpCoeffs(invInterpMat);

      }
 
      /*--- Determining the interpolation error, of the non-control boundary nodes. ---*/
      GetInterpError(geometry, config, type, radius, Derivative, ctrltypes, maxErrorNodeLocal, maxErrorLocal); 
      
      SU2_MPI::Allreduce(&maxErrorLocal, &MaxErrorGlobal, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
      
      if(rank == MASTER_NODE) cout << "Greedy iteration: " << greedyIter << ". Max error: " << MaxErrorGlobal << ". Global nr. of ctrl nodes: "  << nCtrlNodesGlobal << "\n" << endl;
      
      
      if (MaxErrorGlobal < dataReductionTolerance && ctrltypes.size() < nDim - 1) {
        // ensuring error is above tolerance:
        MaxErrorGlobal = 2*dataReductionTolerance;
        // including the surface nodes in the interpolation
        ctrltypes.push_back("surface");
      }
      
      greedyIter++;
    } 
    
    /*--- Once the data reduction tolerance is reached the error for the periodic nodes has to be determined. ---*/
    
    GetPeriodicNodeErrors(geometry, config, type, radius, Derivative);

    
    if (Derivative){  
      SetInternalNodeDerivatives(geometry, config, internalNodes, ForwardProjectionDerivative);
      ComputeSensitivity(geometry, type, radius, invInterpMat, internalNodes); 
    }
  }else{

    /*--- Obtaining the interpolation coefficients. ---*/
    GetInterpCoeffs(geometry, config, type, radius, Derivative, internalNodes, ForwardProjectionDerivative);

    vector<string> types = {"edge", "surface"};

    // iter over nodetype
    for (auto iType : types) {

      // check if type is present among local nodes.  
      // if (node_indices[iType].size() > 0) {
      
      // get markers with specific node types. 
      const auto& markers = node_indices[iType];
      for (const auto& iMarker : markers) {  
        ProjectBoundNodes(geometry, config, type, radius, iType, false, nodes, iMarker);
      }
      
      CtrlTypeVec.push_back(iType);
      Get_nCtrlNodesGlobal(config);
      GetInterpCoeffs(geometry, config, type, radius, Derivative, internalNodes, ForwardProjectionDerivative);
    }
  }
}

void CRadialBasisFunctionInterpolation::GetInterpCoeffs(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, bool Derivative, vector<unsigned long>& internalNodes, bool ForwardProjectionDerivative){
  
  /*--- Obtaining the control nodes coordinates and distributing over all processes. ---*/
  SetCtrlNodeCoords(geometry, config);

  /*--- Obtaining the deformation of the control nodes. ---*/
  if (!Derivative){
    SetBoundaryDisplacements(geometry, config);
  }else{
    SetInternalNodeDerivatives(geometry, config, internalNodes, ForwardProjectionDerivative);
  } 

  /*--- Computation of the (inverse) interpolation matrix. ---*/
  su2passivematrix invInterpMat;
  GetInvInterpMat(geometry, config, type, radius, invInterpMat);

  if(!Derivative){
    /*--- Obtaining the interpolation coefficients. ---*/
    ComputeInterpCoeffs(invInterpMat);
  }else{
    ComputeSensitivity(geometry, type, radius, invInterpMat, internalNodes);
  }
}

void CRadialBasisFunctionInterpolation::ComputeSensitivity(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius,su2passivematrix &invInterpMat, vector<unsigned long>& internalNodes){

  vector<su2double> inter_res(nCtrlNodesGlobal * nDim, 0.0);

  //loop over all global control nodes
  for (auto iNode = 0ul; iNode < nCtrlNodesGlobal; iNode++){
    for (auto iDim =0u; iDim < nDim; iDim++){
      
      // loop over the local internal nodes
      for (auto jNode =0ul; jNode < internalNodes.size(); jNode++){

        su2double dist = GetDistance(CtrlCoords[iNode*nDim], geometry->nodes->GetCoord(internalNodes[jNode]));
        su2double rbf_eval = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));
        
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

void CRadialBasisFunctionInterpolation::SetBoundNodes(CGeometry* geometry, CConfig* config, bool Derivative, bool& surfaceCorrection){
  /*--- Storing of the local node, marker and vertex information of the boundary nodes ---*/

  
  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {
    if (!config->GetMarker_All_SendRecv(iMarker) && !config->GetMarker_All_PerBound(iMarker)) {


      if (config->GetMarker_All_Deform_Mesh_Slide(iMarker)){
        for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          auto iNode = geometry->vertex[iMarker][iVertex]->GetNode(); 
          if (!geometry->nodes->GetDomain(iNode)){
            continue;
          }
          
          // counting how many boundaries the node is part of
          unsigned short nVertex = 0;
          
          for (auto jMarker = 0u; jMarker < config->GetnMarker_All(); jMarker++){
            if (geometry->nodes->GetVertex(iNode, jMarker) != -1 && !config->GetMarker_All_SendRecv(jMarker)){
               nVertex++;
            }
          }

          if (nVertex == 1) {
            nodes.push_back(new CRadialBasisFunctionNode(iNode, iMarker, iVertex));

            if (nDim == 2) {
              nodes.back()->setNodetype("edge");
            } else {
              nodes.back()->setNodetype("surface");            
            }
          }


          else if (nVertex == nDim - 1) {

            if(geometry->nodes->GetPeriodicBoundary(iNode)) {
              unsigned short perMarker;
              // iterate over markers
              for (auto jMarker = 0u; jMarker < config->GetnMarker_All(); jMarker++){

                // if node is on marker and marker is periodic
                if(geometry->nodes->GetVertex(iNode, jMarker) != -1 && config->GetMarker_All_PerBound(jMarker)){

                  // if periodic marker is the smaller of the two periodic markers then its saved.
                  auto idx_PeriodicMarker = config->GetMarker_Periodic(config->GetMarker_All_TagBound(jMarker));
                  if(idx_PeriodicMarker < config->GetMarker_Periodic_Donor2(idx_PeriodicMarker)){
                    nodes.push_back(new CRadialBasisFunctionNode(iNode, iMarker, iVertex));
                    nodes.back()->setNodetype("surface");

                  }else{
                    per_nodes.push_back(new CRadialBasisFunctionNode(iNode, iMarker, iVertex));
                    per_nodes.back()->setNodetype("surface");
                  }
                  break;
                }
              }
            } else {
              
              // add check to ensure it's not part of an moving boundary. 
              bool isMoving = false;
              for (auto jMarker = 0u; jMarker < config->GetnMarker_All(); jMarker++){
                if ( config->GetMarker_All_DV(jMarker) && geometry->nodes->GetVertex(iNode, jMarker) != -1) {
                  surfaceCorrection = true;
                  isMoving = true;
                  break;
                }
              }
              
              if (!isMoving) {
                nodes.push_back(new CRadialBasisFunctionNode(iNode, iMarker, iVertex));
                nodes.back()->setNodetype(Derivative ? "displaced" : "edge");
              }
            }
            
          } 
          
          else {
            // if on periodic boundary
            if(geometry->nodes->GetPeriodicBoundary(iNode)){
              unsigned short perMarker;
              // iterate over markers
              for (auto jMarker = 0u; jMarker < config->GetnMarker_All(); jMarker++){

                // if node is on marker and marker is periodic
                if(geometry->nodes->GetVertex(iNode, jMarker) != -1 && config->GetMarker_All_PerBound(jMarker)){

                  // if periodic marker is the smaller of the two periodic markers then its saved.
                  auto idx_PeriodicMarker = config->GetMarker_Periodic(config->GetMarker_All_TagBound(jMarker));
                  if(idx_PeriodicMarker < config->GetMarker_Periodic_Donor2(idx_PeriodicMarker)){
                    nodes.push_back(new CRadialBasisFunctionNode(iNode, iMarker, iVertex));
                    
                    if (!Derivative){
                      nodes.back()->setNodetype("edge"); 
                    } else {
                      nodes.back()->setNodetype("displaced");
                    }
                    
                  }else{
                    per_nodes.push_back(new CRadialBasisFunctionNode(iNode, iMarker, iVertex));
                    if (!Derivative) {
                      per_nodes.back()->setNodetype("edge");
                    } else {
                      per_nodes.back()->setNodetype("displaced");
                    }
                  }
                  break;
                }
              }
            } else {
              nodes.push_back(new CRadialBasisFunctionNode(iNode, iMarker, iVertex));
              nodes.back()->setNodetype("displaced");
            }
          }
        }
      }

      else {
        for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          auto iNode = geometry->vertex[iMarker][iVertex]->GetNode(); 
        
          // if on periodic boundary
          if(geometry->nodes->GetPeriodicBoundary(iNode)){
            unsigned short perMarker;
            // iterate over markers
            for (auto jMarker = 0u; jMarker < config->GetnMarker_All(); jMarker++){

              // if node is on marker and marker is periodic
              if(geometry->nodes->GetVertex(iNode, jMarker) != -1 && config->GetMarker_All_PerBound(jMarker)){

                // if periodic marker is the smaller of the two periodic markers then its saved.
                auto idx_PeriodicMarker = config->GetMarker_Periodic(config->GetMarker_All_TagBound(jMarker));
                if(idx_PeriodicMarker < config->GetMarker_Periodic_Donor2(idx_PeriodicMarker)){
                  nodes.push_back(new CRadialBasisFunctionNode(iNode, iMarker, iVertex));
                  nodes.back()->setNodetype("displaced");
                }else{
                  per_nodes.push_back(new CRadialBasisFunctionNode(iNode, iMarker, iVertex));
                  per_nodes.back()->setNodetype("displaced");
                }
                break;
              }
            }
          } else {
            if (geometry->nodes->GetDomain(iNode)){
              nodes.push_back(new CRadialBasisFunctionNode(iNode, iMarker, iVertex)); 
              nodes.back()->setNodetype("displaced");
            }
          }
        }
      }
    }    
  }

  for (auto iMarker = 0u; iMarker < config->GetnMarker_CfgFile(); iMarker++) {
    if (config->GetMarker_CfgFile_Deform_Mesh_Slide(config->GetMarker_CfgFile_TagBound(iMarker))) {
      node_indices["edge"][iMarker]; 
    }
  }

  stable_sort(nodes.begin(), nodes.end(), HasSmallerIndex);
  nodes.resize(distance(nodes.begin(), unique(nodes.begin(), nodes.end(), HasEqualIndex)));

  stable_sort(per_nodes.begin(), per_nodes.end(), HasSmallerIndex);
  per_nodes.resize(distance(per_nodes.begin(), unique(per_nodes.begin(), per_nodes.end(), HasEqualIndex)));

  for(unsigned long x = 0ul; x < nodes.size(); x++){
    auto type = nodes[x]->getNodetype();
    auto marker = nodes[x]->GetMarker(); // this is the local marker index
    auto tag = config->GetMarker_All_TagBound(marker);
    auto global_marker = config->GetMarker_CfgFile_TagBound(tag);
    node_indices[type][global_marker].push_back(x);
  }

  

  for (unsigned long x = 0ul; x < per_nodes.size(); x++) {
    auto type = per_nodes[x]->getNodetype();
    auto marker = per_nodes[x]->GetMarker();
    
    per_node_indices[type][marker].push_back(x);
  }
  // TODO -  debug output 
  
  ofstream out("edge"+to_string(rank)+".txt");
  auto idx = node_indices["edge"];
  for (auto& pair : idx) {
      const unsigned short& marker = pair.first;      // marker (key)
      const std::vector<unsigned long>& indices = pair.second;  // indices (value)

      for (auto i : indices) {
          out << geometry->nodes->GetGlobalIndex(nodes[i]->GetIndex())  /*<< "\t" << nodes[i]->GetMarker() << "\t" << marker*/ << endl;
      }
  }
  out.close();

  ofstream out2("surface.txt");
  auto idx2 = node_indices["surface"];
  for (auto& pair : idx2) {
      const unsigned short& marker = pair.first;      // marker (key)
      const std::vector<unsigned long>& indices = pair.second;  // indices (value)

      for (auto i : indices) {
          out2 << nodes[i]->GetIndex() << endl;
      }
  }
  out2.close();


  ofstream out3("disp"+to_string(rank)+".txt");
  auto idx3 = node_indices["displaced"];
  for (auto& pair : idx3) {
      const unsigned short& marker = pair.first;      // marker (key)
      const std::vector<unsigned long>& indices = pair.second;  // indices (value)

      for (auto i : indices) {
          out3 << geometry->nodes->GetGlobalIndex(nodes[i]->GetIndex()) << endl;
      }
  }
  out3.close();

  ofstream out4("per_nodes"+to_string(rank)+".txt");
  for (auto i : per_nodes){
    out4 << geometry->nodes->GetGlobalIndex(i->GetIndex()) << /*"\t" <<  i->getNodetype()  <<*/  endl;
  }
  out4.close();
}

void CRadialBasisFunctionInterpolation::SetCtrlNodes(CConfig* config){
  
  CtrlTypeVec.resize(1);
  CtrlTypeVec[0] = "displaced";

  /*--- Obtaining the total number of control nodes. ---*/
  Get_nCtrlNodesGlobal(config);
};

void CRadialBasisFunctionInterpolation::GetInvInterpMat(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, su2passivematrix& invInterpMat) {
  
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
}

void CRadialBasisFunctionInterpolation::GetInterpMat_parallel(CGeometry* geometry, const RADIAL_BASIS& type, const su2double radius, CSymmetricMatrix& interpMat){
  /*--- Initialization of the interpolation matrix ---*/
  interpMat.Initialize(nCtrlNodesGlobal);

  // Total number of elems in the lower triangular matrix
  su2double N_lowerTriangle = (nCtrlNodesGlobal*(nCtrlNodesGlobal+1))/2; // NOTE this must be a su2double, otherwise the ceil operation does not work correctly.

  // global rbf evaluations are stored in rbf_vals_all
  vector<su2double> rbf_vals_all(N_lowerTriangle);
  // Average number of elements per process
  su2double N_perProcess = ceil(N_lowerTriangle/size);

  // For balancing the number of elements per process, the start and end rows are determined.
  // The number of elements in a lower triangle matrix of size n x n is given by: T_n = n(n+1)/2
  // Starting and ending row is determined by solving this equations for the number of elements (rank*N_perProcess) up to that row:
  // row(row+1)/2 = rank*N_perProcess. Using quadratic formula results in: row = (-1 + sqrt(1+8*N_perProcess*rank))/2.
  // Ceil is used to obtain integer numbers.

  unsigned long start_row = ceil((-1 + sqrt(1+8*N_perProcess*(rank))) / 2);
  unsigned long end_row = ceil((-1 + sqrt(1+8*N_perProcess*(rank+1))) / 2);
  if (start_row > nCtrlNodesGlobal) start_row = nCtrlNodesGlobal;
  if (end_row > nCtrlNodesGlobal)  end_row = nCtrlNodesGlobal;  
  
  // Number of elements to be evaluated
  int nr_elems = (end_row*(end_row+1) - start_row*(start_row+1))/2;
  // Finding RBF evaluations
  vector<su2double> rbf_vals(nr_elems);
  unsigned long cnt = 0;

  for (auto row_i = start_row; row_i < end_row; row_i++){
    for (auto col_i = 0ul; col_i <= row_i; col_i++){
      su2double dist = GetDistance(CtrlCoords[row_i*nDim], CtrlCoords[col_i*nDim]);
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
      su2double dist = GetDistance(CtrlCoords[iNode*nDim], CtrlCoords[jNode*nDim]);
      /*--- Evaluation of RBF ---*/
      interpMat(iNode, jNode) = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));
    }
  }
}


void CRadialBasisFunctionInterpolation::SetBoundaryDisplacements(CGeometry* geometry, CConfig* config){
  
  /* --- Initialization of the deformation vector ---*/
  CtrlNodeDeformation.resize(nCtrlNodesLocal*nDim, 0.0); 

  /*--- If requested (no by default) impose the surface deflections in
    increments and solve the grid deformation with
    successive small deformations. ---*/
  // const su2double VarIncrement = 1.0 / ((su2double)config->GetGridDef_Nonlinear_Iter());

  SU2_COMPONENT Kind_SU2 = config->GetKind_SU2();

  unsigned long idx = 0;

  for (auto type : CtrlTypeVec) {
    auto test = config->GetRBF_DataReduction()
    ? ConvertToVectorMap(ctrl_node_indices[type])
    : node_indices[type];
    
    for (auto mark : test){
      for (auto idx_i : mark.second) {
        auto iMarker = nodes[idx_i]->GetMarker();
        if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_CFD)) ||
          ((config->GetMarker_All_DV(iMarker) == YES || config->GetMarker_All_Deform_Mesh_Slide(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_DEF)) ||
          ((config->GetDirectDiff() == D_DESIGN) && (Kind_SU2 == SU2_COMPONENT::SU2_CFD) &&
          (config->GetMarker_All_DV(iMarker) == YES)) /*NOTE: This feature has not been tested for RBF interpolation*/ ||
          ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_COMPONENT::SU2_DOT))) {

          su2double varCoord[nDim];
          GetNodalDeformation(geometry, config, nodes[idx_i] ,varCoord);

          auto baseIndex = idx * nDim;
          for (auto iDim = 0u; iDim < nDim; iDim++){
            CtrlNodeDeformation[baseIndex + iDim] = SU2_TYPE::GetValue(varCoord[iDim] /** VarIncrement*/);
          }
          
        } else {
          for (auto iDim = 0u; iDim < nDim; iDim++) {
            CtrlNodeDeformation[idx * nDim + iDim] = 0.0;
          }
        }
        idx++;
      }  
    }
  }

}

void CRadialBasisFunctionInterpolation::GetNodalDeformation(CGeometry* geometry, CConfig* config, CRadialBasisFunctionNode* iNode, su2double* varCoord) {

  /*--- If requested (no by default) impose the surface deflections in
    increments and solve the grid deformation with
    successive small deformations. ---*/
  // su2double VarIncrement = 1;
  // if (iNode->getNodetype() == "displaced"){
  const su2double VarIncrement = 1.0 / ((su2double)config->GetGridDef_Nonlinear_Iter());
  // } 
  auto varCoordCartOrig = geometry->vertex[iNode->GetMarker()][iNode->GetVertex()]->GetVarCoord();
  su2double varCoordCart[nDim];
  for (auto iDim = 0u; iDim < nDim; iDim++) {
    varCoordCart[iDim] = varCoordCartOrig[iDim]*VarIncrement;
  }

  if (IsCylindrical) {
    auto coord = geometry->nodes->GetCoord(iNode->GetIndex());
    CartDispToCyl(coord, varCoordCart, varCoord);
  } else{
    for (auto iDim = 0u; iDim < nDim; iDim++){
      varCoord[iDim] = varCoordCart[iDim];
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
    CtrlNodeDeformation.resize(nCtrlNodesLocal * nDim, 0.0);

    unsigned long idx = 0;
    for (auto type : CtrlTypeVec) {
      auto test = config->GetRBF_DataReduction()
      ? ConvertToVectorMap(ctrl_node_indices[type])
      : node_indices[type];

      // loop through markers
      for (auto mark : test){
        // loop through indices of each marker
        for (auto idx_i : mark.second) {
          auto iMarker = nodes[idx_i]->GetMarker();

          for (auto iDim = 0u; iDim < nDim; iDim++) {
            CtrlNodeDeformation[idx *nDim + iDim] = SU2_TYPE::GetValue(geometry->GetSensitivity(nodes[idx_i]->GetIndex(), iDim));
          }
          idx++;
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
    } else if (geometry->nodes->GetPeriodicBoundary(iNode)) {
      if (find_if (nodes.begin(), nodes.end(), [&](CRadialBasisFunctionNode* i){return i->GetIndex() == iNode;}) == nodes.end()) {
        internalNodes.push_back(iNode);
      }  
    }
  }  
  
  /*--- sorting of the local indices ---*/
  sort(internalNodes.begin(), internalNodes.end());

  /*--- Obtaining unique set of internal nodes ---*/
  internalNodes.resize(std::distance(internalNodes.begin(), unique(internalNodes.begin(), internalNodes.end())));

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
          if (find_if (nodes.begin(), nodes.end(), [&](CRadialBasisFunctionNode* i){return i->GetIndex() == iNode;}) == nodes.end()) {
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


  unsigned long nCtrlNode = nCtrlNodesLocal;
  unsigned long nCtrlNodes[size];
  SU2_MPI::Allgather(&nCtrlNode, 1, MPI_UNSIGNED_LONG, nCtrlNodes, 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm()); 
  
  unsigned long start_idx = 0;
  for (auto iProc = 0; iProc < rank; iProc++){
    start_idx += nCtrlNodes[iProc];
  }

  for (auto iNode = 0; iNode < nCtrlNodesLocal; iNode++){
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

void CRadialBasisFunctionInterpolation::UpdateGridCoord(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, const vector<unsigned long>& internalNodes, bool surfaceCorrection){

  if(rank == MASTER_NODE){
    cout << "updating the grid coordinates" << endl;
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
      ? ConvertToVectorMap(ctrl_node_indices[type])
      : node_indices[type];
    
      for (auto mark : test){
        for (auto idx_i : mark.second) {
          auto iMarker = nodes[idx_i]->GetMarker();

          if (config->GetSolid_Wall(iMarker) || (config->GetMarker_All_DV(iMarker) == YES)) {
            auto iPoint = nodes[idx_i]->GetIndex();
            for (auto iDim = 0u; iDim < nDim; iDim++) {
              // summation of current sensitivity and the computed update
              su2double sens_new =  geometry->GetSensitivity(iPoint, iDim) + sensitivity_update[idx * nDim + iDim];
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
    ? ConvertToVectorMap(ctrl_node_indices[type])
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

void CRadialBasisFunctionInterpolation::GetInitMaxErrorNode(CGeometry* geometry, CConfig* config, bool Derivative,  unsigned long& maxErrorNodeLocal, su2double& maxErrorLocal){

  /*--- Set max error to zero ---*/
  maxErrorLocal = 0.0;

  auto indices = node_indices["displaced"]; 
  for (auto idx : indices) {
    /*--- Loop over the nodes ---*/
    for (auto iNode : idx.second){

      su2double normSquaredDeformation = 0;
      /*--- Compute to squared norm of the deformation ---*/
      if (Derivative){
        auto iPoint = nodes[iNode]->GetIndex();
        for (auto iDim = 0u; iDim < nDim; iDim++){
          passivedouble sens_i = SU2_TYPE::GetValue(geometry->GetSensitivity(iPoint, iDim));
          normSquaredDeformation += sens_i * sens_i;
        }
      } else {
        normSquaredDeformation = GeometryToolbox::SquaredNorm(nDim, geometry->vertex[nodes[iNode]->GetMarker()][nodes[iNode]->GetVertex()]->GetVarCoord());    
      }
      
      /*--- In case squared norm deformation is larger than the error, update the error ---*/
      if(normSquaredDeformation > maxErrorLocal){
        maxErrorLocal = normSquaredDeformation;
        maxErrorNodeLocal = iNode;
      }
    }
  }

  /*--- Account for the possibility of applying the deformation in multiple steps ---*/
  maxErrorLocal = sqrt(maxErrorLocal) / ((su2double)config->GetGridDef_Nonlinear_Iter());
}


void CRadialBasisFunctionInterpolation::SetCtrlNodeCoords(CGeometry* geometry, CConfig* config){
  /*--- The coordinates of all control nodes are made available on all processes ---*/

  /*--- resizing the matrix containing the global control node coordinates ---*/
  CtrlCoords.resize(nCtrlNodesGlobal*nDim);
  
  /*--- Array containing the local control node coordinates ---*/ 
  su2double localCoords[nDim*nCtrlNodesLocal];
  ofstream ctrl("ctrl.txt");
  /*--- Storing local control node coordinates ---*/
  unsigned long idx = 0;
  for (auto type : CtrlTypeVec) {

    auto indices = config->GetRBF_DataReduction() ? ConvertToVectorMap(ctrl_node_indices[type]) : node_indices[type];

    for (auto ctrl_idx : indices) {
      for (auto idx_i : ctrl_idx.second) {
        ctrl << nodes[idx_i]->GetIndex() << endl;
        su2double* coord = IsCylindrical ? nodes[idx_i]->GetCylCoord() : geometry->nodes->GetCoord(nodes[idx_i]->GetIndex());

        for (auto iDim = 0u; iDim < nDim; iDim++) {
          localCoords[idx++] = coord[iDim];
        }
      }
    }
  }
  ctrl.close();

  /*--- Gathering local control node coordinate sizes on all processes. ---*/
  int LocalCoordsSizes[size];
  int localCoordsSize = nDim * nCtrlNodesLocal;
  SU2_MPI::Allgather(&localCoordsSize, 1, MPI_INT, LocalCoordsSizes, 1, MPI_INT, SU2_MPI::GetComm()); 

  /*--- Array containing the starting indices for the allgatherv operation */
  int disps[SU2_MPI::GetSize()] = {0};    

  for(auto iProc = 1; iProc < SU2_MPI::GetSize(); iProc++){
    disps[iProc] = disps[iProc-1]+LocalCoordsSizes[iProc-1];
  }
  
  /*--- Distributing global control node coordinates among all processes ---*/
  SU2_MPI::Allgatherv(&localCoords, localCoordsSize, MPI_DOUBLE, CtrlCoords.data(), LocalCoordsSizes, disps, MPI_DOUBLE, SU2_MPI::GetComm()); 
};


void CRadialBasisFunctionInterpolation::GetInterpError(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, bool Derivative, vector<string>& ctrltypes,  unsigned long& maxErrorNodeLocal, su2double& maxErrorLocal){
  /*--- Array containing the local error ---*/
  su2double localError[nDim];

  /*--- Magnitude of the local maximum error ---*/
  maxErrorLocal = 0.0;

  /*--- Loop over non-selected boundary nodes ---*/
  string errortype; 
  auto nodes_disp = node_indices["displaced"];
  for (auto mark : nodes_disp) {
    for (auto idx : mark.second) {
      if(!nodes[idx]->GetControl()){

        /*--- Compute nodal error ---*/
        GetNodalError(geometry, config, type, radius, nodes[idx], Derivative, localError);

        /*--- Setting error ---*/
        nodes[idx]->SetError(localError, nDim);
        
        /*--- Compute error magnitude and update local maximum error if necessary ---*/
        su2double errorMagnitude = GeometryToolbox::Norm(nDim, localError);
        if(errorMagnitude > maxErrorLocal){
          maxErrorLocal = errorMagnitude;
          maxErrorNodeLocal = idx;
          errortype = "displaced";
        }
      }
    }
  }

  
  for (auto i_type : ctrltypes) {
    
    auto nodes_edge = node_indices[i_type];

    for (const auto& mark : nodes_edge) {
      
      ProjectBoundNodes(geometry, config, type, radius, i_type, true, nodes, mark);

      const auto& indices = mark.second;  
      for (const auto idx : indices) {
        if(!nodes[idx]->GetControl()){
          auto localError = nodes[idx]->GetError();

          su2double errorMagnitude = GeometryToolbox::Norm(nDim, localError);
          if(errorMagnitude > maxErrorLocal){
            maxErrorLocal = errorMagnitude;
            maxErrorNodeLocal = idx;
            errortype = nodes[idx]->getNodetype();
          }
        }
      }
    }
  }
}

void CRadialBasisFunctionInterpolation::GetNodalError(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, CRadialBasisFunctionNode* iNode, bool Derivative, su2double* localError){ 
  
  /*--- If requested (no by default) impose the surface deflections in increments ---*/
  const su2double VarIncrement = 1.0 / ((su2double)config->GetGridDef_Nonlinear_Iter());
  
  /*--- If node is part of a moving boundary then the error is defined as the difference
           between the found and prescribed displacements. Thus, here the displacement is substracted from the error ---*/

  if (Derivative) {
    for (auto iDim =0u; iDim < nDim; iDim++) {
      localError[iDim] = -geometry->GetSensitivity(iNode->GetIndex(), iDim) * VarIncrement;
    }
  }else{  
    auto disp_true = geometry->vertex[iNode->GetMarker()][iNode->GetVertex()]->GetVarCoord();
    for(auto iDim = 0u; iDim < nDim; iDim++){
      localError[iDim] = -disp_true[iDim] * VarIncrement;
    }
  }


  su2double disp_interp[nDim] = {0.0};
  /*--- Resulting displacement from the RBF interpolation is added to the error ---*/ 

  auto targetCoord = geometry->nodes->GetCoord(iNode->GetIndex());

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

  if (IsCylindrical) {
    CylDispToCart(geometry->nodes->GetCoord(iNode->GetIndex()),  disp_interp);
  }

  for (auto iDim = 0u; iDim < nDim; iDim++) {
    localError[iDim] += disp_interp[iDim];
  }
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
  auto disp_map = node_indices["displaced"];
  for (auto& pair : disp_map){
     nDisplaced += pair.second.size(); // conservative size;
  }

  vector<su2double> coord_moving(nDisplaced*nDim);
  vector<unsigned long> ID_moving(nDisplaced);

  unsigned long cnt = 0;
  for(auto iVertex = 0ul; iVertex < nodes.size(); iVertex++){
    if (nodes[iVertex]->getNodetype() == "displaced" /*&& !nodes[iVertex]->GetControl()*/) {

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



  auto surf_map = node_indices["surface"];
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

void CRadialBasisFunctionInterpolation::AddControlNode(unsigned long maxErrorNode){ 
  nodes[maxErrorNode]->setControl();

  vector<su2double> zero(nDim, 0.0);
  nodes[maxErrorNode]->SetError(zero.data(), nDim);
  // control_node_indices.insert(maxErrorNode);
  // ctrl_nodes_type[nodes[maxErrorNode]->getNodetype()].insert(maxErrorNode);
  ctrl_node_indices[nodes[maxErrorNode]->getNodetype()][nodes[maxErrorNode]->GetMarker()].insert(maxErrorNode);
}


void CRadialBasisFunctionInterpolation::Get_nCtrlNodesGlobal(CConfig* config){
  /*--- Determining the global number of control nodes ---*/

  /*--- Local number of control nodes ---*/
  nCtrlNodesLocal = 0;
  
  if (config->GetRBF_DataReduction()) {
    for (auto type : CtrlTypeVec) {
      auto markers = ctrl_node_indices[type];
      for (auto mark : markers) {
        auto indices = mark.second;
        nCtrlNodesLocal += indices.size();
      } 
    } 
  } else {
    for (auto type : CtrlTypeVec) { 
      auto markers = node_indices[type];
      for (auto mark : markers) {
        auto indices = mark.second;
        nCtrlNodesLocal += indices.size();
      }
    }
  }

  /*--- Summation of local number of control nodes ---*/
  SU2_MPI::Allreduce(&nCtrlNodesLocal, &nCtrlNodesGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
}


void CRadialBasisFunctionInterpolation::SetPeriodicVars(CConfig* config){

  /*--- Counter number of periodic angles encountered ---*/
  unsigned short rotationalAngleCnt = 0;
  
  /*--- Loop over periodic markers and find rotational periodic parameters ---*/
  for (auto iMarker = 0u; iMarker < config->GetnMarker_Periodic(); iMarker++) {

    auto periodicAngles = config->GetPeriodic_RotAngles(iMarker);
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
          PeriodicAngle = periodicAngles[iDim];
          IsCylindrical = true;
        } else if (RotationalAxis != iDim) {
          SU2_MPI::Error("Only a single rotationally periodic angle can be provided with MARKER_PERIODIC.", CURRENT_FUNCTION);
        } else if (fabs(periodicAngles[iDim]) != fabs(PeriodicAngle)) {
          SU2_MPI::Error("Two different values of periodic angles detected in MARKER_PERIODIC.", CURRENT_FUNCTION);
        }
        rotationalAngleCnt++;
      }
    }
  }

  /*--- Loop over periodic markers to find translational periodic parameters ---*/
  for (auto iMarker = 0u; iMarker < config->GetnMarker_Periodic(); iMarker++) {
    
    auto periodicTranslation = config->GetPeriodic_Translation(iMarker);
    for (auto iDim = 0u; iDim < 3; iDim++) {

        if (periodicTranslation[iDim] != 0.0) {

          /*--- For 2D, periodicity cannot be rotational and translational.
                  For 3D domain with rotational periodicity, the domain can only be translationally peridic along the rotational periodic axis. ---*/
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

void CRadialBasisFunctionInterpolation::TransformBoundNodesToCylindrical(CGeometry* geometry, CConfig* config){
  su2double cyl_coord[nDim];

  auto transform = [&](const vector<CRadialBasisFunctionNode*>& nodes) {
    for (auto* iNode : nodes) {
      CartToCyl(geometry->nodes->GetCoord(iNode->GetIndex()), cyl_coord);
      iNode->SetCylCoord(cyl_coord, nDim);
    }
  };

  transform(nodes);
  transform(per_nodes);
}

void CRadialBasisFunctionInterpolation::CartDispToCyl(const su2double* init_coord_cart, const su2double* var_coord_cart, su2double* var_coord_cyl) {
  // new Carthesian coord
  su2double new_coord[nDim];
  for (auto iDim = 0u; iDim < nDim; iDim++){
    new_coord[iDim] = init_coord_cart[iDim] + var_coord_cart[iDim];
  }

  su2double init_coord_cyl[nDim];
  CartToCyl(init_coord_cart, init_coord_cyl);
  CartToCyl(new_coord, var_coord_cyl);

  for (auto iDim = 0u; iDim < nDim; iDim++){
    var_coord_cyl[iDim] -= init_coord_cyl[iDim]; 
  }
}

void CRadialBasisFunctionInterpolation::CylDispToCart(const su2double* init_coord_cart, su2double* var_coord) {
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

  // subtract initial coordinate position to find the delta
  for (auto iDim = 0u; iDim < nDim; iDim++){
    var_coord[iDim] -= init_coord_cart[iDim]; 
  }
}


unique_ptr<CADTPointsOnlyClass> CRadialBasisFunctionInterpolation::CreateADT(CGeometry* geometry, const string& type, const short marker, bool periodic) {
  
  auto indices = node_indices[type][marker];
  // TODO -  
  cout << rank << " " << marker << " " << indices.size() << endl;

 
    
  unsigned long size = indices.size();
  int nr_markers = periodic ? 3 : 1; // 3: original + `2 periodic images
  
  /*--- Vector storing the coordinates of the boundary nodes ---*/ 
  vector<su2double> Coord_bound(size * nDim * nr_markers);

  /*--- Vector storing the IDs of the boundary nodes ---*/
  vector<unsigned long> PointIDs(size * nr_markers);

   if (marker == -1) {
    cout << "SIZE: " << size << endl;
    return unique_ptr<CADTPointsOnlyClass>( new CADTPointsOnlyClass(nDim, size*3, Coord_bound.data(), PointIDs.data(), true));
  }

  auto type_idx = node_indices[type];

  auto nVertex = geometry->GetnVertex(marker); // it breaks here

  for (auto idx_i = 0ul; idx_i < indices.size(); idx_i++) {

    auto iNode = nodes[indices[idx_i]];

    for (auto i = 0u; i < nr_markers; i++) {
      PointIDs[ idx_i + size* i] = iNode->GetVertex() + i * nVertex;
    }
    
    su2double* coord = IsCylindrical ? iNode->GetCylCoord() : geometry->nodes->GetCoord(iNode->GetIndex());
    // 0 is original image
    // 1 is positive
    // 2 is negative
    
    // loop through dimensions 
    su2double angular_periodic_shift[3] = {0.0, +PeriodicAngle, -PeriodicAngle};
    for (auto iDim = 0u; iDim < nDim; iDim++){
      su2double translation_periodic_shift[3] = {0.0, +PeriodicLength[iDim], -PeriodicLength[iDim]};
      for (auto image = 0; image < nr_markers; image++) {
        // NOTE for angular periodicity the rotational axis is always the theta/second coordinate
        Coord_bound[(idx_i + image * size)*nDim + iDim] = coord[iDim] + translation_periodic_shift[image] + ((iDim == 1) ? angular_periodic_shift[image] : 0.0); //original 
      }
    }
  }

  /*--- Construction of AD tree ---*/
  // TODO -  
  // cout << rank << " grootte: " << size*nr_markers << endl;
  return unique_ptr<CADTPointsOnlyClass>( new CADTPointsOnlyClass(nDim, size*3, Coord_bound.data(), PointIDs.data(), true));
}


template <typename T>
void CRadialBasisFunctionInterpolation::ProjectBoundNodes(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, const string& nodetype, bool SetError, vector<CRadialBasisFunctionNode*>& target_nodes, T marker) {
    
  auto markerIdx = marker.first;
  auto tag = config->GetMarker_CfgFile_TagBound(markerIdx);
  auto localIdx = config->GetMarker_All_TagBound(tag);
  
  auto nodeIndices = ToVector(marker.second); // returns reference to vector or creates vector from set.
  
  auto BoundADT = CreateADT(geometry, nodetype, markerIdx, config->GetnMarker_Periodic() != 0);

  
  SU2_MPI::Barrier(SU2_MPI::GetComm()); // TODO -  still required? 
  
  su2double projection[nDim];

  for (auto x : nodeIndices) {

    // When calculating error skip nodes that are selected as control nodes. 
    if( SetError && target_nodes[x]->GetControl()) continue;
    
    auto node = target_nodes[x]->GetIndex();

    auto coord = IsCylindrical ? target_nodes[x]->GetCylCoord() : geometry->nodes->GetCoord(node); 

    su2double new_coord[nDim];
    ApplyRBF(coord, type, radius, new_coord); // RBF displaced coord is stored in new_coord.
    if (rank == 0 ) {
    cout << node << endl;
    }
    ApplyProjection(geometry, config, localIdx, BoundADT.get(), nodetype, new_coord, projection);  // This should return the projection  

    if (SetError) {      

      if(IsCylindrical){
        su2double temp[nDim];
        CylToCart(new_coord, temp);
        CylDispToCart(temp, projection);
      }

      target_nodes[x]->SetError(projection, nDim);
    } else {
      
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        new_coord[iDim] += - projection[iDim] - coord[iDim]; // contains coordinate variation
      }
      
      if(IsCylindrical){
        CylDispToCart(geometry->nodes->GetCoord(node), new_coord);
      }

      const su2double iNonLinear_iter = ((su2double)config->GetGridDef_Nonlinear_Iter());
      for (auto iDim = 0u; iDim < nDim; iDim++){
        new_coord[iDim] = iNonLinear_iter * new_coord[iDim];
      }
    
      geometry->vertex[localIdx][target_nodes[x]->GetVertex()]->SetVarCoord(new_coord);
    }
  }
}

void CRadialBasisFunctionInterpolation::ApplyRBF(const su2double* coord, const RADIAL_BASIS& type, const su2double radius, su2double* new_coord) {
  std::copy(coord, coord + nDim, new_coord);

  
  for(auto jNode = 0ul; jNode < nCtrlNodesGlobal; jNode++){

    /*--- Determine distance between considered internal and control node ---*/
    su2double dist = GetDistance(CtrlCoords[jNode*nDim], coord);
    
    /*--- Evaluate RBF based on distance ---*/
    auto rbf = SU2_TYPE::GetValue(CRadialBasisFunction::Get_RadialBasisValue(type, radius, dist));

    /*--- Add contribution to new coordinates -- -*/
    for(auto iDim = 0u; iDim < nDim; iDim++){
      new_coord[iDim] += rbf*InterpCoeff[jNode*nDim+iDim];
    } 
  }
}

void CRadialBasisFunctionInterpolation::ApplyProjection(CGeometry* geometry, CConfig* config, unsigned short iMarker, CADTPointsOnlyClass* BoundADT, const string& nodetype, su2double* new_coord, su2double* projection) {
  unsigned long pointID;
  su2double dist;
  int rankID;
  BoundADT->DetermineNearestNode(new_coord, dist, pointID, rankID);
  
  // retrieve nearest node and nearest image from pointID
  auto nVertex = geometry->GetnVertex(iMarker);
  auto nearestNode = pointID%nVertex; 
  auto nearestImg = pointID / nVertex;
  
  
  auto closestNode = geometry->vertex[iMarker][nearestNode]->GetNode();
  
  su2double closest_point_coord[nDim];
  auto cart_coord = geometry->nodes->GetCoord(closestNode);

  if (IsCylindrical) {
    // auto cart_coord = geometry->nodes->GetCoord(closestNode);
    CartToCyl(cart_coord, closest_point_coord);
  } else{
    copy(cart_coord, cart_coord + nDim, closest_point_coord);
  }

  
  if (nearestImg == 1) {
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      closest_point_coord[iDim] += PeriodicLength[iDim];
    }
  }
  else if(nearestImg == 2) {
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      closest_point_coord[iDim] -= PeriodicLength[iDim];
    }
  }
  
  auto test_var = geometry->GetnVertex(iMarker);
  if (test_var < nearestNode){
    cout << "warning" << endl;
  }


  auto cart_normal = geometry->vertex[iMarker][nearestNode]->GetNormal();
  su2double normal[nDim];
  if (IsCylindrical){
    CartDispToCyl(geometry->vertex[iMarker][nearestNode]->GetCoord(), cart_normal, normal);
  } else{
    copy(cart_normal, cart_normal + nDim, normal);
  }


  // magnitude first normal
  auto mag_normal = GeometryToolbox::Norm(nDim, normal);

  for (auto iDim = 0u; iDim < nDim; iDim++){
    if (abs(normal[iDim])/mag_normal < 0.017){
      normal[iDim] = 0;
    }
  }
 // TODO -  I don't think the issue is with the fact that the normal of the periodic nodes is calculated slightly off. 
 // since the node that is placed wrongly, for that node the nearest node is not the periodic one.
 // so what else could cause this issue, since the normal seems fine for the other nodes. Is it part of send/receive markers and updated twice or something. Else weird that when more ranks are used, it seems like it occures more often at that boundary.... 
  if (rank == MASTER_NODE) { cout << "node: " << geometry->nodes->GetGlobalIndex(closestNode) << " normal: " << normal[0] << " " << normal[1] << " " << "nearest rank: " << rankID <<  endl;}

  su2double dist_vec[nDim];
  GeometryToolbox::Distance(nDim, new_coord, closest_point_coord, dist_vec);

  // get dot product
  auto dot_product = GeometryToolbox::DotProduct(nDim, normal, dist_vec);

  // get projection
  for(auto iDim = 0u; iDim < nDim; iDim++){ 
    projection[iDim] = dot_product*normal[iDim]/pow(mag_normal,2);
  }

  // part for the double projection of some edges 
  if (nodetype=="edge" && nDim == 3){

    // find the second marker that is not periodic
    unsigned short mark = 0; 
    for( auto i = 0u; i < config->GetnMarker_All(); i++){
      if(geometry->nodes->GetVertex(closestNode, i) != -1 && i != iMarker && !config->GetMarker_All_PerBound(i)){
        break;
      } 
      mark++;    
    }

    auto iNode = geometry->nodes->GetVertex(closestNode, mark);

    cart_normal = geometry->vertex[mark][iNode]->GetNormal();
    su2double normal2[nDim];
    if (IsCylindrical){
      CartDispToCyl(geometry->vertex[iMarker][nearestNode]->GetCoord(), cart_normal, normal2);
    } else{
      copy(cart_normal, cart_normal + nDim, normal2);
    }


    auto mag_normal_test = GeometryToolbox::Norm(nDim, normal2);

    for (auto iDim = 0u; iDim < nDim; iDim++){
      if (abs(normal2[iDim])/mag_normal_test < 0.017){
        normal2[iDim] = 0;
      }
    }

    // get dot product
    auto dot_product2 = GeometryToolbox::DotProduct(nDim, normal2, dist_vec);

    // magnitude second normal
    auto mag_normal2 = GeometryToolbox::Norm(nDim, normal2);

    // get projection
    for(auto iDim = 0u; iDim < nDim; iDim++){ 
      projection[iDim] += dot_product2*normal2[iDim]/pow(mag_normal2,2);
    }
  }
}

void CRadialBasisFunctionInterpolation::GetPeriodicNodeErrors(CGeometry* geometry, CConfig* config, const RADIAL_BASIS& type, const su2double radius, bool Derivative) {
  
  const auto& markers = per_node_indices["displaced"];
  for (const auto& mark : markers) {
      for (auto idx : mark.second) {
          su2double localError[nDim];
          GetNodalError(geometry, config, type, radius, per_nodes[idx], Derivative, localError);
          per_nodes[idx]->SetError(localError, nDim);
      }
  }

  for (auto i = 1u; i < CtrlTypeVec.size(); i++) {
    const auto& iType = CtrlTypeVec[i];
    const auto& markers = per_node_indices[iType];
    for (const auto& mark : markers) {
        ProjectBoundNodes(geometry, config, type, radius, iType, true, per_nodes, mark);
    }
  }
}