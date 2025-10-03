/*!
 * \file CMeshOutput.cpp
 * \brief Main subroutines for the heat solver output
 * \author R. Sanchez
 * \version 8.0.1 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../include/output/CMeshOutput.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"
// #include "CMeshOutput.hpp"

CMeshOutput::CMeshOutput(CConfig *config, unsigned short nDim) : COutput(config, nDim, false) {

  /*--- Set the default history fields if nothing is set in the config file ---*/

  requestedVolumeFields.emplace_back("COORDINATES");
  nRequestedVolumeFields = requestedVolumeFields.size();

  /*--- Set the volume filename --- */

  volumeFilename = config->GetMesh_Out_FileName();

  /*--- Set the surface filename ---*/

  surfaceFilename = "surface_mesh";

}

CMeshOutput::~CMeshOutput() = default;

void CMeshOutput::SetVolumeOutputFields(CConfig *config){

  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES", "x-component of the coordinate vector");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES", "z-component of the coordinate vector");

  // Mesh quality metrics, computed in CPhysicalGeometry::ComputeMeshQualityStatistics.
  AddVolumeOutput("ORTHOGONALITY", "Orthogonality", "MESH_QUALITY", "Orthogonality Angle (deg.)");
  AddVolumeOutput("ASPECT_RATIO",  "Aspect_Ratio",  "MESH_QUALITY", "CV Face Area Aspect Ratio");
  AddVolumeOutput("VOLUME_RATIO",  "Volume_Ratio",  "MESH_QUALITY", "CV Sub-Volume Ratio");

}

void CMeshOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  CPoint*    Node_Geo  = geometry->nodes;

  SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(iPoint, 0));
  SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(iPoint, 1));
  if (nDim == 3)
    SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(iPoint, 2));

  // Mesh quality metrics
  if (config->GetWrt_MeshQuality()) {
    SetVolumeOutputValue("ORTHOGONALITY", iPoint, geometry->Orthogonality[iPoint]);
    SetVolumeOutputValue("ASPECT_RATIO",  iPoint, geometry->Aspect_Ratio[iPoint]);
    SetVolumeOutputValue("VOLUME_RATIO",  iPoint, geometry->Volume_Ratio[iPoint]);
  }

}


void CMeshOutput::WriteMeshQualityStatistics(CConfig *config, CGeometry *geometry) {
  auto nZone = config->GetnZone();
  auto iZone = config->GetiZone();

  su2double orthoMin = *min_element(geometry->Orthogonality.begin(), geometry->Orthogonality.end()),
  orthoMax = *max_element(geometry->Orthogonality.begin(), geometry->Orthogonality.end()),
  arMin = *min_element(geometry->Aspect_Ratio.begin(), geometry->Aspect_Ratio.end()),
  arMax = *max_element(geometry->Aspect_Ratio.begin(), geometry->Aspect_Ratio.end()),
  vrMin = *min_element(geometry->Volume_Ratio.begin(), geometry->Volume_Ratio.end()),
  vrMax = *max_element(geometry->Volume_Ratio.begin(), geometry->Volume_Ratio.end());

  su2double Global_Ortho_Min, Global_Ortho_Max;
  SU2_MPI::Allreduce(&orthoMin, &Global_Ortho_Min, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&orthoMax, &Global_Ortho_Max, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());

  su2double Global_AR_Min, Global_AR_Max;
  SU2_MPI::Allreduce(&arMin, &Global_AR_Min, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&arMax, &Global_AR_Max, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());

  su2double Global_VR_Min, Global_VR_Max;
  SU2_MPI::Allreduce(&vrMin, &Global_VR_Min, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&vrMax, &Global_VR_Max, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());

  auto fileName = config->GetMesh_Qual_FileName();
  if (nZone > 1) fileName += "_" + PrintingToolbox::to_string(iZone);
  fileName += ".dat";

  ofstream file;
  file.open(fileName);
  file.precision(6);
  file << "Min. Orthogonality [deg]";
  file.width(25);
  file << "Max. Orthogonality [deg]";
  file.width(25);
  file << "Min. Aspect Ratio [deg]";
  file.width(25);
  file << "Max. Aspect Ratio [deg]";
  file.width(25);
  file << "Min. Volume Ratio [deg]";
  file.width(25);
  file << "Max. Volume Ratio [deg]\n";
  
  file << Global_Ortho_Min;
  file.width(25);
  file << Global_Ortho_Max;
  file.width(25);
  file << Global_AR_Min;
  file.width(25);
  file << Global_AR_Max;
  file.width(25);
  file << Global_VR_Min;
  file.width(25);
  file << Global_VR_Max;

  file.close();
}
