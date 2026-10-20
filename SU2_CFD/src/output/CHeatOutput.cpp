/*!
 * \file CHeatOutput.cpp
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


#include "../../include/output/CHeatOutput.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

namespace {
/*--- CHT heat balance outputs: the solid has a CHT interface with the Robin coupling. Uses the config file markers,
      the output fields are set before the geometry. ---*/
bool CHTBalanceSolid(const CConfig* config) {
  if (config->GetKind_CHT_Coupling() != CHT_COUPLING::DIRECT_TEMPERATURE_ROBIN_HEATFLUX) return false;
  for (auto iMarker = 0u; iMarker < config->GetnMarker_CfgFile(); iMarker++)
    if (config->GetMarker_CfgFile_KindBC(config->GetMarker_CfgFile_TagBound(iMarker)) == CHT_WALL_INTERFACE)
      return true;
  return false;
}
}  // namespace

CHeatOutput::CHeatOutput(CConfig *config, unsigned short nDim) : CFVMOutput(config, nDim, false) {

  /*--- Set the default history fields if nothing is set in the config file ---*/

  if (nRequestedHistoryFields == 0){
    RequestCommonHistory(config->GetTime_Domain());
    nRequestedHistoryFields = requestedHistoryFields.size();
  }
  if (nRequestedScreenFields == 0){
    requestedScreenFields.emplace_back("OUTER_ITER");
    requestedScreenFields.emplace_back("INNER_ITER");
    requestedScreenFields.emplace_back("RMS_TEMPERATURE");
    nRequestedScreenFields = requestedScreenFields.size();
  }
  if (nRequestedVolumeFields == 0){
    requestedVolumeFields.emplace_back("COORDINATES");
    requestedVolumeFields.emplace_back("SOLUTION");
    nRequestedVolumeFields = requestedVolumeFields.size();
  }

  stringstream ss;
  ss << "Zone " << config->GetiZone() << " (Solid Heat)";
  multiZoneHeaderString = ss.str();

  /*--- Set the volume filename --- */

  volumeFilename = config->GetVolume_FileName();

  /*--- Set the surface filename --- */

  surfaceFilename = config->GetSurfCoeff_FileName();

  /*--- Set the restart filename --- */

  restartFilename = config->GetRestart_FileName();

  /*--- Set the default convergence field --- */

  if (convFields.empty() ) convFields.emplace_back("RMS_TEMPERATURE");

}

void CHeatOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) {

  CSolver* heat_solver = solver[HEAT_SOL];

  SetHistoryOutputValue("TOTAL_HEATFLUX", heat_solver->GetTotal_HeatFlux());
  SetHistoryOutputValue("MAXIMUM_HEATFLUX", heat_solver->GetTotal_MaxHeatFlux());
  SetHistoryOutputValue("AVG_TEMPERATURE", heat_solver->GetTotal_AvgTemperature());
  SetHistoryOutputValue("RMS_TEMPERATURE", log10(heat_solver->GetRes_RMS(0)));
  SetHistoryOutputValue("MAX_TEMPERATURE", log10(heat_solver->GetRes_Max(0)));
  if (multiZone)
    SetHistoryOutputValue("BGS_TEMPERATURE", log10(heat_solver->GetRes_BGS(0)));

  SetHistoryOutputValue("LINSOL_ITER", heat_solver->GetIterLinSolver());
  SetHistoryOutputValue("LINSOL_RESIDUAL", log10(heat_solver->GetResLinSolver()));
  SetHistoryOutputValue("CFL_NUMBER", config->GetCFL(MESH_0));

  if (CHTBalanceSolid(config)) {
    su2double cht[4];
    heat_solver->GetCHTHeatBalance(geometry, config, cht);
    SetHistoryOutputValue("CHT_Q_IF", cht[0]);
    if (config->GetStreamwise_Periodic_ExpDecay()) {
      SetHistoryOutputValue("CHT_Q_LAMBDA", cht[1]);
      SetHistoryOutputValue("CHT_Q_SINK", cht[2]);
      SetHistoryOutputValue("CHT_Q_DECAY", cht[3]);
      SetHistoryOutputValue("CHT_ZONE_BALANCE", cht[3] + cht[1] + cht[0] - cht[2]);
    }
  }

  ComputeSimpleCustomOutputs(config);

  /*--- Keep this as last, since it uses the history values that were set. ---*/
  SetCustomAndComboObjectives(HEAT_SOL, config, solver);
}


void CHeatOutput::SetHistoryOutputFields(CConfig *config){

  AddHistoryOutput("LINSOL_ITER", "LinSolIter", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver.");
  AddHistoryOutput("LINSOL_RESIDUAL", "LinSolRes", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver.");

  AddHistoryOutput("RMS_TEMPERATURE", "rms[T]", ScreenOutputFormat::FIXED, "RMS_RES", "Root mean square residual of the temperature", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("MAX_TEMPERATURE", "max[T]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the temperature", HistoryFieldType::RESIDUAL);
  AddHistoryOutput("BGS_TEMPERATURE", "bgs[T]", ScreenOutputFormat::FIXED, "BGS_RES", "Block-Gauss-Seidel residual of the temperature", HistoryFieldType::RESIDUAL);

  AddHistoryOutput("TOTAL_HEATFLUX", "HF", ScreenOutputFormat::SCIENTIFIC, "HEAT", "Total heatflux on all surfaces defined in MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
  AddHistoryOutput("MAXIMUM_HEATFLUX", "MaxHF", ScreenOutputFormat::SCIENTIFIC, "HEAT", "Maximum heatflux on all surfaces defined in MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
  AddHistoryOutput("AVG_TEMPERATURE", "AvgTemp", ScreenOutputFormat::SCIENTIFIC, "HEAT", "Average temperature on all surfaces defined in MARKER_MONITORING", HistoryFieldType::COEFFICIENT);
  AddHistoryOutput("CFL_NUMBER", "CFL number", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current value of the CFL number");

  if (CHTBalanceSolid(config)) {
    AddHistoryOutput("CHT_Q_IF", "Q_IF", ScreenOutputFormat::SCIENTIFIC, "CHT_BALANCE", "Heat into the solid through the CHT interfaces by the Robin exchange term", HistoryFieldType::COEFFICIENT);
    if (config->GetStreamwise_Periodic_ExpDecay()) {
      AddHistoryOutput("CHT_Q_LAMBDA", "Q_Lambda", ScreenOutputFormat::SCIENTIFIC, "CHT_BALANCE", "Heat into the solid by the lambda term of the interface flux (decaying formulation)", HistoryFieldType::COEFFICIENT);
      AddHistoryOutput("CHT_Q_SINK", "Q_Sink", ScreenOutputFormat::SCIENTIFIC, "CHT_BALANCE", "Heat out of the solid through isothermal walls (decaying formulation)", HistoryFieldType::COEFFICIENT);
      AddHistoryOutput("CHT_Q_DECAY", "Q_Decay", ScreenOutputFormat::SCIENTIFIC, "CHT_BALANCE", "Volume terms of the decaying formulation in the solid", HistoryFieldType::COEFFICIENT);
      AddHistoryOutput("CHT_ZONE_BALANCE", "ZoneBal", ScreenOutputFormat::SCIENTIFIC, "CHT_BALANCE", "Solid heat balance Q_DECAY + Q_LAMBDA + Q_IF - Q_SINK, zero at convergence", HistoryFieldType::COEFFICIENT);
    }
  }

  AddHistoryOutput("COMBO", "ComboObj", ScreenOutputFormat::SCIENTIFIC, "COMBO", "Combined obj. function value.", HistoryFieldType::COEFFICIENT);
}


void CHeatOutput::SetVolumeOutputFields(CConfig *config){

  // Grid coordinates
  AddCoordinates();

  // SOLUTION
  AddVolumeOutput("TEMPERATURE", "Temperature", "SOLUTION", "Temperature");

  // Primitives
  AddVolumeOutput("HEAT_FLUX", "Heat_Flux", "PRIMITIVE", "Heatflux");

  // Residuals
  AddVolumeOutput("RES_TEMPERATURE", "Residual_Temperature", "RESIDUAL", "Residual of the temperature");

  AddCommonFVMOutputs(config);
}


void CHeatOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

  CVariable* Node_Heat = solver[HEAT_SOL]->GetNodes();
  const auto Node_Geo  = geometry->nodes;

  // Grid coordinates
  LoadCoordinates(Node_Geo->GetCoord(iPoint), iPoint);

  // SOLUTION
  SetVolumeOutputValue("TEMPERATURE", iPoint, Node_Heat->GetSolution(iPoint, 0));

  // Residuals
  SetVolumeOutputValue("RES_TEMPERATURE", iPoint, solver[HEAT_SOL]->LinSysRes(iPoint, 0));

  LoadCommonFVMOutputs(config, geometry, iPoint);
}

void CHeatOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){

  if (!config->GetViscous_Wall(iMarker)) return;

  /* Heat flux value at each surface grid node. */
  SetVolumeOutputValue("HEAT_FLUX", iPoint, solver[HEAT_SOL]->GetHeatFlux(iMarker, iVertex));

}

