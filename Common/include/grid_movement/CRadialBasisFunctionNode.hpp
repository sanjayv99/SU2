/*!
 * \file CRadialBasisFunctionNode.hpp
 * \brief Declaration of the RBF node class that stores nodal information for the RBF interpolation.
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

#include "../../../Common/include/geometry/CGeometry.hpp"

/*!
 * \class CRadialBasisFunctionNode
 * \brief Class for defining a Radial Basis Function Node
 * \author F. van Steen
 */

class CRadialBasisFunctionNode{
 protected:
  unsigned long idx;          /*!< \brief Global index. */
  unsigned short marker_idx;  /*!< \brief Marker index. */
  unsigned long vertex_idx;   /*!< \brief Vertex index. */
    
  su2double error[3];         /*!< \brief Nodal data reduction error; */ // TODO - set initially to zero? 
  bool control = false; // TODO add description
  bool periodic = false;
  bool DomainVertex = false;
  
  su2double cylindrical_coord[3];

  su2double new_coord[3];
  unsigned long nearest_node;
  int nearest_rank;
  su2double nearest_normal[3];
  su2double nearest_normal2[3];
  su2double nearest_coord[3];
  su2double var_coord[3];
  // su2double init_coord[3];
  su2double il_height;

  

  
  public:

  enum class NODETYPE {
    DISPLACED,
    EDGE,
    SURFACE,
    IL_WALL,
    IL_EDGE
  };

  NODETYPE NodeType;

  /*!
  * \brief Constructor of the class.
  * \param[in] idx_val - Local node index.
  * \param[in] marker_val - Local marker index.
  * \param[in] vertex_val - Local vertex index.
  */
  CRadialBasisFunctionNode(unsigned long idx_val, unsigned short marker_val, unsigned long vertex_val);

  /*!
  * \brief Returns local global index.
  * \return Local node index.
  */
  inline unsigned long GetIndex() const /* // TODO -  const for more functions?*/ {return idx;}

  /*!
  * \brief Returns local vertex index.
  * \return Local vertex index.
  */
  inline unsigned long GetVertex() const {return vertex_idx;}

  /*!
  * \brief Returns local marker index.
  * \return Local marker index.
  */
  inline unsigned short GetMarker() const {return marker_idx;}

  /*!
  * \brief Set error of the RBF node.
  * \param val_error - Nodal error.
  * \param nDim - Number of dimensions.
  */

  inline void SetError(const su2double* val_error, unsigned short nDim) {
    for (auto iDim = 0u; iDim < nDim; iDim++) error[iDim] = val_error[iDim];
  }

  inline void AddError(const su2double val_error, unsigned short iDim) {
    error[iDim] += val_error;
  }

  /*!
  * \brief Get nodal error.
  * \return Nodal error.
  */
  inline su2double* GetError(){ return error;}

  // inline void setNodetype(const string& type) { nodetype = type;}
  inline void SetNodeType(NODETYPE type) {NodeType = type;}

  // inline const string& getNodetype(){return nodetype;}
  inline NODETYPE GetNodeType() const { return NodeType; }

  inline void setControl(){control = true;}

  inline bool GetControl(){return control;}
  inline void resetControl(){control = false;}
  
  inline bool GetPeriodic(){return periodic;}
  
  inline void SetCylCoord(const su2double* cyl_coord, unsigned short nDim) {
    for (auto iDim = 0u; iDim < nDim; iDim++) cylindrical_coord[iDim] = cyl_coord[iDim];
  }

  inline const su2double* GetCylCoord() const { return cylindrical_coord;}


  inline void SetNewCoord(const su2double* newCoord, unsigned short nDim) {
    for (auto iDim = 0u; iDim < nDim; iDim++) new_coord[iDim] = newCoord[iDim];
  }

  inline const su2double* GetNewCoord() const { return new_coord;}
  
  inline void AddNewCoord(const su2double val_coord, unsigned short iDim) {
    new_coord[iDim] += val_coord;
  }

  inline void SetNearestNode(const unsigned long ID) {
    nearest_node = ID;
  }

  inline unsigned long GetNearestNode() const {
    return nearest_node;
  }

  inline void SetNearestRank( const int ID) {
    nearest_rank = ID;
  }
  inline int GetNearestRank() const {
    return nearest_rank;
  }

  inline void SetNearestNormal(const su2double* val_normal, unsigned short nDim) {
    for (auto iDim = 0u; iDim < nDim; iDim++) nearest_normal[iDim] = val_normal[iDim];
  }

  inline const su2double* GetNearestNormal() const { return nearest_normal;}

  inline void SetNearestNormal2(const su2double* val_normal, unsigned short nDim) {
    for (auto iDim = 0u; iDim < nDim; iDim++) nearest_normal2[iDim] = val_normal[iDim];
  }

  inline const su2double* GetNearestNormal2() const { return nearest_normal2;}

  inline void SetNearestCoord(const su2double* val_coord, unsigned short nDim) {
    for (auto iDim = 0u; iDim < nDim; iDim++) nearest_coord[iDim] = val_coord[iDim];
  }

  inline const su2double* GetNearestCoord() const { return nearest_coord;}

  inline void SetVarCoord(const su2double* val_coord, unsigned short nDim) {
    for (auto iDim = 0u; iDim < nDim; iDim++) var_coord[iDim] = val_coord[iDim];
  }

  inline const su2double* GetVarCoord() const { return var_coord;}

  
  // inline void SetInitCoord(const su2double* initCoord, unsigned short nDim) {
  //   for (auto iDim = 0u; iDim < nDim; iDim++) init_coord[iDim] = initCoord[iDim];
  // }

  // inline const su2double* GetInitCoord() const { return init_coord;}

  inline void SetILHeight(su2double val) {
    il_height = val;
  }

  inline su2double GetILHeight() const { return il_height;}
};