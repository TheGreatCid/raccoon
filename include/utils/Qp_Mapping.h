//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "MooseError.h"

#include <string>
#include <unordered_map>

/**
 * Quadrature point numbering and variable naming of the recovery files, in which material
 * properties are stored as one elemental variable per quadrature point. The quadrature points of
 * a file are numbered from 1 in the SIERRA order, which differs from the MOOSE order.
 */
namespace QpMapping
{
/// Element type and quadrature order
enum class Element
{
  TET4_2nd,
  TET4_4th,
  TET10_4th,
  HEX8_3rd
};

inline constexpr const char * ELEMENT_ENUM_DEFINITION = "TET4_2nd TET4_4th TET10_4th HEX8_3rd";

/**
 * The quadrature point map of an element type, with 1-based indices on both sides: from the
 * file's numbering to the MOOSE numbering, or (reversed) from the MOOSE numbering to the file's.
 * Also gives the number of quadrature points.
 */
const std::unordered_map<int, int> *
getLookup(Element elem, unsigned int & qpnum, bool reversed = false);

/// Look up a (1-based) quadrature point in a map from getLookup()
unsigned int getQP(unsigned int qp, const std::unordered_map<int, int> * lookup);

/**
 * Name of the variable holding a scalar property at a quadrature point of the file (1-based).
 * The quadrature point number has two digits when the element has ten or more of them.
 */
std::string variableName(const std::string & property, unsigned int file_qp, unsigned int n_qps);

/**
 * Name of the variable holding component (i, j) of a tensor property at a quadrature point of the
 * file (1-based). The property name is truncated so that the variable name fits the 32 characters
 * of an Exodus name.
 */
std::string variableName(const std::string & property,
                         unsigned int i,
                         unsigned int j,
                         unsigned int file_qp,
                         unsigned int n_qps);
}
