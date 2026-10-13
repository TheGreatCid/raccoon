//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "Qp_Mapping.h"

#include <array>

namespace QpMapping
{
namespace
{
typedef std::unordered_map<int, int> Lookup;

Lookup
inverse(const Lookup & lookup)
{
  Lookup inv;
  for (const auto & [from, to] : lookup)
    inv[to] = from;
  return inv;
}

// File (SIERRA) quadrature point -> MOOSE quadrature point
const Lookup TET4_2nd_lookup = {{1, 4}, {2, 1}, {3, 2}, {4, 3}};
const Lookup TET4_4th_lookup = {{1, 1}, {2, 5}, {3, 4}, {4, 3}, {5, 2}};
const Lookup TET10_4th_lookup = {
    {1, 1}, {2, 5}, {3, 4}, {4, 3}, {5, 2}, {6, 6}, {7, 11}, {8, 10}, {9, 7}, {10, 9}, {11, 8}};
const Lookup HEX8_3rd_lookup = {{1, 8}, {2, 7}, {3, 6}, {4, 5}, {5, 4}, {6, 3}, {7, 2}, {8, 1}};

// MOOSE quadrature point -> file (SIERRA) quadrature point
const Lookup TET4_2nd_lookup_rev = inverse(TET4_2nd_lookup);
const Lookup TET4_4th_lookup_rev = inverse(TET4_4th_lookup);
const Lookup TET10_4th_lookup_rev = inverse(TET10_4th_lookup);
const Lookup HEX8_3rd_lookup_rev = inverse(HEX8_3rd_lookup);

std::string
qpSuffix(const unsigned int file_qp, const unsigned int n_qps)
{
  return "_" + std::string(n_qps >= 10 && file_qp < 10 ? "0" : "") + std::to_string(file_qp);
}
}

unsigned int
getQP(unsigned int qp, const std::unordered_map<int, int> * lookup)
{
  auto it = lookup->find(qp);
  if (it == lookup->end())
    mooseError("Quadrature point ", qp, " is not in the quadrature point map of the element type");
  return it->second;
}

const std::unordered_map<int, int> *
getLookup(Element elem, unsigned int & qpnum, bool reversed)
{
  switch (elem)
  {
    case Element::TET4_2nd:
      qpnum = 4;
      return reversed ? &TET4_2nd_lookup_rev : &TET4_2nd_lookup;
    case Element::TET4_4th:
      qpnum = 5;
      return reversed ? &TET4_4th_lookup_rev : &TET4_4th_lookup;
    case Element::TET10_4th:
      qpnum = 11;
      return reversed ? &TET10_4th_lookup_rev : &TET10_4th_lookup;
    case Element::HEX8_3rd:
      qpnum = 8;
      return reversed ? &HEX8_3rd_lookup_rev : &HEX8_3rd_lookup;
  }
  mooseError("Unknown element type in QpMapping");
}

std::string
variableName(const std::string & property, const unsigned int file_qp, const unsigned int n_qps)
{
  return property + qpSuffix(file_qp, n_qps);
}

std::string
variableName(const std::string & property,
             const unsigned int i,
             const unsigned int j,
             const unsigned int file_qp,
             const unsigned int n_qps)
{
  static const std::array<std::string, 3> component = {{"x", "y", "z"}};
  // 26 characters of the name + "_xx" + "_NN" = 32, the length of an Exodus name
  return property.substr(0, 26) + "_" + component[i] + component[j] + qpSuffix(file_qp, n_qps);
}
}
