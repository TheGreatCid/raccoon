#pragma once
#include <unordered_map>

class Qp_Mapping
{
public:
  // Static lookup tables accessible without constructing Qp_Mapping
  static const std::unordered_map<int, int> TET4_2nd_lookup;
  static const std::unordered_map<int, int> TET4_4th_lookup;
  static const std::unordered_map<int, int> TET10_4th_lookup;
  static const std::unordered_map<int, int> HEX8_3rd_lookup;
};