#include "Qp_Mapping.h"

// Moose to SIERRA
const std::unordered_map<int, int> Qp_Mapping::TET4_2nd_lookup = {{1, 4}, {2, 1}, {3, 2}, {4, 3}};

const std::unordered_map<int, int> Qp_Mapping::TET4_4th_lookup = {
    {1, 1}, {2, 5}, {3, 4}, {4, 3}, {5, 2}};

const std::unordered_map<int, int> Qp_Mapping::TET10_4th_lookup = {
    {1, 1}, {2, 5}, {3, 4}, {4, 3}, {5, 2}, {6, 6}, {7, 11}, {8, 10}, {9, 7}, {10, 9}, {11, 8}};

const std::unordered_map<int, int> Qp_Mapping::HEX8_3rd_lookup = {
    {1, 8}, {2, 7}, {3, 6}, {4, 5}, {5, 4}, {6, 3}, {7, 2}, {8, 1}};

// SIERRA to MOOSE

const std::unordered_map<int, int> Qp_Mapping::TET4_2nd_lookup_rev = {
    {4, 1}, {1, 2}, {2, 3}, {3, 4}};

const std::unordered_map<int, int> Qp_Mapping::TET4_4th_lookup_rev = {
    {1, 1}, {5, 2}, {4, 3}, {3, 4}, {2, 5}};

const std::unordered_map<int, int> Qp_Mapping::TET10_4th_lookup_rev = {
    {1, 1}, {5, 2}, {4, 3}, {3, 4}, {2, 5}, {6, 6}, {11, 7}, {10, 8}, {7, 9}, {9, 10}, {8, 11}};

const std::unordered_map<int, int> Qp_Mapping::HEX8_3rd_lookup_rev = {
    {8, 1}, {7, 2}, {6, 3}, {5, 4}, {4, 5}, {3, 6}, {2, 7}, {1, 8}};
