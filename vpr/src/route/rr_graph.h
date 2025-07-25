#ifndef RR_GRAPH_H
#define RR_GRAPH_H

/* Include track buffers or not. Track buffers isolate the tracks from the
 * input connection block. However, they are difficult to lay out in practice,
 * and so are not currently used in commercial architectures. */
#define INCLUDE_TRACK_BUFFERS false

#include "device_grid.h"
#include "vpr_types.h"
#include "rr_graph_type.h"
#include "describe_rr_node.h"

/* Warnings about the routing graph that can be returned.
 * This is to avoid output messages during a value sweep */
enum {
    RR_GRAPH_NO_WARN = 0x00,
    RR_GRAPH_WARN_FC_CLIPPED = 0x01,
    RR_GRAPH_WARN_CHAN_X_WIDTH_CHANGED = 0x02,
    RR_GRAPH_WARN_CHAN_Y_WIDTH_CHANGED = 0x03
};
// (PARSA) Luka, 2025: This was moved to route_common.h
// std::pair<ClusterNetId, int> get_netid_and_sinkid(std::string connection_id);
void create_rr_graph(const t_graph_type graph_type,
                     const std::vector<t_physical_tile_type>& block_types,
                     const DeviceGrid& grid,
                     t_chan_width nodes_per_chan,
                     const int num_arch_switches,
                     t_det_routing_arch* det_routing_arch,
                     const std::vector<t_segment_inf>& segment_inf,
                     const t_router_opts& router_opts,
                     const t_direct_inf* directs,
                     const int num_directs,
                     int* Warnings,
                     bool is_flat = false);

void free_rr_graph();

t_rr_switch_inf create_rr_switch_from_arch_switch(int arch_switch_idx,
                                                  const float R_minW_nmos,
                                                  const float R_minW_pmos);
// Sets the spec for the rr_switch based on the arch switch
void load_rr_switch_from_arch_switch(int arch_switch_idx,
                                     int rr_switch_idx,
                                     int fanin,
                                     const float R_minW_nmos,
                                     const float R_minW_pmos);

t_non_configurable_rr_sets identify_non_configurable_rr_sets();

#endif