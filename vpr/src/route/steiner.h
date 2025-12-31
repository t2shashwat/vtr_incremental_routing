#pragma once

/** @file Functions for Steiner minimum tree guided PathFinder routing. */

#include <string>
#include <vector>
#include <map>
#include <utility>
#include <ranges>

#include "netlist.h"
#include "vtr_strong_id.h"
#include "vpr_context.h"
#include "globals.h"

#include "flute.h"
// UTILS
int calculate_manhattan_distance(std::string sb1_id, std::string sb2_id);

std::pair<int, int> make_pair_from_str_id(const std::string& str_id);

std::string make_str_id(int x, int y);

std::string make_str_id(int x, int y, int pin_id);
// end of UTILS


struct CLB {
    int x;
    int y;
    int pin_id; // -1 if not a pin, 1 if driver, 2-n if sink
    std::string str_id;

    CLB() : x(-1), y(-1), pin_id(-1) {
        this->str_id = std::to_string(x)+"_"+std::to_string(y)+"_"+std::to_string(pin_id);
    }
    CLB(int x, int y) : x(x), y(y), pin_id(-1) {
        this->str_id = std::to_string(x)+"_"+std::to_string(y)+"_"+std::to_string(pin_id);
    }
    CLB(int x, int y, int pin_id) : x(x), y(y), pin_id(pin_id) {
        this->str_id = std::to_string(x)+"_"+std::to_string(y)+"_"+std::to_string(pin_id);
    }
};

struct SB {
    int x;
    int y;
    int pin_id; // -1 if not a pin, 1 if driver, 2-n if sink
    std::string str_id;

    SB() : x(-1), y(-1), pin_id(-1) {
        this->str_id = std::to_string(x)+"_"+std::to_string(y)+"_"+std::to_string(pin_id);
    }
    SB(int x, int y) : x(x), y(y), pin_id(-1) {
        this->str_id = std::to_string(x)+"_"+std::to_string(y)+"_"+std::to_string(pin_id);
    }
    SB(int x, int y, int pin_id) : x(x), y(y), pin_id(pin_id) {
        this->str_id = std::to_string(x)+"_"+std::to_string(y)+"_"+std::to_string(pin_id);
    }
};

struct RSMTEdgeSet {
    std::unordered_map<std::string, std::unordered_map<std::string, bool>> edge_map;

    void add_edge(int x1, int y1, int x2, int y2) {
        // When adding edges, we need to make them bidirectional, because we do not know
        // which direction the signal will be coming from while constructing tree's edges
        edge_map[make_str_id(x1, y1)][make_str_id(x2, y2)] = true;
        edge_map[make_str_id(x2, y2)][make_str_id(x1, y1)] = true;
    }
    bool has_edge(std::string str_id_src, std::string str_id_sink) {
        bool has_edge = edge_map[str_id_src][str_id_sink];
        if (!has_edge)
            edge_map[str_id_src].erase(str_id_sink);
        return has_edge;
    }
};


class Steiner {
private:
    std::unordered_map<std::string, CLB> clb_map;

    std::unordered_map<std::string, SB> sb_map;

    // UTILS
    bool has_node_at_loc(int x, int y) {
        for (const auto& sb : this->sb_map) {
            if (sb.second.x == x && sb.second.y == y)
                return true;
        }
        return false;
    }

    int loc_to_pin_id(int x, int y) const {
        int pin_id = -1;
        for (const auto& sb : this->sb_map) {
            if (sb.second.x == x && sb.second.y == y) {
                pin_id = sb.second.pin_id;
            }
        }
        return pin_id;
    }

    std::vector<int> loc_to_sink_ids(int x, int y) {
        std::vector<int> sink_ids;
        for (const auto& sb : this->sb_map) {
            if (sb.second.pin_id > 0 && sb.second.x == x && sb.second.y == y)
                sink_ids.push_back(sb.second.pin_id);
        }

        return sink_ids;
    }

public:
    ClusterNetId net_id;
    
    int source_x, source_y;
    
    RSMTEdgeSet sb_edges;    
    
    Steiner(const Netlist<ClusterBlockId, ClusterPortId, ClusterPinId, ClusterNetId>& net_list,
        ClusterNetId net_id,
        const vtr::vector_map<ClusterBlockId, t_block_loc>& block_locs,
        //Flute::FluteState *flute1, 
	std::ofstream& fluteOutfile, bool dump_raw_flute_trees);

    void add_clb(int x, int y) {
        clb_map[make_str_id(x, y, -1)] = CLB(x, y);
    }
    void add_pin_clb(int x, int y, int pin_id) {
        clb_map[make_str_id(x, y, pin_id)] = CLB(x, y, pin_id);
    }
    void add_sb(int x, int y) {
        sb_map[make_str_id(x, y, -1)] = SB(x, y);
    }
    void add_pin_sb(int x, int y, int pin_id) {
        sb_map[make_str_id(x, y, pin_id)] = SB(x, y, pin_id);
    }

    void get_clb_coordinates(const Netlist<ClusterBlockId, ClusterPortId, ClusterPinId, ClusterNetId>& net_list,
        ClusterNetId net_id,
        const vtr::vector_map<ClusterBlockId, t_block_loc>& blocks_locs);

    void map_clbs_to_sbs();

    //bool detect_diagonal(Flute::FluteState *flute1);
    //void build_sb_rsmt_flute(Flute::FluteState *flute1, std::ofstream& fluteOutfile);
    
    void build_sb_rsmt_post_process_flute();
    
    void make_tree_graph_directed(std::string source_sb_id, std::set<std::string>& visited, std::set<std::string>& edges_to_remove, const std::string& parent_sb_id = "");
    
    void prune_subtree_from(const std::string& node_id, std::set<std::string>& edges_to_remove);

    void create_coarsened_regions(std::string source_sb_id);
    
    std::vector<int> traverse_rsmt_build_paths(std::string source_sb_id);
    
    void map_net_sink_to_dnodes(int source_x, int source_y, int sink_x, int sink_y, std::vector<int>& sink_ids);

    std::tuple<std::unordered_map<int, std::vector<Corridor>>, std::unordered_map<int, std::vector<unsigned short>>, std::unordered_map<int, bool>> build_corridor_list_per_connection(std::string source_sb_id) const;

    void compute_dependency_graph_sink_order(std::string source_sb_id, std::unordered_map<size_t, std::unordered_map<int, int>>& sink_order);
    
    void create_dependency_graph(std::string source_sb_id, std::unordered_map<std::string, std::vector<std::string>>& dependency_graph);

    std::pair<int, std::string> create_dependency_subgraph(std::string source_sb_id, std::unordered_map<std::string, std::vector<std::string>>& dependency_graph);

    void compute_sink_order(std::string source_sb_id, std::unordered_map<std::string, std::vector<std::string>>& dependency_graph, std::unordered_map<size_t, std::unordered_map<int, int>>& sink_order);
    
    int compute_sink_order(std::string source_sb_id, std::unordered_map<std::string, std::vector<std::string>>& dependency_graph, std::unordered_map<size_t, std::unordered_map<int, int>>& sink_order, int i);
};

std::set<size_t> load_nets_to_skip();

void load_gnode_maps(SteinerContext& steiner_ctx);

void create_connections_per_dnode_file(SteinerContext& steiner_ctx);

void create_sink_order_file(const ClusteringContext& cluster_ctx, const vtr::vector_map<ClusterBlockId, t_block_loc>& blocks_locs);

void create_branch_node_map_file(SteinerContext& steiner_ctx);

void steiner_pre_processing(bool create_steiner_constraints, bool compute_dependency_graph_sink_orders, bool dump_raw_flute_trees, std::string global_routing_algorithm, int global_channel_capacity);
