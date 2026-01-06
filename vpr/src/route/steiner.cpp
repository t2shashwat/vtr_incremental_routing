#include "steiner.h"
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <vector>
#include <regex>
#include <unordered_map>
#include <string>

#include "vpr_types.h"
#include "place_macro.h"
#include "vpr_context.h"
#include <queue>
// For purposes of creating RSMTs
#include "flute.h"
#include "FastRoute.h"
#include <chrono>
using namespace FastRoute;

// For parallelization attempts
// #include <omp.h>

Steiner::Steiner(
    const Netlist<ClusterBlockId, ClusterPortId, ClusterPinId, ClusterNetId>& net_list,
        ClusterNetId net_id,
    const vtr::vector_map<ClusterBlockId, t_block_loc>& blocks_locs,
    //Flute::FluteState *flute1, 
    std::ofstream& fluteOutfile, bool dump_raw_flute_trees, std::string global_router_algorithm) {
    // VTR_LOG("----------------------------------------\nSTEINER CONSTRUCTOR CALLED FOR %zu\n---------------------------------------\n", size_t(net_id));
    
    // Get the coordinates of the net's CLBs
    this->get_clb_coordinates(net_list, net_id, blocks_locs);

    // Select SB coordinates
    this->map_clbs_to_sbs();
    
    // Find the source to start the traversal from
    SB source_sb;
    for (const auto& sb : sb_map) {
        if (sb.second.pin_id == 0) {
            source_sb = SB(sb.second.x, sb.second.y, sb.second.pin_id); // Source SB
            break;
        }
    }
    this->source_x = source_sb.x;
    this->source_y = source_sb.y;
    
    // Build SB RSMT
    bool diagonal_detected = true;
    //diagonal_detected = this->detect_diagonal(flute1);
    if (this->sb_map.size()>1 && diagonal_detected == true && global_router_algorithm == "FLUTE" && dump_raw_flute_trees == true) {
        this->build_sb_rsmt_flute(fluteOutfile);
    }
    //
    
    // commenting this out only temporarily as the flute used by FastRoute does not support Flute_type datatype
    /*if (dump_raw_flute_trees == true) {
    	this->build_sb_rsmt_flute(flute1, fluteOutfile);
    }*/

    if (dump_raw_flute_trees == false) {
        this->build_sb_rsmt_post_process_flute();
    }

    // Make tree directed

    /*std::set<std::string> visited;
    std::set<std::string> edges_to_remove;
    // Traverse the tree in the correct direction 
    this->make_tree_graph_directed(make_str_id(source_sb.x, source_sb.y), visited, edges_to_remove);
    visited.clear();
    
    for (const auto& edge : edges_to_remove) {
        std::stringstream ss(edge);
        std::string str_id_src, str_id_sink;
        std::getline(ss, str_id_src, '-');
        std::getline(ss, str_id_sink, '-');
        this->sb_edges.edge_map[str_id_src].erase(str_id_sink);
    }
    edges_to_remove.clear();*/

// FOR DEBUGGING PURPOSES: Print all edges in the RSMT so that you can plot it.
    // for (const auto& edge : this->sb_edges.edge_map) {
    //     for (const auto& sink_edge : edge.second) {
    //         if (sink_edge.second) { // Only print true edges
    //             VTR_LOG("Edge: %s -> %s\n", edge.first.c_str(), sink_edge.first.c_str());
    //         }
    //     }
    // }

}

void Steiner::get_clb_coordinates(
    const Netlist<ClusterBlockId, ClusterPortId, ClusterPinId, ClusterNetId>& net_list,
    ClusterNetId net_id,
    const vtr::vector_map<ClusterBlockId, t_block_loc>& blocks_locs) {
    
    this->net_id = net_id;

    auto net_pins = net_list.net_pins(net_id);

    int pin_count = 0;

    // Inspired by load_net_rr_terminals()
    for (auto pin_id : net_pins) {
        auto block_id = net_list.pin_block(pin_id);
        int x = blocks_locs[block_id].loc.x;
        int y = blocks_locs[block_id].loc.y;

        auto type = physical_tile_type(block_id);

        int phys_pin = tile_pin_index(pin_id);

        int iclass = type->pin_class[phys_pin];
        
        this->add_pin_clb(x, y, pin_count);
        pin_count++;
    }
    return;
}

void Steiner::map_clbs_to_sbs() {
    for (auto& clb : this->clb_map) {
        this->add_pin_sb(clb.second.x-1, clb.second.y, clb.second.pin_id);
    }
}
/*bool Steiner::detect_diagonal(Flute::FluteState *flute1) {
    Flute::Tree flutetree;
    int d=0;
    int n = this->sb_map.size();
    int x[n], y[n];

    for (const auto& sb : this->sb_map) {
        x[d] = sb.second.x;
        y[d] = sb.second.y;

        d++;
    }
    
    flutetree = Flute::flute(flute1, d, x, y, FLUTE_ACCURACY);
    //flutetree = Flute::flute(d, x, y, FLUTE_ACCURACY);
    for (int i = 0; i<2*flutetree.deg-2; i++) {
        int o_x, o_y, i_x, i_y;
        o_x = flutetree.branch[i].x;
        o_y = flutetree.branch[i].y;
        i_x = flutetree.branch[flutetree.branch[i].n].x;
        i_y = flutetree.branch[flutetree.branch[i].n].y;

        if (i_x != o_x && i_y != o_y) {
    		Flute::free_tree(flute1, flutetree);
		return true;
        } 
    }

    Flute::free_tree(flute1, flutetree);
    return false;
}*/

void Steiner::build_sb_rsmt_flute(std::ofstream& fluteOutfile) {
    //Flute::Tree flutetree;
    int d=0;
    int n = this->sb_map.size();
    int x[n], y[n];

    fluteOutfile << "Net id: " << size_t(this->net_id) << "\n";
    fluteOutfile << "Source at: (" << this->source_x << ", " << this->source_y << ")\n";

    for (const auto& sb : this->sb_map) {
        x[d] = sb.second.x;
        y[d] = sb.second.y;
	if (sb.second.pin_id != 0){
            //VTR_LOG(" Sink %d: (%d, %d)\n", d, x[d], y[d]);
	    fluteOutfile << " Sink " << d << ": (" << x[d] << ", " << y[d] << ")\n";
	}
        d++;
    }
    
    //flutetree = Flute::flute(flute1, d, x, y, FLUTE_ACCURACY);
    Flute::Tree flutetree = Flute::flute(d, x, y, FLUTE_ACCURACY);

    //flutetree = Flute::flute(d, x, y, FLUTE_ACCURACY);
    for (int i = 0; i<2*flutetree.deg-2; i++) {
        int o_x, o_y, i_x, i_y;
        o_x = flutetree.branch[i].x;
        o_y = flutetree.branch[i].y;
        i_x = flutetree.branch[flutetree.branch[i].n].x;
        i_y = flutetree.branch[flutetree.branch[i].n].y;

        // Add SBs (with pin_id == -1) if they are not already present
        if (!this->has_node_at_loc(o_x, o_y)) {
            this->add_sb(o_x, o_y);
        }
        if (!has_node_at_loc(i_x, i_y)) {
            this->add_sb(i_x, i_y);
        }
        // Add corner SBs if FLUTE output is not a horizontal/vertical line
        // Add edge/s
        if (i_x != o_x && i_y != o_y) {
            // Use the essentially uniformly distributed x coordinates to serve as a random function to resolve "diagonal edges" FLUTE outputs for certain terminals
	    //if (size_t(this->net_id) == 415){
	     //VTR_LOG(" FLUTE: (%d, %d) -> (%d, %d)\n", o_x, o_y, i_x, i_y);
	     fluteOutfile << " FLUTE: (" << o_x << ", " << o_y << ") -> (" << i_x << ", " << i_y << ")\n";
	    //}
            this->sb_edges.add_edge(o_x, o_y, i_x, i_y);
            //if (i_x % 2 == 0) {
            //    this->add_sb(i_x, o_y);
            //    this->sb_edges.add_edge(o_x, o_y, i_x, o_y);
            //    this->sb_edges.add_edge(i_x, o_y, i_x, i_y);
            //}
            //else {
            //    this->add_sb(o_x, i_y);
            //    this->sb_edges.add_edge(o_x, o_y, o_x, i_y);
            //    this->sb_edges.add_edge(o_x, i_y, i_x, i_y);
            //}
        } 
        else if (i_x != o_x || i_y != o_y) {
	    //if (size_t(this->net_id) == 415){
	     //VTR_LOG(" FLUTE: (%d, %d) -> (%d, %d)\n", o_x, o_y, i_x, i_y);
	     fluteOutfile << " FLUTE: (" << o_x << ", " << o_y << ") -> (" << i_x << ", " << i_y << ")\n";
	    //}
            this->sb_edges.add_edge(o_x, o_y, i_x, i_y);
        }
    }
    //int total_wirelength = Flute::flute_wirelength(flutetree);
    int total_wirelength = Flute::wirelength(flutetree);

    //VTR_LOG("Total wirelength: %d\n", total_wirelength);
    fluteOutfile << "Total wirelength: " << total_wirelength << "\n";

    //Flute::free_tree(flute1, flutetree);
    Flute::free_tree(flutetree);

    //Flute::free_tree(flutetree);
}

void Steiner::build_sb_rsmt_post_process_flute() {
    
    SteinerContext& steiner_ctx = g_vpr_ctx.mutable_steiner(); 
    int net = size_t(this->net_id);
    VTR_LOG("Net id: %d\n", net);

    // Check if net_id exists in net_edges
    // here the net_edges are either populated from an external file or is computed from fastroute
    auto it = steiner_ctx.net_edges.find(net);
    if (it == steiner_ctx.net_edges.end() || it->second.empty()) {
        VTR_LOG("  No FLUTE edges for net_id: %d\n", net);
        return;
    }

    const auto& edge_list = it->second; // vector<Edge>
    for (const auto& edge : edge_list) {
        int o_x = edge.a.x;
        int o_y = edge.a.y;
        int i_x = edge.b.x;
        int i_y = edge.b.y;

        // Add SBs at endpoints if missing
        if (!this->has_node_at_loc(o_x, o_y)) {
            this->add_sb(o_x, o_y);
        }
        if (!this->has_node_at_loc(i_x, i_y)) {
            this->add_sb(i_x, i_y);
        }

        // Check for diagonal (invalid in RSMT)
        if (i_x != o_x && i_y != o_y) {
            VTR_ASSERT_MSG(false, "Diagonal FLUTE edge encountered (shouldn't happen)");
            VTR_LOG("  Diagonal FLUTE edge: (%d, %d) -> (%d, %d)\n", o_x, o_y, i_x, i_y);
        } else if (i_x != o_x || i_y != o_y) {
            // Valid edge, add to SB graph
            this->sb_edges.add_edge(o_x, o_y, i_x, i_y);
        }
    }
}

void Steiner::prune_subtree_from(const std::string& node_id, std::set<std::string>& edges_to_remove) {
    if (this->sb_edges.edge_map.find(node_id) == this->sb_edges.edge_map.end()) return;

    for (auto& [child_id, is_valid] : this->sb_edges.edge_map[node_id]) {
        if (is_valid) {
            is_valid = false;
            this->sb_edges.edge_map[child_id][node_id] = false;

            edges_to_remove.insert(node_id + "-" + child_id);
            prune_subtree_from(child_id, edges_to_remove);
        }
    }
}
void Steiner::make_tree_graph_directed(std::string source_sb_id, std::set<std::string>& visited, std::set<std::string>& edges_to_remove, const std::string& parent_sb_id) {    

    visited.insert(source_sb_id);
    if (this->sb_edges.edge_map.find(source_sb_id) == this->sb_edges.edge_map.end()) {
        return; // No edges to process
    }
    for (auto& sink_sb : this->sb_edges.edge_map[source_sb_id]) {
        const std::string& sink_id = sink_sb.first;
        bool& is_valid = sink_sb.second;

        if (!is_valid) {
            edges_to_remove.insert(source_sb_id + "-" + sink_id);
            continue;
        }
	auto parse_coords = [](const std::string& sb_id) -> std::pair<int, int> {
	    auto underscore_pos = sb_id.find('_');
	    int x = std::stoi(sb_id.substr(0, underscore_pos));
	    int y = std::stoi(sb_id.substr(underscore_pos + 1));
	    return {x, y};
	};

	if (!parent_sb_id.empty()) {
            auto [ax, ay] = parse_coords(parent_sb_id);
            auto [bx, by] = parse_coords(source_sb_id);
            auto [cx, cy] = parse_coords(sink_id);

    	    // Detect reversal in direction — a U-turn
    	    bool is_horizontal_uturn = (ay == by && by == cy) &&
    	                               ((bx < ax && cx > bx) || (bx > ax && cx < bx));
    	    bool is_vertical_uturn = (ax == bx && bx == cx) &&
    	                             ((by < ay && cy > by) || (by > ay && cy < by));
            if (is_horizontal_uturn || is_vertical_uturn) {
                is_valid = false;
                this->sb_edges.edge_map[sink_id][source_sb_id] = false;
                edges_to_remove.insert(source_sb_id + "-" + sink_id);

                //prune_subtree_from(sink_id, edges_to_remove);
                continue;
            }
	}

        if (visited.find(sink_id) != visited.end()) {
            edges_to_remove.insert(source_sb_id + "-" + sink_id);
            is_valid = false;
            this->sb_edges.edge_map[sink_id][source_sb_id] = false;
            continue;
        }

        this->sb_edges.edge_map[sink_id][source_sb_id] = false;

        this->make_tree_graph_directed(sink_id, visited, edges_to_remove, source_sb_id);
    }

    /*for (auto& sink_sb : this->sb_edges.edge_map[source_sb_id]) {
        if (!sink_sb.second) {
            edges_to_remove.insert(source_sb_id + "-" + sink_sb.first);
            continue; // Skip false edges
        }
        
        this->sb_edges.edge_map[sink_sb.first][source_sb_id] = false;
        if (visited.find(sink_sb.first) != visited.end()) {
            edges_to_remove.insert(source_sb_id + "-" + sink_sb.first);
        }
        else {
            this->make_tree_graph_directed(sink_sb.first, visited, edges_to_remove);
        }
    } */
}

void Steiner::create_coarsened_regions(std::string source_sb_id) {
    for (auto& source_sb : this->sb_edges.edge_map) {
        std::vector<std::string> to_remove;
        for (auto& sink_sb : this->sb_edges.edge_map[source_sb.first]) {
            if (!sink_sb.second) {
                
                to_remove.push_back(sink_sb.first);
            }
        }
        for (const auto& sb_id : to_remove) {
            this->sb_edges.edge_map[source_sb.first].erase(sb_id);
        }
    }
    // Traverse the RSMT and build paths for connections
    std::vector<int> all_sinks = this->traverse_rsmt_build_paths(source_sb_id);
    int i_net_id = size_t(this->net_id);
    for (auto sink_id : all_sinks) {
        std::pair<int,int> loc_pin(this->sb_map[source_sb_id+"_0"].x+1, this->sb_map[source_sb_id+"_0"].y);
        std::pair<int,int> loc_term(this->sb_map[source_sb_id+"_0"].x+1, this->sb_map[source_sb_id+"_0"].y);
        
        SteinerContext& steiner_ctx = g_vpr_ctx.mutable_steiner();
        for (const auto gnode : steiner_ctx.loc_type_to_gnode[loc_term]) {
            if (gnode.first == "SINK")
                continue;
            for (auto dnode : steiner_ctx.gnode_to_dnodes[gnode.second]) {
                steiner_ctx.connections_per_dnode[dnode].push_back(std::make_pair(i_net_id, sink_id));
            }
        }
    }

    // Define the order in which the sinks will be routed (nearest neighbour policy)
    SteinerContext& steiner_ctx = g_vpr_ctx.mutable_steiner();
    std::unordered_map<int, int>& steiner_sink_order = steiner_ctx.steiner_sink_orders[size_t(this->net_id)];
    int counter = 1;
    for (int sink : all_sinks) {
        // The sink ID is the pin_id of the CLB
        steiner_sink_order[sink] = counter;
        counter++;
    }
}

std::vector<int> Steiner::traverse_rsmt_build_paths(std::string source_sb_id) {
    std::vector<int> all_sink_ids;
    int o_x, o_y, i_x, i_y;
    char sep;
    
    std::istringstream iss(source_sb_id);
    iss >> o_x >> sep >> o_y;
    for (auto& sink_sb : this->sb_edges.edge_map[source_sb_id]) {
        if (!sink_sb.second) {
            continue; // Skip false edges if some remained after making the tree directed
        }

        std::vector<int> sink_ids = this->traverse_rsmt_build_paths(sink_sb.first); // Recursive call to traverse further

        iss = std::istringstream(sink_sb.first);
        iss >> i_x >> sep >> i_y;
        
        this->map_net_sink_to_dnodes(o_x, o_y, i_x, i_y, sink_ids);

        for (auto sink_id : sink_ids) {
            all_sink_ids.push_back(sink_id); 
        }
    }

    // Add all net_pin_ids whoose connections end at this location's grid tile (this Steiner node)
    for (int sink_id : loc_to_sink_ids(o_x, o_y)) {
        all_sink_ids.insert(all_sink_ids.begin(), sink_id);
        std::pair<int,int> loc_pin(o_x+1, o_y);
        std::pair<int,int> loc_term(o_x+1, o_y);
        std::pair<int,int> loc(o_x+1, o_y);
        int i_net_id = size_t(this->net_id);
        SteinerContext& steiner_ctx = g_vpr_ctx.mutable_steiner();
        for (auto gnode : steiner_ctx.loc_type_to_gnode[loc_term]) {
            std::stringstream ss(gnode.first);
            std::string type;
            std::getline(ss, type, '_');
            if (type == "SOURCE") {
                continue;
            }
            for (auto dnode : steiner_ctx.gnode_to_dnodes[gnode.second]) {
                steiner_ctx.connections_per_dnode[dnode].push_back(std::make_pair(i_net_id, sink_id));
            }
        }
    }

    return all_sink_ids; // Return the list of sink IDs whoose connections branch off from this SB
}

void Steiner::map_net_sink_to_dnodes(int source_x, int source_y, int sink_x, int sink_y, std::vector<int>& sink_ids) {
    SteinerContext& steiner_ctx = g_vpr_ctx.mutable_steiner();
    int i_net_id = size_t(this->net_id);
    std::string dir;
    // Determine the dircetion of the wire
    if (source_x == sink_x) {
        if (source_y > sink_y) { // Southwards wires
            dir = "S";
        }
        else if (source_y < sink_y) { // Northwards wires
            dir = "N";
        }
    }
    else if (source_x > sink_x) // Westwards wires
        dir = "W";
    else // Eastwards wires
        dir = "E";

    int x_low = std::min(source_x, sink_x);
    int y_low = std::min(source_y, sink_y);

    // Determine edge length 
    int e_len = std::max(std::abs(source_y - sink_y), std::abs(source_x - sink_x));

    for (auto sink_id : sink_ids) {
        for (int i = 0; i < e_len; i++) {
            int w_len = e_len - i;
            std::pair<int, int> loc;
            if (dir == "W")
                loc = std::make_pair(x_low+i+1, y_low);
            else if (dir == "E")
                loc = std::make_pair(x_low+i+1, y_low);
            else if (dir == "S")
                loc = std::make_pair(x_low, y_low+i+1);
            else if (dir == "N")
                loc = std::make_pair(x_low, y_low+i+1);

            if (w_len >= 12) {
                int g_id_12, g_id_4, g_id_2, g_id_1;
                // Get gnode IDs for available wirelengths
                g_id_12 = steiner_ctx.loc_dir_len_to_gnode[loc][dir+std::to_string(12)];
                g_id_4 = steiner_ctx.loc_dir_len_to_gnode[loc][dir+std::to_string(4)];
                g_id_2 = steiner_ctx.loc_dir_len_to_gnode[loc][dir+std::to_string(2)];
                g_id_1 = steiner_ctx.loc_dir_len_to_gnode[loc][dir+std::to_string(1)];
                // Iterate over dnodes for each gnode ID, and map them to net+sink pairs
                for (auto d_id_12 : steiner_ctx.gnode_to_dnodes[g_id_12]) {
                    steiner_ctx.connections_per_dnode[d_id_12].push_back(std::make_pair(i_net_id, sink_id));
                }
                for (auto d_id_4 : steiner_ctx.gnode_to_dnodes[g_id_4]) {
                    steiner_ctx.connections_per_dnode[d_id_4].push_back(std::make_pair(i_net_id, sink_id));
                }
                for (auto d_id_2 : steiner_ctx.gnode_to_dnodes[g_id_2]) {
                    steiner_ctx.connections_per_dnode[d_id_2].push_back(std::make_pair(i_net_id, sink_id));
                }
                for (auto d_id_1 : steiner_ctx.gnode_to_dnodes[g_id_1]) {
                    steiner_ctx.connections_per_dnode[d_id_1].push_back(std::make_pair(i_net_id, sink_id));
                }
            }
            else if (w_len >= 4) {
                int g_id_4, g_id_2, g_id_1;
                // Get gnode IDs for available wirelengths
                g_id_4 = steiner_ctx.loc_dir_len_to_gnode[loc][dir+std::to_string(4)];
                g_id_2 = steiner_ctx.loc_dir_len_to_gnode[loc][dir+std::to_string(2)];
                g_id_1 = steiner_ctx.loc_dir_len_to_gnode[loc][dir+std::to_string(1)];
                // Iterate over dnodes for each gnode ID, and map them to net+sink pairs
                for (auto d_id_4 : steiner_ctx.gnode_to_dnodes[g_id_4]) {
                    steiner_ctx.connections_per_dnode[d_id_4].push_back(std::make_pair(i_net_id, sink_id));
                }
                for (auto d_id_2 : steiner_ctx.gnode_to_dnodes[g_id_2]) {
                    steiner_ctx.connections_per_dnode[d_id_2].push_back(std::make_pair(i_net_id, sink_id));
                }
                for (auto d_id_1 : steiner_ctx.gnode_to_dnodes[g_id_1]) {
                    steiner_ctx.connections_per_dnode[d_id_1].push_back(std::make_pair(i_net_id, sink_id));
                }
            }
            else if (w_len >= 2) {
                int g_id_2, g_id_1;
                // Get gnode IDs for available wirelengths
                g_id_2 = steiner_ctx.loc_dir_len_to_gnode[loc][dir+std::to_string(2)];
                g_id_1 = steiner_ctx.loc_dir_len_to_gnode[loc][dir+std::to_string(1)];
                // Iterate over dnodes for each gnode ID, and map them to net+sink pairs
                for (auto d_id_2 : steiner_ctx.gnode_to_dnodes[g_id_2]) {
                    steiner_ctx.connections_per_dnode[d_id_2].push_back(std::make_pair(i_net_id, sink_id));
                }
                for (auto d_id_1 : steiner_ctx.gnode_to_dnodes[g_id_1]) {
                    steiner_ctx.connections_per_dnode[d_id_1].push_back(std::make_pair(i_net_id, sink_id));
                }
            }
            else {
                int g_id_1;
                // Get gnode ID for available wirelength
                g_id_1 = steiner_ctx.loc_dir_len_to_gnode[loc][dir+std::to_string(1)];
                // Iterate over dnodes for the gnode ID, and map them to net+sink pairs
                for (auto d_id_1 : steiner_ctx.gnode_to_dnodes[g_id_1]) {
                    steiner_ctx.connections_per_dnode[d_id_1].push_back(std::make_pair(i_net_id, sink_id));
                }
            }   
        }
    }
}

void Steiner::compute_dependency_graph_sink_order(std::string source_sb_id, std::unordered_map<size_t, std::unordered_map<int, int>>& sink_order) {
    // Create dependency graph
    std::unordered_map<std::string, std::vector<std::string>> dependency_graph;
    this->create_dependency_graph(source_sb_id, dependency_graph);
    VTR_LOG("Dependency graph created for source SB: %s\n", source_sb_id.c_str());
    // Compute sink order
    this->compute_sink_order(source_sb_id, dependency_graph, sink_order);
}

void Steiner::create_dependency_graph(std::string source_sb_id, std::unordered_map<std::string, std::vector<std::string>>& dependency_graph) {
    std::vector<std::pair<int, std::string>> results;
    for (const auto& sink_sb : this->sb_edges.edge_map[source_sb_id]) {
        if (!sink_sb.second) {
            continue; // Skip false edges
        }
        std::pair<int, std::string> dist_sink = this->create_dependency_subgraph(sink_sb.first, dependency_graph);
        dependency_graph[source_sb_id].push_back(dist_sink.second);
    }
}

std::pair<int, std::string> Steiner::create_dependency_subgraph(std::string source_sb_id, std::unordered_map<std::string, std::vector<std::string>>& dependency_graph) {
    std::vector<std::pair<int, std::string>> results;
    std::pair<int, std::string> maxdist_sink(0, source_sb_id); // If the current SB is a terminal, it is the furthest SB reachable from itself
    VTR_LOG("Creating dependency subgraph for source SB: %s\n", source_sb_id.c_str());
    for (const auto& sink_sb : this->sb_edges.edge_map[source_sb_id]) {
        if (!sink_sb.second) {
            continue; // Skip false edges
        }
        std::pair<int, std::string> dist_sink = this->create_dependency_subgraph(sink_sb.first, dependency_graph);
        dist_sink.first += calculate_manhattan_distance(source_sb_id, sink_sb.first);
        results.push_back(dist_sink);
    }

    if (!results.empty()) {
        // Sort results in descending order by distance
        std::stable_sort(results.begin(), results.end(), [](const std::pair<int, std::string>& a, const std::pair<int, std::string>& b) {
            return a.first > b.first; 
        });
        // The subcluster's source SB is dependant on the furthest sink SB within it
        maxdist_sink = results[0];
        dependency_graph[results[0].second].push_back(source_sb_id);
        // Other sinks are dependant on the source SB
        for (int i = 1; i < results.size(); i++) {
            dependency_graph[source_sb_id].push_back(results[i].second);
        }
    }
    return maxdist_sink;
}

void Steiner::compute_sink_order(std::string source_sb_id, std::unordered_map<std::string, std::vector<std::string>>& dependency_graph, std::unordered_map<size_t, std::unordered_map<int, int>>& sink_order) {
    int i_net_id = size_t(this->net_id), num_sinks = 0;
    std::pair<int,int> source_sb_loc = make_pair_from_str_id(source_sb_id);
    for (auto sink_id : loc_to_sink_ids(source_sb_loc.first, source_sb_loc.second)) {
        sink_order[i_net_id][sink_id] = num_sinks+1;
        num_sinks++;
    }
    
    num_sinks = this->compute_sink_order(source_sb_id, dependency_graph, sink_order, num_sinks);
    VTR_LOG("Number of sinks: %d\nNumber of orders %d\n", g_vpr_ctx.clustering().clb_nlist.net_sinks(this->net_id).size(), num_sinks);
    VTR_ASSERT(num_sinks == g_vpr_ctx.clustering().clb_nlist.net_sinks(this->net_id).size());
}

int Steiner::compute_sink_order(std::string source_sb_id, std::unordered_map<std::string, std::vector<std::string>>& dependency_graph, std::unordered_map<size_t, std::unordered_map<int, int>>& sink_order, int i) {
    size_t i_net_id = size_t(this->net_id);
    VTR_LOG("Computing sink order for source SB: %s\n", source_sb_id.c_str());
    for (auto sink_sb_id : dependency_graph[source_sb_id]) {
        std::pair<int, int> sink_sb_loc = make_pair_from_str_id(sink_sb_id);
        std::vector<int> sink_ids = this->loc_to_sink_ids(sink_sb_loc.first, sink_sb_loc.second);
        for (auto sink_id : sink_ids) {
            sink_order[i_net_id][sink_id] = i+1;
            i++;
        }
    }
    for (auto sink_sb_id : dependency_graph[source_sb_id]) {
        i = this->compute_sink_order(sink_sb_id, dependency_graph, sink_order, i);
    }
    return i;
}

// This was used to avoid preprocessing of "global connect?" nets, as some of them were too large (and none of them are routed by the detailed router)
std::set<size_t> load_nets_to_skip() {
    std::set<size_t> nets_to_skip;
    // Load nets to skip from a file
    const std::string filename = "nets_to_skip_st.txt";
    std::ifstream infile(filename);
    // Each line contains a single NetId to skip
    std::string line;
    while (std::getline(infile, line)) {
        ClusterNetId net_id = ClusterNetId(std::stoul(line));
        nets_to_skip.insert(size_t(net_id));
    }
    infile.close();
    return nets_to_skip;
}

void load_gnode_maps(SteinerContext& steiner_ctx) {
    //const std::string to_gnode_filename = "../../../scripts/node_dict_gr_1L_reg_v4.map";
    //const std::string to_dnode_filename = "../../../scripts/gr_dr_map_reg_v4_full_iib.map";
    const std::string post_processed_flute_filename = "./post_processed_flute_trees.txt";

    VTR_LOG("Loading post porcessed flute trees");
    std::ifstream post_flute_file(post_processed_flute_filename);
    if (!post_flute_file.is_open()) {
        std::cerr << "Error: Could not open file " << post_processed_flute_filename << std::endl;
        return;
    }


    std::string line;
    int current_net = -1;
    std::regex net_re(R"(^\s*Net id:\s*(\d+))");
    std::regex edge_re(R"(^\s*FLUTE:\s*\((\-?\d+),\s*(\-?\d+)\)\s*->\s*\((\-?\d+),\s*(\-?\d+)\))");

    while (std::getline(post_flute_file, line)) {
        std::smatch match;

        // Check for Net id
        if (std::regex_search(line, match, net_re)) {
            current_net = std::stoi(match[1].str());
            continue;
        }

        // Check for FLUTE edge
	// Note: for nets with no edges, nets which are intra-cluster, will not have an entry in net_edges
	// in route_timing.cpp ensure not to query for such a net
        if (std::regex_search(line, match, edge_re)) {
            int x1 = std::stoi(match[1].str());
            int y1 = std::stoi(match[2].str());
            int x2 = std::stoi(match[3].str());
            int y2 = std::stoi(match[4].str());

            steiner_ctx.net_edges[current_net].emplace_back(Point{x1, y1}, Point{x2, y2});
        }

    }

    post_flute_file.close();


    /*VTR_LOG("Loading loc_to_gnode map:\n");
    std::ifstream toGfile(to_gnode_filename);
    if (!toGfile.is_open()) {
        std::cerr << "Error: Could not open file " << to_gnode_filename << std::endl;
        return;
    }

    std::unordered_set<std::string> wires = {
        "N1", "N2", "N4", "N12",
        "E1", "E2", "E4", "E12",
        "S1", "S2", "S4", "S12",
        "W1", "W2", "W4", "W12"
    };
    // std::unordered_set<std::string> tile_nodes = {
    //     "SOURCE", "SINK", "OPIN", "IPIN"
    // };

    while (std::getline(toGfile, line)) {
        std::istringstream iss(line);
        int x, y, gnode_id, throwaway;
        std::string dir_len;

        if (!(iss >> x >> y >> dir_len >> gnode_id >> throwaway)) {
            std::cerr << "Warning: Malformed line skipped: " << line << std::endl;
            continue;
        }

        std::pair<int, int> loc(x, y);
        if (wires.find(dir_len) == wires.end()) {
            // Other RRNodes
            // e.g. (147.46), A_o5_out
            steiner_ctx.loc_type_to_gnode[loc][dir_len] = gnode_id; // Store the gnode ID, dir_len contains the non-wire type of the RRNode here
        }
        else {
            // Wire RRNodes
            steiner_ctx.loc_dir_len_to_gnode[loc][dir_len] = gnode_id; // Store the gnode ID
        }
    }

    toGfile.close();*/

    /*VTR_LOG("Loading gnode_to_dnode map:\n");
    std::ifstream toDfile(to_dnode_filename);
    if (!toDfile.is_open()) {
        std::cerr << "Error: Could not open file " << to_dnode_filename << std::endl;
        return;
    }

    while (std::getline(toDfile, line)) {
        std::istringstream iss(line);
        int gnode_id;
        if (!(iss >> gnode_id)) {
            std::cerr << "Warning: Malformed line skipped: " << line << std::endl;
            continue;
        }

        int dnode_id;
        std::vector<int> dnode_ids;

        while (iss >> dnode_id) {
            dnode_ids.push_back(dnode_id);
        }

        steiner_ctx.gnode_to_dnodes[gnode_id] = std::move(dnode_ids);
    }

    toDfile.close();*/
}

void create_connections_per_dnode_file(SteinerContext& steiner_ctx) {
    const std::string filename = "connections_per_dnode.txt";
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "Error: Cannot open " << filename << " for writing.\n";
        return;
    }

    for (const auto& [dnode_id, net_sink_list] : steiner_ctx.connections_per_dnode) {
        outfile << dnode_id;
        for (const auto& [net_id, sink_id] : net_sink_list) {
            int hop = 0;  // fixed hop
            outfile << " " << net_id << "_" << sink_id << "_" << hop;
        }
        outfile << "\n";
    }

    outfile.close();
}

void create_sink_order_file(const ClusteringContext& cluster_ctx, const vtr::vector_map<ClusterBlockId, t_block_loc>& blocks_locs) {
    const std::string filename = "sink_order.txt";
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "Error: Cannot open " << filename << " for writing.\n";
        return;
    }

    for (ClusterNetId net_id : cluster_ctx.clb_nlist.nets()) {
        // Enter sink order for [net_id]
        outfile << size_t(net_id);

        int num_sinks = cluster_ctx.clb_nlist.net_sinks(net_id).size();

        for (int i = 1; i <= num_sinks; i++) {
            outfile << " " << i;
        }

        outfile << "\n";
    }

    outfile.close();
}

void create_branch_node_map_file(SteinerContext& steiner_ctx) {
    std::unordered_map<std::pair<int, int>, std::vector<int>, vtr::hash_pair> branch_node_map;
    for (auto dnode_claimers : steiner_ctx.connections_per_dnode) {
        for (auto pair : dnode_claimers.second) {
            branch_node_map[pair].push_back(dnode_claimers.first);
        }
    }

    const std::string filename = "branch_node_map.txt";
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Cannot open " << filename << " for writing.\n";
        return;
    }
    for (auto& pair : branch_node_map) {
        outfile << pair.first.first << "_" << pair.first.second;
        for (auto dnode : pair.second) {
            outfile << " " << dnode;
        }
        outfile << "\n";
    }
    outfile.close();
    branch_node_map.clear();
}

/*
    (1) Load the mapping from gnodes to dnodes, as well as the 
    mapping from location, [direction_length]/[type] to gnodes.
    (2) For every net in the clustered netlist, build an RSMT 
    with net's pins as terminals.
        (3) For every net in the clustered netlist, create RSMT
        (3.s) For every edge in the net's RSMT, find the corresponding 
        routing resource nodes in the RR graph and map it to all the 
        nets and sinks that use it. Assign the net a sink order that
        routes the sinks in the order they appear along the RSMT.
        (3.o) Use the net's RSMT to compute the sink order using Francois' 
        dependency graph algorithm
    (4) Write the mapping to the "connections_per_dnode.txt" 
    file.
    (5) Remove the gnode maps from the global Steiner context to
    free up memory.
*/

std::tuple<std::unordered_map<int, std::vector<Corridor>>, std::unordered_map<int, std::vector<unsigned short>>, std::unordered_map<int, bool>> Steiner::build_corridor_list_per_connection_from_flute(std::string source_sb_id) const {
   
    std::unordered_map<int, std::vector<Corridor>> corridor_list_per_connection;
    std::unordered_map<int, std::vector<unsigned short>> corridors_lookahead;
    std::unordered_map<int, bool> connection_intra_tile;

    std::unordered_map<std::string, std::string> child_to_parent;
    // 1. Build reverse edge map: child → parent

    if (source_sb_id.empty()) {
        // Handle case with no source pin
	VTR_LOG(" [Building corridors: No source for this net\n");
        return {corridor_list_per_connection, corridors_lookahead, connection_intra_tile};
    }

    auto parse_coords = [](const std::string& sb_id) -> std::pair<int, int> {
            auto underscore_pos = sb_id.find('_');
            int x = std::stoi(sb_id.substr(0, underscore_pos));
            int y = std::stoi(sb_id.substr(underscore_pos + 1));
            return {x, y};
    };


    auto direction_from_to = [&parse_coords](const std::string& from, const std::string& to) -> std::string {
        auto [x1, y1] = parse_coords(from);
        auto [x2, y2] = parse_coords(to);
        if (x1 == x2) return (y2 > y1) ? "N" : "S";
        if (y1 == y2) return (x2 > x1) ? "E" : "W";
        return ""; // invalid
    };


    auto is_opposite = [](const std::string& dir1, const std::string& dir2) -> bool {
        return (dir1 == "N" && dir2 == "S") || (dir1 == "S" && dir2 == "N") ||
               (dir1 == "E" && dir2 == "W") || (dir1 == "W" && dir2 == "E");
    };


    std::set<std::string> visited;
    std::queue<std::pair<std::string, std::string>> q;
    q.push({source_sb_id, ""}); // Source has no parent
    visited.insert(source_sb_id); 

    while (!q.empty()) {
        auto [current_node, parent_node] = q.front();
        q.pop();

        // Check for edges from the current node
        if (this->sb_edges.edge_map.count(current_node)) {
            for (const auto& [neighbor_node, is_valid] : this->sb_edges.edge_map.at(current_node)) {

		// all edges are marked true in both directions
                if (visited.count(neighbor_node) || !is_valid) {
                    continue;
                }

                // If it passes all checks, add to tree
                visited.insert(neighbor_node);
                child_to_parent[neighbor_node] = current_node;
                q.push({neighbor_node, current_node});
            }
        }
    }
    for (const auto& [str_id, sb] : this->sb_map) {
        if (sb.pin_id <= 0) {
            continue; // Skip non-sinks and the source
        }

        int sink_id = sb.pin_id;
        std::string current = std::to_string(sb.x) + "_" + std::to_string(sb.y);;
	// if sink and source in the same tile, set to true
	connection_intra_tile[sink_id] = (source_sb_id == current) ? true : false;
        std::vector<Corridor> corridor_path;
        std::vector<unsigned short> corridor_lookahead;

	bool path_found = false;//(child_to_parent.count(current) || current == source_sb_id);
	unsigned short previous_distance = 0;
        while (child_to_parent.count(current)) {
            const std::string& parent = child_to_parent.at(current);
            auto [from_x, from_y] = parse_coords(parent);
            auto [to_x, to_y] = parse_coords(current);
	    // assumes that corridor are not horizontal, they are either horizontal or vertical
	    unsigned short corridor_len = std::abs(from_x - to_x) + std::abs(from_y - to_y);

            corridor_path.emplace_back(from_x, from_y, to_x, to_y);

	    corridor_lookahead.emplace_back(previous_distance);

 	    previous_distance += corridor_len;

            if (parent == source_sb_id) {
                path_found = true;
                break;
            }
            current = parent;
        }

        if (path_found) {
            std::reverse(corridor_path.begin(), corridor_path.end());
            std::reverse(corridor_lookahead.begin(), corridor_lookahead.end());
            corridor_list_per_connection[sink_id] = std::move(corridor_path);
	    corridors_lookahead[sink_id] = std::move(corridor_lookahead);
        } else if (current != source_sb_id) {
            VTR_LOG("Sink pin_id %d at %d is disconnected from the source. No path found.\n", sink_id, this->net_id);
        }
	else if (current == source_sb_id){ // intra-CLB connections having no edges, but we need to initialize the connection
	    corridor_list_per_connection[sink_id] = std::move(corridor_path);
	    corridors_lookahead[sink_id] = std::move(corridor_lookahead);
	}
    }

    return {corridor_list_per_connection, corridors_lookahead, connection_intra_tile};
}

std::tuple<std::unordered_map<int, std::vector<Corridor>>, std::unordered_map<int, std::vector<unsigned short>>, std::unordered_map<int, bool>> Steiner::build_corridor_list_per_connection(std::string source_sb_id) const {

    std::unordered_map<int, std::vector<Corridor>> corridor_list_per_connection;
    std::unordered_map<int, std::vector<unsigned short>> corridors_lookahead;
    std::unordered_map<int, bool> connection_intra_tile;

    std::unordered_map<std::string, std::pair<std::string, bool>> child_to_parent;
    // 1. Build reverse edge map: child → parent

    if (source_sb_id.empty()) {
        // Handle case with no source pin
	VTR_LOG(" [Building corridors: No source for this net\n");
        return {corridor_list_per_connection, corridors_lookahead, connection_intra_tile};
    }
   
    auto parse_coords = [](const std::string& sb_id) -> std::pair<int, int> {
            auto underscore_pos = sb_id.find('_');
            int x = std::stoi(sb_id.substr(0, underscore_pos));
            int y = std::stoi(sb_id.substr(underscore_pos + 1));
            return {x, y};
    };
    

    auto direction_from_to = [&parse_coords](const std::string& from, const std::string& to) -> std::string {
        auto [x1, y1] = parse_coords(from);
        auto [x2, y2] = parse_coords(to);
        if (x1 == x2) return (y2 > y1) ? "N" : "S";
        if (y1 == y2) return (x2 > x1) ? "E" : "W";
        return ""; // invalid
    };
    
    auto orientation_from_to = [&parse_coords](const std::string& from, const std::string& to) -> std::string {
        auto [x1, y1] = parse_coords(from);
        auto [x2, y2] = parse_coords(to);
        if (x1 == x2) return "V";
        if (y1 == y2) return "H";
        return ""; // invalid
    };

    
    auto is_opposite = [](const std::string& dir1, const std::string& dir2) -> bool {
        return (dir1 == "N" && dir2 == "S") || (dir1 == "S" && dir2 == "N") ||
               (dir1 == "E" && dir2 == "W") || (dir1 == "W" && dir2 == "E");
    };
    
    std::set<std::string> visited;
    std::queue<std::pair<std::string, std::string>> q;
    q.push({source_sb_id, ""}); // Source has no parent
    visited.insert(source_sb_id); 

    bool corridor_end_point = false;
    std::string edge_orientation;
    while (!q.empty()) {
        auto [current_node, parent_node] = q.front();
        q.pop();
	corridor_end_point = false;
	if (parent_node != "") {
	    edge_orientation = orientation_from_to(parent_node, current_node);
	}
        // Check for edges from the current node
        if (this->sb_edges.edge_map.count(current_node)) {
	    // 1. first check if current node is a terminal, if yes, mark as end point
	    // this does a for loop over all sbs // can be optimized for runtime later
	    auto [x, y] = parse_coords(current_node);
	    int pin_id = this->loc_to_pin_id(x, y);
	    if (pin_id > 0) {
	        corridor_end_point = true;
	    }
	    // 2. if total neighbours more than two, mark as corridor end point
	    int total_neighbours = this->sb_edges.edge_map.at(current_node).size();
 	    if (total_neighbours > 2 || (current_node == source_sb_id)) {
		// always mark source as SP
	    	corridor_end_point = true;
	    }

	    // this code for now marks source as a SP if it has more than two neighbours
            for (const auto& [neighbor_node, is_valid] : this->sb_edges.edge_map.at(current_node)) {
                
		// all edges are marked true in both directions
                if (visited.count(neighbor_node) || !is_valid) {
                    continue;
                }
		// 3. check if the direction changed, if yes, mark as end of corridor
		// only check change of direction for neighbour 2 because if more neighbours then 
		// it is already marked as steiner point
		if (total_neighbours == 2 && parent_node != ""){
		    std::string next_edge_orientation = orientation_from_to(current_node, neighbor_node);
		    if (next_edge_orientation != edge_orientation) {
		        corridor_end_point = true;
		    }
		}
                
                // If it passes all checks, add to tree
                visited.insert(neighbor_node);
                child_to_parent[neighbor_node] = {current_node, corridor_end_point};
                q.push({neighbor_node, current_node});
            }
        }
    }

    for (const auto& [str_id, sb] : this->sb_map) {
        if (sb.pin_id <= 0) {
            continue; // Skip non-sinks and the source
        }
	
        int sink_id = sb.pin_id;
        std::string current = std::to_string(sb.x) + "_" + std::to_string(sb.y);;
	// if sink and source in the same tile, set to true
	connection_intra_tile[sink_id] = (source_sb_id == current) ? true : false;
        std::vector<Corridor> corridor_path;
        std::vector<unsigned short> corridor_lookahead;

	bool path_found = false;
	unsigned short previous_distance = 0;

	std::string stored_to = current;  //ipin
        while (child_to_parent.count(current)) {
	    const auto& [parent, parent_is_steiner_point] = child_to_parent.at(current);
            // const std::string& parent = child_to_parent.at(current);
	    // now, if parent is SP
	    // we have hit end of corridor
	    // record it in corridors like the current state of code
	    // if parent is not SP
	    // record only the "to" coordinate for the corridor
	    // keep iterating until you hit a SP and for that record "from" coordinate
	    // i should always mark source as SP above 
            
	    if (parent_is_steiner_point == true ) {
                auto [to_x, to_y] = parse_coords(stored_to);
	        auto [from_x, from_y] = parse_coords(parent);

	    	// assumes that corridor are not diagonal, they are either horizontal or vertical
		VTR_ASSERT(from_x == to_x || from_y == to_y);
	    	unsigned short corridor_len = std::abs(from_x - to_x) + std::abs(from_y - to_y);

            	corridor_path.emplace_back(from_x, from_y, to_x, to_y);

	    	corridor_lookahead.emplace_back(previous_distance);

 	    	previous_distance += corridor_len;
            	if (parent == source_sb_id) {
            	    path_found = true;
            	    break;
            	}
		stored_to = parent; 
	    }

            current = parent;
        }

        if (path_found) {
            std::reverse(corridor_path.begin(), corridor_path.end());
            std::reverse(corridor_lookahead.begin(), corridor_lookahead.end());
            corridor_list_per_connection[sink_id] = std::move(corridor_path);
	    corridors_lookahead[sink_id] = std::move(corridor_lookahead);
        } else if (current != source_sb_id) {
            VTR_LOG("Sink pin_id %d at %d is disconnected from the source. No path found.\n", sink_id, this->net_id);
        }
	else if (current == source_sb_id){ // intra-CLB connections having no edges, but we need to initialize the connection
	    corridor_list_per_connection[sink_id] = std::move(corridor_path);
	    corridors_lookahead[sink_id] = std::move(corridor_lookahead);
	}
    }
    
    return {corridor_list_per_connection, corridors_lookahead, connection_intra_tile};
}


void steiner_pre_processing(bool create_steiner_constraints, bool compute_dependency_graph_sink_orders, bool dump_raw_flute_trees, std::string global_router_algorithm, int global_channel_capacity) {
    // Initialize context refereces
    const ClusteringContext& cluster_ctx = g_vpr_ctx.clustering();
    const PlacementContext& placement_ctx = g_vpr_ctx.placement();
    auto& device_ctx  = g_vpr_ctx.device();

    SteinerContext& steiner_ctx = g_vpr_ctx.mutable_steiner();

    // Nets belonging to the "global connecting nets" are not routed by the global router
    // net is global should replace this
    //std::set<size_t> nets_to_skip = load_nets_to_skip();

    // Create RSMT for each net in the clustered netlist (except for
    // "global connecting nets") and update the global Steiner context 
    // with information about the corresponding constrained regions
    const vtr::vector_map<ClusterBlockId, t_block_loc>& blocks_locs = g_vpr_ctx.placement().block_locs;
    
    const std::string fluteOutfilename = "nets_with_and_without_diagonals_and_uturns.txt";
    std::ofstream fluteOutfile(fluteOutfilename);

    if (!fluteOutfile.is_open()) {
        std::cerr << "Error: Cannot open " << fluteOutfilename << " for writing.\n";
        return;
    }


    if (global_router_algorithm == "FLUTE") {
	// go in the if only to load the post processed flute trees
    	if (create_steiner_constraints && dump_raw_flute_trees == false)
       	    load_gnode_maps(steiner_ctx);
    	// Initialize FLUTE
    	//Flute::FluteState *flute1 = Flute::flute_init(FLUTE_POWVFILE, FLUTE_PORTFILE);
	//flute initialization for the version FastRoute uses
	static bool flute_initialized = false;
	if (!flute_initialized) {
		Flute::readLUT();   // or readLUT(FLUTE_POWVFILE, FLUTE_PORTFILE)
    		flute_initialized = true;
	}
    	VTR_LOG("FLUTE initialized.\n");

    	for (ClusterNetId net_id : cluster_ctx.clb_nlist.nets()) {
    	    // Skip "global connecting nets"; this might also be done by checking clb_nlist.net_is_global(net_id), but I haven't been able to get that to work
    	    //unsigned int num_sinks = cluster_ctx.clb_nlist.net_sinks(net_id).size()
    	    //VTR_LOG("Net id: %d (num_sinks: %d)\n", size_t(net_id), num_sinks);
    	    if (cluster_ctx.clb_nlist.net_is_ignored(net_id)) {
    	        printf("Net ID: %d\n", size_t(net_id));
    	        continue;
    	    }

    	    // Construct the RSMT for the net
    	    //Steiner steiner(cluster_ctx.clb_nlist, net_id, blocks_locs, flute1, fluteOutfile, dump_raw_flute_trees);
    	    Steiner steiner(cluster_ctx.clb_nlist, net_id, blocks_locs, fluteOutfile, dump_raw_flute_trees, global_router_algorithm);

    	    if (create_steiner_constraints) {
    	        // LUKA's implementation that fits in the FCCM implementation; but this implementation leads to high runtime due to repeated memory accesses
    	        // Create RSMT constrained regions for the net and
    	        // update the global Steiner context with the information 
    	        //steiner.create_coarsened_regions(make_str_id(steiner.source_x, steiner.source_y));

    	        // The RSMT is used to compute the nearest neighbour sink order;
    	        // by doing this, we reduce heap pops and redundant wires
    	        //
    	        // Below is the new implemenation driven to optimize runtime by reducing memory accesses
    	        auto [net_corridors, net_corridors_lookahead, connection_intra_tile] = steiner.build_corridor_list_per_connection_from_flute(make_str_id(steiner.source_x, steiner.source_y));
    	        steiner_ctx.all_corridors[net_id] = std::move(net_corridors);
    	        steiner_ctx.all_corridors_lookahead[net_id] = std::move(net_corridors_lookahead);
    	        steiner_ctx.net_connection_intra_tile[net_id] = std::move(connection_intra_tile);
    	        // Debug print: Dump corridors for this net
    	        /*if (size_t(net_id) == 415 || size_t(net_id) == 0 || size_t(net_id) == 1) {
    	        VTR_LOG("Corridors for Net %zu:\n", size_t(net_id));
    	        for (const auto& sink_entry : steiner_ctx.all_corridors[net_id]) {
    	            int sink_id = sink_entry.first;
    	            const std::vector<Corridor>& corridor_list = sink_entry.second;
    	        
    	            VTR_LOG("  Sink %d has %zu corridors:\n", sink_id, corridor_list.size());
    	            for (size_t i = 0; i < corridor_list.size(); ++i) {
    	                const Corridor& c = corridor_list[i];
    	                VTR_LOG("    Corridor %zu: (%d, %d) → (%d, %d)\n", i, c.from_x, c.from_y, c.to_x, c.to_y);
    	            }
    	        }
    	        VTR_LOG("\n");
    	        }*/
    	    }
    	    // Calculate Francois' dependency graph sink order
    	    if (compute_dependency_graph_sink_orders) {
    	        steiner.compute_dependency_graph_sink_order(make_str_id(steiner.source_x, steiner.source_y), steiner_ctx.steiner_sink_orders);
    	        VTR_LOG("Sink order for net ID: %zu\n", size_t(net_id));
    	        for (const auto& [sink_id, order] : steiner_ctx.steiner_sink_orders[size_t(net_id)]) {
    	            VTR_LOG("Sink ID: %d, Order: %d\n", sink_id, order);
    	        }       
    	    }
    	}
	if (dump_raw_flute_trees == true){
	    VTR_LOG("Finished running flute\n");
	    VTR_LOG("Finished running flute\n");
    	    fluteOutfile.close();
	    VTR_ASSERT(false);
	}
    }
    else if (global_router_algorithm == "FastRoute4.0") {
	struct PinLess {
   	    bool operator()(const PIN& a, const PIN& b) const {
                if (a.x != b.x) return a.x < b.x;
                if (a.y != b.y) return a.y < b.y;
                return a.layer < b.layer;
    	    }
	};
	auto start_time_fr = std::chrono::high_resolution_clock::now();
    	// 1. Define FastRoute grid
    	int Gx = device_ctx.grid.width() - 2; // 168
    	int Gy = device_ctx.grid.height() - 2; // 480
	VTR_LOG("Width = %d height = %d\n", Gx, Gy);
    	int nLayers = 1; // FPGA-style: single layer

	FastRoute::FT fr_router;

    	fr_router.setGridsAndLayers(Gx, Gy, nLayers);

    	int layer = 1;
    	int capacity = global_channel_capacity; // choose something; or derive from channel width
	VTR_LOG("Global channel Capacity = %d \n", capacity);
    	fr_router.addVCapacity(capacity, layer);
    	fr_router.addHCapacity(capacity, layer); 

	fr_router.addMinWidth(1, layer);
    	fr_router.addMinSpacing(1, layer);
    	fr_router.addViaSpacing(1, layer);

    	fr_router.setLowerLeft(1, 1);
    	fr_router.setTileSize(1, 1);
   
	// 2. Count nets and set number
    	int num_nets = 0;
	int valid_global_nets = 0;
	int max_net_degree = 0;
    	for (ClusterNetId net_id : cluster_ctx.clb_nlist.nets()) {
    	    if (cluster_ctx.clb_nlist.net_is_ignored(net_id)) {
		continue;
    	    }
            bool is_global = false;
            int first_x = -1, first_y = -1;
            int pin_count = 0;

            // Quick check of pins to determine global status
            std::set<PIN, PinLess> unique_pins;	
            for (ClusterPinId pin_id : cluster_ctx.clb_nlist.net_pins(net_id)) {
                auto block_id = cluster_ctx.clb_nlist.pin_block(pin_id);
                int x = blocks_locs[block_id].loc.x - 1;
                int y = blocks_locs[block_id].loc.y;

                if (pin_count == 0) {
                    first_x = x; first_y = y;
                } 
		else {
                    if (x != first_x || y != first_y) is_global = true;
                }
                pin_count++;
		
		PIN p;
                p.x = x;
                p.y = y;
                p.layer = 1;
		unique_pins.insert(PIN{(long)x, (long)y, 1});
            }
    	    num_nets++;
	    int unique_pins_fanout = static_cast<int>(unique_pins.size());
	    max_net_degree = std::max(unique_pins_fanout, max_net_degree);

            // Only count if it's a valid global net
            if (is_global) {
                valid_global_nets++;
            }
    	}

    	fr_router.setNumberNets(valid_global_nets);
	fr_router.setMaxNetDegree(max_net_degree);
	VTR_LOG("valid net = %d max_net_deg = %d\n",valid_global_nets, max_net_degree);
   
	std::vector<ClusterNetId> fr_to_vpr_net;  
	fr_to_vpr_net.reserve(valid_global_nets);
	// For FastRoute: Keep pin storage alive for all nets and also the net names
	std::vector<std::vector<PIN>> fr_pins(valid_global_nets);
	std::vector<std::string> stored_net_names;
	stored_net_names.reserve(valid_global_nets);
	
	std::vector<ClusterNetId> intra_tile_nets;  
	int intra_tile_nets_count = num_nets - valid_global_nets;
	intra_tile_nets.reserve(intra_tile_nets_count);
      
       	// 3. Add nets to FastRoute
    	int netIdx = 0;
    	for (ClusterNetId net_id : cluster_ctx.clb_nlist.nets()) {
    	    if (cluster_ctx.clb_nlist.net_is_ignored(net_id)) continue;

    	    std::string net_name = cluster_ctx.clb_nlist.net_name(net_id);

    	    int nPins = cluster_ctx.clb_nlist.net_pins(net_id).size();

            std::set<PIN, PinLess> unique_pins;	
    	    int pin_i = 0;
	    // C. Extract Pins & Detect Global vs. Local
            bool is_global_net = false;
            int first_x = -1;
            int first_y = -1;

    	    for (ClusterPinId pin_id : cluster_ctx.clb_nlist.net_pins(net_id)) {
		auto block_id = cluster_ctx.clb_nlist.pin_block(pin_id);
		//TODO: for now assuming, loc.x will not be 0, which will result in x  = -1, out of grid problem
        	int x = blocks_locs[block_id].loc.x - 1;
        	int y = blocks_locs[block_id].loc.y;
		
		// pin_i is the source coordinate
		if (pin_i == 0) {
                    first_x = x;
                    first_y = y;
                } else {
                    if (x != first_x || y != first_y) {
                        is_global_net = true;
                    }
                }

		if (x < 0 || y < 0 || x > Gx || y > Gy) {
    			VTR_LOG_ERROR("Pin out of bounds! x=%d, y=%d, grid=(%d,%d)\n",
                  	x, y, Gx, Gy);
		}

		PIN p;
                p.x = x;
                p.y = y;
                p.layer = 1;
		unique_pins.insert(PIN{(long)x, (long)y, 1});
    	        ++pin_i;
    	    }
	    
	    if (!is_global_net) {
            	// We SKIP this net entirely.
            	// We do NOT increment netIdx.
            	// We do NOT add to fr_to_vpr_net.
		intra_tile_nets.push_back(net_id);
            	continue;
            }
	
	    stored_net_names.push_back(net_name);
	    std::vector<PIN> temp_pins(unique_pins.begin(), unique_pins.end());
	    fr_pins[netIdx] = temp_pins;
	    int nPins_fr = (int)temp_pins.size();
	    // atleast two unique pins, then it is not an intra-tile net
	    VTR_ASSERT(nPins_fr > 1);

    	    int minWidth = 1; // not really important for FPGA-style use
	    /*VTR_LOG("FastRoute pins for net %zu (%s), nPins=%d\n",
            size_t(net_id), net_name.c_str(), nPins_fr);
            for (int i = 0; i < nPins_fr; ++i) {
        	VTR_LOG("  pin %d: (%ld,%ld,l%d)\n",
                i, temp_pins[i].x, temp_pins[i].y, temp_pins[i].layer);
    	    }*/
	    auto& s = stored_net_names.back();
	    if (s.size() >= 200) {
    		VTR_LOG_ERROR("FastRoute net name too long (%zu chars): %s\n",
                  s.size(), s.c_str());
	    }	

    	    fr_router.addNet(const_cast<char*>(stored_net_names.back().c_str()), netIdx, nPins_fr, minWidth, fr_pins[netIdx].data());
	    fr_to_vpr_net.push_back(net_id);
    	    ++netIdx;
    	}
	VTR_ASSERT(netIdx == valid_global_nets);

	// 4. Run FastRoute
    	fr_router.initEdges();
    	fr_router.setNumAdjustments(0); // no regional capacity tweaks for now
    	fr_router.initAuxVar();

	VTR_LOG("Calling FastRoute::run on %d nets, but total_nets including intra_tile nets %d\n", valid_global_nets, num_nets);
    	std::vector<NET> routed_nets;
    	int fr_status = fr_router.run(routed_nets);
	VTR_LOG("FastRoute::run returned %d, routed_nets.size() = %zu\n",
        	fr_status, routed_nets.size());	
	
	double fr_runtime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - start_time_fr).count();
	auto start_time_global_tree_processing = std::chrono::high_resolution_clock::now();	
	// 5. Extract Steiner trees and populate sb_edges
	ClusterNetId debug_net_id = ClusterNetId(0); 
	for (int fr_idx = 0; fr_idx < routed_nets.size(); fr_idx++) {
    		ClusterNetId net_id = fr_to_vpr_net[fr_idx];
    		const NET& fr_net = routed_nets[fr_idx];


    		bool debug_this_net = (net_id == debug_net_id);

    		if (dump_raw_flute_trees) {
    		    auto net_name = cluster_ctx.clb_nlist.net_name(net_id);
    		    VTR_LOG("==== (%d) FastRoute tree for net %zu (%s) ====\n",
    		            fr_idx, size_t(net_id), net_name.c_str());

    		    VTR_LOG("Number of FR edges: %zu\n", fr_net.route.size());
    		}
    		for (const ROUTE& r : fr_net.route) {
			int x1 = r.initX;
			int y1 = r.initY;
			int x2 = r.finalX;
			int y2 = r.finalY;
			steiner_ctx.net_edges[size_t(net_id)].emplace_back(Point{x1, y1}, Point{x2, y2});

			if (dump_raw_flute_trees) {
                			VTR_LOG("  edge: (%ld,%ld,l%d) -> (%ld,%ld,l%d)\n",
                        		r.initX, r.initY, r.initLayer,
                        		r.finalX, r.finalY, r.finalLayer);
        		}
    		}
    	    	
		// this class may not be usable as it is
		Steiner steiner(cluster_ctx.clb_nlist, net_id, blocks_locs, fluteOutfile, dump_raw_flute_trees, global_router_algorithm);
    	        
		auto [net_corridors, net_corridors_lookahead, connection_intra_tile] = steiner.build_corridor_list_per_connection(make_str_id(steiner.source_x, steiner.source_y));
    	        steiner_ctx.all_corridors[net_id] = std::move(net_corridors);
    	        steiner_ctx.all_corridors_lookahead[net_id] = std::move(net_corridors_lookahead);
    	        steiner_ctx.net_connection_intra_tile[net_id] = std::move(connection_intra_tile);
	}

	for (auto net_id: intra_tile_nets) {
	    
	    Steiner steiner(cluster_ctx.clb_nlist, net_id, blocks_locs, fluteOutfile, dump_raw_flute_trees, global_router_algorithm);
	    auto [net_corridors, net_corridors_lookahead, connection_intra_tile] = steiner.build_corridor_list_per_connection(make_str_id(steiner.source_x, steiner.source_y));
    	    steiner_ctx.all_corridors[net_id] = std::move(net_corridors);
    	    steiner_ctx.all_corridors_lookahead[net_id] = std::move(net_corridors_lookahead);
    	    steiner_ctx.net_connection_intra_tile[net_id] = std::move(connection_intra_tile);
	
	}
	double global_tree_processing_runtime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - start_time_global_tree_processing).count();
        VTR_LOG("[FastRoute] Runtime: %lf ms\n", fr_runtime); 
        VTR_LOG("[FastRoute post processing] Runtime: %lf ms\n", global_tree_processing_runtime); 

    }

    fluteOutfile.close();
    if (create_steiner_constraints) {
        // Create "connections_per_dnode.txt"
        //create_connections_per_dnode_file(steiner_ctx);

        // Create "sink_order.txt" with default net pin orders
        // dummy file, not dumping dependacy graph algorithm
        //create_sink_order_file(cluster_ctx, blocks_locs);

        // Create "branch_node_map.txt"
        //create_branch_node_map_file(steiner_ctx);

        // Remove the gnode maps from the global Steiner context to free up memory
        steiner_ctx.loc_dir_len_to_gnode.clear();
        steiner_ctx.loc_type_to_gnode.clear();
        steiner_ctx.gnode_to_dnodes.clear();
    }
}

// UTILS

int calculate_manhattan_distance(std::string sb1_id, std::string sb2_id) {
    int x1, y1, x2, y2;
    char sep;

    std::istringstream iss(sb1_id);
    iss >> x1 >> sep >> y1;

    iss = std::istringstream(sb2_id);
    iss >> x2 >> sep >> y2;

    return std::abs(x1 - x2) + std::abs(y1 - y2);
}
std::pair<int, int> make_pair_from_str_id(const std::string& str_id) {
    int x, y;
    char sep;
    std::istringstream iss(str_id);
    iss >> x >> sep >> y;
    return std::make_pair(x, y);
}

std::string make_str_id(int x, int y) {
    return std::to_string(x)+"_"+std::to_string(y);
}

std::string make_str_id(int x, int y, int pin_id) {
    return std::to_string(x)+"_"+std::to_string(y)+"_"+std::to_string(pin_id);
}
