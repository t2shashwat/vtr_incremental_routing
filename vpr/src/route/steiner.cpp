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

// For parallelization attempts
// #include <omp.h>

Steiner::Steiner(
    const Netlist<ClusterBlockId, ClusterPortId, ClusterPinId, ClusterNetId>& net_list,
        ClusterNetId net_id,
    const vtr::vector_map<ClusterBlockId, t_block_loc>& blocks_locs,
    Flute::FluteState *flute1) {
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
    //if (this->sb_map.size()>1 && diagonal_detected == true) {
    //    this->build_sb_rsmt_flute(flute1);
    //}
    //
    this->build_sb_rsmt_post_process_flute();

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
bool Steiner::detect_diagonal(Flute::FluteState *flute1) {
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
}

void Steiner::build_sb_rsmt_flute(Flute::FluteState *flute1) {
    Flute::Tree flutetree;
    int d=0;
    int n = this->sb_map.size();
    int x[n], y[n];

    VTR_LOG("Net id: %d\n", size_t(this->net_id));
    VTR_LOG("Source at: (%d, %d)\n", this->source_x, this->source_y);

    for (const auto& sb : this->sb_map) {
        x[d] = sb.second.x;
        y[d] = sb.second.y;
	if (sb.second.pin_id != 0){
            VTR_LOG(" Sink %d: (%d, %d)\n", d, x[d], y[d]);
	}
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
	     VTR_LOG(" FLUTE: (%d, %d) -> (%d, %d)\n", o_x, o_y, i_x, i_y);
	    //}
            this->sb_edges.add_edge(o_x, o_y, i_x, i_y);
            /*if (i_x % 2 == 0) {
                this->add_sb(i_x, o_y);
                this->sb_edges.add_edge(o_x, o_y, i_x, o_y);
                this->sb_edges.add_edge(i_x, o_y, i_x, i_y);
            }
            else {
                this->add_sb(o_x, i_y);
                this->sb_edges.add_edge(o_x, o_y, o_x, i_y);
                this->sb_edges.add_edge(o_x, i_y, i_x, i_y);
            }*/
        } 
        else if (i_x != o_x || i_y != o_y) {
	    //if (size_t(this->net_id) == 415){
	     VTR_LOG(" FLUTE: (%d, %d) -> (%d, %d)\n", o_x, o_y, i_x, i_y);
	    //}
            this->sb_edges.add_edge(o_x, o_y, i_x, i_y);
        }
    }
    int total_wirelength = Flute::flute_wirelength(flutetree);
    VTR_LOG("Total wirelength: %d\n", total_wirelength);

    Flute::free_tree(flute1, flutetree);
    //Flute::free_tree(flutetree);
}

void Steiner::build_sb_rsmt_post_process_flute() {
    
    SteinerContext& steiner_ctx = g_vpr_ctx.mutable_steiner(); 
    int net = size_t(this->net_id);
    VTR_LOG("Net id: %d\n", net);

    // Check if net_id exists in net_edges
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




/*void Steiner::decompose_diagonal() {
   // create adjacency lis for a undirected grapht
   // for each nondiagonal edge get the neighbours
   // the neighbours would be from both end points
   // four categories of neighbours: 1. short-acute 2. short-obtuse 3. long-acute 4. long-obtuse (short is if less than the length of current edge)
   // obtuse diagonals can be left to be processed with edges which form acute angle with them, and if non form an acute angle with them, then they can be processed at the end by checking which diagonals still remain (mark diagonals (?)) 
   // also when the diagonal is processed on an edge, do not remove it from the adjaceny list, as we will decompose it along other edges too, and later we figure out the shortest path to it when we create the tree
   // now, for acute diagonals which are long and short
   //  
   //
   // each node and its neighbours (does not work as we do not which diagonals are acute and obtuse, as it depends on the incoming edge to the node under processing)
   // for each nodes, get the neighbours, if a diagonal is detected (can be long and short), a short diagonal can be decomposed on the edge in the neighbourhood of current node, but for a long diagonal
   // a diagonal can also be 
   // maintain a list of terminals and steiner points, before decomposing a diagonal to new steiner points, check if the point exists, if not then create a new one, append to the existing list and add edges between the points, both horizontal and vertical direction
   // In case of short diagonals, one of the direction already would exist, so that will have to be broken until the new steiner point, which will entail updating the edge map too
   // does this mean updating the adjaceny list? with the new neighbours and adding new steiner points so that their neighbours are also updated (?)
   // if decomposing the diagonal on two edges, then if the other new steiner ends up not being a steiner point, then this point given it is not a termincal should be removed.
   //
   // traverse over all the  
   // add further check to compare total wirelngth of FLUTE tree and post processed tree, to ensure not going too far from optimiality of FLUTE



}*/

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
    const std::string to_gnode_filename = "../../../scripts/node_dict_gr_1L_reg_v4.map";
    const std::string to_dnode_filename = "../../../scripts/gr_dr_map_reg_v4_full_iib.map";
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


    VTR_LOG("Loading loc_to_gnode map:\n");
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

    toGfile.close();

    VTR_LOG("Loading gnode_to_dnode map:\n");
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

    toDfile.close();
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


/*
std::unordered_map<int, std::vector<Corridor>> Steiner::build_corridor_list_per_connection() const {

    std::unordered_map<int, std::vector<Corridor>> corridor_list_per_connection;

    // 1. Build reverse edge map: child → parent
    // cannot use this because the graph has not being directed, and thus no relationship of parent and children is establisehd yet
    std::unordered_map<std::string, std::string> child_to_parent;
    for (const auto& [parent, children] : this->sb_edges.edge_map) {
        for (const auto& [child, is_valid] : children) {
            if (is_valid) {
                child_to_parent[child] = parent;
            }
        }
    }

    // 2. Traverse from each sink SB up to the source
    for (const auto& [str_id, sb] : this->sb_map) {
        if (sb.pin_id <= 0) continue;  // Only terminal sinks

        int sink_id = sb.pin_id;
        std::string current = str_id;
        std::vector<Corridor> corridor_path;

        while (child_to_parent.count(current)) {
            const std::string& parent = child_to_parent.at(current);
            const SB& from_sb = this->sb_map.at(parent);
            const SB& to_sb = this->sb_map.at(current);

            corridor_path.emplace_back(from_sb.x, from_sb.y, to_sb.x, to_sb.y);
            current = parent;
        }

        std::reverse(corridor_path.begin(), corridor_path.end());
        corridor_list_per_connection[sink_id] = std::move(corridor_path);
    }

    return corridor_list_per_connection;
}
*/
std::tuple<std::unordered_map<int, std::vector<Corridor>>, std::unordered_map<int, bool>> Steiner::build_corridor_list_per_connection(std::string source_sb_id) const {

    std::unordered_map<int, std::vector<Corridor>> corridor_list_per_connection;
    std::unordered_map<int, bool> connection_intra_tile;

    std::unordered_map<std::string, std::string> child_to_parent;
    // 1. Build reverse edge map: child → parent

    if (source_sb_id.empty()) {
        // Handle case with no source pin
	VTR_LOG(" [Building corridors: No source for this net\n");
        return {corridor_list_per_connection, connection_intra_tile};
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
    
    // ------------------------------------------------------------------
    
    /*std::set<std::pair<std::string, std::string>> visited; // (node, arrival_dir)
    std::queue<std::pair<std::string, std::string>> q;     // (current_node, parent_node)
    
    q.push({source_sb_id, ""});
    visited.insert({source_sb_id, ""});
    
    while (!q.empty()) {
        auto [current_node, parent_node] = q.front();
        q.pop();
    
        std::string incoming_dir = parent_node.empty() ? "" : direction_from_to(parent_node, current_node);
    
        if (this->sb_edges.edge_map.count(current_node)) {
            for (const auto& [neighbor_node, is_valid] : this->sb_edges.edge_map.at(current_node)) {
    
                if (!is_valid) continue;
    
                std::string outgoing_dir = direction_from_to(current_node, neighbor_node);
    
                // skip invalid (diagonal) edges
                if (outgoing_dir.empty()) continue;
    
                // U-turn check: prevent going back in the opposite direction
                if (!incoming_dir.empty() && is_opposite(incoming_dir, outgoing_dir)) {
                    if (size_t(this->net_id) == 415) {
                        VTR_LOG("Detected u-turn: %s → %s → %s (dir=%s → %s)\n",
                                parent_node.c_str(), current_node.c_str(), neighbor_node.c_str(),
                                incoming_dir.c_str(), outgoing_dir.c_str());
                    }
                    continue;
                }
    
                std::pair<std::string, std::string> state = {neighbor_node, outgoing_dir};
                if (visited.count(state)) continue; // already visited with this direction
    
                // mark visited with arrival direction
                visited.insert(state);
                child_to_parent[neighbor_node] = current_node;
                q.push({neighbor_node, current_node});
            }
        }
    }*/
 
    /*std::unordered_map<std::string, int> cost; // Stores the minimum cost (distance) from source
    
    // Priority queue to find the shortest path.
    // Stores tuples of (cost, current_node, parent_node)
    // std::priority_queue is a max-heap by default, so we store negative costs.
    std::priority_queue<std::tuple<int, std::string, std::string>> pq;

    cost[source_sb_id] = 0;
    pq.push({0, source_sb_id, ""}); // Source has no parent and a cost of 0
    int uturn_detected = 0;
    while (!pq.empty()) {
        auto [current_cost, current_node, parent_node] = pq.top();
        pq.pop();
        
        // Skip if we've found a better path already
        if (-current_cost > cost[current_node]) {
            continue;
        }

        if (this->sb_edges.edge_map.count(current_node)) {
            for (const auto& [neighbor_node, is_valid] : this->sb_edges.edge_map.at(current_node)) {
                
                if (!is_valid) {
                    continue;
                }
		uturn_detected = 0;
                
                // Perform geometric u-turn check  but not for cycles A->B B->A
                if (!parent_node.empty() && neighbor_node != parent_node) {
                    auto [parent_x, parent_y] = parse_coords(parent_node);
                    auto [current_x, current_y] = parse_coords(current_node);
                    auto [neighbor_x, neighbor_y] = parse_coords(neighbor_node);
                    
                    bool is_horizontal_uturn = (parent_y == current_y && current_y == neighbor_y) &&
                                               ((current_x < parent_x && neighbor_x > current_x) || 
                                                (current_x > parent_x && neighbor_x < current_x));
                    
                    bool is_vertical_uturn = (parent_x == current_x && current_x == neighbor_x) &&
                                             ((current_y < parent_y && neighbor_y > current_y) || 
                                              (current_y > parent_y && neighbor_y < current_y));

                    if (is_horizontal_uturn || is_vertical_uturn) {
                        VTR_LOG_WARN("Detected u-turn: %s -> %s -> %s. Skipping this path.\n", 
                            parent_node.c_str(), current_node.c_str(), neighbor_node.c_str());
                        //uturn_detected = 1;
			continue;
                    }
                }

                // Calculate cost to neighbor
                auto [current_x, current_y] = parse_coords(current_node);
                auto [neighbor_x, neighbor_y] = parse_coords(neighbor_node);
		if (size_t(this->net_id) == 415){
		     VTR_LOG(" [DFS] parent: (%d, %d) -> (%d, %d) \n", current_x, current_y, neighbor_x, neighbor_y);
		}

                //int edge_cost = (uturn_detected == 1) ? std::numeric_limits<int>::max() / 4 : abs(current_x - neighbor_x) + abs(current_y - neighbor_y);
                int edge_cost = abs(current_x - neighbor_x) + abs(current_y - neighbor_y);
                int new_cost = -current_cost + edge_cost;
                uturn_detected = 0;
                // If we found a shorter path, update it
                if (!cost.count(neighbor_node) || new_cost < cost[neighbor_node]) {
                    cost[neighbor_node] = new_cost;
                    child_to_parent[neighbor_node] = current_node;
                    pq.push({-new_cost, neighbor_node, current_node});
                }
            }
        }
    }*/




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
                
                // Perform geometric u-turn check
                /*if (!parent_node.empty()) {
                    auto [parent_x, parent_y] = parse_coords(parent_node);
                    auto [current_x, current_y] = parse_coords(current_node);
                    auto [neighbor_x, neighbor_y] = parse_coords(neighbor_node);

                    // Check for horizontal u-turn
                    bool is_horizontal_uturn = (parent_y == current_y && current_y == neighbor_y) &&
                                               ((current_x < parent_x && neighbor_x > current_x) || 
                                                (current_x > parent_x && neighbor_x < current_x));
                    
                    // Check for vertical u-turn
                    bool is_vertical_uturn = (parent_x == current_x && current_x == neighbor_x) &&
                                             ((current_y < parent_y && neighbor_y > current_y) || 
                                              (current_y > parent_y && neighbor_y < current_y));

                    if (is_horizontal_uturn || is_vertical_uturn) {
			if (size_t(this->net_id) == 415)
			    //VTR_LOG("Detected horizontal u-turn: %s -> %s -> %s. Path segment ignored.\n", parent_node.c_str(), current_node.c_str(), neighbor_node.c_str());
                        continue; // Skip this edge to avoid a u-turn
                    }
                }*/
                
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

	bool path_found = false;//(child_to_parent.count(current) || current == source_sb_id);
        //if (path_found) {
        while (child_to_parent.count(current)) {
            const std::string& parent = child_to_parent.at(current);
            auto [from_x, from_y] = parse_coords(parent);
            auto [to_x, to_y] = parse_coords(current);

            corridor_path.emplace_back(from_x, from_y, to_x, to_y);

            if (parent == source_sb_id) {
                path_found = true;
                break;
            }
            current = parent;
        }

        if (path_found) {
            std::reverse(corridor_path.begin(), corridor_path.end());
            corridor_list_per_connection[sink_id] = std::move(corridor_path);
        } else if (current != source_sb_id) {
            VTR_LOG("Sink pin_id %d at %d is disconnected from the source. No path found.\n", sink_id, this->net_id);
        }
	else if (current == source_sb_id){ // intra-CLB connections having no edges, but we need to initialize the connection
	    corridor_list_per_connection[sink_id] = std::move(corridor_path);
	}
    }
    
    return {corridor_list_per_connection, connection_intra_tile};
}


//    VTR_LOG("Building child_to_parent map from sb_edges:\n");
    /*for (const auto& [parent, children] : this->sb_edges.edge_map) {
        for (const auto& [child, is_valid] : children) {
            //VTR_LOG("  Edge: %s -> %s = %d\n", parent.c_str(), child.c_str(), is_valid);
            if (is_valid) {
                child_to_parent[child] = parent;
            }
        }
    }*/
//    VTR_LOG("  Total valid edges in child_to_parent: %zu\n", child_to_parent.size());

    // 2. Traverse from each sink SB up to the source
    /*for (const auto& [str_id, sb] : this->sb_map) {
//        VTR_LOG("Checking SB: %s, pin_id = %d\n", str_id.c_str(), sb.pin_id);
        if (sb.pin_id <= 0) {
//            VTR_LOG("  Skipping non-terminal SB.\n");
            continue;
        }

        int sink_id = sb.pin_id;
	std::string current = std::to_string(sb.x) + "_" + std::to_string(sb.y);
        std::vector<Corridor> corridor_path;

//        VTR_LOG("  Traversing back from sink %d (%s):\n", sink_id, current.c_str());

        while (child_to_parent.count(current)) {
            const std::string& parent = child_to_parent.at(current);
//            VTR_LOG("    Parent of %s is %s\n", current.c_str(), parent.c_str());
    	    // Unpack parent: "x_y"
    	    int from_x, from_y;
    	    char sep;
    	    std::istringstream iss(parent);
    	    iss >> from_x >> sep >> from_y;
    	    
	    int to_x, to_y;
    	    std::istringstream iss2(current);
    	    iss2 >> to_x >> sep >> to_y;

            //const SB& from_sb = this->sb_map.at(parent);
            //const SB& to_sb = this->sb_map.at(current);
//            VTR_LOG("      from_sb: (%d, %d), to_sb: (%d, %d)\n", from_x, from_y, to_x, to_y);

            corridor_path.emplace_back(from_x, from_y, to_x, to_y);
            current = parent;
        }

//        if (corridor_path.empty()) {
//            VTR_LOG("    Warning: No path found for sink_id %d\n", sink_id);
//        } else {
            std::reverse(corridor_path.begin(), corridor_path.end());
            corridor_list_per_connection[sink_id] = std::move(corridor_path);
//            VTR_LOG("    Added path of %zu corridors for sink_id %d\n", corridor_list_per_connection[sink_id].size(), sink_id);
//        }
          if(corridor_list_per_connection.size() == 0){
		  VTR_LOG("Empty corridor list for: net: %d sink: %d \n",this->net_id, sink_id);
	  }
    }

//    VTR_LOG("Finished building corridor_list_per_connection. Total sinks: %zu\n", corridor_list_per_connection.size());
    
    return corridor_list_per_connection;
}*/



void steiner_pre_processing(bool create_steiner_constraints, bool compute_dependency_graph_sink_orders) {
    // Initialize context refereces
    const ClusteringContext& cluster_ctx = g_vpr_ctx.clustering();
    const PlacementContext& placement_ctx = g_vpr_ctx.placement();
    SteinerContext& steiner_ctx = g_vpr_ctx.mutable_steiner();

    // Initialize FLUTE
    Flute::FluteState *flute1 = Flute::flute_init(FLUTE_POWVFILE, FLUTE_PORTFILE);
    VTR_LOG("FLUTE initialized.\n");

    // Load gnode maps onto global Steiner context
    if (create_steiner_constraints)
        load_gnode_maps(steiner_ctx);

    // Nets belonging to the "global connecting nets" are not routed by the global router
    // net is global should replace this
    //std::set<size_t> nets_to_skip = load_nets_to_skip();

    // Create RSMT for each net in the clustered netlist (except for
    // "global connecting nets") and update the global Steiner context 
    // with information about the corresponding constrained regions
    const vtr::vector_map<ClusterBlockId, t_block_loc>& blocks_locs = g_vpr_ctx.placement().block_locs;

    for (ClusterNetId net_id : cluster_ctx.clb_nlist.nets()) {
        // Skip "global connecting nets"; this might also be done by checking clb_nlist.net_is_global(net_id), but I haven't been able to get that to work
	//unsigned int num_sinks = cluster_ctx.clb_nlist.net_sinks(net_id).size()
	//VTR_LOG("Net id: %d (num_sinks: %d)\n", size_t(net_id), num_sinks);
        if (cluster_ctx.clb_nlist.net_is_ignored(net_id)) {
            printf("Net ID: %d\n", size_t(net_id));
            continue;
        }

        // Construct the RSMT for the net
        Steiner steiner(cluster_ctx.clb_nlist, net_id, blocks_locs, flute1);

        if (create_steiner_constraints) {
	    // LUKA's implementation that fits in the FCCM implementation; but this implementation leads to high runtime due to repeated memory accesses
            // Create RSMT constrained regions for the net and
            // update the global Steiner context with the information 
            //steiner.create_coarsened_regions(make_str_id(steiner.source_x, steiner.source_y));

            // The RSMT is used to compute the nearest neighbour sink order;
            // by doing this, we reduce heap pops and redundant wires
	    //
	    // Below is the new implemenation driven to optimize runtime by reducing memory accesses
	    auto [net_corridors, connection_intra_tile] = steiner.build_corridor_list_per_connection(make_str_id(steiner.source_x, steiner.source_y));
            steiner_ctx.all_corridors[net_id] = std::move(net_corridors);
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
