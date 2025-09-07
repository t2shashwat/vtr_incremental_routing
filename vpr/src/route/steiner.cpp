#include "steiner.h"
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <algorithm>

#include "vpr_types.h"
#include "place_macro.h"
#include "vpr_context.h"

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
    
    // Build SB RSMT
    if (this->sb_map.size()>1) {
        this->build_sb_rsmt_flute(flute1);
    }

    // Make tree directed

    // Find the source to start the traversal from
    SB source_sb;
    for (const auto& sb : sb_map) {
        if (sb.second.pin_id == 0) {
            source_sb = SB(sb.second.x, sb.second.y, sb.second.pin_id); // Source SB
            break;
        }
    }
    std::set<std::string> visited;
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
    edges_to_remove.clear();

// FOR DEBUGGING PURPOSES: Print all edges in the RSMT so that you can plot it.
    // for (const auto& edge : this->sb_edges.edge_map) {
    //     for (const auto& sink_edge : edge.second) {
    //         if (sink_edge.second) { // Only print true edges
    //             VTR_LOG("Edge: %s -> %s\n", edge.first.c_str(), sink_edge.first.c_str());
    //         }
    //     }
    // }

    this->source_x = source_sb.x;
    this->source_y = source_sb.y;
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

void Steiner::build_sb_rsmt_flute(Flute::FluteState *flute1) {
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
            if (i_x % 2 == 0) {
                this->add_sb(i_x, o_y);
                this->sb_edges.add_edge(o_x, o_y, i_x, o_y);
                this->sb_edges.add_edge(i_x, o_y, i_x, i_y);
            }
            else {
                this->add_sb(o_x, i_y);
                this->sb_edges.add_edge(o_x, o_y, o_x, i_y);
                this->sb_edges.add_edge(o_x, i_y, i_x, i_y);
            }
        } 
        else if (i_x != o_x || i_y != o_y) {
            this->sb_edges.add_edge(o_x, o_y, i_x, i_y);
        }
    }

    Flute::free_tree(flute1, flutetree);
}

void Steiner::make_tree_graph_directed(std::string source_sb_id, std::set<std::string>& visited, std::set<std::string>& edges_to_remove) {    
    visited.insert(source_sb_id);
    if (this->sb_edges.edge_map.find(source_sb_id) == this->sb_edges.edge_map.end()) {
        return; // No edges to process
    }
    for (auto& sink_sb : this->sb_edges.edge_map[source_sb_id]) {
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
    } 
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

    std::string line;
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
std::unordered_map<int, std::vector<Corridor>> Steiner::build_corridor_list_per_connection() const {

    std::unordered_map<int, std::vector<Corridor>> corridor_list_per_connection;

    // 1. Build reverse edge map: child → parent
    std::unordered_map<std::string, std::string> child_to_parent;

//    VTR_LOG("Building child_to_parent map from sb_edges:\n");
    for (const auto& [parent, children] : this->sb_edges.edge_map) {
        for (const auto& [child, is_valid] : children) {
//            VTR_LOG("  Edge: %s -> %s = %d\n", parent.c_str(), child.c_str(), is_valid);
            if (is_valid) {
                child_to_parent[child] = parent;
            }
        }
    }
//    VTR_LOG("  Total valid edges in child_to_parent: %zu\n", child_to_parent.size());

    // 2. Traverse from each sink SB up to the source
    for (const auto& [str_id, sb] : this->sb_map) {
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
    }

//    VTR_LOG("Finished building corridor_list_per_connection. Total sinks: %zu\n", corridor_list_per_connection.size());

    return corridor_list_per_connection;
}



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
	    auto net_corridors = steiner.build_corridor_list_per_connection();
            steiner_ctx.all_corridors[net_id] = std::move(net_corridors);
	    // Debug print: Dump corridors for this net
	    if (size_t(net_id) == 415 || size_t(net_id) == 0 || size_t(net_id) == 1) {
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
	    }
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
