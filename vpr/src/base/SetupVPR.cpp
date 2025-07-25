#include <cstring>
#include <vector>
#include <sstream>

#include "vtr_assert.h"
#include "vtr_util.h"
#include "vtr_random.h"
#include "vtr_log.h"
#include "vtr_memory.h"
#include "vtr_time.h"

#include "vpr_types.h"
#include "vpr_error.h"
#include "vpr_utils.h"

#include "globals.h"
#include "read_xml_arch_file.h"
#include "read_fpga_interchange_arch.h"
#include "SetupVPR.h"
#include "pb_type_graph.h"
#include "pack_types.h"
#include "lb_type_rr_graph.h"
#include "rr_graph_area.h"
#include "echo_arch.h"
#include "read_options.h"
#include "echo_files.h"
#include "clock_modeling.h"
#include "ShowSetup.h"

static void SetupNetlistOpts(const t_options& Options, t_netlist_opts& NetlistOpts);
static void SetupPackerOpts(const t_options& Options,
                            t_packer_opts* PackerOpts);
static void SetupPlacerOpts(const t_options& Options,
                            t_placer_opts* PlacerOpts);
static void SetupAnnealSched(const t_options& Options,
                             t_annealing_sched* AnnealSched);
static void SetupRouterOpts(const t_options& Options, t_router_opts* RouterOpts);
static void SetupNocOpts(const t_options& Options,
                         t_noc_opts* NocOpts);
static void SetupRoutingArch(const t_arch& Arch, t_det_routing_arch* RoutingArch);
static void SetupTiming(const t_options& Options, const bool TimingEnabled, t_timing_inf* Timing);
static void SetupSwitches(const t_arch& Arch,
                          t_det_routing_arch* RoutingArch,
                          const t_arch_switch_inf* ArchSwitches,
                          int NumArchSwitches);
static void SetupAnalysisOpts(const t_options& Options, t_analysis_opts& analysis_opts);
static void SetupPowerOpts(const t_options& Options, t_power_opts* power_opts, t_arch* Arch);
static int find_ipin_cblock_switch_index(const t_arch& Arch);
static void alloc_and_load_intra_cluster_resources();

/**
 * @brief Sets VPR parameters and defaults.
 *
 * Does not do any error checking as this should have been done by
 * the various input checkers
 */
void SetupVPR(const t_options* Options,
              const bool TimingEnabled,
              const bool readArchFile,
              t_file_name_opts* FileNameOpts,
              t_arch* Arch,
              t_model** user_models,
              t_model** library_models,
              t_netlist_opts* NetlistOpts,
              t_packer_opts* PackerOpts,
              t_placer_opts* PlacerOpts,
              t_annealing_sched* AnnealSched,
              t_router_opts* RouterOpts,
              t_analysis_opts* AnalysisOpts,
              t_noc_opts* NocOpts,
              t_det_routing_arch* RoutingArch,
              std::vector<t_lb_type_rr_node>** PackerRRGraphs,
              std::vector<t_segment_inf>& Segments,
              t_timing_inf* Timing,
              bool* ShowGraphics,
              int* GraphPause,
              bool* SaveGraphics,
              std::string* GraphicsCommands,
              t_power_opts* PowerOpts,
              t_vpr_setup* vpr_setup) {
    using argparse::Provenance;

    auto& device_ctx = g_vpr_ctx.mutable_device();

    if (Options->CircuitName.value() == "") {
        VPR_FATAL_ERROR(VPR_ERROR_BLIF_F,
                        "No blif file found in arguments (did you specify an architecture file?)\n");
    }

    alloc_and_load_output_file_names(Options->CircuitName);

    //TODO: Move FileNameOpts setup into separate function
    FileNameOpts->CircuitName = Options->CircuitName;
    FileNameOpts->ArchFile = Options->ArchFile;
    FileNameOpts->CircuitFile = Options->CircuitFile;
    FileNameOpts->NetFile = Options->NetFile;
    FileNameOpts->PlaceFile = Options->PlaceFile;
    FileNameOpts->RouteFile = Options->RouteFile;
    FileNameOpts->ActFile = Options->ActFile;
    FileNameOpts->PowerFile = Options->PowerFile;
    FileNameOpts->CmosTechFile = Options->CmosTechFile;
    FileNameOpts->out_file_prefix = Options->out_file_prefix;
    FileNameOpts->read_vpr_constraints_file = Options->read_vpr_constraints_file;
    FileNameOpts->write_vpr_constraints_file = Options->write_vpr_constraints_file;
    FileNameOpts->write_block_usage = Options->write_block_usage;

    FileNameOpts->verify_file_digests = Options->verify_file_digests;

    SetupNetlistOpts(*Options, *NetlistOpts);
    SetupPlacerOpts(*Options, PlacerOpts);
    SetupAnnealSched(*Options, AnnealSched);
    SetupRouterOpts(*Options, RouterOpts);
    SetupAnalysisOpts(*Options, *AnalysisOpts);
    SetupPowerOpts(*Options, PowerOpts, Arch);
    SetupNocOpts(*Options, NocOpts);

    if (readArchFile == true) {
        vtr::ScopedStartFinishTimer t("Loading Architecture Description");
        switch (Options->arch_format) {
            case e_arch_format::VTR:
                XmlReadArch(Options->ArchFile.value().c_str(),
                            TimingEnabled,
                            Arch,
                            device_ctx.physical_tile_types,
                            device_ctx.logical_block_types);
                break;
            case e_arch_format::FPGAInterchange:
                VTR_LOG("Use FPGA Interchange device\n");
                FPGAInterchangeReadArch(Options->ArchFile.value().c_str(),
                                        TimingEnabled,
                                        Arch,
                                        device_ctx.physical_tile_types,
                                        device_ctx.logical_block_types);
                break;
            default:
                VPR_FATAL_ERROR(VPR_ERROR_ARCH, "Invalid architecture format!");
        }
    }
    VTR_LOG("\n");

    *user_models = Arch->models;
    *library_models = Arch->model_library;

    device_ctx.EMPTY_PHYSICAL_TILE_TYPE = nullptr;
    int num_inputs = 0;
    int num_outputs = 0;
    for (auto& type : device_ctx.physical_tile_types) {
        if (type.is_empty()) {
            VTR_ASSERT(device_ctx.EMPTY_PHYSICAL_TILE_TYPE == nullptr);
            VTR_ASSERT(type.num_pins == 0);
            device_ctx.EMPTY_PHYSICAL_TILE_TYPE = &type;
        }

        if (type.is_input_type) {
            num_inputs += 1;
        }

        if (type.is_output_type) {
            num_outputs += 1;
        }
    }

    device_ctx.EMPTY_LOGICAL_BLOCK_TYPE = nullptr;
    int max_equivalent_tiles = 0;
    for (const auto& type : device_ctx.logical_block_types) {
        if (type.is_empty()) {
            VTR_ASSERT(device_ctx.EMPTY_LOGICAL_BLOCK_TYPE == nullptr);
            VTR_ASSERT(type.pb_type == nullptr);
            device_ctx.EMPTY_LOGICAL_BLOCK_TYPE = &type;
        }

        max_equivalent_tiles = std::max(max_equivalent_tiles, (int)type.equivalent_tiles.size());
    }

    VTR_ASSERT(max_equivalent_tiles > 0);
    device_ctx.has_multiple_equivalent_tiles = max_equivalent_tiles > 1;

    VTR_ASSERT(device_ctx.EMPTY_PHYSICAL_TILE_TYPE != nullptr);
    VTR_ASSERT(device_ctx.EMPTY_LOGICAL_BLOCK_TYPE != nullptr);

    if (num_inputs == 0) {
        VPR_ERROR(VPR_ERROR_ARCH,
                  "Architecture contains no top-level block type containing '.input' models");
    }

    if (num_outputs == 0) {
        VPR_ERROR(VPR_ERROR_ARCH,
                  "Architecture contains no top-level block type containing '.output' models");
    }

    Segments = Arch->Segments;

    SetupSwitches(*Arch, RoutingArch, Arch->Switches, Arch->num_switches);
    SetupRoutingArch(*Arch, RoutingArch);
    SetupTiming(*Options, TimingEnabled, Timing);
    SetupPackerOpts(*Options, PackerOpts);
    RoutingArch->write_rr_graph_filename = Options->write_rr_graph_file;
    RoutingArch->read_rr_graph_filename = Options->read_rr_graph_file;

    //Setup the default flow, if no specific stages specified
    //do all
    if (!Options->do_packing
        && !Options->do_placement
        && !Options->do_routing
        && !Options->do_analysis) {
        //run all stages if none specified
        PackerOpts->doPacking = STAGE_DO;
        PlacerOpts->doPlacement = STAGE_DO;
        RouterOpts->doRouting = STAGE_DO;
        AnalysisOpts->doAnalysis = STAGE_AUTO; //Deferred until implementation status known
    } else {
        //We run all stages up to the specified stage
        //Note that by checking in reverse order (i.e. analysis to packing)
        //we ensure that earlier stages override the default 'LOAD' action
        //set by later stages

        if (Options->do_analysis) {
            PackerOpts->doPacking = STAGE_LOAD;
            PlacerOpts->doPlacement = STAGE_LOAD;
            RouterOpts->doRouting = STAGE_LOAD;
            AnalysisOpts->doAnalysis = STAGE_DO;
        }

        if (Options->do_routing) {
            PackerOpts->doPacking = STAGE_LOAD;
            PlacerOpts->doPlacement = STAGE_LOAD;
            RouterOpts->doRouting = STAGE_DO;
            AnalysisOpts->doAnalysis = ((Options->do_analysis) ? STAGE_DO : STAGE_AUTO); //Always run analysis after routing
        }

        if (Options->do_placement) {
            PackerOpts->doPacking = STAGE_LOAD;
            PlacerOpts->doPlacement = STAGE_DO;
        }

        if (Options->do_packing) {
            PackerOpts->doPacking = STAGE_DO;
        }
    }

    ShowSetup(*vpr_setup);

    /* init global variables */
    vtr::out_file_prefix = Options->out_file_prefix;

    /* Set seed for pseudo-random placement, default seed to 1 */
    PlacerOpts->seed = Options->Seed;
    vtr::srandom(PlacerOpts->seed);

    {
        vtr::ScopedStartFinishTimer t("Building complex block graph");
        alloc_and_load_all_pb_graphs(PowerOpts->do_power, RouterOpts->flat_routing);
        *PackerRRGraphs = alloc_and_load_all_lb_type_rr_graph();
    }

    if (RouterOpts->flat_routing) {
        alloc_and_load_intra_cluster_resources();
    }

    if ((Options->clock_modeling == ROUTED_CLOCK) || (Options->clock_modeling == DEDICATED_NETWORK)) {
        ClockModeling::treat_clock_pins_as_non_globals();
    }

    if (getEchoEnabled() && isEchoFileEnabled(E_ECHO_LB_TYPE_RR_GRAPH)) {
        echo_lb_type_rr_graphs(getEchoFileName(E_ECHO_LB_TYPE_RR_GRAPH), *PackerRRGraphs);
    }

    if (getEchoEnabled() && isEchoFileEnabled(E_ECHO_PB_GRAPH)) {
        echo_pb_graph(getEchoFileName(E_ECHO_PB_GRAPH));
    }

    *GraphPause = Options->GraphPause;

    *ShowGraphics = Options->show_graphics;

    *SaveGraphics = Options->save_graphics;
    *GraphicsCommands = Options->graphics_commands;

    if (getEchoEnabled() && isEchoFileEnabled(E_ECHO_ARCH)) {
        EchoArch(getEchoFileName(E_ECHO_ARCH), device_ctx.physical_tile_types, device_ctx.logical_block_types, Arch);
    }
}

static void SetupTiming(const t_options& Options, const bool TimingEnabled, t_timing_inf* Timing) {
    /* Don't do anything if they don't want timing */
    if (false == TimingEnabled) {
        Timing->timing_analysis_enabled = false;
        return;
    }

    Timing->timing_analysis_enabled = TimingEnabled;
    Timing->SDCFile = Options.SDCFile;
}

/**
 * @brief This loads up VPR's arch_switch_inf data by combining the switches
 *        from the arch file with the special switches that VPR needs.
 */
static void SetupSwitches(const t_arch& Arch,
                          t_det_routing_arch* RoutingArch,
                          const t_arch_switch_inf* ArchSwitches,
                          int NumArchSwitches) {
    auto& device_ctx = g_vpr_ctx.mutable_device();

    int switches_to_copy = NumArchSwitches;
    device_ctx.num_arch_switches = NumArchSwitches;

    RoutingArch->wire_to_arch_ipin_switch = find_ipin_cblock_switch_index(Arch);

    /* Depends on device_ctx.num_arch_switches */
    RoutingArch->delayless_switch = device_ctx.num_arch_switches++;

    /* Alloc the list now that we know the final num_arch_switches value */
    device_ctx.arch_switch_inf = new t_arch_switch_inf[device_ctx.num_arch_switches];
    for (int iswitch = 0; iswitch < switches_to_copy; iswitch++) {
        device_ctx.arch_switch_inf[iswitch] = ArchSwitches[iswitch];
    }

    /* Delayless switch for connecting sinks and sources with their pins. */
    device_ctx.arch_switch_inf[RoutingArch->delayless_switch].set_type(SwitchType::MUX);
    device_ctx.arch_switch_inf[RoutingArch->delayless_switch].name = vtr::strdup(VPR_DELAYLESS_SWITCH_NAME);
    device_ctx.arch_switch_inf[RoutingArch->delayless_switch].R = 0.;
    device_ctx.arch_switch_inf[RoutingArch->delayless_switch].Cin = 0.;
    device_ctx.arch_switch_inf[RoutingArch->delayless_switch].Cout = 0.;
    device_ctx.arch_switch_inf[RoutingArch->delayless_switch].set_Tdel(t_arch_switch_inf::UNDEFINED_FANIN, 0.);
    device_ctx.arch_switch_inf[RoutingArch->delayless_switch].power_buffer_type = POWER_BUFFER_TYPE_NONE;
    device_ctx.arch_switch_inf[RoutingArch->delayless_switch].mux_trans_size = 0.;
    device_ctx.arch_switch_inf[RoutingArch->delayless_switch].buf_size_type = BufferSize::ABSOLUTE;
    device_ctx.arch_switch_inf[RoutingArch->delayless_switch].buf_size = 0.;
    VTR_ASSERT_MSG(device_ctx.arch_switch_inf[RoutingArch->delayless_switch].buffered(), "Delayless switch expected to be buffered (isolating)");
    VTR_ASSERT_MSG(device_ctx.arch_switch_inf[RoutingArch->delayless_switch].configurable(), "Delayless switch expected to be configurable");

    RoutingArch->global_route_switch = RoutingArch->delayless_switch;

    //Warn about non-zero Cout values for the ipin switch, since these values have no effect.
    //VPR do not model the R/C's of block internal routing connectsion.
    //
    //Note that we don't warn about the R value as it may be used to size the buffer (if buf_size_type is AUTO)
    if (device_ctx.arch_switch_inf[RoutingArch->wire_to_arch_ipin_switch].Cout != 0.) {
        VTR_LOG_WARN("Non-zero switch output capacitance (%g) has no effect when switch '%s' is used for connection block inputs\n",
                     device_ctx.arch_switch_inf[RoutingArch->wire_to_arch_ipin_switch].Cout, Arch.ipin_cblock_switch_name.c_str());
    }
}

/**
 * @brief Sets up routing structures.
 *
 * Since checks are already done, this just copies values across
 */
static void SetupRoutingArch(const t_arch& Arch,
                             t_det_routing_arch* RoutingArch) {
    RoutingArch->switch_block_type = Arch.SBType;
    RoutingArch->R_minW_nmos = Arch.R_minW_nmos;
    RoutingArch->R_minW_pmos = Arch.R_minW_pmos;
    RoutingArch->Fs = Arch.Fs;
    RoutingArch->directionality = BI_DIRECTIONAL;
    if (Arch.Segments.size()) {
        RoutingArch->directionality = Arch.Segments[0].directionality;
    }

    /* copy over the switch block information */
    RoutingArch->switchblocks = Arch.switchblocks;
}

static void SetupRouterOpts(const t_options& Options, t_router_opts* RouterOpts) {
    RouterOpts->do_check_rr_graph = Options.check_rr_graph;
    RouterOpts->astar_fac = Options.astar_fac;
    RouterOpts->router_profiler_astar_fac = Options.router_profiler_astar_fac;
    RouterOpts->bb_factor = Options.bb_factor;
    RouterOpts->criticality_exp = Options.criticality_exp;
    RouterOpts->max_criticality = Options.max_criticality;
    RouterOpts->max_router_iterations = Options.max_router_iterations;
    RouterOpts->init_wirelength_abort_threshold = Options.router_init_wirelength_abort_threshold;
    RouterOpts->min_incremental_reroute_fanout = Options.min_incremental_reroute_fanout;
    RouterOpts->incr_reroute_delay_ripup = Options.incr_reroute_delay_ripup;
    RouterOpts->pres_fac_mult = Options.pres_fac_mult;
    RouterOpts->route_type = Options.RouteType;

    RouterOpts->full_stats = Options.full_stats;
    //mycode==================================
    RouterOpts->incr_route = Options.incr_route;
    RouterOpts->icr_iter = Options.icr_iter;
    RouterOpts->sbNode_lookahead_factor = Options.sbNode_lookahead_factor;
    RouterOpts->detailed_router = Options.detailed_router;
    RouterOpts->tree_type= Options.tree_type;
    RouterOpts->shuffle1 = Options.shuffle1;
    RouterOpts->shuffle2 = Options.shuffle2;
    RouterOpts->save_history_cost_per_iteration= Options.save_history_cost_per_iteration;
    RouterOpts->ripup_all_nets = Options.ripup_all_nets;
    RouterOpts->shuffle_net_order = Options.shuffle_net_order;
    RouterOpts->nets_to_skip = Options.nets_to_skip;
    RouterOpts->preorder_sink_order = Options.preorder_sink_order;
    RouterOpts->relax_hop_order = Options.relax_hop_order;
    RouterOpts->global_occ_factor = Options.global_occ_factor;
    RouterOpts->load_gr_history = Options.load_gr_history;
    RouterOpts->offpath_penalty = Options.offpath_penalty;
    //========================================
    // (PARSA) Luka, 2025
    RouterOpts->steiner_constraints = Options.steiner_constraints;
    RouterOpts->dependency_graph_sink_order = Options.dependency_graph_sink_order;
    RouterOpts->shuffle_first_iteration = Options.shuffle_first_iteration;
    RouterOpts->target_bracket = Options.target_bracket;
    //==========================================

    //TODO document these?
    RouterOpts->congestion_analysis = Options.full_stats;
    RouterOpts->fanout_analysis = Options.full_stats;
    RouterOpts->switch_usage_analysis = Options.full_stats;

    RouterOpts->verify_binary_search = Options.verify_binary_search;
    RouterOpts->router_algorithm = Options.RouterAlgorithm;
    RouterOpts->fixed_channel_width = Options.RouteChanWidth;
    RouterOpts->min_channel_width_hint = Options.min_route_chan_width_hint;
    RouterOpts->read_rr_edge_metadata = Options.read_rr_edge_metadata;
    RouterOpts->reorder_rr_graph_nodes_algorithm = Options.reorder_rr_graph_nodes_algorithm;
    RouterOpts->reorder_rr_graph_nodes_threshold = Options.reorder_rr_graph_nodes_threshold;
    RouterOpts->reorder_rr_graph_nodes_seed = Options.reorder_rr_graph_nodes_seed;

    RouterOpts->initial_pres_fac = Options.initial_pres_fac;
    RouterOpts->base_cost_type = Options.base_cost_type;
    RouterOpts->first_iter_pres_fac = Options.first_iter_pres_fac;
    RouterOpts->acc_fac = Options.acc_fac;
    RouterOpts->bend_cost = Options.bend_cost;
    if (Options.do_routing) {
        RouterOpts->doRouting = STAGE_DO;
    }
    RouterOpts->routing_failure_predictor = Options.routing_failure_predictor;
    RouterOpts->routing_budgets_algorithm = Options.routing_budgets_algorithm;
    RouterOpts->save_routing_per_iteration = Options.save_routing_per_iteration;
    RouterOpts->congested_routing_iteration_threshold_frac = Options.congested_routing_iteration_threshold_frac;
    RouterOpts->route_bb_update = Options.route_bb_update;
    RouterOpts->clock_modeling = Options.clock_modeling;
    RouterOpts->two_stage_clock_routing = Options.two_stage_clock_routing;
    RouterOpts->high_fanout_threshold = Options.router_high_fanout_threshold;
    RouterOpts->high_fanout_max_slope = Options.router_high_fanout_max_slope;
    RouterOpts->router_debug_net = Options.router_debug_net;
    RouterOpts->router_debug_sink_rr = Options.router_debug_sink_rr;
    RouterOpts->router_debug_iteration = Options.router_debug_iteration;
    RouterOpts->lookahead_type = Options.router_lookahead_type;
    RouterOpts->max_convergence_count = Options.router_max_convergence_count;
    RouterOpts->reconvergence_cpd_threshold = Options.router_reconvergence_cpd_threshold;
    RouterOpts->initial_timing = Options.router_initial_timing;
    RouterOpts->update_lower_bound_delays = Options.router_update_lower_bound_delays;
    RouterOpts->first_iteration_timing_report_file = Options.router_first_iteration_timing_report_file;
    RouterOpts->strict_checks = Options.strict_checks;

    RouterOpts->write_router_lookahead = Options.write_router_lookahead;
    RouterOpts->read_router_lookahead = Options.read_router_lookahead;

    RouterOpts->router_heap = Options.router_heap;
    RouterOpts->exit_after_first_routing_iteration = Options.exit_after_first_routing_iteration;

    RouterOpts->check_route = Options.check_route;
    RouterOpts->timing_update_type = Options.timing_update_type;

    RouterOpts->max_logged_overused_rr_nodes = Options.max_logged_overused_rr_nodes;
    RouterOpts->generate_rr_node_overuse_report = Options.generate_rr_node_overuse_report;
    RouterOpts->flat_routing = Options.flat_routing;
}

static void SetupAnnealSched(const t_options& Options,
                             t_annealing_sched* AnnealSched) {
    AnnealSched->alpha_t = Options.PlaceAlphaT;
    if (AnnealSched->alpha_t >= 1 || AnnealSched->alpha_t <= 0) {
        VPR_FATAL_ERROR(VPR_ERROR_OTHER, "alpha_t must be between 0 and 1 exclusive.\n");
    }

    AnnealSched->exit_t = Options.PlaceExitT;
    if (AnnealSched->exit_t <= 0) {
        VPR_FATAL_ERROR(VPR_ERROR_OTHER, "exit_t must be greater than 0.\n");
    }

    AnnealSched->init_t = Options.PlaceInitT;
    if (AnnealSched->init_t <= 0) {
        VPR_FATAL_ERROR(VPR_ERROR_OTHER, "init_t must be greater than 0.\n");
    }

    if (AnnealSched->init_t < AnnealSched->exit_t) {
        VPR_FATAL_ERROR(VPR_ERROR_OTHER, "init_t must be greater or equal to than exit_t.\n");
    }

    AnnealSched->inner_num = Options.PlaceInnerNum;
    if (AnnealSched->inner_num <= 0) {
        VPR_FATAL_ERROR(VPR_ERROR_OTHER, "inner_num must be greater than 0.\n");
    }

    AnnealSched->alpha_min = Options.PlaceAlphaMin;
    if (AnnealSched->alpha_min >= 1 || AnnealSched->alpha_min <= 0) {
        VPR_FATAL_ERROR(VPR_ERROR_OTHER, "alpha_min must be between 0 and 1 exclusive.\n");
    }

    AnnealSched->alpha_max = Options.PlaceAlphaMax;
    if (AnnealSched->alpha_max >= 1 || AnnealSched->alpha_max <= AnnealSched->alpha_min) {
        VPR_FATAL_ERROR(VPR_ERROR_OTHER, "alpha_max must be between alpha_min and 1 exclusive.\n");
    }

    AnnealSched->alpha_decay = Options.PlaceAlphaDecay;
    if (AnnealSched->alpha_decay >= 1 || AnnealSched->alpha_decay <= 0) {
        VPR_FATAL_ERROR(VPR_ERROR_OTHER, "alpha_decay must be between 0 and 1 exclusive.\n");
    }

    AnnealSched->success_min = Options.PlaceSuccessMin;
    if (AnnealSched->success_min >= 1 || AnnealSched->success_min <= 0) {
        VPR_FATAL_ERROR(VPR_ERROR_OTHER, "success_min must be between 0 and 1 exclusive.\n");
    }

    AnnealSched->success_target = Options.PlaceSuccessTarget;
    if (AnnealSched->success_target >= 1 || AnnealSched->success_target <= 0) {
        VPR_FATAL_ERROR(VPR_ERROR_OTHER, "success_target must be between 0 and 1 exclusive.\n");
    }

    AnnealSched->type = Options.anneal_sched_type;
}

/**
 * @brief Sets up the s_packer_opts structure baesd on users inputs and
 *        on the architecture specified.
 *
 * Error checking, such as checking for conflicting params is assumed
 * to be done beforehand
 */
void SetupPackerOpts(const t_options& Options,
                     t_packer_opts* PackerOpts) {
    PackerOpts->output_file = Options.NetFile;

    PackerOpts->circuit_file_name = Options.CircuitFile;

    if (Options.do_packing) {
        PackerOpts->doPacking = STAGE_DO;
    }

    //TODO: document?
    PackerOpts->global_clocks = true;       /* DEFAULT */
    PackerOpts->hill_climbing_flag = false; /* DEFAULT */

    PackerOpts->allow_unrelated_clustering = Options.allow_unrelated_clustering;
    PackerOpts->connection_driven = Options.connection_driven_clustering;
    PackerOpts->timing_driven = Options.timing_driven_clustering;
    PackerOpts->cluster_seed_type = Options.cluster_seed_type;
    PackerOpts->alpha = Options.alpha_clustering;
    PackerOpts->beta = Options.beta_clustering;
    PackerOpts->pack_verbosity = Options.pack_verbosity;
    PackerOpts->enable_pin_feasibility_filter = Options.enable_clustering_pin_feasibility_filter;
    PackerOpts->balance_block_type_utilization = Options.balance_block_type_utilization;
    PackerOpts->target_external_pin_util = Options.target_external_pin_util;
    PackerOpts->target_device_utilization = Options.target_device_utilization;
    PackerOpts->prioritize_transitive_connectivity = Options.pack_prioritize_transitive_connectivity;
    PackerOpts->high_fanout_threshold = Options.pack_high_fanout_threshold;
    PackerOpts->transitive_fanout_threshold = Options.pack_transitive_fanout_threshold;
    PackerOpts->feasible_block_array_size = Options.pack_feasible_block_array_size;
    PackerOpts->use_attraction_groups = Options.use_attraction_groups;

    //TODO: document?
    PackerOpts->inter_cluster_net_delay = 1.0; /* DEFAULT */
    PackerOpts->auto_compute_inter_cluster_net_delay = true;
    PackerOpts->packer_algorithm = PACK_GREEDY; /* DEFAULT */

    PackerOpts->device_layout = Options.device_layout;

    PackerOpts->timing_update_type = Options.timing_update_type;
}

static void SetupNetlistOpts(const t_options& Options, t_netlist_opts& NetlistOpts) {
    NetlistOpts.const_gen_inference = Options.const_gen_inference;
    NetlistOpts.absorb_buffer_luts = Options.absorb_buffer_luts;
    NetlistOpts.sweep_dangling_primary_ios = Options.sweep_dangling_primary_ios;
    NetlistOpts.sweep_dangling_nets = Options.sweep_dangling_nets;
    NetlistOpts.sweep_dangling_blocks = Options.sweep_dangling_blocks;
    NetlistOpts.sweep_constant_primary_outputs = Options.sweep_constant_primary_outputs;
    NetlistOpts.netlist_verbosity = Options.netlist_verbosity;
}

/**
 * @brief Sets up the s_placer_opts structure based on users input.
 *
 * Error checking, such as checking for conflicting params
 * is assumed to be done beforehand
 */
static void SetupPlacerOpts(const t_options& Options, t_placer_opts* PlacerOpts) {
    if (Options.do_placement) {
        PlacerOpts->doPlacement = STAGE_DO;
    }

    PlacerOpts->inner_loop_recompute_divider = Options.inner_loop_recompute_divider;
    PlacerOpts->quench_recompute_divider = Options.quench_recompute_divider;

    //TODO: document?
    PlacerOpts->place_cost_exp = 1;

    PlacerOpts->td_place_exp_first = Options.place_exp_first;

    PlacerOpts->td_place_exp_last = Options.place_exp_last;

    PlacerOpts->place_algorithm = Options.PlaceAlgorithm;
    PlacerOpts->place_quench_algorithm = Options.PlaceQuenchAlgorithm;

    PlacerOpts->constraints_file = Options.constraints_file;

    PlacerOpts->pad_loc_type = Options.pad_loc_type;

    PlacerOpts->place_chan_width = Options.PlaceChanWidth;

    PlacerOpts->recompute_crit_iter = Options.RecomputeCritIter;

    PlacerOpts->timing_tradeoff = Options.PlaceTimingTradeoff;

    /* Depends on PlacerOpts->place_algorithm */
    PlacerOpts->delay_offset = Options.place_delay_offset;
    PlacerOpts->delay_ramp_delta_threshold = Options.place_delay_ramp_delta_threshold;
    PlacerOpts->delay_ramp_slope = Options.place_delay_ramp_slope;
    PlacerOpts->tsu_rel_margin = Options.place_tsu_rel_margin;
    PlacerOpts->tsu_abs_margin = Options.place_tsu_abs_margin;
    PlacerOpts->delay_model_type = Options.place_delay_model;
    PlacerOpts->delay_model_reducer = Options.place_delay_model_reducer;

    //TODO: document?
    PlacerOpts->place_freq = PLACE_ONCE; /* DEFAULT */

    PlacerOpts->post_place_timing_report_file = Options.post_place_timing_report_file;

    PlacerOpts->rlim_escape_fraction = Options.place_rlim_escape_fraction;
    PlacerOpts->move_stats_file = Options.place_move_stats_file;
    PlacerOpts->placement_saves_per_temperature = Options.placement_saves_per_temperature;
    PlacerOpts->place_delta_delay_matrix_calculation_method = Options.place_delta_delay_matrix_calculation_method;

    PlacerOpts->strict_checks = Options.strict_checks;

    PlacerOpts->write_placement_delay_lookup = Options.write_placement_delay_lookup;
    PlacerOpts->read_placement_delay_lookup = Options.read_placement_delay_lookup;

    PlacerOpts->allowed_tiles_for_delay_model = Options.allowed_tiles_for_delay_model;

    PlacerOpts->effort_scaling = Options.place_effort_scaling;
    PlacerOpts->timing_update_type = Options.timing_update_type;
    PlacerOpts->enable_analytic_placer = Options.enable_analytic_placer;
    PlacerOpts->place_static_move_prob = Options.place_static_move_prob;
    PlacerOpts->place_static_notiming_move_prob = Options.place_static_notiming_move_prob;
    PlacerOpts->place_high_fanout_net = Options.place_high_fanout_net;
    PlacerOpts->RL_agent_placement = Options.RL_agent_placement;
    PlacerOpts->place_agent_multistate = Options.place_agent_multistate;
    PlacerOpts->place_checkpointing = Options.place_checkpointing;
    PlacerOpts->place_agent_epsilon = Options.place_agent_epsilon;
    PlacerOpts->place_agent_gamma = Options.place_agent_gamma;
    PlacerOpts->place_dm_rlim = Options.place_dm_rlim;
    PlacerOpts->place_reward_fun = Options.place_reward_fun;
    PlacerOpts->place_crit_limit = Options.place_crit_limit;
    PlacerOpts->place_agent_algorithm = Options.place_agent_algorithm;
    PlacerOpts->place_constraint_expand = Options.place_constraint_expand;
    PlacerOpts->place_constraint_subtile = Options.place_constraint_subtile;
    PlacerOpts->floorplan_num_horizontal_partitions = Options.floorplan_num_horizontal_partitions;
    PlacerOpts->floorplan_num_vertical_partitions = Options.floorplan_num_vertical_partitions;
}

static void SetupAnalysisOpts(const t_options& Options, t_analysis_opts& analysis_opts) {
    if (Options.do_analysis) {
        analysis_opts.doAnalysis = STAGE_DO;
    }

    analysis_opts.gen_post_synthesis_netlist = Options.Generate_Post_Synthesis_Netlist;
    analysis_opts.gen_post_implementation_merged_netlist = Options.Generate_Post_Implementation_Merged_Netlist;

    analysis_opts.timing_report_npaths = Options.timing_report_npaths;
    analysis_opts.timing_report_detail = Options.timing_report_detail;
    analysis_opts.timing_report_skew = Options.timing_report_skew;
    analysis_opts.echo_dot_timing_graph_node = Options.echo_dot_timing_graph_node;

    analysis_opts.post_synth_netlist_unconn_input_handling = Options.post_synth_netlist_unconn_input_handling;
    analysis_opts.post_synth_netlist_unconn_output_handling = Options.post_synth_netlist_unconn_output_handling;

    analysis_opts.timing_update_type = Options.timing_update_type;
    analysis_opts.write_timing_summary = Options.write_timing_summary;
}

static void SetupPowerOpts(const t_options& Options, t_power_opts* power_opts, t_arch* Arch) {
    auto& device_ctx = g_vpr_ctx.mutable_device();

    power_opts->do_power = Options.do_power;

    if (power_opts->do_power) {
        if (!Arch->power)
            Arch->power = new t_power_arch();

        if (!Arch->clocks)
            Arch->clocks = new t_clock_arch();

        device_ctx.clock_arch = Arch->clocks;
    } else {
        Arch->power = nullptr;
        Arch->clocks = nullptr;
        device_ctx.clock_arch = nullptr;
    }
}

/*
 * Go through all the NoC options supplied by the user and store them internally.
 */
static void SetupNocOpts(const t_options& Options, t_noc_opts* NocOpts) {
    // assign the noc specific options from the command line
    NocOpts->noc = Options.noc;
    NocOpts->noc_flows_file = Options.noc_flows_file;
    NocOpts->noc_routing_algorithm = Options.noc_routing_algorithm;

    return;
}

static int find_ipin_cblock_switch_index(const t_arch& Arch) {
    int ipin_cblock_switch_index = UNDEFINED;
    for (int i = 0; i < Arch.num_switches; ++i) {
        if (Arch.Switches[i].name == Arch.ipin_cblock_switch_name) {
            if (ipin_cblock_switch_index != UNDEFINED) {
                VPR_FATAL_ERROR(VPR_ERROR_ARCH, "Found duplicate switches named '%s'\n", Arch.ipin_cblock_switch_name.c_str());
            } else {
                ipin_cblock_switch_index = i;
            }
        }
    }

    if (ipin_cblock_switch_index == UNDEFINED) {
        VPR_FATAL_ERROR(VPR_ERROR_ARCH, "Failed to find connection block input pin switch named '%s'\n", Arch.ipin_cblock_switch_name.c_str());
    }
    return ipin_cblock_switch_index;
}

static void alloc_and_load_intra_cluster_resources() {
    auto& device_ctx = g_vpr_ctx.mutable_device();

    for (auto& physical_type : device_ctx.physical_tile_types) {
        int physical_pin_offset = physical_type.num_pins;
        int physical_class_offset = (int)physical_type.class_inf.size();
        for (auto& sub_tile : physical_type.sub_tiles) {
            sub_tile.starting_internal_class_idx.resize(sub_tile.capacity.total());
            sub_tile.starting_internal_pin_idx.resize(sub_tile.capacity.total());
            for (int sub_tile_inst = 0; sub_tile_inst < sub_tile.capacity.total(); sub_tile_inst++) {
                for (auto logic_block_ptr : sub_tile.equivalent_sites) {
                    sub_tile.starting_internal_class_idx[sub_tile_inst].insert(
                        std::make_pair(logic_block_ptr, physical_class_offset));
                    sub_tile.starting_internal_pin_idx[sub_tile_inst].insert(
                        std::make_pair(logic_block_ptr, physical_pin_offset));

                    auto logical_classes = logic_block_ptr->logical_class_inf;
                    std::for_each(logical_classes.begin(), logical_classes.end(),
                                  [&physical_pin_offset](t_class& l_class) { for(auto &pin : l_class.pinlist) {
                        pin += physical_pin_offset;} });

                    int physical_class_num = physical_class_offset;
                    for (auto& logic_class : logical_classes) {
                        auto result = physical_type.internal_class_inf.insert(std::make_pair(physical_class_num, logic_class));
                        VTR_ASSERT(result.second);
                        physical_class_num++;
                    }

                    vtr::bimap<t_logical_pin, t_physical_pin> directs_map;
                    for (auto pin_to_pb_pin_map : logic_block_ptr->pin_logical_num_to_pb_pin_mapping) {
                        int physical_pin_num = pin_to_pb_pin_map.first + physical_pin_offset;
                        int class_logical_num = logic_block_ptr->pb_pin_to_class_logical_num_mapping.at(pin_to_pb_pin_map.second);

                        auto result = physical_type.internal_pin_class.insert(
                            std::make_pair(physical_pin_num, class_logical_num + physical_class_offset));
                        VTR_ASSERT(result.second);
                    }

                    physical_pin_offset += (int)logic_block_ptr->pin_logical_num_to_pb_pin_mapping.size();
                    physical_class_offset += (int)logic_block_ptr->logical_class_inf.size();
                }
            }
        }
    }
}
