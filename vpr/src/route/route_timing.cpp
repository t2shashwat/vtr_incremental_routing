#include <cstdio>
#include <ctime>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include "vtr_assert.h"
#include "vtr_log.h"
#include "vtr_time.h"

#include "vpr_utils.h"
#include "vpr_types.h"
#include "vpr_error.h"

#include "globals.h"
#include "route_export.h"
#include "route_common.h"
#include "route_tree_timing.h"
#include "route_timing.h"
#include "net_delay.h"
#include "stats.h"
#include "echo_files.h"
#include "draw.h"
#include "breakpoint.h"
#include "move_utils.h"
#include "rr_graph.h"
#include "routing_predictor.h"
#include "VprTimingGraphResolver.h"

// all functions in profiling:: namespace, which are only activated if PROFILE is defined
#include "route_profiling.h"

#include "timing_info.h"
#include "timing_util.h"
#include "route_budgets.h"
#include "binary_heap.h"
#include "bucket.h"
#include "connection_router.h"

#include "tatum/TimingReporter.hpp"
#include "overuse_report.h"

//shashwat's includes=============
#include <fstream>
#include<string>
#include "read_route.h"
#include "check_route.h"
#include "route_breadth_first.h"

//=================================
#define CONGESTED_SLOPE_VAL -0.04

enum class RouterCongestionMode {
    NORMAL,
    CONFLICTED
};

//identifies the two breakpoint types in routing
typedef enum router_breakpoint_type {
    BP_ROUTE_ITER,
    BP_NET_ID
} bp_router_type;

struct RoutingMetrics {
    size_t used_wirelength = 0;

    float sWNS = std::numeric_limits<float>::quiet_NaN();
    float sTNS = std::numeric_limits<float>::quiet_NaN();
    float hWNS = std::numeric_limits<float>::quiet_NaN();
    float hTNS = std::numeric_limits<float>::quiet_NaN();
    tatum::TimingPathInfo critical_path;
};

/*
 * File-scope variables
 */

/**
 * @brief Run-time flag to control when router debug information is printed
 * Note only enables debug output if compiled with VTR_ENABLE_DEBUG_LOGGING defined
 * f_router_debug is used to stop the router when a breakpoint is reached. When a breakpoint is reached, this flag is set to true.
 *
 * In addition f_router_debug is used to print additional debug information during routing, for instance lookahead expected costs
 * information.
 */
bool f_router_debug = false;

//Count the number of times the router has failed
static int num_routing_failed = 0;

/******************** Subroutines local to route_timing.c ********************/

template<typename ConnectionRouter>
static bool timing_driven_route_sink(
    ConnectionRouter& router,
    ClusterNetId net_id,
    unsigned itarget,
    int target_pin,
    const t_conn_cost_params cost_params,
    const t_router_opts& router_opts,
    t_rt_node* rt_root,
    t_rt_node** rt_node_of_sink,
    SpatialRouteTreeLookup& spatial_rt_lookup,
    RouterStats& router_stats,
    route_budgets& budgeting_inf,
    const RoutingPredictor& routing_predictor,
    bool is_flat,
    std::set<int> branch_nodes,
    int itry);

template<typename ConnectionRouter>
static bool timing_driven_pre_route_to_clock_root(
    ConnectionRouter& router,
    ClusterNetId net_id,
    int sink_node,
    const t_conn_cost_params cost_params,
    int high_fanout_threshold,
    t_rt_node* rt_root,
    SpatialRouteTreeLookup& spatial_rt_lookup,
    RouterStats& router_stats,
    bool is_flat);

void disable_expansion_and_remove_sink_from_route_tree_nodes(t_rt_node* node);

static t_rt_node* setup_routing_resources(int itry,
                                          ClusterNetId net_id,
                                          unsigned num_sinks,
                                          int min_incremental_reroute_fanout,
                                          CBRR& incremental_rerouting_res,
                                          t_rt_node** rt_node_of_sink,
                                          const t_router_opts& router_opts,
                                          bool ripup_high_fanout_nets);
//mycode =========================
static t_rt_node* setup_routing_resources_incr_route(const t_file_name_opts& filename_opts,
                                          int itry,
                                          ClusterNetId net_id,
                                          unsigned num_sinks,
                                          int min_incremental_reroute_fanout,
                                          CBRR& incremental_rerouting_res,
                                          t_rt_node** rt_node_of_sink,
                                          const t_router_opts& router_opts,
                                          bool ripup_high_fanout_nets);

//=================================
static bool timing_driven_check_net_delays(ClbNetPinsMatrix<float>& net_delay);

void increase_short_path_crit_if_congested(std::vector<ClusterNetId>& rerouted_nets,
                                           route_budgets& budgeting_inf,
                                           int itry);

static bool should_route_net(ClusterNetId net_id, CBRR& connections_inf, bool if_force_reroute);
static bool early_exit_heuristic(const t_router_opts& router_opts, const WirelengthInfo& wirelength_info);

static bool check_hold(const t_router_opts& router_opts, float worst_neg_slack);

struct more_sinks_than {
    inline bool operator()(const ClusterNetId net_index1, const ClusterNetId net_index2) {
        auto& cluster_ctx = g_vpr_ctx.clustering();
        return cluster_ctx.clb_nlist.net_sinks(net_index1).size() > cluster_ctx.clb_nlist.net_sinks(net_index2).size();
    }
};

struct less_sinks_than {
    inline bool operator()(const ClusterNetId net_index1, const ClusterNetId net_index2) {
        auto& cluster_ctx = g_vpr_ctx.clustering();
        return cluster_ctx.clb_nlist.net_sinks(net_index1).size() < cluster_ctx.clb_nlist.net_sinks(net_index2).size();
    }
};

static size_t calculate_wirelength_available();
static WirelengthInfo calculate_wirelength_info(size_t available_wirelength);

static void print_route_status_header();
static void print_route_status(int itry,
                               double elapsed_sec,
                               float pres_fac,
                               int num_bb_updated,
                               const RouterStats& router_stats,
                               const OveruseInfo& overuse_info,
                               const WirelengthInfo& wirelength_info,
                               std::shared_ptr<const SetupHoldTimingInfo> timing_info,
                               float est_success_iteration);

static void print_overused_nodes_status(const t_router_opts& router_opts, const OveruseInfo& overuse_info);

static void print_router_criticality_histogram(const SetupTimingInfo& timing_info,
                                               const ClusteredPinAtomPinsLookup& netlist_pin_lookup);

static bool is_high_fanout(int fanout, int fanout_threshold);

static size_t dynamic_update_bounding_boxes(const std::vector<ClusterNetId>& nets, int high_fanout_threshold);
static t_bb calc_current_bb(const t_trace* head);

static bool is_better_quality_routing(const vtr::vector<ClusterNetId, t_traceback>& best_routing,
                                      const RoutingMetrics& best_routing_metrics,
                                      const WirelengthInfo& wirelength_info,
                                      std::shared_ptr<const SetupHoldTimingInfo> timing_info);

static bool early_reconvergence_exit_heuristic(const t_router_opts& router_opts,
                                               int itry_since_last_convergence,
                                               std::shared_ptr<const SetupHoldTimingInfo> timing_info,
                                               const RoutingMetrics& best_routing_metrics);

static void generate_route_timing_reports(const t_router_opts& router_opts,
                                          const t_analysis_opts& analysis_opts,
                                          const SetupTimingInfo& timing_info,
                                          const RoutingDelayCalculator& delay_calc);

static void prune_unused_non_configurable_nets(CBRR& connections_inf);

static void init_net_delay_from_lookahead(const RouterLookahead& router_lookahead,
                                          ClbNetPinsMatrix<float>& net_delay, const float bias);

#ifndef NO_GRAPHICS
void update_router_info_and_check_bp(bp_router_type type, int net_id);
#endif

// The reason that try_timing_driven_route_tmpl (and descendents) are being
// templated over is because using a virtual interface instead fully templating
// the router results in a 5% runtime increase.
//
// The reason to template over the router in general is to enable runtime
// selection of core router algorithm's, specifically the router heap.
template<typename ConnectionRouter>
static bool try_timing_driven_route_tmpl(const t_router_opts& router_opts,
                                         const t_analysis_opts& analysis_opts,
                                         const std::vector<t_segment_inf>& segment_inf,
                                         ClbNetPinsMatrix<float>& net_delay,
                                         const ClusteredPinAtomPinsLookup& netlist_pin_lookup,
                                         std::shared_ptr<SetupHoldTimingInfo> timing_info,
                                         std::shared_ptr<RoutingDelayCalculator> delay_calc,
                                         ScreenUpdatePriority first_iteration_priority,
                                         bool is_flat);
//mycode ====================================
template<typename ConnectionRouter>
static bool try_timing_driven_route_tmpl_incr_route(const t_file_name_opts& filename_opts,
                                         const t_router_opts& router_opts,
                                         const t_analysis_opts& analysis_opts,
                                         const std::vector<t_segment_inf>& segment_inf,
                                         ClbNetPinsMatrix<float>& net_delay,
                                         const ClusteredPinAtomPinsLookup& netlist_pin_lookup,
                                         std::shared_ptr<SetupHoldTimingInfo> timing_info,
                                         std::shared_ptr<RoutingDelayCalculator> delay_calc,
                                         ScreenUpdatePriority first_iteration_priority,
                                         bool is_flat);
//===========================================
/************************ Subroutine definitions *****************************/
bool try_timing_driven_route(const t_router_opts& router_opts,
                             const t_analysis_opts& analysis_opts,
                             const std::vector<t_segment_inf>& segment_inf,
                             ClbNetPinsMatrix<float>& net_delay,
                             const ClusteredPinAtomPinsLookup& netlist_pin_lookup,
                             std::shared_ptr<SetupHoldTimingInfo> timing_info,
                             std::shared_ptr<RoutingDelayCalculator> delay_calc,
                             ScreenUpdatePriority first_iteration_priority,
                             bool is_flat) {
    switch (router_opts.router_heap) {
        case e_heap_type::BINARY_HEAP:
            return try_timing_driven_route_tmpl<ConnectionRouter<BinaryHeap>>(
                router_opts,
                analysis_opts,
                segment_inf,
                net_delay,
                netlist_pin_lookup,
                timing_info,
                delay_calc,
                first_iteration_priority,
                is_flat);
            break;
        case e_heap_type::BUCKET_HEAP_APPROXIMATION:
            return try_timing_driven_route_tmpl<ConnectionRouter<Bucket>>(
                router_opts,
                analysis_opts,
                segment_inf,
                net_delay,
                netlist_pin_lookup,
                timing_info,
                delay_calc,
                first_iteration_priority,
                is_flat);
        default:
            VPR_FATAL_ERROR(VPR_ERROR_ROUTE, "Unknown heap type %d", router_opts.router_heap);
    }
}

template<typename ConnectionRouter>
bool try_timing_driven_route_tmpl(const t_router_opts& router_opts,
                                  const t_analysis_opts& analysis_opts,
                                  const std::vector<t_segment_inf>& segment_inf,
                                  ClbNetPinsMatrix<float>& net_delay,
                                  const ClusteredPinAtomPinsLookup& netlist_pin_lookup,
                                  std::shared_ptr<SetupHoldTimingInfo> timing_info,
                                  std::shared_ptr<RoutingDelayCalculator> delay_calc,
                                  ScreenUpdatePriority first_iteration_priority,
                                  bool is_flat) {
    /* Timing-driven routing algorithm.  The timing graph (includes slack)   *
     * must have already been allocated, and net_delay must have been allocated. *
     * Returns true if the routing succeeds, false otherwise.                    */

    const auto& device_ctx = g_vpr_ctx.device();
    const auto& atom_ctx = g_vpr_ctx.atom();
    const auto& cluster_ctx = g_vpr_ctx.clustering();
    auto& route_ctx = g_vpr_ctx.mutable_routing();

    //Initially, the router runs normally trying to reduce congestion while
    //balancing other metrics (timing, wirelength, run-time etc.)
    RouterCongestionMode router_congestion_mode = RouterCongestionMode::NORMAL;

    //Initialize and properly size the lookups for profiling
    profiling::profiling_initialization(get_max_pins_per_net());

    //sort so net with most sinks is routed first.
    auto sorted_nets = std::vector<ClusterNetId>(cluster_ctx.clb_nlist.nets().begin(), cluster_ctx.clb_nlist.nets().end());
    // DEFAULT
    std::sort(sorted_nets.begin(), sorted_nets.end(), more_sinks_than());
    /*
     * Configure the routing predictor
     */
    RoutingPredictor routing_predictor;
    float abort_iteration_threshold = std::numeric_limits<float>::infinity(); //Default no early abort
    if (router_opts.routing_failure_predictor == SAFE) {
        abort_iteration_threshold = ROUTING_PREDICTOR_ITERATION_ABORT_FACTOR_SAFE * router_opts.max_router_iterations;
    } else if (router_opts.routing_failure_predictor == AGGRESSIVE) {
        abort_iteration_threshold = ROUTING_PREDICTOR_ITERATION_ABORT_FACTOR_AGGRESSIVE * router_opts.max_router_iterations;
    } else {
        VTR_ASSERT_MSG(router_opts.routing_failure_predictor == OFF, "Unrecognized routing failure predictor setting");
    }

    float high_effort_congestion_mode_iteration_threshold = router_opts.congested_routing_iteration_threshold_frac * router_opts.max_router_iterations;

    /* Set delay of ignored signals to zero. Non-ignored net delays are set by
     * update_net_delays_from_route_tree() inside timing_driven_route_net(),
     * which is only called for non-ignored nets. */
    for (auto net_id : cluster_ctx.clb_nlist.nets()) {
        if (cluster_ctx.clb_nlist.net_is_ignored(net_id)) {
            for (unsigned int ipin = 1; ipin < cluster_ctx.clb_nlist.net_pins(net_id).size(); ++ipin) {
                net_delay[net_id][ipin] = 0.;
            }
        }
    }

    CBRR connections_inf{};

    route_budgets budgeting_inf;

    const auto* router_lookahead = get_cached_router_lookahead(
        router_opts.lookahead_type,
        router_opts.sbNode_lookahead_factor,
        router_opts.write_router_lookahead,
        router_opts.read_router_lookahead,
        segment_inf,
        is_flat);

    /*
     * Routing parameters
     */
    float pres_fac = update_pres_fac(router_opts.first_iter_pres_fac); /* Typically 0 -> ignore cong. */
    int bb_fac = router_opts.bb_factor;

    //When routing conflicts are detected the bounding boxes are scaled
    //by BB_SCALE_FACTOR every BB_SCALE_ITER_COUNT iterations
    constexpr float BB_SCALE_FACTOR = 2;
    constexpr int BB_SCALE_ITER_COUNT = 5;

    size_t available_wirelength = calculate_wirelength_available();

    /*
     * Routing status and metrics
     */
    bool routing_is_successful = false;
    WirelengthInfo wirelength_info;
    OveruseInfo overuse_info(device_ctx.rr_graph.num_nodes());
    tatum::TimingPathInfo critical_path;
    int itry; //Routing iteration number
    int itry_conflicted_mode = 0;

    /*
     * Best result so far
     */
    vtr::vector<ClusterNetId, t_traceback> best_routing;
    t_clb_opins_used best_clb_opins_used_locally;
    RoutingMetrics best_routing_metrics;
    int legal_convergence_count = 0;
    std::vector<int> scratch;

    ConnectionRouter router(
        device_ctx.grid,
        *router_lookahead,
        device_ctx.rr_graph.rr_nodes(),
        &device_ctx.rr_graph,
        device_ctx.rr_rc_data,
        device_ctx.rr_graph.rr_switch(),
        route_ctx.rr_node_route_inf,
        is_flat);

    // Make sure template type ConnectionRouter is a ConnectionRouterInterface.
    static_assert(std::is_base_of<ConnectionRouterInterface, ConnectionRouter>::value, "ConnectionRouter must implement the ConnectionRouterInterface");

    /*
     * On the first routing iteration ignore congestion to get reasonable net
     * delay estimates. Set criticalities to 1 when timing analysis is on to
     * optimize timing, and to 0 when timing analysis is off to optimize routability.
     *
     * Subsequent iterations use the net delays from the previous iteration.
     */
    std::shared_ptr<SetupHoldTimingInfo> route_timing_info;
    {
        vtr::ScopedStartFinishTimer init_timing_timer("Initializing router criticalities");
        if (timing_info) {
            if (router_opts.initial_timing == e_router_initial_timing::ALL_CRITICAL) {
                //First routing iteration, make all nets critical for a min-delay routing
                route_timing_info = make_constant_timing_info(1.);
            } else {
                VTR_ASSERT(router_opts.initial_timing == e_router_initial_timing::LOOKAHEAD);

                {
                    //Estimate initial connection delays from the router lookahead
                    init_net_delay_from_lookahead(*router_lookahead, net_delay, router_opts.sbNode_lookahead_factor);

                    //Run STA to get estimated criticalities
                    timing_info->update();
                }
                route_timing_info = timing_info;
            }
        } else {
            //Not timing driven, force criticality to zero for a routability-driven routing
            route_timing_info = make_constant_timing_info(0.);
        }
        VTR_LOG("Initial Net Connection Criticality Histogram:\n");
        print_router_criticality_histogram(*route_timing_info, netlist_pin_lookup);
    }

    std::unique_ptr<ClusteredPinTimingInvalidator> pin_timing_invalidator;
    if (timing_info) {
        pin_timing_invalidator = std::make_unique<ClusteredPinTimingInvalidator>(cluster_ctx.clb_nlist,
                                                                                 netlist_pin_lookup,
                                                                                 atom_ctx.nlist,
                                                                                 atom_ctx.lookup,
                                                                                 *timing_info->timing_graph());
    }

    RouterStats router_stats;
    timing_driven_route_structs route_structs;
    float prev_iter_cumm_time = 0;
    vtr::Timer iteration_timer;
    int num_net_bounding_boxes_updated = 0;
    int itry_since_last_convergence = -1;

    // This heap is used for reserve_locally_used_opins.
    BinaryHeap small_heap;
    small_heap.init_heap(device_ctx.grid);

    // When RCV is enabled the router will not stop unless negative hold slack is 0
    // In some cases this isn't doable, due to global nets or intracluster routing issues
    // In these cases RCV will finish early if it goes RCV_FINISH_EARLY_COUNTDOWN iterations without detecting resolvable negative hold slack
    // Increasing this will make the router fail occasionally, decreasing will sometimes not let all hold violations be resolved
    constexpr int RCV_FINISH_EARLY_COUNTDOWN = 15;

    int rcv_finished_count = RCV_FINISH_EARLY_COUNTDOWN;

    print_route_status_header();
    //For each node, load nets allowed to use that node
    for (itry = 1; itry <= router_opts.max_router_iterations; ++itry) {
        RouterStats router_iteration_stats;
        std::vector<ClusterNetId> rerouted_nets;

        /* Reset "is_routed" and "is_fixed" flags to indicate nets not pre-routed (yet) */
        for (auto net_id : cluster_ctx.clb_nlist.nets()) {
            route_ctx.net_status.set_is_routed(net_id, false);
            route_ctx.net_status.set_is_fixed(net_id, false);
        }

        if (itry_since_last_convergence >= 0) {
            ++itry_since_last_convergence;
        }

        // Calculate this once and pass it into net routing to check if should reroute for hold
        float worst_negative_slack = 0;
        if (budgeting_inf.if_set()) {
            worst_negative_slack = timing_info->hold_total_negative_slack();
        }

        /*
         * Route each net
         */
	
        for (auto net_id : sorted_nets) {
            bool was_rerouted = false;
            bool is_routable = try_timing_driven_route_net(router,
                                                           net_id,
                                                           itry,
                                                           pres_fac,
                                                           router_opts,
                                                           connections_inf,
                                                           router_iteration_stats,
                                                           route_structs.pin_criticality,
                                                           route_structs.rt_node_of_sink,
                                                           net_delay,
                                                           netlist_pin_lookup,
                                                           route_timing_info,
                                                           pin_timing_invalidator.get(),
                                                           budgeting_inf,
                                                           was_rerouted,
                                                           worst_negative_slack,
                                                           routing_predictor,
                                                           is_flat);

            if (!is_routable) {
                return (false); //Impossible to route
            }

            if (was_rerouted) {
                rerouted_nets.push_back(net_id);
#ifndef NO_GRAPHICS
                update_router_info_and_check_bp(BP_NET_ID, size_t(net_id));
#endif
            }
        }

        // Make sure any CLB OPINs used up by subblocks being hooked directly to them are reserved for that purpose
        bool rip_up_local_opins = (itry == 1 ? false : true);
        reserve_locally_used_opins(&small_heap, pres_fac,
                                   router_opts.acc_fac, router_opts.global_occ_factor, rip_up_local_opins);

        /*
         * Calculate metrics for the current routing
         */
        bool routing_is_feasible = feasible_routing();
        float est_success_iteration = routing_predictor.estimate_success_iteration();

        //Update resource costs and overuse info
        if (itry == 1) {
            pathfinder_update_acc_cost_and_overuse_info(0., overuse_info); /* Acc_fac=0 for first iter. */
        } else {
            pathfinder_update_acc_cost_and_overuse_info(router_opts.acc_fac, overuse_info);
        }

        wirelength_info = calculate_wirelength_info(available_wirelength);
        routing_predictor.add_iteration_overuse(itry, overuse_info.overused_nodes);

        if (timing_info) {
            //Update timing based on the new routing
            //Note that the net delays have already been updated by timing_driven_route_net
            timing_info->update();
            timing_info->set_warn_unconstrained(false); //Don't warn again about unconstrained nodes again during routing
            pin_timing_invalidator->reset();

            //Use the real timing analysis criticalities for subsequent routing iterations
            //  'route_timing_info' is what is actually passed into the net/connection routers,
            //  and for the 1st iteration may not be the actual STA results (e.g. all criticalities set to 1)
            route_timing_info = timing_info;

            critical_path = timing_info->least_slack_critical_path();

            VTR_ASSERT_SAFE(timing_driven_check_net_delays(net_delay));

            if (itry == 1) {
                generate_route_timing_reports(router_opts, analysis_opts, *timing_info, *delay_calc);
            }
        }

        float iter_cumm_time = iteration_timer.elapsed_sec();
        float iter_elapsed_time = iter_cumm_time - prev_iter_cumm_time;

        //Output progress
        print_route_status(itry, iter_elapsed_time, pres_fac, num_net_bounding_boxes_updated, router_iteration_stats, overuse_info, wirelength_info, timing_info, est_success_iteration);

        prev_iter_cumm_time = iter_cumm_time;

        //Update graphics
        if (itry == 1) {
            update_screen(first_iteration_priority, "Routing...", ROUTING, timing_info);
        } else {
            update_screen(ScreenUpdatePriority::MINOR, "Routing...", ROUTING, timing_info);
        }

        if (router_opts.save_routing_per_iteration) {
            std::string filename = vtr::string_fmt("iteration_%03d.route", itry);
            print_route(nullptr, filename.c_str());
        }

        //Update router stats (total)
        router_stats.connections_routed += router_iteration_stats.connections_routed;
        router_stats.nets_routed += router_iteration_stats.nets_routed;
        router_stats.heap_pushes += router_iteration_stats.heap_pushes;
        router_stats.heap_pops += router_iteration_stats.heap_pops;

        /*
         * Are we finished?
         */
        if (is_iteration_complete(routing_is_feasible, router_opts, itry, timing_info, rcv_finished_count == 0)) {
            auto& router_ctx = g_vpr_ctx.routing();

            if (is_better_quality_routing(best_routing, best_routing_metrics, wirelength_info, timing_info)) {
                //Save routing
                best_routing = router_ctx.trace;
                best_clb_opins_used_locally = router_ctx.clb_opins_used_locally;

                routing_is_successful = true;

                //Update best metrics
                if (timing_info) {
                    timing_driven_check_net_delays(net_delay);

                    best_routing_metrics.sTNS = timing_info->setup_total_negative_slack();
                    best_routing_metrics.sWNS = timing_info->setup_worst_negative_slack();
                    best_routing_metrics.hTNS = timing_info->hold_total_negative_slack();
                    best_routing_metrics.hWNS = timing_info->hold_worst_negative_slack();
                    best_routing_metrics.critical_path = critical_path;
                }
                best_routing_metrics.used_wirelength = wirelength_info.used_wirelength();
            }

            //Decrease pres_fac so that critical connections will take more direct routes
            //Note that we use first_iter_pres_fac here (typically zero), and switch to
            //use initial_pres_fac on the next iteration.
            pres_fac = update_pres_fac(router_opts.first_iter_pres_fac);

            //Reduce timing tolerances to re-route more delay-suboptimal signals
            connections_inf.set_connection_criticality_tolerance(0.7);
            connections_inf.set_connection_delay_tolerance(1.01);

            ++legal_convergence_count;
            itry_since_last_convergence = 0;

            VTR_ASSERT(routing_is_successful);
        }

        if (itry_since_last_convergence == 1) {
            //We used first_iter_pres_fac when we started routing again
            //after the first routing convergence. Since that is often zero,
            //we want to set pres_fac to a reasonable (i.e. typically non-zero)
            //value afterwards -- so it grows when multiplied by pres_fac_mult
            pres_fac = update_pres_fac(router_opts.initial_pres_fac);
        }

        //Have we converged the maximum number of times, did not make any changes, or does it seem
        //unlikely additional convergences will improve QoR?
        if (legal_convergence_count >= router_opts.max_convergence_count
            || router_iteration_stats.connections_routed == 0
            || early_reconvergence_exit_heuristic(router_opts, itry_since_last_convergence, timing_info, best_routing_metrics)) {
#ifndef NO_GRAPHICS
            update_router_info_and_check_bp(BP_ROUTE_ITER, -1);
#endif
            break; //Done routing
        }

        /*
         * Abort checks: Should we give-up because this routing problem is unlikely to converge to a legal routing?
         */
        if (itry == 1 && early_exit_heuristic(router_opts, wirelength_info)) {
#ifndef NO_GRAPHICS
            update_router_info_and_check_bp(BP_ROUTE_ITER, -1);
#endif
            //Abort
            break;
        }

        //Estimate at what iteration we will converge to a legal routing
        if (overuse_info.overused_nodes > ROUTING_PREDICTOR_MIN_ABSOLUTE_OVERUSE_THRESHOLD) {
            //Only consider aborting if we have a significant number of overused resources

            if (!std::isnan(est_success_iteration) && est_success_iteration > abort_iteration_threshold && router_opts.routing_budgets_algorithm != YOYO) {
                VTR_LOG("Routing aborted, the predicted iteration for a successful route (%.1f) is too high.\n", est_success_iteration);
#ifndef NO_GRAPHICS
                update_router_info_and_check_bp(BP_ROUTE_ITER, -1);
#endif
                break; //Abort
            }
        }

        if (itry == 1 && router_opts.exit_after_first_routing_iteration) {
            VTR_LOG("Exiting after first routing iteration as requested\n");
#ifndef NO_GRAPHICS
            update_router_info_and_check_bp(BP_ROUTE_ITER, -1);
#endif
            break;
        }

        /*
         * Prepare for the next iteration
         */

        if (router_opts.route_bb_update == e_route_bb_update::DYNAMIC) {
            num_net_bounding_boxes_updated = dynamic_update_bounding_boxes(rerouted_nets, router_opts.high_fanout_threshold);
        }

        if (itry >= high_effort_congestion_mode_iteration_threshold) {
            //We are approaching the maximum number of routing iterations,
            //and still do not have a legal routing. Switch to a mode which
            //focuses more on attempting to resolve routing conflicts.
            router_congestion_mode = RouterCongestionMode::CONFLICTED;
        }

        //Update pres_fac
        if (itry == 1) {
            pres_fac = update_pres_fac(router_opts.initial_pres_fac);
        } else {
            pres_fac *= router_opts.pres_fac_mult;

            /* Avoid overflow for high iteration counts, even if acc_cost is big */
            pres_fac = update_pres_fac(std::min(pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5)));

            // Increase short path criticality if it's having a hard time resolving hold violations due to congestion
            if (budgeting_inf.if_set()) {
                bool rcv_finished = false;

                /* This constant represents how much extra delay the budget increaser adds to the minimum and maximum delay budgets
                 * Experimentally this value delivers fast hold slack resolution, while not overwhelming the router 
                 * Increasing this will make it resolve hold faster, but could result in lower circuit quality */
                constexpr float budget_increase_factor = 300e-12;

                if (itry > 5 && worst_negative_slack != 0) rcv_finished = budgeting_inf.increase_min_budgets_if_struggling(budget_increase_factor, timing_info, worst_negative_slack, netlist_pin_lookup);
                if (rcv_finished)
                    rcv_finished_count--;
                else
                    rcv_finished_count = RCV_FINISH_EARLY_COUNTDOWN;
            }
        }

        if (router_congestion_mode == RouterCongestionMode::CONFLICTED) {
            //The design appears to have routing conflicts which are difficult to resolve:
            //  1) Don't re-route legal connections due to delay. This allows
            //     the router to focus on the actual conflicts
            //  2) Increase the net bounding boxes. This potentially allows
            //     the router to route around otherwise congested regions
            //     (at the cost of high run-time).

            //Increase the size of the net bounding boxes to give the router more
            //freedom to find alternate paths.
            //
            //In the case of routing conflicts there are multiple connections competing
            //for the same resources which can not resolve the congestion themselves.
            //In normal routing mode we try to keep the bounding boxes small to minimize
            //run-time, but this can limits how far signals can detour (i.e. they can't
            //route outside the bounding box), which can cause conflicts to oscillate back
            //and forth without resolving.
            //
            //By scaling the bounding boxes here, we slowly increase the router's search
            //space in hopes of it allowing signals to move further out of the way to
            //alleviate the conflicts.
            if (itry_conflicted_mode % BB_SCALE_ITER_COUNT == 0) {
                //We scale the bounding boxes by BB_SCALE_FACTOR,
                //every BB_SCALE_ITER_COUNT iterations. This ensures
                //that we give the router some time (BB_SCALE_ITER_COUNT) to try
                //resolve/negotiate congestion at the new BB factor.
                //
                //Note that we increase the BB factor slowly to try and minimize
                //the bounding box size (since larger bounding boxes slow the router down).
                auto& grid = g_vpr_ctx.device().grid;
                int max_grid_dim = std::max(grid.width(), grid.height());

                //Scale by BB_SCALE_FACTOR but clip to grid size to avoid overflow
                bb_fac = std::min<int>(max_grid_dim, bb_fac * BB_SCALE_FACTOR);

                route_ctx.route_bb = load_route_bb(bb_fac);
            }

            ++itry_conflicted_mode;
        }

        if (timing_info) {
            if (should_setup_lower_bound_connection_delays(itry, router_opts)) {
                // first iteration sets up the lower bound connection delays since only timing is optimized for
                connections_inf.set_stable_critical_path_delay(critical_path.delay());
                connections_inf.set_lower_bound_connection_delays(net_delay);

                //load budgets using information from uncongested delay information
                budgeting_inf.load_route_budgets(net_delay, timing_info, netlist_pin_lookup, router_opts);
                /*for debugging purposes*/
                // if (budgeting_inf.if_set()) {
                //     budgeting_inf.print_route_budget(std::string("route_budgets_") + std::to_string(itry) + ".txt", net_delay);
                // }

                if (router_opts.routing_budgets_algorithm == YOYO) router.set_rcv_enabled(true);

            } else {
                bool stable_routing_configuration = true;

                /*
                 * Determine if any connection need to be forcibly re-routed due to timing
                 */

                //Yes, if explicitly enabled
                bool should_ripup_for_delay = (router_opts.incr_reroute_delay_ripup == e_incr_reroute_delay_ripup::ON);

                //Or, if things are not too congested
                should_ripup_for_delay |= (router_opts.incr_reroute_delay_ripup == e_incr_reroute_delay_ripup::AUTO
                                           && router_congestion_mode == RouterCongestionMode::NORMAL);

                if (should_ripup_for_delay) {
                    if (connections_inf.critical_path_delay_grew_significantly(critical_path.delay())) {
                        // only need to forcibly reroute if critical path grew significantly
                        stable_routing_configuration = connections_inf.forcibly_reroute_connections(router_opts.max_criticality,
                                                                                                    timing_info,
                                                                                                    netlist_pin_lookup,
                                                                                                    net_delay);
                    }
                }

                // not stable if any connection needs to be forcibly rerouted
                if (stable_routing_configuration) {
                    connections_inf.set_stable_critical_path_delay(critical_path.delay());
                }
            }
        } else {
            /* If timing analysis is not enabled, make sure that the criticalities and the
             * net_delays stay as 0 so that wirelength can be optimized. */

            for (auto net_id : cluster_ctx.clb_nlist.nets()) {
                for (unsigned int ipin = 1; ipin < cluster_ctx.clb_nlist.net_pins(net_id).size(); ++ipin) {
                    net_delay[net_id][ipin] = 0.;
                }
            }
        }

        if (router_opts.congestion_analysis) profiling::congestion_analysis();
        if (router_opts.fanout_analysis) profiling::time_on_fanout_analysis();
        // profiling::time_on_criticality_analysis();
    }

    if (routing_is_successful) {
        VTR_LOG("Restoring best routing\n");

        auto& router_ctx = g_vpr_ctx.mutable_routing();

        /* Restore congestion from best route */
        for (auto net_id : cluster_ctx.clb_nlist.nets()) {
            pathfinder_update_path_occupancy(route_ctx.trace[net_id].head, -1);
            pathfinder_update_path_occupancy(best_routing[net_id].head, 1);
        }
        router_ctx.trace = best_routing;
        router_ctx.clb_opins_used_locally = best_clb_opins_used_locally;

        prune_unused_non_configurable_nets(connections_inf);

        if (timing_info) {
            VTR_LOG("Critical path: %g ns\n", 1e9 * best_routing_metrics.critical_path.delay());
        }

        VTR_LOG("Successfully routed after %d routing iterations.\n", itry);
    } else {
        VTR_LOG("Routing failed.\n");

        //If the routing fails, print the overused info
        print_overused_nodes_status(router_opts, overuse_info);

        ++num_routing_failed;

#ifdef VTR_ENABLE_DEBUG_LOGGING
        if (f_router_debug) print_invalid_routing_info(is_flat);
#endif
    }

    VTR_LOG("Final Net Connection Criticality Histogram:\n");
    print_router_criticality_histogram(*route_timing_info, netlist_pin_lookup);

    VTR_LOG("Router Stats: total_nets_routed: %zu total_connections_routed: %zu total_heap_pushes: %zu total_heap_pops: %zu\n",
            router_stats.nets_routed, router_stats.connections_routed, router_stats.heap_pushes, router_stats.heap_pops);

    return routing_is_successful;
}

template<typename ConnectionRouter>
bool try_timing_driven_route_net(ConnectionRouter& router,
                                 ClusterNetId net_id,
                                 int itry,
                                 float pres_fac,
                                 const t_router_opts& router_opts,
                                 CBRR& connections_inf,
                                 RouterStats& router_stats,
                                 float* pin_criticality,
                                 t_rt_node** rt_node_of_sink,
                                 ClbNetPinsMatrix<float>& net_delay,
                                 const ClusteredPinAtomPinsLookup& netlist_pin_lookup,
                                 std::shared_ptr<SetupHoldTimingInfo> timing_info,
                                 ClusteredPinTimingInvalidator* pin_timing_invalidator,
                                 route_budgets& budgeting_inf,
                                 bool& was_rerouted,
                                 float worst_negative_slack,
                                 const RoutingPredictor& routing_predictor,
                                 bool is_flat) {
    auto& cluster_ctx = g_vpr_ctx.clustering();
    auto& route_ctx = g_vpr_ctx.mutable_routing();

    bool is_routed = false;

    connections_inf.prepare_routing_for_net(net_id);

    bool reroute_for_hold = false;
    if (budgeting_inf.if_set()) {
        reroute_for_hold = (budgeting_inf.get_should_reroute(net_id));
        reroute_for_hold &= worst_negative_slack != 0;
    }

    if (route_ctx.net_status.is_fixed(net_id)) { /* Skip pre-routed nets. */
        is_routed = true;
    } else if (cluster_ctx.clb_nlist.net_is_ignored(net_id)) { /* Skip ignored nets. */
        is_routed = true;
    } else if (!(reroute_for_hold) && should_route_net(net_id, connections_inf, true) == false) {
        is_routed = true;
    } else {
        // track time spent vs fanout
        profiling::net_fanout_start();

        is_routed = timing_driven_route_net(router,
                                            net_id,
                                            itry,
                                            pres_fac,
                                            router_opts,
                                            connections_inf,
                                            router_stats,
                                            pin_criticality,
                                            rt_node_of_sink,
                                            net_delay[net_id].data(),
                                            netlist_pin_lookup,
                                            timing_info,
                                            pin_timing_invalidator,
                                            budgeting_inf,
                                            worst_negative_slack,
                                            routing_predictor,
                                            is_flat);

        profiling::net_fanout_end(cluster_ctx.clb_nlist.net_sinks(net_id).size());

        /* Impossible to route? (disconnected rr_graph) */
        if (is_routed) {
            route_ctx.net_status.set_is_routed(net_id, true);
        } else {
            VTR_LOG("Routing failed for net %d\n", net_id);
        }

        was_rerouted = true; //Flag to record whether routing was actually changed
    }
    return (is_routed);
}

/*
 * NOTE:
 * Suggest using a timing_driven_route_structs struct. Memory is managed for you
 */
void alloc_timing_driven_route_structs(float** pin_criticality_ptr,
                                       int** sink_order_ptr,
                                       t_rt_node*** rt_node_of_sink_ptr) {
    /* Allocates all the structures needed only by the timing-driven router.   */

    int max_sinks = std::max(get_max_pins_per_net() - 1, 0);

    *pin_criticality_ptr = new float[max_sinks + 1]; /* First sink is pin #1.*/
    *sink_order_ptr = new int[max_sinks + 1];
    *rt_node_of_sink_ptr = new t_rt_node*[max_sinks + 1];

    /* Element 0 should be an invalid value so we are likely to crash if we accidentally use it. */
    (*pin_criticality_ptr)[0] = -1;
    (*sink_order_ptr)[0] = -1;
    (*rt_node_of_sink_ptr)[0] = nullptr;

    alloc_route_tree_timing_structs(false);
}

/*
 * NOTE:
 * Suggest using a timing_driven_route_structs struct. Memory is managed for you
 */
void free_timing_driven_route_structs(float* pin_criticality, int* sink_order, t_rt_node** rt_node_of_sink) {
    /* Frees all the structures needed only by the timing-driven router.        */

    // coverity[offset_free : Intentional]
    delete[](pin_criticality);
    // coverity[offset_free : Intentional]
    delete[](sink_order);
    // coverity[offset_free : Intentional]
    delete[](rt_node_of_sink);

    free_route_tree_timing_structs();
}

timing_driven_route_structs::timing_driven_route_structs() {
    alloc_timing_driven_route_structs(&pin_criticality,
                                      &sink_order,
                                      &rt_node_of_sink);
}

timing_driven_route_structs::~timing_driven_route_structs() {
    free_timing_driven_route_structs(pin_criticality,
                                     sink_order,
                                     rt_node_of_sink);
}

void increase_short_path_crit_if_congested(std::vector<ClusterNetId>& rerouted_nets,
                                           route_budgets& budgeting_inf,
                                           int itry) {
    if (budgeting_inf.if_set() && itry > 9) {
        for (auto net_id : rerouted_nets) {
            if (budgeting_inf.get_should_reroute(net_id)) {
                budgeting_inf.update_congestion_times(net_id);
            } else {
                budgeting_inf.not_congested_this_iteration(net_id);
            }
            budgeting_inf.increase_short_crit(net_id, 4);
        }
    }
}

int get_max_pins_per_net() {
    int max_pins_per_net = 0;

    auto& cluster_ctx = g_vpr_ctx.clustering();

    for (auto net_id : cluster_ctx.clb_nlist.nets()) {
        if (!cluster_ctx.clb_nlist.net_is_ignored(net_id))
            max_pins_per_net = std::max(max_pins_per_net, (int)cluster_ctx.clb_nlist.net_pins(net_id).size());
    }

    return (max_pins_per_net);
}

struct Criticality_comp {
    const float* criticality;

    Criticality_comp(const float* calculated_criticalities)
        : criticality{calculated_criticalities} {
    }

    bool operator()(int a, int b) const {
        return criticality[a] > criticality[b];
    }
};

template<typename ConnectionRouter>
bool timing_driven_route_net(ConnectionRouter& router,
                             ClusterNetId net_id,
                             int itry,
                             float pres_fac,
                             const t_router_opts& router_opts,
                             CBRR& connections_inf,
                             RouterStats& router_stats,
                             float* pin_criticality,
                             t_rt_node** rt_node_of_sink,
                             float* net_delay,
                             const ClusteredPinAtomPinsLookup& netlist_pin_lookup,
                             std::shared_ptr<SetupHoldTimingInfo> timing_info,
                             ClusteredPinTimingInvalidator* pin_timing_invalidator,
                             route_budgets& budgeting_inf,
                             float worst_neg_slack,
                             const RoutingPredictor& routing_predictor,
                             bool is_flat) {
    /* Returns true as long as found some way to hook up this net, even if that *
     * way resulted in overuse of resources (congestion).  If there is no way   *
     * to route this net, even ignoring congestion, it returns false.  In this  *
     * case the rr_graph is disconnected and you can give up.                   */
    auto& cluster_ctx = g_vpr_ctx.clustering();
    auto& device_ctx = g_vpr_ctx.device();
    const auto& rr_graph = device_ctx.rr_graph;
    auto& route_ctx = g_vpr_ctx.routing();

    unsigned int num_sinks = cluster_ctx.clb_nlist.net_sinks(net_id).size();

    VTR_LOGV_DEBUG(f_router_debug, "Routing Net %zu (%zu sinks)\n", size_t(net_id), num_sinks);

    t_rt_node* rt_root;
    rt_root = setup_routing_resources(itry,
                                      net_id,
                                      num_sinks,
                                      router_opts.min_incremental_reroute_fanout,
                                      connections_inf,
                                      rt_node_of_sink,
                                      router_opts,
                                      check_hold(router_opts, worst_neg_slack));

    bool high_fanout = is_high_fanout(num_sinks, router_opts.high_fanout_threshold);

    SpatialRouteTreeLookup spatial_route_tree_lookup;
    if (high_fanout) {
        spatial_route_tree_lookup = build_route_tree_spatial_lookup(net_id, rt_root);
    }

    // after this point the route tree is correct
    // remaining_targets from this point on are the **pin indices** that have yet to be routed
    auto& remaining_targets = connections_inf.get_remaining_targets();

    // calculate criticality of remaining target pins
    for (int ipin : remaining_targets) {
        if (timing_info) {
            auto clb_pin = cluster_ctx.clb_nlist.net_pin(net_id, ipin);
            if (!route_ctx.is_clock_net[net_id]) {
                pin_criticality[ipin] = calculate_clb_net_pin_criticality(*timing_info, netlist_pin_lookup, clb_pin);
            } else {
                // Use max_criticality for clock nets.
                // calculate_clb_net_pin_criticality likely doesn't generate
                // good values for clock nets.
                //
                // This will cause them to use min delay paths rather than
                // avoid congestion. As a future enchancement, the clock nets
                // should likely route for min slew, but that is a larger
                // change.
                pin_criticality[ipin] = router_opts.max_criticality;
            }

            /* Pin criticality is between 0 and 1.
             * Shift it downwards by 1 - max_criticality (max_criticality is 0.99 by default,
             * so shift down by 0.01) and cut off at 0.  This means that all pins with small
             * criticalities (<0.01) get criticality 0 and are ignored entirely, and everything
             * else becomes a bit less critical. This effect becomes more pronounced if
             * max_criticality is set lower. */
            // VTR_ASSERT(pin_criticality[ipin] > -0.01 && pin_criticality[ipin] < 1.01);
            pin_criticality[ipin] = std::max(pin_criticality[ipin] - (1.0 - router_opts.max_criticality), 0.0);

            /* Take pin criticality to some power (1 by default). */
            pin_criticality[ipin] = std::pow(pin_criticality[ipin], router_opts.criticality_exp);

            /* Cut off pin criticality at max_criticality. */
            pin_criticality[ipin] = std::min(pin_criticality[ipin], router_opts.max_criticality);
        } else {
            //No timing info, implies we want a min delay routing, so use criticality of 1.
            pin_criticality[ipin] = 1.;
        }
    }

    // compare the criticality of different sink nodes
    sort(begin(remaining_targets), end(remaining_targets), Criticality_comp{pin_criticality});

    /* Update base costs according to fanout and criticality rules */
    update_rr_base_costs(num_sinks);

    t_conn_delay_budget conn_delay_budget;
    t_conn_cost_params cost_params;
    cost_params.astar_fac = router_opts.astar_fac;
    cost_params.bend_cost = router_opts.bend_cost;
    cost_params.pres_fac = pres_fac;
    cost_params.delay_budget = ((budgeting_inf.if_set()) ? &conn_delay_budget : nullptr);
    cost_params.bias = router_opts.sbNode_lookahead_factor;
    cost_params.offpath_penalty = router_opts.offpath_penalty;
    cost_params.detailed_router = router_opts.detailed_router;

    // Pre-route to clock source for clock nets (marked as global nets)
    if (cluster_ctx.clb_nlist.net_is_global(net_id) && router_opts.two_stage_clock_routing) {
        //VTR_ASSERT(router_opts.clock_modeling == DEDICATED_NETWORK);
        int sink_node = device_ctx.virtual_clock_network_root_idx;

        enable_router_debug(router_opts, net_id, sink_node, itry, &router);

        VTR_LOGV_DEBUG(f_router_debug, "Pre-routing global net %zu\n", size_t(net_id));

        // Set to the max timing criticality which should intern minimize clock insertion
        // delay by selecting a direct route from the clock source to the virtual sink
        cost_params.criticality = router_opts.max_criticality;
    	cost_params.bias = router_opts.sbNode_lookahead_factor;
        if (!timing_driven_pre_route_to_clock_root(
                router,
                net_id,
                sink_node,
                cost_params,
                router_opts.high_fanout_threshold,
                rt_root,
                spatial_route_tree_lookup,
                router_stats,
                is_flat)) {
            return false;
        }
    }

    if (budgeting_inf.if_set()) {
        budgeting_inf.set_should_reroute(net_id, false);
    }

    // explore in order of decreasing criticality (no longer need sink_order array)
    for (unsigned itarget = 0; itarget < remaining_targets.size(); ++itarget) {
        int target_pin = remaining_targets[itarget];

        int sink_rr = route_ctx.net_rr_terminals[net_id][target_pin];

        enable_router_debug(router_opts, net_id, sink_rr, itry, &router);

        VTR_LOGV_DEBUG(f_router_debug, "Routing Net %zu (%zu sinks)\n", size_t(net_id), num_sinks);

        cost_params.criticality = pin_criticality[target_pin];
    	cost_params.bias = router_opts.sbNode_lookahead_factor;

        if (budgeting_inf.if_set()) {
            conn_delay_budget.max_delay = budgeting_inf.get_max_delay_budget(net_id, target_pin);
            conn_delay_budget.target_delay = budgeting_inf.get_delay_target(net_id, target_pin);
            conn_delay_budget.min_delay = budgeting_inf.get_min_delay_budget(net_id, target_pin);
            conn_delay_budget.short_path_criticality = budgeting_inf.get_crit_short_path(net_id, target_pin);
            conn_delay_budget.routing_budgets_algorithm = router_opts.routing_budgets_algorithm;
        }

        profiling::conn_start();
	std::set<int> branch_nodes;
        // build a branch in the route tree to the target
        if (!timing_driven_route_sink(router,
                                      net_id,
                                      itarget,
                                      target_pin,
                                      cost_params,
                                      router_opts,
                                      rt_root, rt_node_of_sink,
                                      spatial_route_tree_lookup,
                                      router_stats,
                                      budgeting_inf,
                                      routing_predictor,
                                      is_flat,
				      branch_nodes,
				      itry))
            return false;

        profiling::conn_finish(route_ctx.net_rr_terminals[net_id][0],
                               sink_rr,
                               pin_criticality[target_pin]);

        ++router_stats.connections_routed;
    } // finished all sinks

    ++router_stats.nets_routed;
    profiling::net_finish();

    /* For later timing analysis. */

    // may have to update timing delay of the previously legally reached sinks since downstream capacitance could be changed
    update_net_delays_from_route_tree(net_delay, rt_node_of_sink, net_id, timing_info.get(), pin_timing_invalidator);

    if (router_opts.update_lower_bound_delays) {
        for (int ipin : remaining_targets) {
            connections_inf.update_lower_bound_connection_delay(net_id, ipin, net_delay[ipin]);
        }
    }

    if (!cluster_ctx.clb_nlist.net_is_ignored(net_id)) {
        for (unsigned ipin = 1; ipin < cluster_ctx.clb_nlist.net_pins(net_id).size(); ++ipin) {
            if (net_delay[ipin] == 0) { // should be SOURCE->OPIN->IPIN->SINK
                VTR_ASSERT(rr_graph.node_type(RRNodeId(rt_node_of_sink[ipin]->parent_node->parent_node->inode)) == OPIN);
            }
        }
    }
    VTR_ASSERT_MSG(route_ctx.rr_node_route_inf[rt_root->inode].occ() <= rr_graph.node_capacity(RRNodeId(rt_root->inode)), "SOURCE should never be congested");

    // route tree is not kept persistent since building it from the traceback the next iteration takes almost 0 time
    VTR_LOGV_DEBUG(f_router_debug, "Routed Net %zu (%zu sinks)\n", size_t(net_id), num_sinks);

    free_route_tree(rt_root);
    router.empty_rcv_route_tree_set();
    return (true);
}

template<typename ConnectionRouter>
static bool timing_driven_pre_route_to_clock_root(
    ConnectionRouter& router,
    ClusterNetId net_id,
    int sink_node,
    const t_conn_cost_params cost_params,
    int high_fanout_threshold,
    t_rt_node* rt_root,
    SpatialRouteTreeLookup& spatial_rt_lookup,
    RouterStats& router_stats,
    bool is_flat) {
    const auto& device_ctx = g_vpr_ctx.device();
    auto& route_ctx = g_vpr_ctx.mutable_routing();
    auto& cluster_ctx = g_vpr_ctx.clustering();
    auto& m_route_ctx = g_vpr_ctx.mutable_routing();

    bool high_fanout = is_high_fanout(cluster_ctx.clb_nlist.net_sinks(net_id).size(), high_fanout_threshold);
    VTR_LOG("[SHA] Routing pre route to clock root\n");
    VTR_LOGV_DEBUG(f_router_debug, "Net %zu pre-route to (%s)\n", size_t(net_id), describe_rr_node(device_ctx.rr_graph, device_ctx.grid, device_ctx.rr_indexed_data, sink_node, is_flat).c_str());

    profiling::sink_criticality_start();

    VTR_ASSERT_DEBUG(verify_traceback_route_tree_equivalent(route_ctx.trace[net_id].head, rt_root));

    t_bb bounding_box = route_ctx.route_bb[net_id];

    router.clear_modified_rr_node_info();

    bool found_path;
    t_heap cheapest;
    //std::string conn_id = std::to_string(static_cast<std::size_t>(net_id));
    int sink_id = 0;
    std::set<int> branch_nodes;
    int itry = 1;
    std::tie(found_path, cheapest) = router.timing_driven_route_connection_from_route_tree(
        rt_root,
        sink_node,
        cost_params,
        bounding_box,
        router_stats,
        net_id, sink_id, branch_nodes, itry);

    // TODO: Parts of the rest of this function are repetitive to code in timing_driven_route_sink. Should refactor.
    if (!found_path) {
        ClusterBlockId src_block = cluster_ctx.clb_nlist.net_driver_block(net_id);
        VTR_LOG("Failed to route connection from '%s' to '%s' for net '%s' (#%zu)\n",
                cluster_ctx.clb_nlist.block_name(src_block).c_str(),
                describe_rr_node(device_ctx.rr_graph, device_ctx.grid, device_ctx.rr_indexed_data, sink_node, is_flat).c_str(),
                cluster_ctx.clb_nlist.net_name(net_id).c_str(),
                size_t(net_id));
        if (f_router_debug) {
            update_screen(ScreenUpdatePriority::MAJOR, "Unable to route connection.", ROUTING, nullptr);
        }
        return false;
    }

    profiling::sink_criticality_end(cost_params.criticality);

    /* NB:  In the code below I keep two records of the partial routing:  the   *
     * traceback and the route_tree.  The route_tree enables fast recomputation *
     * of the Elmore delay to each node in the partial routing.  The traceback  *
     * lets me reuse all the routines written for breadth-first routing, which  *
     * all take a traceback structure as input.                                 */

    /* This is a special pre-route to a sink that does not correspond to any    *
     * netlist pin, but which can be reached from the global clock root drive   *
     * points. Therefore, we can set the net pin index of the sink node to      *
     * OPEN (meaning illegal) as it is not meaningful for this sink.            */

    t_trace* new_route_start_tptr = update_traceback(&cheapest, OPEN, net_id);
    VTR_ASSERT_DEBUG(validate_traceback(route_ctx.trace[net_id].head));
    update_route_tree(&cheapest, OPEN, ((high_fanout) ? &spatial_rt_lookup : nullptr));
    VTR_ASSERT_DEBUG(verify_route_tree(rt_root));
    VTR_ASSERT_DEBUG(verify_traceback_route_tree_equivalent(route_ctx.trace[net_id].head, rt_root));
    VTR_ASSERT_DEBUG(!high_fanout || validate_route_tree_spatial_lookup(rt_root, spatial_rt_lookup));
    if (f_router_debug) {
        std::string msg = vtr::string_fmt("Routed Net %zu connection to RR node %d successfully", size_t(net_id), sink_node);
        update_screen(ScreenUpdatePriority::MAJOR, msg.c_str(), ROUTING, nullptr);
    }
    pathfinder_update_path_occupancy(new_route_start_tptr, 1);

    // need to guarantee ALL nodes' path costs are HUGE_POSITIVE_FLOAT at the start of routing to a sink
    // do this by resetting all the path_costs that have been touched while routing to the current sink
    router.reset_path_costs();

    // Post route trace back and route tree clean up:
    // - remove sink from trace back and route tree
    // - fix routing for all nodes leading to the sink
    // - free up virtual sink occupancy
    disable_expansion_and_remove_sink_from_route_tree_nodes(rt_root);
    VTR_LOGV_DEBUG(f_router_debug, "Traceback tail before update %d \n",
                   route_ctx.trace[net_id].tail->index);
    drop_traceback_tail(net_id);
    VTR_LOGV_DEBUG(f_router_debug, "Updated traceback ptrs: %d %d \n",
                   route_ctx.trace[net_id].head->index, route_ctx.trace[net_id].tail->index);
    m_route_ctx.rr_node_route_inf[sink_node].set_occ(0);

    // routed to a sink successfully
    return true;
}

template<typename ConnectionRouter>
static bool timing_driven_route_sink(
    ConnectionRouter& router,
    ClusterNetId net_id,
    unsigned itarget,
    int target_pin,
    const t_conn_cost_params cost_params,
    const t_router_opts& router_opts,
    t_rt_node* rt_root,
    t_rt_node** rt_node_of_sink,
    SpatialRouteTreeLookup& spatial_rt_lookup,
    RouterStats& router_stats,
    route_budgets& budgeting_inf,
    const RoutingPredictor& routing_predictor,
    bool is_flat,
    std::set<int> branch_nodes,
    int itry) {
    /* Build a path from the existing route tree rooted at rt_root to the target_node
     * add this branch to the existing route tree and update pathfinder costs and rr_node_route_inf to reflect this */
    const auto& device_ctx = g_vpr_ctx.device();
    auto& route_ctx = g_vpr_ctx.mutable_routing();
    auto& cluster_ctx = g_vpr_ctx.clustering();

    profiling::sink_criticality_start();

    int sink_node = route_ctx.net_rr_terminals[net_id][target_pin];
    VTR_LOGV_DEBUG(f_router_debug, "Net %zu Target %d (%s)\n", size_t(net_id), itarget, describe_rr_node(device_ctx.rr_graph, device_ctx.grid, device_ctx.rr_indexed_data, sink_node, is_flat).c_str());

    VTR_ASSERT_DEBUG(verify_traceback_route_tree_equivalent(route_ctx.trace[net_id].head, rt_root));

    router.clear_modified_rr_node_info();

    bool found_path;
    t_heap cheapest;
    t_bb bounding_box = route_ctx.route_bb[net_id];

    bool net_is_global = cluster_ctx.clb_nlist.net_is_global(net_id);
    bool high_fanout = is_high_fanout(cluster_ctx.clb_nlist.net_sinks(net_id).size(), router_opts.high_fanout_threshold);
    constexpr float HIGH_FANOUT_CRITICALITY_THRESHOLD = 0.9;
    bool sink_critical = (cost_params.criticality > HIGH_FANOUT_CRITICALITY_THRESHOLD);
    bool net_is_clock = route_ctx.is_clock_net[net_id] != 0;
    //std::string conn_id = std::to_string(static_cast<std::size_t>(net_id)) + "_" + std::to_string(target_pin);
    
    //VTR_LOG("connection_id: %s %zu %d\n", conn_id.c_str(), size_t(net_id), target_pin);
    //std::cout << "This is a string: " << conn_id << std::endl;
    //ClusterNetId net_sink_id = ClusterNetId(net_id_with_sink_id);
    //We normally route high fanout nets by only adding spatially close-by routing to the heap (reduces run-time).
    //However, if the current sink is 'critical' from a timing perspective, we put the entire route tree back onto
    //the heap to ensure it has more flexibility to find the best path.
    if (high_fanout && !sink_critical && !net_is_global && !net_is_clock && -routing_predictor.get_slope() > router_opts.high_fanout_max_slope) {
        std::tie(found_path, cheapest) = router.timing_driven_route_connection_from_route_tree_high_fanout(rt_root,
                                                                                                           sink_node,
                                                                                                           cost_params,
                                                                                                           bounding_box,
                                                                                                           spatial_rt_lookup,
                                                                                                           router_stats, net_id, target_pin, branch_nodes, itry);
    } else {
        std::tie(found_path, cheapest) = router.timing_driven_route_connection_from_route_tree(rt_root,
                                                                                               sink_node,
                                                                                               cost_params,
                                                                                               bounding_box,
                                                                                               router_stats, net_id, target_pin, branch_nodes, itry);
    }

    if (!found_path) {
        ClusterBlockId src_block = cluster_ctx.clb_nlist.net_driver_block(net_id);
        ClusterBlockId sink_block = cluster_ctx.clb_nlist.pin_block(*(cluster_ctx.clb_nlist.net_pins(net_id).begin() + target_pin));
        VTR_LOG("Failed to route connection from '%s' to '%s' for net '%s' (#%zu)\n",
                cluster_ctx.clb_nlist.block_name(src_block).c_str(),
                cluster_ctx.clb_nlist.block_name(sink_block).c_str(),
                cluster_ctx.clb_nlist.net_name(net_id).c_str(),
                size_t(net_id));
        if (f_router_debug) {
            update_screen(ScreenUpdatePriority::MAJOR, "Unable to route connection.", ROUTING, nullptr);
        }
        return false;
    }

    profiling::sink_criticality_end(cost_params.criticality);

    /* NB:  In the code below I keep two records of the partial routing:  the   *
     * traceback and the route_tree.  The route_tree enables fast recomputation *
     * of the Elmore delay to each node in the partial routing.  The traceback  *
     * lets me reuse all the routines written for breadth-first routing, which  *
     * all take a traceback structure as input.                                 */

    int inode = cheapest.index;
    route_ctx.rr_node_route_inf[inode].target_flag--; /* Connected to this SINK. */
    t_trace* new_route_start_tptr = update_traceback(&cheapest, target_pin, net_id);

    VTR_ASSERT_DEBUG(validate_traceback(route_ctx.trace[net_id].head));

    rt_node_of_sink[target_pin] = update_route_tree(&cheapest, target_pin, ((high_fanout) ? &spatial_rt_lookup : nullptr));
    VTR_ASSERT_DEBUG(verify_route_tree(rt_root));
    VTR_ASSERT_DEBUG(verify_traceback_route_tree_equivalent(route_ctx.trace[net_id].head, rt_root));
    VTR_ASSERT_DEBUG(!high_fanout || validate_route_tree_spatial_lookup(rt_root, spatial_rt_lookup));
    if (f_router_debug) {
        std::string msg = vtr::string_fmt("Routed Net %zu connection %d to RR node %d successfully", size_t(net_id), itarget, sink_node);
        update_screen(ScreenUpdatePriority::MAJOR, msg.c_str(), ROUTING, nullptr);
    }

    if (budgeting_inf.if_set() && cheapest.path_data != nullptr && cost_params.delay_budget) {
        if (cheapest.path_data->backward_delay < cost_params.delay_budget->min_delay) {
            budgeting_inf.set_should_reroute(net_id, true);
        }
    }

    pathfinder_update_path_occupancy(new_route_start_tptr, 1);

    // need to guarantee ALL nodes' path costs are HUGE_POSITIVE_FLOAT at the start of routing to a sink
    // do this by resetting all the path_costs that have been touched while routing to the current sink
    router.reset_path_costs();

    // routed to a sink successfully
    return true;
}

static t_rt_node* setup_routing_resources(int itry,
                                          ClusterNetId net_id,
                                          unsigned num_sinks,
                                          int min_incremental_reroute_fanout,
                                          CBRR& connections_inf,
                                          t_rt_node** rt_node_of_sink,
                                          const t_router_opts& router_opts,
                                          bool ripup_high_fanout_nets) {
    /* Build and return a partial route tree from the legal connections from last iteration.
     * along the way do:
     * 	update pathfinder costs to be accurate to the partial route tree
     *	update the net's traceback to be accurate to the partial route tree
     * 	find and store the pins that still need to be reached in incremental_rerouting_resources.remaining_targets
     * 	find and store the rt nodes that have been reached in incremental_rerouting_resources.reached_rt_sinks
     *	mark the rr_node sinks as targets to be reached */

    auto& route_ctx = g_vpr_ctx.routing();

    t_rt_node* rt_root;

    // for nets below a certain size (min_incremental_reroute_fanout), rip up any old routing
    // otherwise, we incrementally reroute by reusing legal parts of the previous iteration
    // convert the previous iteration's traceback into the starting route tree for this iteration
    if ((int)num_sinks < min_incremental_reroute_fanout || itry == 1 || ripup_high_fanout_nets) {
        profiling::net_rerouted();

        // rip up the whole net
        pathfinder_update_path_occupancy(route_ctx.trace[net_id].head, -1);
        free_traceback(net_id);

        rt_root = init_route_tree_to_source(net_id);
        for (unsigned int sink_pin = 1; sink_pin <= num_sinks; ++sink_pin)
            connections_inf.toreach_rr_sink(sink_pin);
        // since all connections will be rerouted for this net, clear all of net's forced reroute flags
        connections_inf.clear_force_reroute_for_net();

        // when we don't prune the tree, we also don't know the sink node indices
        // thus we'll use functions that act on pin indices like mark_ends instead
        // of their versions that act on node indices directly like mark_remaining_ends
        mark_ends(net_id);
    } else {
        auto& reached_rt_sinks = connections_inf.get_reached_rt_sinks();
        auto& remaining_targets = connections_inf.get_remaining_targets();

        profiling::net_rebuild_start();

        // convert the previous iteration's traceback into a route tree
        rt_root = traceback_to_route_tree(net_id);

        //Sanity check that route tree and traceback are equivalent before pruning
        VTR_ASSERT_DEBUG(verify_traceback_route_tree_equivalent(route_ctx.trace[net_id].head, rt_root));

        // check for edge correctness
        VTR_ASSERT_SAFE(is_valid_skeleton_tree(rt_root));

        // Skip this check if RCV is enabled, as RCV can use another method to cause reroutes
        VTR_ASSERT_SAFE(should_route_net(net_id, connections_inf, true) || router_opts.routing_budgets_algorithm == YOYO);

        //Prune the branches of the tree that don't legally lead to sinks
        rt_root = prune_route_tree(rt_root, connections_inf);

        //Now that the tree has been pruned, we can free the old traceback
        // NOTE: this must happen *after* pruning since it changes the
        //       recorded congestion
        pathfinder_update_path_occupancy(route_ctx.trace[net_id].head, -1);
        free_traceback(net_id);

        if (rt_root) { //Partially pruned
            profiling::route_tree_preserved();

            //Since we have a valid partial routing (to at least one SINK)
            //we need to make sure the traceback is synchronized to the route tree
            traceback_from_route_tree(net_id, rt_root, reached_rt_sinks.size());

            //Sanity check the traceback for self-consistency
            VTR_ASSERT_DEBUG(validate_traceback(route_ctx.trace[net_id].head));

            //Sanity check that route tree and traceback are equivalent after pruning
            VTR_ASSERT_DEBUG(verify_traceback_route_tree_equivalent(route_ctx.trace[net_id].head, rt_root));

            // put the updated occupancies of the route tree nodes back into pathfinder
            pathfinder_update_path_occupancy(route_ctx.trace[net_id].head, 1);

        } else { //Fully destroyed
            profiling::route_tree_pruned();

            //Initialize only to source
            rt_root = init_route_tree_to_source(net_id);

            //NOTE: We leave the traceback uninitialized, so update_traceback()
            //      will correctly add the SOURCE node when the branch to
            //      the first SINK is found.
            VTR_ASSERT(route_ctx.trace[net_id].head == nullptr);
            VTR_ASSERT(route_ctx.trace[net_id].tail == nullptr);
            VTR_ASSERT(route_ctx.trace_nodes[net_id].empty());
        }

        //Update R/C
        load_new_subtree_R_upstream(rt_root);
        load_new_subtree_C_downstream(rt_root);

        VTR_ASSERT(reached_rt_sinks.size() + remaining_targets.size() == num_sinks);

        //Record current routing
        add_route_tree_to_rr_node_lookup(rt_root);

        // give lookup on the reached sinks
        for (t_rt_node* sink_node : reached_rt_sinks) {
            rt_node_of_sink[sink_node->net_pin_index] = sink_node;
        }

        profiling::net_rebuild_end(num_sinks, remaining_targets.size());

        // check for R_upstream C_downstream and edge correctness
        VTR_ASSERT_SAFE(is_valid_route_tree(rt_root));
        // congestion should've been pruned away
        VTR_ASSERT_SAFE(is_uncongested_route_tree(rt_root));

        // mark remaining ends
        mark_remaining_ends(net_id, remaining_targets);

        // still need to calculate the tree's time delay (0 Tarrival means from SOURCE)
        load_route_tree_Tdel(rt_root, 0);

        // mark the lookup (rr_node_route_inf) for existing tree elements as NO_PREVIOUS so add_to_path stops when it reaches one of them
        load_route_tree_rr_route_inf(rt_root);
    }

    // completed constructing the partial route tree and updated all other data structures to match
    return rt_root;
}

void disable_expansion_and_remove_sink_from_route_tree_nodes(t_rt_node* rt_node) {
    /* Remove sink in route tree and mark all nodes
     * leading to the sink as unexpandable.
     */
    auto& device_ctx = g_vpr_ctx.device();
    const auto& rr_graph = device_ctx.rr_graph;
    t_rt_node* child_node;
    t_linked_rt_edge* linked_rt_edge;
    linked_rt_edge = rt_node->u.child_list;

    while (linked_rt_edge != nullptr) {
        child_node = linked_rt_edge->child;
        if (rr_graph.node_type(RRNodeId(child_node->inode)) == SINK) {
            VTR_LOGV_DEBUG(f_router_debug,
                           "Removing sink %d from route tree\n", child_node->inode);
            rt_node->u.child_list = nullptr;
            rt_node->u.next = nullptr;
            free(child_node);
            break;
        } else {
            rt_node->re_expand = false;
            VTR_LOGV_DEBUG(f_router_debug,
                           "unexpanding: %d in route tree\n", rt_node->inode);
        }
        disable_expansion_and_remove_sink_from_route_tree_nodes(child_node);
        linked_rt_edge = linked_rt_edge->next;
    }
}

void update_rr_base_costs(int fanout) {
    /* Changes the base costs of different types of rr_nodes according to the  *
     * criticality, fanout, etc. of the current net being routed (net_id).       */
    auto& device_ctx = g_vpr_ctx.mutable_device();

    float factor;
    size_t index;

    /* Other reasonable values for factor include fanout and 1 */
    factor = sqrt(fanout);

    for (index = CHANX_COST_INDEX_START; index < device_ctx.rr_indexed_data.size(); index++) {
        if (device_ctx.rr_indexed_data[RRIndexedDataId(index)].T_quadratic > 0.) { /* pass transistor */
            device_ctx.rr_indexed_data[RRIndexedDataId(index)].base_cost = device_ctx.rr_indexed_data[RRIndexedDataId(index)].saved_base_cost * factor;
        } else {
            device_ctx.rr_indexed_data[RRIndexedDataId(index)].base_cost = device_ctx.rr_indexed_data[RRIndexedDataId(index)].saved_base_cost;
        }
    }
}

static bool timing_driven_check_net_delays(ClbNetPinsMatrix<float>& net_delay) {
    constexpr float ERROR_TOL = 0.0001;

    /* Checks that the net delays computed incrementally during timing driven    *
     * routing match those computed from scratch by the net_delay.c module.      */
    auto& cluster_ctx = g_vpr_ctx.clustering();

    unsigned int ipin;
    ClbNetPinsMatrix<float> net_delay_check = make_net_pins_matrix<float>(cluster_ctx.clb_nlist);

    load_net_delay_from_routing(net_delay_check);

    for (auto net_id : cluster_ctx.clb_nlist.nets()) {
        for (ipin = 1; ipin < cluster_ctx.clb_nlist.net_pins(net_id).size(); ipin++) {
            if (net_delay_check[net_id][ipin] == 0.) { /* Should be only GLOBAL nets */
                if (fabs(net_delay[net_id][ipin]) > ERROR_TOL) {
                    VPR_ERROR(VPR_ERROR_ROUTE,
                              "in timing_driven_check_net_delays: net %lu pin %d.\n"
                              "\tIncremental calc. net_delay is %g, but from scratch net delay is %g.\n",
                              size_t(net_id), ipin, net_delay[net_id][ipin], net_delay_check[net_id][ipin]);
                }
            } else {
                float error = fabs(1.0 - net_delay[net_id][ipin] / net_delay_check[net_id][ipin]);
                if (error > ERROR_TOL) {
                    VPR_ERROR(VPR_ERROR_ROUTE,
                              "in timing_driven_check_net_delays: net %d pin %lu.\n"
                              "\tIncremental calc. net_delay is %g, but from scratch net delay is %g.\n",
                              size_t(net_id), ipin, net_delay[net_id][ipin], net_delay_check[net_id][ipin]);
                }
            }
        }
    }

    return true;
}

/* Detect if net should be routed or not */
static bool should_route_net(ClusterNetId net_id, CBRR& connections_inf, bool if_force_reroute) {
    auto& route_ctx = g_vpr_ctx.routing();
    auto& device_ctx = g_vpr_ctx.device();
    const auto& rr_graph = device_ctx.rr_graph;

    t_trace* tptr = route_ctx.trace[net_id].head;

    if (tptr == nullptr) {
        /* No routing yet. */
        return true;
    }

    for (;;) {
        int inode = tptr->index;
        int occ = route_ctx.rr_node_route_inf[inode].occ();
        int capacity = rr_graph.node_capacity(RRNodeId(inode));

        if (occ > capacity) {
            return true; /* overuse detected */
        }

        if (tptr->iswitch == OPEN) { //End of a branch
            // even if net is fully routed, not complete if parts of it should get ripped up (EXPERIMENTAL)
            if (if_force_reroute) {
                if (connections_inf.should_force_reroute_connection(inode)) {
                    return true;
                }
            }
            tptr = tptr->next; /* Skip next segment (duplicate of original branch node). */
            if (tptr == nullptr)
                break;
        }

        tptr = tptr->next;

    } /* End while loop -- did an entire traceback. */

    VTR_ASSERT(connections_inf.get_remaining_targets().empty());

    return false; /* Current route has no overuse */
}

static bool early_exit_heuristic(const t_router_opts& router_opts, const WirelengthInfo& wirelength_info) {
    /* Early exit code for cases where it is obvious that a successful route will not be found
     * Heuristic: If total wirelength used in first routing iteration is X% of total available wirelength, exit */

    if (wirelength_info.used_wirelength_ratio() > router_opts.init_wirelength_abort_threshold) {
        VTR_LOG("Wire length usage ratio %g exceeds limit of %g, fail routing.\n",
                wirelength_info.used_wirelength_ratio(),
                router_opts.init_wirelength_abort_threshold);
        return true;
    }
    return false;
}

static bool check_hold(const t_router_opts& router_opts, float worst_neg_slack) {
    /* When RCV is enabled, it's necessary to be able to completely ripup high fanout nets if there is still negative hold slack
     * Normally the router will prune the illegal branches of high fanout nets, this will bypass this */

    if (router_opts.routing_budgets_algorithm != YOYO) {
        return false;
    } else if (worst_neg_slack != 0) {
        return true;
    }
    return false;
}

static size_t calculate_wirelength_available() {
    auto& device_ctx = g_vpr_ctx.device();
    const auto& rr_graph = device_ctx.rr_graph;

    size_t available_wirelength = 0;
    // But really what's happening is that this for loop iterates over every node and determines the available wirelength
    VTR_LOG("============ Analysis: Calculating wirelength without considering IIB =========");
    for (const RRNodeId& rr_id : device_ctx.rr_graph.nodes()) {
        const t_rr_type channel_type = rr_graph.node_type(rr_id);
	int ptc_val = rr_graph.node_ptc_num(rr_id);//getting to check if in the IIB
	if (ptc_val < 304){
        	if (channel_type == CHANX || channel_type == CHANY) {
            		available_wirelength += rr_graph.node_capacity(rr_id) * rr_graph.node_length(rr_id);
        	}
	}
    }
    return available_wirelength;
}

static WirelengthInfo calculate_wirelength_info(size_t available_wirelength) {
    auto& cluster_ctx = g_vpr_ctx.clustering();

    size_t used_wirelength = 0;
    VTR_ASSERT(available_wirelength > 0);

    for (auto net_id : cluster_ctx.clb_nlist.nets()) {
        if (!cluster_ctx.clb_nlist.net_is_ignored(net_id)
            && cluster_ctx.clb_nlist.net_sinks(net_id).size() != 0) { /* Globals don't count. */
            int bends, wirelength, segments;
            get_num_bends_and_length(net_id, &bends, &wirelength, &segments);
            used_wirelength += wirelength;
        }
    }

    return WirelengthInfo(available_wirelength, used_wirelength);
}

static void print_route_status_header() {
    VTR_LOG("---- ------ ------- ---- ------- ------- ------- ----------------- --------------- -------- ---------- ---------- ---------- ---------- --------\n");
    VTR_LOG("Iter   Time    pres  BBs    Heap  Re-Rtd  Re-Rtd Overused RR Nodes      Wirelength      CPD       sTNS       sWNS       hTNS       hWNS Est Succ\n");
    VTR_LOG("      (sec)     fac Updt    push    Nets   Conns                                       (ns)       (ns)       (ns)       (ns)       (ns)     Iter\n");
    VTR_LOG("---- ------ ------- ---- ------- ------- ------- ----------------- --------------- -------- ---------- ---------- ---------- ---------- --------\n");
}

static void print_route_status(int itry, double elapsed_sec, float pres_fac, int num_bb_updated, const RouterStats& router_stats, const OveruseInfo& overuse_info, const WirelengthInfo& wirelength_info, std::shared_ptr<const SetupHoldTimingInfo> timing_info, float est_success_iteration) {
    //Iteration
    VTR_LOG("%4d", itry);

    //Elapsed Time
    VTR_LOG(" %6.1f", elapsed_sec);

    //pres_fac
    constexpr int PRES_FAC_DIGITS = 7;
    constexpr int PRES_FAC_SCI_PRECISION = 1;
    pretty_print_float(" ", pres_fac, PRES_FAC_DIGITS, PRES_FAC_SCI_PRECISION);
    //VTR_LOG(" %5.1f", pres_fac);

    //Number of bounding boxes updated
    VTR_LOG(" %4d", num_bb_updated);

    //Heap push/pop
    constexpr int HEAP_OP_DIGITS = 7;
    constexpr int HEAP_OP_SCI_PRECISION = 2;
    pretty_print_uint(" ", router_stats.heap_pushes, HEAP_OP_DIGITS, HEAP_OP_SCI_PRECISION);
    VTR_ASSERT(router_stats.heap_pops <= router_stats.heap_pushes);

    //Rerouted nets
    constexpr int NET_ROUTED_DIGITS = 7;
    constexpr int NET_ROUTED_SCI_PRECISION = 2;
    pretty_print_uint(" ", router_stats.nets_routed, NET_ROUTED_DIGITS, NET_ROUTED_SCI_PRECISION);

    //Rerouted connections
    constexpr int CONN_ROUTED_DIGITS = 7;
    constexpr int CONN_ROUTED_SCI_PRECISION = 2;
    pretty_print_uint(" ", router_stats.connections_routed, CONN_ROUTED_DIGITS, CONN_ROUTED_SCI_PRECISION);

    //Overused RR nodes
    constexpr int OVERUSE_DIGITS = 7;
    constexpr int OVERUSE_SCI_PRECISION = 2;
    pretty_print_uint(" ", overuse_info.overused_nodes, OVERUSE_DIGITS, OVERUSE_SCI_PRECISION);
    VTR_LOG(" (%6.3f%%)", overuse_info.overused_node_ratio() * 100);

    //Wirelength
    constexpr int WL_DIGITS = 7;
    constexpr int WL_SCI_PRECISION = 2;
    pretty_print_uint(" ", wirelength_info.used_wirelength(), WL_DIGITS, WL_SCI_PRECISION);
    VTR_LOG(" (%4.1f%%)", wirelength_info.used_wirelength_ratio() * 100);

    //CPD
    if (timing_info) {
        float cpd = timing_info->least_slack_critical_path().delay();
        VTR_LOG(" %#8.3f", 1e9 * cpd);
    } else {
        VTR_LOG(" %8s", "N/A");
    }

    //sTNS
    if (timing_info) {
        float sTNS = timing_info->setup_total_negative_slack();
        VTR_LOG(" % #10.4g", 1e9 * sTNS);
    } else {
        VTR_LOG(" %10s", "N/A");
    }

    //sWNS
    if (timing_info) {
        float sWNS = timing_info->setup_worst_negative_slack();
        VTR_LOG(" % #10.3f", 1e9 * sWNS);
    } else {
        VTR_LOG(" %10s", "N/A");
    }

    //hTNS
    if (timing_info) {
        float hTNS = timing_info->hold_total_negative_slack();
        VTR_LOG(" % #10.4g", 1e9 * hTNS);
    } else {
        VTR_LOG(" %10s", "N/A");
    }

    //hWNS
    if (timing_info) {
        float hWNS = timing_info->hold_worst_negative_slack();
        VTR_LOG(" % #10.3f", 1e9 * hWNS);
    } else {
        VTR_LOG(" %10s", "N/A");
    }

    //Estimated success iteration
    if (std::isnan(est_success_iteration)) {
        VTR_LOG(" %8s", "N/A");
    } else {
        VTR_LOG(" %8.0f", est_success_iteration);
    }

    VTR_LOG("\n");

    fflush(stdout);
}

static void print_overused_nodes_status(const t_router_opts& router_opts, const OveruseInfo& overuse_info) {
    //Print the index of this routing failure
    VTR_LOG("\nFailed routing attempt #%d\n", num_routing_failed);

    size_t num_overused = overuse_info.overused_nodes;
    size_t max_logged_overused_rr_nodes = router_opts.max_logged_overused_rr_nodes;

    //Overused nodes info logging upper limit
    VTR_LOG("Total number of overused nodes: %d\n", num_overused);
    if (num_overused > max_logged_overused_rr_nodes) {
        VTR_LOG("Total number of overused nodes is larger than the logging limit (%d).\n", max_logged_overused_rr_nodes);
        VTR_LOG("Displaying the first %d entries.\n", max_logged_overused_rr_nodes);
    }

    log_overused_nodes_status(max_logged_overused_rr_nodes);
    VTR_LOG("\n");
}

static void print_router_criticality_histogram(const SetupTimingInfo& timing_info, const ClusteredPinAtomPinsLookup& netlist_pin_lookup) {
    print_histogram(create_criticality_histogram(timing_info, netlist_pin_lookup, 10));
}

//Returns true if the specified net fanout is classified as high fanout
static bool is_high_fanout(int fanout, int fanout_threshold) {
    if (fanout_threshold < 0 || fanout < fanout_threshold) return false;
    return true;
}

//In heavily congested designs a static bounding box (BB) can
//become problematic for routability (it effectively enforces a
//hard blockage restricting where a net can route).
//
//For instance, the router will try to route non-critical connections
//away from congested regions, but may end up hitting the edge of the
//bounding box. Limiting how far out-of-the-way it can be routed, and
//preventing congestion from resolving.
//
//To alleviate this, we dynamically expand net bounding boxes if the net's
//*current* routing uses RR nodes 'close' to the edge of it's bounding box.
//
//The result is that connections trying to move out of the way and hitting
//their BB will have their bounding boxes will expand slowly in that direction.
//This helps spread out regions of heavy congestion (over several routing
//iterations).
//
//By growing the BBs slowly and only as needed we minimize the size of the BBs.
//This helps keep the router's graph search fast.
//
//Typically, only a small minority of nets (typically > 10%) have their BBs updated
//each routing iteration.
static size_t dynamic_update_bounding_boxes(const std::vector<ClusterNetId>& updated_nets, int high_fanout_threshold) {
    auto& device_ctx = g_vpr_ctx.device();
    auto& cluster_ctx = g_vpr_ctx.clustering();
    auto& route_ctx = g_vpr_ctx.mutable_routing();

    auto& clb_nlist = cluster_ctx.clb_nlist;
    auto& grid = device_ctx.grid;

    //Controls how close a net's routing needs to be to it's bounding box
    //before the bounding box is expanded.
    //
    //A value of zero indicates that the routing needs to be at the bounding box
    //edge
    constexpr int DYNAMIC_BB_DELTA_THRESHOLD = 0;

    //Walk through each net, calculating the bounding box of its current routing,
    //and then increase the router's bounding box if the two are close together

    int grid_xmax = grid.width() - 1;
    int grid_ymax = grid.height() - 1;

    size_t num_bb_updated = 0;

    for (ClusterNetId net : updated_nets) {
        t_trace* routing_head = route_ctx.trace[net].head;

        if (routing_head == nullptr) continue; //Skip if no routing

        //We do not adjust the bounding boxes of high fanout nets, since they
        //use different bounding boxes based on the target location.
        //
        //This ensures that the delta values calculated below are always non-negative
        if (is_high_fanout(clb_nlist.net_sinks(net).size(), high_fanout_threshold)) continue;

        t_bb curr_bb = calc_current_bb(routing_head);

        t_bb& router_bb = route_ctx.route_bb[net];

        //Calculate the distances between the net's used RR nodes and
        //the router's bounding box
        int delta_xmin = curr_bb.xmin - router_bb.xmin;
        int delta_xmax = router_bb.xmax - curr_bb.xmax;
        int delta_ymin = curr_bb.ymin - router_bb.ymin;
        int delta_ymax = router_bb.ymax - curr_bb.ymax;

        //Note that if the net uses non-configurable switches it's routing
        //may end-up outside the bounding boxes, so the delta values may be
        //negative. The code below will expand the bounding box in those
        //cases.

        //Expand each dimension by one if within DYNAMIC_BB_DELTA_THRESHOLD threshold
        bool updated_bb = false;
        if (delta_xmin <= DYNAMIC_BB_DELTA_THRESHOLD && router_bb.xmin > 0) {
            --router_bb.xmin;
            updated_bb = true;
        }

        if (delta_ymin <= DYNAMIC_BB_DELTA_THRESHOLD && router_bb.ymin > 0) {
            --router_bb.ymin;
            updated_bb = true;
        }

        if (delta_xmax <= DYNAMIC_BB_DELTA_THRESHOLD && router_bb.xmax < grid_xmax) {
            ++router_bb.xmax;
            updated_bb = true;
        }

        if (delta_ymax <= DYNAMIC_BB_DELTA_THRESHOLD && router_bb.ymax < grid_ymax) {
            ++router_bb.ymax;
            updated_bb = true;
        }

        if (updated_bb) {
            ++num_bb_updated;
            //VTR_LOG("Expanded net %6zu router BB to (%d,%d)x(%d,%d) based on net RR node BB (%d,%d)x(%d,%d)\n", size_t(net),
            //router_bb.xmin, router_bb.ymin, router_bb.xmax, router_bb.ymax,
            //curr_bb.xmin, curr_bb.ymin, curr_bb.xmax, curr_bb.ymax);
        }
    }
    return num_bb_updated;
}

//Returns the bounding box of a net's used routing resources
static t_bb calc_current_bb(const t_trace* head) {
    auto& device_ctx = g_vpr_ctx.device();
    const auto& rr_graph = device_ctx.rr_graph;
    auto& grid = device_ctx.grid;

    t_bb bb;
    bb.xmin = grid.width() - 1;
    bb.ymin = grid.height() - 1;
    bb.xmax = 0;
    bb.ymax = 0;

    for (const t_trace* elem = head; elem != nullptr; elem = elem->next) {
        const t_rr_node& node = device_ctx.rr_graph.rr_nodes()[elem->index];
        //The router interprets RR nodes which cross the boundary as being
        //'within' of the BB. Only those which are *strictly* out side the
        //box are excluded, hence we use the nodes xhigh/yhigh for xmin/xmax,
        //and xlow/ylow for xmax/ymax calculations
        bb.xmin = std::min<int>(bb.xmin, rr_graph.node_xhigh(node.id()));
        bb.ymin = std::min<int>(bb.ymin, rr_graph.node_yhigh(node.id()));
        bb.xmax = std::max<int>(bb.xmax, rr_graph.node_xlow(node.id()));
        bb.ymax = std::max<int>(bb.ymax, rr_graph.node_ylow(node.id()));
    }

    VTR_ASSERT(bb.xmin <= bb.xmax);
    VTR_ASSERT(bb.ymin <= bb.ymax);

    return bb;
}

void enable_router_debug(
    const t_router_opts& router_opts,
    ClusterNetId net,
    int sink_rr,
    int router_iteration,
    ConnectionRouterInterface* router) {
    bool active_net_debug = (router_opts.router_debug_net >= -1);
    bool active_sink_debug = (router_opts.router_debug_sink_rr >= 0);
    bool active_iteration_debug = (router_opts.router_debug_iteration >= 0);

    bool match_net = (ClusterNetId(router_opts.router_debug_net) == net || router_opts.router_debug_net == -1);
    bool match_sink = (router_opts.router_debug_sink_rr == sink_rr || router_opts.router_debug_sink_rr < 0);
    bool match_iteration = (router_opts.router_debug_iteration == router_iteration || router_opts.router_debug_iteration < 0);

    f_router_debug = active_net_debug || active_sink_debug || active_iteration_debug;

    if (active_net_debug) f_router_debug &= match_net;
    if (active_sink_debug) f_router_debug &= match_sink;
    if (active_iteration_debug) f_router_debug &= match_iteration;
    
    router->set_router_debug(f_router_debug);

#ifndef VTR_ENABLE_DEBUG_LOGGING
    VTR_LOGV_WARN(f_router_debug, "Limited router debug output provided since compiled without VTR_ENABLE_DEBUG_LOGGING defined\n");
#endif
}

bool is_iteration_complete(bool routing_is_feasible, const t_router_opts& router_opts, int itry, std::shared_ptr<const SetupHoldTimingInfo> timing_info, bool rcv_finished) {
    //This function checks if a routing iteration has completed.
    //When VPR is run normally, we check if routing_budgets_algorithm is disabled, and if the routing is legal
    //With the introduction of yoyo budgeting algorithm, we must check if there are no hold violations
    //in addition to routing being legal and the correct budgeting algorithm being set.

    if (routing_is_feasible) {
        if (router_opts.routing_budgets_algorithm != YOYO) {
            return true;
        } else if (router_opts.routing_budgets_algorithm == YOYO && (timing_info->hold_worst_negative_slack() == 0 || rcv_finished) && itry != 1) {
            return true;
        }
    }
    return false;
}

bool should_setup_lower_bound_connection_delays(int itry, const t_router_opts& /*router_opts*/) {
    /* Checks to see if router should (re)calculate route budgets
     * It's currently set to only calculate after the first routing iteration */

    if (itry == 1) return true;
    return false;
}

static bool is_better_quality_routing(const vtr::vector<ClusterNetId, t_traceback>& best_routing,
                                      const RoutingMetrics& best_routing_metrics,
                                      const WirelengthInfo& wirelength_info,
                                      std::shared_ptr<const SetupHoldTimingInfo> timing_info) {
    if (best_routing.empty()) {
        return true; //First legal routing
    }

    //Rank first based on sWNS, followed by other timing metrics
    if (timing_info) {
        if (timing_info->setup_worst_negative_slack() > best_routing_metrics.sWNS) {
            return true;
        } else if (timing_info->setup_worst_negative_slack() < best_routing_metrics.sWNS) {
            return false;
        }

        if (timing_info->setup_total_negative_slack() > best_routing_metrics.sTNS) {
            return true;
        } else if (timing_info->setup_total_negative_slack() < best_routing_metrics.sTNS) {
            return false;
        }

        if (timing_info->hold_worst_negative_slack() > best_routing_metrics.hWNS) {
            return true;
        } else if (timing_info->hold_worst_negative_slack() > best_routing_metrics.hWNS) {
            return false;
        }

        if (timing_info->hold_total_negative_slack() > best_routing_metrics.hTNS) {
            return true;
        } else if (timing_info->hold_total_negative_slack() > best_routing_metrics.hTNS) {
            return false;
        }
    }

    //Finally, wirelength tie breaker
    return wirelength_info.used_wirelength() < best_routing_metrics.used_wirelength;
}

static bool early_reconvergence_exit_heuristic(const t_router_opts& router_opts,
                                               int itry_since_last_convergence,
                                               std::shared_ptr<const SetupHoldTimingInfo> timing_info,
                                               const RoutingMetrics& best_routing_metrics) {
    //Give-up on reconvergent routing if the CPD improvement after the
    //first iteration since convergence is small, compared to the best
    //CPD seen so far
    if (itry_since_last_convergence == 1) {
        float cpd_ratio = timing_info->setup_worst_negative_slack() / best_routing_metrics.sWNS;

        //Give up if we see less than a 1% CPD improvement,
        //after reducing pres_fac. Typically larger initial
        //improvements are needed to see an actual improvement
        //in final legal routing quality.
        if (cpd_ratio >= router_opts.reconvergence_cpd_threshold) {
            VTR_LOG("Giving up routing since additional routing convergences seem unlikely to improve quality (CPD ratio: %g)\n", cpd_ratio);
            return true; //Potential CPD improvement is small, don't spend run-time trying to improve it
        }
    }

    return false; //Don't give up
}

static void generate_route_timing_reports(const t_router_opts& router_opts,
                                          const t_analysis_opts& analysis_opts,
                                          const SetupTimingInfo& timing_info,
                                          const RoutingDelayCalculator& delay_calc) {
    auto& timing_ctx = g_vpr_ctx.timing();
    auto& atom_ctx = g_vpr_ctx.atom();

    VprTimingGraphResolver resolver(atom_ctx.nlist, atom_ctx.lookup, *timing_ctx.graph, delay_calc);
    resolver.set_detail_level(analysis_opts.timing_report_detail);

    tatum::TimingReporter timing_reporter(resolver, *timing_ctx.graph, *timing_ctx.constraints);

    timing_reporter.report_timing_setup(router_opts.first_iteration_timing_report_file, *timing_info.setup_analyzer(), analysis_opts.timing_report_npaths);
}

// If a route is ripped up during routing, non-configurable sets are left
// behind.  As a result, the final routing may have stubs at
// non-configurable sets.  This function tracks non-configurable set usage,
// and if the sets are unused, prunes them.
static void prune_unused_non_configurable_nets(CBRR& connections_inf) {
    auto& device_ctx = g_vpr_ctx.device();
    auto& cluster_ctx = g_vpr_ctx.clustering();
    auto& route_ctx = g_vpr_ctx.routing();

    std::vector<int> non_config_node_set_usage(device_ctx.rr_non_config_node_sets.size(), 0);
    for (auto net_id : cluster_ctx.clb_nlist.nets()) {
        connections_inf.prepare_routing_for_net(net_id);
        connections_inf.clear_force_reroute_for_net();

        std::fill(non_config_node_set_usage.begin(), non_config_node_set_usage.end(), 0);
        t_rt_node* rt_root = traceback_to_route_tree(net_id, &non_config_node_set_usage);
        if (rt_root == nullptr) {
            continue;
        }

        //Sanity check that route tree and traceback are equivalent before pruning
        VTR_ASSERT(verify_traceback_route_tree_equivalent(
            route_ctx.trace[net_id].head, rt_root));

        // check for edge correctness
        VTR_ASSERT_SAFE(is_valid_skeleton_tree(rt_root));

        //Prune the branches of the tree that don't legally lead to sinks
        rt_root = prune_route_tree(rt_root, connections_inf,
                                   &non_config_node_set_usage);

        // Free old traceback.
        free_traceback(net_id);

        // Update traceback with pruned tree.
        auto& reached_rt_sinks = connections_inf.get_reached_rt_sinks();
        traceback_from_route_tree(net_id, rt_root, reached_rt_sinks.size());
        VTR_ASSERT(verify_traceback_route_tree_equivalent(route_ctx.trace[net_id].head, rt_root));

        free_route_tree(rt_root);
    }
}

//Initializes net_delay based on best-case delay estimates from the router lookahead
static void init_net_delay_from_lookahead(const RouterLookahead& router_lookahead,
                                          ClbNetPinsMatrix<float>& net_delay, const float bias) {
    auto& cluster_ctx = g_vpr_ctx.clustering();
    auto& route_ctx = g_vpr_ctx.routing();

    t_conn_cost_params cost_params;
    cost_params.criticality = 1.; //Ensures lookahead returns delay value
    cost_params.bias = bias;
    //cost_params.offpath_penalty = router_opts.offpath_penalty;
    //cost_params.detailed_router = router_opts.detailed_router;

    for (auto net_id : cluster_ctx.clb_nlist.nets()) {
        if (cluster_ctx.clb_nlist.net_is_ignored(net_id)) continue;

        int source_rr = route_ctx.net_rr_terminals[net_id][0];

        for (size_t ipin = 1; ipin < cluster_ctx.clb_nlist.net_pins(net_id).size(); ++ipin) {
            int sink_rr = route_ctx.net_rr_terminals[net_id][ipin];

            float est_delay = router_lookahead.get_expected_cost(RRNodeId(source_rr), RRNodeId(sink_rr), cost_params, /*R_upstream=*/0.);
            VTR_ASSERT(std::isfinite(est_delay) && est_delay < std::numeric_limits<float>::max());

            net_delay[net_id][ipin] = est_delay;
        }
    }
}
//mycode==========================
bool try_timing_driven_route_incr_route(const t_router_opts& router_opts,
                             const t_file_name_opts& filename_opts,
                             const t_analysis_opts& analysis_opts,
                             const std::vector<t_segment_inf>& segment_inf,
                             ClbNetPinsMatrix<float>& net_delay,
                             const ClusteredPinAtomPinsLookup& netlist_pin_lookup,
                             std::shared_ptr<SetupHoldTimingInfo> timing_info,
                             std::shared_ptr<RoutingDelayCalculator> delay_calc,
                             ScreenUpdatePriority first_iteration_priority,
                             bool is_flat) {
    //printf("******** Entering try_timing_driven_route_incr_route **************\n");
    switch (router_opts.router_heap) {
        case e_heap_type::BINARY_HEAP:
            return try_timing_driven_route_tmpl_incr_route<ConnectionRouter<BinaryHeap>>(
                filename_opts,
                router_opts,
                analysis_opts,
                segment_inf,
                net_delay,
                netlist_pin_lookup,
                timing_info,
                delay_calc,
                first_iteration_priority,
                is_flat);
            break;
        case e_heap_type::BUCKET_HEAP_APPROXIMATION:
            return try_timing_driven_route_tmpl_incr_route<ConnectionRouter<Bucket>>(
                filename_opts,
                router_opts,
                analysis_opts,
                segment_inf,
                net_delay,
                netlist_pin_lookup,
                timing_info,
                delay_calc,
                first_iteration_priority,
                is_flat);
        default:
            VPR_FATAL_ERROR(VPR_ERROR_ROUTE, "Unknown heap type %d", router_opts.router_heap);
    }
}

std::pair<ClusterNetId, int> get_netid_and_sinkid(std::string connection_id) {
    /* configure connection id format in this function */
    // assume 'netid_sinkid' format
    int loc = 0;
    for (char c: connection_id) {
        loc++;
        if (c == '_') break;
    }
    if (loc == 1 or loc == connection_id.length()) {
        VTR_LOG("connection id '%s' does not conform to the format.", connection_id);
        return std::pair<ClusterNetId, int>(ClusterNetId(-1), int(-1));
    }
    return std::pair<ClusterNetId, int>(
        ClusterNetId(std::stoi(connection_id.substr(0, loc))),
        int(std::stoi(connection_id.substr(loc))));
};

int complete_rip_up_net_count;
int partial_rip_up_net_count;
template<typename ConnectionRouter>
bool try_timing_driven_route_tmpl_incr_route(const t_file_name_opts& filename_opts,
                                  const t_router_opts& router_opts,
                                  const t_analysis_opts& analysis_opts,
                                  const std::vector<t_segment_inf>& segment_inf,
                                  ClbNetPinsMatrix<float>& net_delay,
                                  const ClusteredPinAtomPinsLookup& netlist_pin_lookup,
                                  std::shared_ptr<SetupHoldTimingInfo> timing_info,
                                  std::shared_ptr<RoutingDelayCalculator> delay_calc,
                                  ScreenUpdatePriority first_iteration_priority,
                                  bool is_flat) {
    /* Timing-driven routing algorithm.  The timing graph (includes slack)   *
     * must have already been allocated, and net_delay must have been allocated. *
     * Returns true if the routing succeeds, false otherwise.                    */

    const auto& device_ctx = g_vpr_ctx.device();
    const auto& atom_ctx = g_vpr_ctx.atom();
    const auto& cluster_ctx = g_vpr_ctx.clustering();
    auto& route_ctx = g_vpr_ctx.mutable_routing();
    //Initially, the router runs normally trying to reduce congestion while
    //balancing other metrics (timing, wirelength, run-time etc.)
    RouterCongestionMode router_congestion_mode = RouterCongestionMode::NORMAL;

    //Initialize and properly size the lookups for profiling
    profiling::profiling_initialization(get_max_pins_per_net());

    //DEFAULT
    //sort so net with most sinks is routed first.
    auto sorted_nets = std::vector<ClusterNetId>(cluster_ctx.clb_nlist.nets().begin(), cluster_ctx.clb_nlist.nets().end());
    std::sort(sorted_nets.begin(), sorted_nets.end(), more_sinks_than());
    
    // // (PARSA) Luka, 2025: debug
    // auto sorted_nets = std::vector<ClusterNetId>(cluster_ctx.clb_nlist.nets().begin(), cluster_ctx.clb_nlist.nets().end());
    // std::sort(sorted_nets.begin(), sorted_nets.end(), less_sinks_than());
    /*
     * Configure the routing predictor
     */
    RoutingPredictor routing_predictor;
    float abort_iteration_threshold = std::numeric_limits<float>::infinity(); //Default no early abort
    if (router_opts.routing_failure_predictor == SAFE) {
        abort_iteration_threshold = ROUTING_PREDICTOR_ITERATION_ABORT_FACTOR_SAFE * router_opts.max_router_iterations;
    } else if (router_opts.routing_failure_predictor == AGGRESSIVE) {
        abort_iteration_threshold = ROUTING_PREDICTOR_ITERATION_ABORT_FACTOR_AGGRESSIVE * router_opts.max_router_iterations;
    } else {
        VTR_ASSERT_MSG(router_opts.routing_failure_predictor == OFF, "Unrecognized routing failure predictor setting");
    }

    float high_effort_congestion_mode_iteration_threshold = router_opts.congested_routing_iteration_threshold_frac * router_opts.max_router_iterations;

    /* Set delay of ignored signals to zero. Non-ignored net delays are set by
     * update_net_delays_from_route_tree() inside timing_driven_route_net(),
     * which is only called for non-ignored nets. */
    for (auto net_id : cluster_ctx.clb_nlist.nets()) {
        if (cluster_ctx.clb_nlist.net_is_ignored(net_id)) {
            for (unsigned int ipin = 1; ipin < cluster_ctx.clb_nlist.net_pins(net_id).size(); ++ipin) {
                net_delay[net_id][ipin] = 0.;
            }
        }
    }

    CBRR connections_inf{};

    route_budgets budgeting_inf;

    const auto* router_lookahead = get_cached_router_lookahead(
        router_opts.lookahead_type,
        router_opts.sbNode_lookahead_factor,
        router_opts.write_router_lookahead,
        router_opts.read_router_lookahead,
        segment_inf,
        is_flat);

    /*
     * Routing parameters
     */
    float pres_fac = update_pres_fac(router_opts.first_iter_pres_fac); /* Typically 0 -> ignore cong. */
    int bb_fac = router_opts.bb_factor;

    //When routing conflicts are detected the bounding boxes are scaled
    //by BB_SCALE_FACTOR every BB_SCALE_ITER_COUNT iterations
    constexpr float BB_SCALE_FACTOR = 2;
    constexpr int BB_SCALE_ITER_COUNT = 5;

    size_t available_wirelength = calculate_wirelength_available();

    /*
     * Routing status and metrics
     */
    bool routing_is_successful = false;
    WirelengthInfo wirelength_info;
    OveruseInfo overuse_info(device_ctx.rr_graph.num_nodes());
    tatum::TimingPathInfo critical_path;
    int itry; //Routing iteration number
    int itry_conflicted_mode = 0;

    /*
     * Best result so far
     */
    vtr::vector<ClusterNetId, t_traceback> best_routing;
    t_clb_opins_used best_clb_opins_used_locally;
    RoutingMetrics best_routing_metrics;
    int legal_convergence_count = 0;
    std::vector<int> scratch;

    ConnectionRouter router(
        device_ctx.grid,
        *router_lookahead,
        device_ctx.rr_graph.rr_nodes(),
        &device_ctx.rr_graph,
        device_ctx.rr_rc_data,
        device_ctx.rr_graph.rr_switch(),
        route_ctx.rr_node_route_inf,
        is_flat);
   
    //loading files needed for detailed router
    std::unordered_map<size_t, float> history_cost_map;
    std::unordered_map<size_t, float> iib_history_cost_map;
    std::unordered_map<size_t, size_t> node_id_map;
    std::unordered_map<size_t, std::set<size_t>> get_detailed_nodes;
    std::unordered_map<size_t, std::unordered_map<int, int>> net_id_to_sink_order_map;
    std::unordered_map<ClusterNetId, std::unordered_map<int, std::set<int>>> branch_node_map;
    std::set<size_t> nets_to_skip;
    std::vector<ClusterNetId> golden_net_order;
    std::set<size_t> congested_nets;
    if(router_opts.detailed_router == 0 && router_opts.nets_to_skip == 1) {
    	//reading file with nets to skip
    	std::ifstream nets_skip_fp;
    	std::string nets_skip_filename = "reconvergent_nets.txt";
    	nets_skip_fp.open(nets_skip_filename);
    	int lineno = 0;
    	if (!nets_skip_fp.is_open()) {
    	    vpr_throw(VPR_ERROR_ROUTE, get_arch_file_name(), lineno,
    	        "Cannot open reconvergent_nets.txt file");
    	}
    	int net_id;
    	while (nets_skip_fp >> net_id)
    	{
    	     nets_to_skip.insert(net_id);
    	}
    }
    std::string net_order_filename = "net_order_per_iteration.txt";
    std::ofstream net_order_file(net_order_filename);
    if(router_opts.detailed_router == 1) {
    	//reading file with net order
    	std::ifstream net_order_fp;
    	std::string golden_net_order_filename = "golden_net_order.txt";
    	net_order_fp.open(golden_net_order_filename);
    	int lineno = 0;
    	if (!net_order_fp.is_open()) {
    	    vpr_throw(VPR_ERROR_ROUTE, get_arch_file_name(), lineno,
    	        "Cannot open golden_net_order.txt file");
    	}
    	int net_id1;
    	while (net_order_fp >> net_id1)
    	{
    	     golden_net_order.push_back(ClusterNetId(net_id1));
    	}
    	net_order_fp.close();
    	//sorted_nets = golden_net_order;
    	//std::sort(sorted_nets.begin(), sorted_nets.end(), more_sinks_than());
    	//std::reverse(sorted_nets.begin(), sorted_nets.end());
        std::ifstream sink_order_fp;
        std::string sink_order_filename = "sink_order.txt";
        sink_order_fp.open(sink_order_filename);
        lineno = 0;
        VTR_LOG("[SHA] Reading sink_order.txt file\n");
        if (!sink_order_fp.is_open()) {
            vpr_throw(VPR_ERROR_ROUTE, get_arch_file_name(), lineno,
                "Cannot open sink order file");
        }

        std::string line;
        while (getline(sink_order_fp, line)) {
            if (line.empty()) {
                    // If the line is empty, skip processing
                    continue;
                }
                std::istringstream iss(line);
            size_t net_id;
            int net_pin_index;
            std::vector<int> sink_order;
            std::unordered_map<int, int> sink_order_index;
                
            iss >> net_id;  // First read the net ID
            //VTR_LOG("Net ID read: %d  ", node_id);
                while (iss >> net_pin_index) {  // Then read all the following net IDs
                    sink_order.push_back(net_pin_index);
                //VTR_LOG("%s ", net_id.c_str());
                }
                for (int i = 0; i < sink_order.size(); ++i) {
                sink_order_index[sink_order[i]] = i;
                }
            if (net_id == 48){
            VTR_LOG("%zu", net_id);
                    for (int i = 0; i < sink_order.size(); ++i) {
                VTR_LOG(" %d", sink_order[i]);
                    }
                VTR_LOG("\n");
            VTR_LOG("%zu", net_id);
                    for (int i = 0; i < sink_order_index.size(); ++i) {
                VTR_LOG(" %d", sink_order_index[i]);
                    }
                VTR_LOG("\n");
            }
            net_id_to_sink_order_map[net_id] = sink_order_index;
        }
        // reading branch node 
        // net, net_pin_index --> all branch dnode id
        std::ifstream branch_dnode_map_fp;
        std::string branch_dnode_map_filename = "branch_node_map.txt";
        branch_dnode_map_fp.open(branch_dnode_map_filename);
        lineno = 0;
        VTR_LOG("[SHA] Reading branch_node_map.txt file\n");
        if (!branch_dnode_map_fp.is_open()) {
            vpr_throw(VPR_ERROR_ROUTE, get_arch_file_name(), lineno,
                "Cannot open branch node map file");
        }

        while (getline(branch_dnode_map_fp, line)) {
            if (line.empty()) {	
                continue;	    
            }
            std::istringstream iss(line);
            int dnode_id;
            std::set<int> dnode_ids;

            ClusterNetId net_id;
            int sink_id;
            std::string connection_id;
            iss >> connection_id;  // First read the node ID
            
            std::tie(net_id, sink_id) = get_netid_and_sinkid(connection_id);
                
            while (iss >> dnode_id) {  // Then read all the following net IDs
                dnode_ids.insert(dnode_id);
            }
            
            branch_node_map[net_id][sink_id] = dnode_ids;
        }
        branch_dnode_map_fp.close();
    	//reading file with nets to skip
    	std::ifstream nets_skip_fp;
    	std::string nets_skip_filename = "reconvergent_nets.txt";
    	nets_skip_fp.open(nets_skip_filename);
    	lineno = 0;
    	if (!nets_skip_fp.is_open()) {
    	    vpr_throw(VPR_ERROR_ROUTE, get_arch_file_name(), lineno,
    	        "Cannot open reconvergent_nets.txt file");
    	}
    	int net_id;
    	while (nets_skip_fp >> net_id)
    	{
    	     nets_to_skip.insert(net_id);
    	}
    	nets_skip_fp.close();
    	//reading file with congested nets
    	std::ifstream congested_nets_fp;
    	std::string congested_nets_filename = "congested_nets.txt";
    	congested_nets_fp.open(congested_nets_filename);
    	lineno = 0;
    	if (!congested_nets_fp.is_open()) {
    	    vpr_throw(VPR_ERROR_ROUTE, get_arch_file_name(), lineno,
    	        "Cannot open congested_nets.txt file");
    	}
    	while (congested_nets_fp >> net_id)
    	{
    	     congested_nets.insert(net_id);
    	}
    	congested_nets_fp.close();
        //reading node map file
        std::ifstream node_map_fp;
        std::string node_map_filename = "../gr_dr_map.map";
        node_map_fp.open(node_map_filename);
        lineno = 0;
	if (!node_map_fp.is_open()) {
            vpr_throw(VPR_ERROR_ROUTE, get_arch_file_name(), lineno,
                "Cannot open node map file");
        }
        while (getline(node_map_fp, line)) {
            std::istringstream iss(line);
            size_t gr_node_id;
            size_t dr_node_id;
            iss >> gr_node_id;  // First read the node ID
            while (iss >> dr_node_id) {  // Then read all the following net IDs
                node_id_map[dr_node_id] = gr_node_id;
		get_detailed_nodes[gr_node_id].insert(dr_node_id);
            }
        }
        node_map_fp.close();
	// setting the initial cost of nodes equal to occupancy
        for (const RRNodeId& rr_id : device_ctx.rr_graph.nodes()) {
	    int g_occupancy = device_ctx.rr_graph.get_global_occupancy(rr_id); 
            if (g_occupancy > 0) {
		//route_ctx.rr_node_route_inf[size_t(rr_id)].acc_cost = g_occupancy;
            	route_ctx.rr_node_route_inf[size_t(rr_id)].set_g_occ(g_occupancy);
		//VTR_LOG("Loading occupancy as hist of node: %zu gocc: %d\n", size_t(rr_id), g_occupancy);
	    }
	    else {
		//route_ctx.rr_node_route_inf[size_t(rr_id)].acc_cost = 1.;
            	route_ctx.rr_node_route_inf[size_t(rr_id)].set_g_occ(1);
	    
	    }
	}
    }
    if(router_opts.detailed_router == 1 && router_opts.load_gr_history == 1){
        printf("####### BEGINING loading history files for detailed router: %d ######################\n", itry);
        
        //std::string suffix = "_to_legalise.route";
        //std::string prefix = filename_opts.RouteFile.substr(0,filename_opts.RouteFile.size()-6);
        //std::string route_file_to_legalise =  prefix + suffix;
        //printf("******* route_file_to_legalise after concatenation %s\n", route_file_to_legalise);

        //read_route_incr_route(temp_net_id, route_file_to_legalise.c_str(), router_opts, filename_opts.verify_file_digests);
        //VTR_LOG("******* Successfully loaded partial route file\n");
        
        //reading hist file
        std::ifstream hist_fp;
        std::string hist_filename = "history_cost_file_gr.txt";
        hist_fp.open(hist_filename);
        int lineno = 0;
        if (!hist_fp.is_open()) {
            vpr_throw(VPR_ERROR_ROUTE, get_arch_file_name(), lineno,
                "Cannot open history cost file");
        }
        int node_id, history_congestion_cost;
        while (hist_fp >> node_id >> history_congestion_cost)
        {
            history_cost_map[node_id] = float(history_congestion_cost);
        }
        hist_fp.close();
        //======================   
        //reading node map file
        std::ifstream node_map_fp;
        std::string node_map_filename = "../gr_dr_map.map";
        node_map_fp.open(node_map_filename);
        if (!node_map_fp.is_open()) {
            vpr_throw(VPR_ERROR_ROUTE, get_arch_file_name(), lineno,
                "Cannot open node map file");
        }
        std::string line;
        while (getline(node_map_fp, line)) {
            std::istringstream iss(line);
            size_t gr_node_id;
            size_t dr_node_id;
            iss >> gr_node_id;  // First read the node ID
            while (iss >> dr_node_id) {  // Then read all the following net IDs
                node_id_map[dr_node_id] = gr_node_id;
            }
        }
        node_map_fp.close();
        for (const RRNodeId& rr_id : device_ctx.rr_graph.nodes()) {
            size_t dr_node_id = (size_t)rr_id;
            route_ctx.rr_node_route_inf[dr_node_id].acc_cost = history_cost_map[dr_node_id];//1 + (history_cost_map[node_id_from]-1) * 0.4;
             
            //size_t gr_node_id = node_id_map[dr_node_id];
            //if (history_cost_map.find(gr_node_id) != history_cost_map.end()){//GI
            //	route_ctx.rr_node_route_inf[dr_node_id].acc_cost = history_cost_map[gr_node_id];//1 + (history_cost_map[node_id_from]-1) * 0.4;
            //}
         } 
    }

    // Make sure template type ConnectionRouter is a ConnectionRouterInterface.
    static_assert(std::is_base_of<ConnectionRouterInterface, ConnectionRouter>::value, "ConnectionRouter must implement the ConnectionRouterInterface");

    /*
     * On the first routing iteration ignore congestion to get reasonable net
     * delay estimates. Set criticalities to 1 when timing analysis is on to
     * optimize timing, and to 0 when timing analysis is off to optimize routability.
     *
     * Subsequent iterations use the net delays from the previous iteration.
     */
    std::shared_ptr<SetupHoldTimingInfo> route_timing_info;
    {
        vtr::ScopedStartFinishTimer init_timing_timer("Initializing router criticalities");
        if (timing_info) {
            if (router_opts.initial_timing == e_router_initial_timing::ALL_CRITICAL) {
                //First routing iteration, make all nets critical for a min-delay routing
                route_timing_info = make_constant_timing_info(1.);
            } else {
                VTR_ASSERT(router_opts.initial_timing == e_router_initial_timing::LOOKAHEAD);

                {
                    //Estimate initial connection delays from the router lookahead
                    init_net_delay_from_lookahead(*router_lookahead, net_delay, router_opts.sbNode_lookahead_factor);

                    //Run STA to get estimated criticalities
                    timing_info->update();
                }
                route_timing_info = timing_info;
            }
        } else {
            //Not timing driven, force criticality to zero for a routability-driven routing
            route_timing_info = make_constant_timing_info(0.);
        }
        VTR_LOG("Initial Net Connection Criticality Histogram:\n");
        print_router_criticality_histogram(*route_timing_info, netlist_pin_lookup);
    }

    std::unique_ptr<ClusteredPinTimingInvalidator> pin_timing_invalidator;
    if (timing_info) {
        pin_timing_invalidator = std::make_unique<ClusteredPinTimingInvalidator>(cluster_ctx.clb_nlist,
                                                                                 netlist_pin_lookup,
                                                                                 atom_ctx.nlist,
                                                                                 atom_ctx.lookup,
                                                                                 *timing_info->timing_graph());
    }

    RouterStats router_stats;
    timing_driven_route_structs route_structs;
    float prev_iter_cumm_time = 0;
    vtr::Timer iteration_timer;
    int num_net_bounding_boxes_updated = 0;
    int itry_since_last_convergence = -1;

    // This heap is used for reserve_locally_used_opins.
    BinaryHeap small_heap;
    small_heap.init_heap(device_ctx.grid);

    // When RCV is enabled the router will not stop unless negative hold slack is 0
    // In some cases this isn't doable, due to global nets or intracluster routing issues
    // In these cases RCV will finish early if it goes RCV_FINISH_EARLY_COUNTDOWN iterations without detecting resolvable negative hold slack
    // Increasing this will make the router fail occasionally, decreasing will sometimes not let all hold violations be resolved
    constexpr int RCV_FINISH_EARLY_COUNTDOWN = 15;

    int rcv_finished_count = RCV_FINISH_EARLY_COUNTDOWN;
    //mycode
    float final_pres_fac;
    float pres_fac_new;
    size_t connections_routed_first_iteration, nets_routed_first_iteration, heap_pushes_first_iteration, heap_pops_first_iteration; 
    print_route_status_header();
    for (itry = 1; itry <= router_opts.max_router_iterations; ++itry) {
    	RouterStats router_iteration_stats;
        std::vector<ClusterNetId> rerouted_nets;

        /* Reset "is_routed" and "is_fixed" flags to indicate nets not pre-routed (yet) */
        for (auto net_id : cluster_ctx.clb_nlist.nets()) {
            route_ctx.net_status.set_is_routed(net_id, false);
            route_ctx.net_status.set_is_fixed(net_id, false);
        }

        if (itry_since_last_convergence >= 0) {
            ++itry_since_last_convergence;
        }

        // Calculate this once and pass it into net routing to check if should reroute for hold
        float worst_negative_slack = 0;
        if (budgeting_inf.if_set()) {
            worst_negative_slack = timing_info->hold_total_negative_slack();
        }

        /*
         * Route each net
         */
        ClusterNetId temp_net_id;
        complete_rip_up_net_count = 0;
        partial_rip_up_net_count = 0;
	//std::vector<ClusterNetId> first_100_nets;
	//size_t num_elements_to_copy = 100;
	//first_100_nets.reserve(num_elements_to_copy);
	//std::copy_n(sorted_nets.begin(), num_elements_to_copy, std::back_inserter(first_100_nets));
	if(router_opts.shuffle_net_order == 1){
	  //std::random_device rd;  // Seed for random number generator
          //std::mt19937 g(rd());   // Standard Mersenne Twister engine
          //std::shuffle(begin(sorted_nets), end(sorted_nets), g);
          std::sort(sorted_nets.begin(), sorted_nets.end(), less_sinks_than());
	}
	
	for (auto net_id : sorted_nets) {
	    if (nets_to_skip.find(size_t(net_id)) != nets_to_skip.end()){
	    	continue;
	    }
	    /*if (congested_nets.find(size_t(net_id)) == congested_nets.end()){
	        continue;
	    }*/

            temp_net_id = net_id;
	    std::unordered_map<int, int> sink_order_index;
	    if (router_opts.detailed_router == 1){
	    	sink_order_index = net_id_to_sink_order_map[size_t(net_id)];
	    }

            bool was_rerouted = false;
            bool is_routable;
	    is_routable = try_timing_driven_route_net_incr_route(filename_opts,
                                                           router,
                                                           net_id,
                                                           itry,
                                                           pres_fac,
                                                           sink_order_index,
                                                           branch_node_map,
                                                           router_opts,
                                                           connections_inf,
                                                           router_iteration_stats,
                                                           route_structs.pin_criticality,
                                                           route_structs.rt_node_of_sink,
                                                           net_delay,
                                                           netlist_pin_lookup,
                                                           route_timing_info,
                                                           pin_timing_invalidator.get(),
                                                           budgeting_inf,
                                                           was_rerouted,
                                                           worst_negative_slack,
                                                           routing_predictor,
                                                           is_flat,
							   net_order_file);
            if (!is_routable) {
                return (false); //Impossible to route
            }

            if (was_rerouted) {
                rerouted_nets.push_back(net_id);
#ifndef NO_GRAPHICS
                update_router_info_and_check_bp(BP_NET_ID, size_t(net_id));
#endif
            }
        }

        /*
            (PARSA) Luka, 2025: This has no use in the RSMT constrained setting, so && !steiner_constraints was added.
            A similar thing should be done for the rest of the file loading that the old constrained setting router does.
        */ 
        // For locking and loading branch changing nets
        if (router_opts.detailed_router == 1 && itry == 1 && !router_opts.steiner_constraints) {
            std::string suffix = "_ssr_to_load.route";
            std::string prefix = filename_opts.RouteFile.substr(0,filename_opts.RouteFile.size()-6);
            std::string route_file_to_legalise =  prefix + suffix;
            std::cout << prefix << " " << suffix << " " << route_file_to_legalise << std::endl;
            VTR_LOG("******* route_file_to load after concatenation %s\n", route_file_to_legalise);
            read_route_incr_route(temp_net_id, route_file_to_legalise.c_str(), router_opts, filename_opts.verify_file_digests);
            VTR_LOG("******* Successfully loaded route file\n");
	    }
	//}
        //std::unordered_map<size_t, float> history_cost_map;
        //std::unordered_map<size_t, float> iib_history_cost_map;
        //std::unordered_map<size_t, size_t> node_id_map;
        if(router_opts.incr_route == 1 && router_opts.icr_iter >= 1 && itry == 1){
            std::string suffix = "_to_legalise.route";
            std::string prefix = filename_opts.RouteFile.substr(0,filename_opts.RouteFile.size()-6);
            std::string route_file_to_legalise =  prefix + suffix;
            printf("******* route_file_to_legalise after concatenation %s\n", route_file_to_legalise);

            read_route_incr_route(temp_net_id, route_file_to_legalise.c_str(), router_opts, filename_opts.verify_file_digests);
            VTR_LOG("******* Successfully loaded partial route file\n");
            
            //reading hist file
            std::ifstream hist_fp;
            std::string prev_icr_iter = std::to_string(router_opts.icr_iter - 1);
            std::string hist_filename = "history_cost_file_"+prev_icr_iter+".txt";
            hist_fp.open(hist_filename);
            int lineno = 0;
            if (!hist_fp.is_open()) {
                vpr_throw(VPR_ERROR_ROUTE, get_arch_file_name(), lineno,
                    "Cannot open history cost file");
            }
            int node_id, history_congestion_cost;
            while (hist_fp >> node_id >> history_congestion_cost)
            {
                history_cost_map[node_id] = history_congestion_cost;
            }
            hist_fp.close();
            //======================
            //reading presfac file
            std::ifstream pres_fac_fp;
            std::string pres_fac_filename = "pres_fac_file_"+prev_icr_iter+".txt";
            pres_fac_fp.open(pres_fac_filename);
            if (!pres_fac_fp.is_open()) {
                vpr_throw(VPR_ERROR_ROUTE, get_arch_file_name(), lineno,
                    "Cannot open pres fac file");
            }
            float from_prev_icr_iter_pres_fac;
            while (pres_fac_fp >> from_prev_icr_iter_pres_fac)
            {
                pres_fac_new = from_prev_icr_iter_pres_fac;
                printf("New loaded pres fac from prev iteration is: %f", pres_fac_new);
            }
            pres_fac_fp.close();            
            //======================   
            //reading node map file
            std::ifstream node_map_fp;
            std::string node_map_filename = "node_map_file_"+prev_icr_iter+"_to_"+std::to_string(router_opts.icr_iter)+".txt";
            node_map_fp.open(node_map_filename);
            if (!node_map_fp.is_open()) {
                vpr_throw(VPR_ERROR_ROUTE, get_arch_file_name(), lineno,
                    "Cannot open node map file");
            }
            size_t node_id_from, node_id_to;
            while (node_map_fp >> node_id_from >> node_id_to)
            {
                node_id_map[node_id_to] = node_id_from;
            }
            node_map_fp.close();
            //======================
            //reading hist file
            std::ifstream hist_iib_fp;
            std::string iib_hist_filename = "iib_history_cost_file_"+prev_icr_iter+".txt";
            hist_iib_fp.open(iib_hist_filename);
            if (!hist_iib_fp.is_open()) {
                vpr_throw(VPR_ERROR_ROUTE, get_arch_file_name(), lineno,
                    "Cannot open history cost file");
            }
            int node_id_, history_congestion_cost_;
            while (hist_iib_fp >> node_id_ >> history_congestion_cost_)
            {
                iib_history_cost_map[node_id_] = history_congestion_cost_;
            }
            hist_iib_fp.close();
        }
          
        // Make sure any CLB OPINs used up by subblocks being hooked directly to them are reserved for that purpose
        bool rip_up_local_opins = (itry == 1 ? false : true);
        reserve_locally_used_opins(&small_heap, pres_fac,
                                   router_opts.acc_fac, router_opts.global_occ_factor, rip_up_local_opins);
        if (itry == 1) {
	    recompute_occupancy_from_scratch();            
	}

        if(router_opts.incr_route == 1 && router_opts.icr_iter >= 1 && itry == 1){
            recompute_occupancy_from_scratch();            
            pathfinder_update_acc_cost_and_overuse_info(0., overuse_info); 
            reserve_locally_used_opins(&small_heap, pres_fac,
                                   router_opts.acc_fac, router_opts.global_occ_factor, true);
            recompute_occupancy_from_scratch();            
            //loading history costs
            for (const RRNodeId& rr_id : device_ctx.rr_graph.nodes()) {
                size_t node_id_to = (size_t)rr_id;
		auto iter = node_id_map.find(node_id_to);
		if(iter != node_id_map.end()){
                	size_t node_id_from = node_id_map[node_id_to];
			if (history_cost_map.find(node_id_from) != history_cost_map.end()){//GI
                    		route_ctx.rr_node_route_inf[node_id_to].acc_cost = history_cost_map[node_id_from];//1 + (history_cost_map[node_id_from]-1) * 0.4;
                	}
			else{
				VTR_LOG("SHOULD NOT BE HERE\n");
			}
			
		}
		else if (iib_history_cost_map.find(node_id_to) != iib_history_cost_map.end()){//IIB
		    if (route_ctx.rr_node_route_inf[node_id_to].occ() == device_ctx.rr_graph.node_capacity(RRNodeId(node_id_to))){//legal iib
                    	route_ctx.rr_node_route_inf[node_id_to].acc_cost = iib_history_cost_map[node_id_to];// * 0.06 + 1.0;
		    }
		    else {//illegal iibs
                    	route_ctx.rr_node_route_inf[node_id_to].acc_cost = 1.0;//std::floor(iib_history_cost_map[node_id_to]*0.1) + 1.0;
		    }
                }
                else{// will enter here but will not affect routing
		    route_ctx.rr_node_route_inf[node_id_to].acc_cost = 1.0;
                }
            }
        }
        else if (router_opts.incr_route == 1 && itry > 1){
            pathfinder_update_acc_cost_and_overuse_info(router_opts.acc_fac, overuse_info);
        }
        
        /*
         * Calculate metrics for the current routing
         */
        bool routing_is_feasible = feasible_routing();
        float est_success_iteration = routing_predictor.estimate_success_iteration();
        printf("Is routing feasible: %d\n", routing_is_feasible);
        //Update resource costs and overuse info
        if (router_opts.incr_route == 0 && itry == 1) {
	    if (router_opts.first_iter_pres_fac != 0) {
            	pathfinder_update_acc_cost_and_overuse_info(router_opts.acc_fac, overuse_info); /* Acc_fac=0 for first iter. */
	    }
	    else {
            	pathfinder_update_acc_cost_and_overuse_info(0., overuse_info); /* Acc_fac=0 for first iter. */
	    }

        } else if (router_opts.incr_route == 0 && itry > 1){
            	pathfinder_update_acc_cost_and_overuse_info(router_opts.acc_fac, overuse_info);
        }

        wirelength_info = calculate_wirelength_info(available_wirelength);
        routing_predictor.add_iteration_overuse(itry, overuse_info.overused_nodes);

        if (timing_info) {
            //Update timing based on the new routing
            //Note that the net delays have already been updated by timing_driven_route_net
            timing_info->update();
            timing_info->set_warn_unconstrained(false); //Don't warn again about unconstrained nodes again during routing
            pin_timing_invalidator->reset();

            //Use the real timing analysis criticalities for subsequent routing iterations
            //  'route_timing_info' is what is actually passed into the net/connection routers,
            //  and for the 1st iteration may not be the actual STA results (e.g. all criticalities set to 1)
            route_timing_info = timing_info;

            critical_path = timing_info->least_slack_critical_path();

            VTR_ASSERT_SAFE(timing_driven_check_net_delays(net_delay));

            if (itry == 1) {
                generate_route_timing_reports(router_opts, analysis_opts, *timing_info, *delay_calc);
            }
        }

        float iter_cumm_time = iteration_timer.elapsed_sec();
        float iter_elapsed_time = iter_cumm_time - prev_iter_cumm_time;

        //Output progress
        print_route_status(itry, iter_elapsed_time, pres_fac, num_net_bounding_boxes_updated, router_iteration_stats, overuse_info, wirelength_info, timing_info, est_success_iteration);
        VTR_LOG("NETS COMPLETELY RIPPED UP: %d\n", complete_rip_up_net_count);
        VTR_LOG("NETS PARTIALLY RIPPED UP: %d\n", partial_rip_up_net_count);
	if (timing_info) {
            VTR_LOG("Critical path: %g ns  sWNS: %e sTNS: %e\n", 1e9 * critical_path.delay(), timing_info->setup_worst_negative_slack(), timing_info->setup_total_negative_slack());
        }
        prev_iter_cumm_time = iter_cumm_time;

        //Update graphics
        if (itry == 1) {
            update_screen(first_iteration_priority, "Routing...", ROUTING, timing_info);
        } else {
            update_screen(ScreenUpdatePriority::MINOR, "Routing...", ROUTING, timing_info);
        }

        if (router_opts.save_routing_per_iteration) {
            std::string filename = vtr::string_fmt("iteration_%03d.route", itry);
            print_route(nullptr, filename.c_str());
        }
	if (router_opts.save_history_cost_per_iteration){
            std::string hist_filename = vtr::string_fmt("history_cost_iteration_%03d.hcost", itry);
            std::ofstream hist_file(hist_filename);
            for (const RRNodeId& rr_id : device_ctx.rr_graph.nodes()) {
                auto& node_inf = route_ctx.rr_node_route_inf[(size_t)rr_id];
                float history_cost = node_inf.acc_cost;
                hist_file << (size_t)rr_id << " " << history_cost << "\n";
            }
            hist_file.close();
	}
        if (itry == 1){
            connections_routed_first_iteration = router_iteration_stats.connections_routed;
            nets_routed_first_iteration = router_iteration_stats.nets_routed;
            heap_pushes_first_iteration = router_iteration_stats.heap_pushes;
            heap_pops_first_iteration = router_iteration_stats.heap_pops;
        }

        //Update router stats (total)
        router_stats.connections_routed += router_iteration_stats.connections_routed;
        router_stats.nets_routed += router_iteration_stats.nets_routed;
        router_stats.heap_pushes += router_iteration_stats.heap_pushes;
        router_stats.heap_pops += router_iteration_stats.heap_pops;
        
	router_stats.wire_heap_pushes += router_iteration_stats.wire_heap_pushes;
	router_stats.wire_heap_pops += router_iteration_stats.wire_heap_pops;

	router_stats.min_heap_pushes += router_iteration_stats.min_heap_pushes;
	router_stats.min_heap_pops += router_iteration_stats.min_heap_pops;

	router_stats.max_heap_pushes += router_iteration_stats.max_heap_pushes;
	router_stats.max_heap_pops += router_iteration_stats.max_heap_pops;
        /*
         * Are we finished?
         */
	// For locking and loading branch changing nets
        if (itry > 1 && is_iteration_complete(routing_is_feasible, router_opts, itry, timing_info, rcv_finished_count == 0)) {
        //if (is_iteration_complete(routing_is_feasible, router_opts, itry, timing_info, rcv_finished_count == 0)) {
            auto& router_ctx = g_vpr_ctx.routing();
            if (is_better_quality_routing(best_routing, best_routing_metrics, wirelength_info, timing_info)) {
                //Save routing
                best_routing = router_ctx.trace;
                best_clb_opins_used_locally = router_ctx.clb_opins_used_locally;

                routing_is_successful = true;

                //Update best metrics
                if (timing_info) {
                    timing_driven_check_net_delays(net_delay);

                    best_routing_metrics.sTNS = timing_info->setup_total_negative_slack();
                    best_routing_metrics.sWNS = timing_info->setup_worst_negative_slack();
                    best_routing_metrics.hTNS = timing_info->hold_total_negative_slack();
                    best_routing_metrics.hWNS = timing_info->hold_worst_negative_slack();
                    best_routing_metrics.critical_path = critical_path;
                }
                best_routing_metrics.used_wirelength = wirelength_info.used_wirelength();
            }

            //Decrease pres_fac so that critical connections will take more direct routes
            //Note that we use first_iter_pres_fac here (typically zero), and switch to
            //use initial_pres_fac on the next iteration.
            pres_fac = update_pres_fac(router_opts.first_iter_pres_fac);

            //Reduce timing tolerances to re-route more delay-suboptimal signals
            connections_inf.set_connection_criticality_tolerance(0.7);
            connections_inf.set_connection_delay_tolerance(1.01);

            ++legal_convergence_count;
            itry_since_last_convergence = 0;

            VTR_ASSERT(routing_is_successful);
        }

        if (itry_since_last_convergence == 1) {
            //We used first_iter_pres_fac when we started routing again
            //after the first routing convergence. Since that is often zero,
            //we want to set pres_fac to a reasonable (i.e. typically non-zero)
            //value afterwards -- so it grows when multiplied by pres_fac_mult
            pres_fac = update_pres_fac(router_opts.initial_pres_fac);
        }

        //Have we converged the maximum number of times, did not make any changes, or does it seem
        //unlikely additional convergences will improve QoR?
        if (legal_convergence_count >= router_opts.max_convergence_count
            || router_iteration_stats.connections_routed == 0
            || early_reconvergence_exit_heuristic(router_opts, itry_since_last_convergence, timing_info, best_routing_metrics)) {
#ifndef NO_GRAPHICS
            update_router_info_and_check_bp(BP_ROUTE_ITER, -1);
#endif
            break; //Done routing
        }

        /*
         * Abort checks: Should we give-up because this routing problem is unlikely to converge to a legal routing?
         */
        if (itry == 1 && early_exit_heuristic(router_opts, wirelength_info)) {
#ifndef NO_GRAPHICS
            update_router_info_and_check_bp(BP_ROUTE_ITER, -1);
#endif
            //Abort
            break;
        }

        //Estimate at what iteration we will converge to a legal routing
        if (overuse_info.overused_nodes > ROUTING_PREDICTOR_MIN_ABSOLUTE_OVERUSE_THRESHOLD) {
            //Only consider aborting if we have a significant number of overused resources

            if (!std::isnan(est_success_iteration) && est_success_iteration > abort_iteration_threshold && router_opts.routing_budgets_algorithm != YOYO) {
                VTR_LOG("Routing aborted, the predicted iteration for a successful route (%.1f) is too high.\n", est_success_iteration);
#ifndef NO_GRAPHICS
                update_router_info_and_check_bp(BP_ROUTE_ITER, -1);
#endif
                break; //Abort
            }
        }

        if (itry == 1 && router_opts.exit_after_first_routing_iteration) {
            VTR_LOG("Exiting after first routing iteration as requested\n");
#ifndef NO_GRAPHICS
            update_router_info_and_check_bp(BP_ROUTE_ITER, -1);
#endif
            break;
        }

        /*
         * Prepare for the next iteration
         */

        if (router_opts.route_bb_update == e_route_bb_update::DYNAMIC) {
            num_net_bounding_boxes_updated = dynamic_update_bounding_boxes(rerouted_nets, router_opts.high_fanout_threshold);
        }

        if (itry >= high_effort_congestion_mode_iteration_threshold) {
            //We are approaching the maximum number of routing iterations,
            //and still do not have a legal routing. Switch to a mode which
            //focuses more on attempting to resolve routing conflicts.
            router_congestion_mode = RouterCongestionMode::CONFLICTED;
        }

        //Update pres_fac
        if (itry == 1) {
            pres_fac = update_pres_fac(router_opts.initial_pres_fac);
        }
        else {
            pres_fac *= router_opts.pres_fac_mult;

            /* Avoid overflow for high iteration counts, even if acc_cost is big */
            //pres_fac = update_pres_fac(std::min(pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5)));
            pres_fac = update_pres_fac(std::min(pres_fac, static_cast<float>(1000)));
            VTR_LOG("pres fac: %f\n", pres_fac);
            // Increase short path criticality if it's having a hard time resolving hold violations due to congestion
            if (budgeting_inf.if_set()) {
                bool rcv_finished = false;

                /* This constant represents how much extra delay the budget increaser adds to the minimum and maximum delay budgets
                 * Experimentally this value delivers fast hold slack resolution, while not overwhelming the router 
                 * Increasing this will make it resolve hold faster, but could result in lower circuit quality */
                constexpr float budget_increase_factor = 300e-12;

                if (itry > 5 && worst_negative_slack != 0) rcv_finished = budgeting_inf.increase_min_budgets_if_struggling(budget_increase_factor, timing_info, worst_negative_slack, netlist_pin_lookup);
                if (rcv_finished)
                    rcv_finished_count--;
                else
                    rcv_finished_count = RCV_FINISH_EARLY_COUNTDOWN;
            }
        }

        if (router_congestion_mode == RouterCongestionMode::CONFLICTED) {
            //The design appears to have routing conflicts which are difficult to resolve:
            //  1) Don't re-route legal connections due to delay. This allows
            //     the router to focus on the actual conflicts
            //  2) Increase the net bounding boxes. This potentially allows
            //     the router to route around otherwise congested regions
            //     (at the cost of high run-time).

            //Increase the size of the net bounding boxes to give the router more
            //freedom to find alternate paths.
            //
            //In the case of routing conflicts there are multiple connections competing
            //for the same resources which can not resolve the congestion themselves.
            //In normal routing mode we try to keep the bounding boxes small to minimize
            //run-time, but this can limits how far signals can detour (i.e. they can't
            //route outside the bounding box), which can cause conflicts to oscillate back
            //and forth without resolving.
            //
            //By scaling the bounding boxes here, we slowly increase the router's search
            //space in hopes of it allowing signals to move further out of the way to
            //alleviate the conflicts.
            if (itry_conflicted_mode % BB_SCALE_ITER_COUNT == 0) {
                //We scale the bounding boxes by BB_SCALE_FACTOR,
                //every BB_SCALE_ITER_COUNT iterations. This ensures
                //that we give the router some time (BB_SCALE_ITER_COUNT) to try
                //resolve/negotiate congestion at the new BB factor.
                //
                //Note that we increase the BB factor slowly to try and minimize
                //the bounding box size (since larger bounding boxes slow the router down).
                auto& grid = g_vpr_ctx.device().grid;
                int max_grid_dim = std::max(grid.width(), grid.height());

                //Scale by BB_SCALE_FACTOR but clip to grid size to avoid overflow
                bb_fac = std::min<int>(max_grid_dim, bb_fac * BB_SCALE_FACTOR);

                route_ctx.route_bb = load_route_bb(bb_fac);
            }

            ++itry_conflicted_mode;
        }

        if (timing_info) {
            if (should_setup_lower_bound_connection_delays(itry, router_opts)) {
                // first iteration sets up the lower bound connection delays since only timing is optimized for
                connections_inf.set_stable_critical_path_delay(critical_path.delay());
                connections_inf.set_lower_bound_connection_delays(net_delay);

                //load budgets using information from uncongested delay information
                budgeting_inf.load_route_budgets(net_delay, timing_info, netlist_pin_lookup, router_opts);
                /*for debugging purposes*/
                // if (budgeting_inf.if_set()) {
                //     budgeting_inf.print_route_budget(std::string("route_budgets_") + std::to_string(itry) + ".txt", net_delay);
                // }

                if (router_opts.routing_budgets_algorithm == YOYO) router.set_rcv_enabled(true);

            } else {
                bool stable_routing_configuration = true;

                /*
                 * Determine if any connection need to be forcibly re-routed due to timing
                 */

                //Yes, if explicitly enabled
                bool should_ripup_for_delay = (router_opts.incr_reroute_delay_ripup == e_incr_reroute_delay_ripup::ON);

                //Or, if things are not too congested
                should_ripup_for_delay |= (router_opts.incr_reroute_delay_ripup == e_incr_reroute_delay_ripup::AUTO
                                           && router_congestion_mode == RouterCongestionMode::NORMAL);

                if (should_ripup_for_delay) {
                    if (connections_inf.critical_path_delay_grew_significantly(critical_path.delay())) {
                        // only need to forcibly reroute if critical path grew significantly
                        stable_routing_configuration = connections_inf.forcibly_reroute_connections(router_opts.max_criticality,
                                                                                                    timing_info,
                                                                                                    netlist_pin_lookup,
                                                                                                    net_delay);
                    }
                }

                // not stable if any connection needs to be forcibly rerouted
                if (stable_routing_configuration) {
                    connections_inf.set_stable_critical_path_delay(critical_path.delay());
                }
            }
        } else {
            /* If timing analysis is not enabled, make sure that the criticalities and the
             * net_delays stay as 0 so that wirelength can be optimized. */

            for (auto net_id : cluster_ctx.clb_nlist.nets()) {
                for (unsigned int ipin = 1; ipin < cluster_ctx.clb_nlist.net_pins(net_id).size(); ++ipin) {
                    net_delay[net_id][ipin] = 0.;
                }
            }
        }

        if (router_opts.congestion_analysis) profiling::congestion_analysis();
        if (router_opts.fanout_analysis) profiling::time_on_fanout_analysis();
        // profiling::time_on_criticality_analysis();
	
    	if (router_opts.detailed_router == 0) {
            net_order_file << "\n";
    	}
    }
    
    net_order_file.close();
    
    if (routing_is_successful) {
        VTR_LOG("Restoring best routing\n");

        auto& router_ctx = g_vpr_ctx.mutable_routing();

        /* Restore congestion from best route */
        for (auto net_id : cluster_ctx.clb_nlist.nets()) {
            pathfinder_update_path_occupancy(route_ctx.trace[net_id].head, -1);
            pathfinder_update_path_occupancy(best_routing[net_id].head, 1);
        }
        router_ctx.trace = best_routing;
        router_ctx.clb_opins_used_locally = best_clb_opins_used_locally;

        prune_unused_non_configurable_nets(connections_inf);
        //mycode
        std::string write_hist_filename = "history_cost_file_"+std::to_string(router_opts.icr_iter)+".txt";
        std::ofstream hist_file(write_hist_filename);
        for (const RRNodeId& rr_id : device_ctx.rr_graph.nodes()) {
            auto& node_inf = route_ctx.rr_node_route_inf[(size_t)rr_id];
            float history_cost = node_inf.acc_cost;
            //printf("%d %f\n", rr_id, history_cost);
            hist_file << (size_t)rr_id << " " << history_cost << "\n";
        }
        hist_file.close();

        std::string write_pres_fac_filename = "pres_fac_file_"+std::to_string(router_opts.icr_iter)+".txt";
        std::ofstream pres_fac_file(write_pres_fac_filename);
        pres_fac_file << final_pres_fac << "\n";
        pres_fac_file.close();




        if (timing_info) {
            VTR_LOG("Critical path: %g ns\n", 1e9 * best_routing_metrics.critical_path.delay());
        }

        VTR_LOG("Successfully routed after %d routing iterations.\n", itry);
    } else {
        VTR_LOG("Routing failed.\n");

        //If the routing fails, print the overused info
        print_overused_nodes_status(router_opts, overuse_info);

        ++num_routing_failed;
        //mycode
        std::string write_hist_filename = "history_cost_file_"+std::to_string(router_opts.icr_iter)+".txt";
        std::ofstream hist_file(write_hist_filename);
        for (const RRNodeId& rr_id : device_ctx.rr_graph.nodes()) {
            auto& node_inf = route_ctx.rr_node_route_inf[(size_t)rr_id];
            float history_cost = node_inf.acc_cost;
            //printf("%d %f\n", rr_id, history_cost);
            hist_file << (size_t)rr_id << " " << history_cost << "\n";
        }
        hist_file.close();

        std::string write_pres_fac_filename = "pres_fac_file_"+std::to_string(router_opts.icr_iter)+".txt";
        std::ofstream pres_fac_file(write_pres_fac_filename);
        pres_fac_file << final_pres_fac << "\n";
        pres_fac_file.close();


#ifdef VTR_ENABLE_DEBUG_LOGGING
        if (f_router_debug) print_invalid_routing_info(is_flat);
#endif
    }

    VTR_LOG("Final Net Connection Criticality Histogram:\n");
    print_router_criticality_histogram(*route_timing_info, netlist_pin_lookup);

    if (router_opts.incr_route == 1){
        size_t nets_routed_wo_first_iteration = router_stats.nets_routed - nets_routed_first_iteration;
        size_t connections_routed_wo_first_iteration = router_stats.connections_routed - connections_routed_first_iteration;
        size_t heap_pushes_wo_first_iteration = router_stats.heap_pushes - heap_pushes_first_iteration;
        size_t heap_pops_wo_first_iteration = router_stats.heap_pops - heap_pops_first_iteration;

        VTR_LOG("Stats for incremental routing, excluding first itertion\n");
        VTR_LOG("Router Stats: total_nets_routed: %zu total_connections_routed: %zu total_heap_pushes: %zu total_heap_pops: %zu\n",
                nets_routed_wo_first_iteration, connections_routed_wo_first_iteration, heap_pushes_wo_first_iteration, heap_pops_wo_first_iteration);
        
        VTR_LOG("=============================================================\n");
        VTR_LOG("Router Stats for first iteration: total_nets_routed: %zu total_connections_routed: %zu total_heap_pushes: %zu total_heap_pops: %zu\n",
                nets_routed_first_iteration, connections_routed_first_iteration, heap_pushes_first_iteration, heap_pops_first_iteration);
        VTR_LOG("Stats for incremental routing including all iterations\n");
        VTR_LOG("xxxxxxxxxxxxxxxxxxx Router Stats: total_nets_routed: %zu total_connections_routed: %zu total_heap_pushes: %zu total_heap_pops: %zu\n",
                router_stats.nets_routed, router_stats.connections_routed, router_stats.heap_pushes, router_stats.heap_pops);
        VTR_LOG("=============================================================\n");
    }
    else {
        VTR_LOG("Stats for one-stage routing\n");
        VTR_LOG("Router Stats: total_nets_routed: %zu total_connections_routed: %zu total_heap_pushes: %zu total_heap_pops: %zu total_wire_heap_pushes: %zu total_wire_heap_pops: %zu min_heap_pushes: %zu min_heap_pops: %zu max_heap_pushes: %zu max_heap_pops: %zu\n",
                router_stats.nets_routed, router_stats.connections_routed, router_stats.heap_pushes, router_stats.heap_pops, router_stats.wire_heap_pushes, router_stats.wire_heap_pops, router_stats.min_heap_pushes, router_stats.min_heap_pops, router_stats.max_heap_pushes, router_stats.max_heap_pops);
    }

    return routing_is_successful;
}

template<typename ConnectionRouter>
bool try_timing_driven_route_net_incr_route(const t_file_name_opts& filename_opts,
                                 ConnectionRouter& router,
                                 ClusterNetId net_id,
                                 int itry,
                                 float pres_fac,
				 std::unordered_map<int, int>& sink_order_index,
				 std::unordered_map<ClusterNetId, std::unordered_map<int, std::set<int>>>& branch_node_map,
                                 const t_router_opts& router_opts,
                                 CBRR& connections_inf,
                                 RouterStats& router_stats,
                                 float* pin_criticality,
                                 t_rt_node** rt_node_of_sink,
                                 ClbNetPinsMatrix<float>& net_delay,
                                 const ClusteredPinAtomPinsLookup& netlist_pin_lookup,
                                 std::shared_ptr<SetupHoldTimingInfo> timing_info,
                                 ClusteredPinTimingInvalidator* pin_timing_invalidator,
                                 route_budgets& budgeting_inf,
                                 bool& was_rerouted,
                                 float worst_negative_slack,
                                 const RoutingPredictor& routing_predictor,
                                 bool is_flat,
				 std::ofstream& net_order_file) {
    auto& cluster_ctx = g_vpr_ctx.clustering();
    auto& route_ctx = g_vpr_ctx.mutable_routing();
    
    // (PARSA) Luka, 2025: FOr debug purposes, don't route nets except 38908
    // if (size_t(net_id) != 38908) {
    //     route_ctx.net_status.set_is_routed(net_id, true);
    //     return true; // Skip routing this net
    // }
    bool is_routed = false;
    connections_inf.prepare_routing_for_net(net_id);
    
    bool reroute_for_hold = false;
    if (budgeting_inf.if_set()) {
        reroute_for_hold = (budgeting_inf.get_should_reroute(net_id));
        reroute_for_hold &= worst_negative_slack != 0;
    }

    if (route_ctx.net_status.is_fixed(net_id)) { /* Skip pre-routed nets. */
        is_routed = true;
    } else if (cluster_ctx.clb_nlist.net_is_ignored(net_id)) { /* Skip ignored nets. */
        is_routed = true;
    } 
    // For locking and loading branch changing nets
    //else if (!(reroute_for_hold) && should_route_net(net_id, connections_inf, true) == false && router_opts.ripup_all_nets == 0) {
    else if (!(reroute_for_hold) && should_route_net(net_id, connections_inf, true) == false && router_opts.ripup_all_nets == 0 && (itry > 2 || router_opts.detailed_router == 0)) {
        is_routed = true;
    } 
    else {
        // track time spent vs fanout
        profiling::net_fanout_start();

        is_routed = timing_driven_route_net_incr_route(filename_opts,
                                            router,
                                            net_id,
                                            itry,
                                            pres_fac,
					    sink_order_index,
					    branch_node_map,
                                            router_opts,
                                            connections_inf,
                                            router_stats,
                                            pin_criticality,
                                            rt_node_of_sink,
                                            net_delay[net_id].data(),
                                            netlist_pin_lookup,
                                            timing_info,
                                            pin_timing_invalidator,
                                            budgeting_inf,
                                            worst_negative_slack,
                                            routing_predictor,
                                            is_flat,
					    net_order_file);

        profiling::net_fanout_end(cluster_ctx.clb_nlist.net_sinks(net_id).size());

        /* Impossible to route? (disconnected rr_graph) */
        if (is_routed) {
            route_ctx.net_status.set_is_routed(net_id, true);
        } else {
            VTR_LOG("Routing failed for net %d\n", net_id);
        }

        was_rerouted = true; //Flag to record whether routing was actually changed
    }
    return (is_routed);
}
template<typename ConnectionRouter>
bool timing_driven_route_net_incr_route(const t_file_name_opts& filename_opts,
                             ConnectionRouter& router,
                             ClusterNetId net_id,
                             int itry,
                             float pres_fac,
                             std::unordered_map<int, int>& sink_order_index,
                             std::unordered_map<ClusterNetId, std::unordered_map<int, std::set<int>>>& branch_node_map,
                             const t_router_opts& router_opts,
                             CBRR& connections_inf,
                             RouterStats& router_stats,
                             float* pin_criticality,
                             t_rt_node** rt_node_of_sink,
                             float* net_delay,
                             const ClusteredPinAtomPinsLookup& netlist_pin_lookup,
                             std::shared_ptr<SetupHoldTimingInfo> timing_info,
                             ClusteredPinTimingInvalidator* pin_timing_invalidator,
                             route_budgets& budgeting_inf,
                             float worst_neg_slack,
                             const RoutingPredictor& routing_predictor,
                             bool is_flat,
			     std::ofstream& net_order_file) {
    /* Returns true as long as found some way to hook up this net, even if that *
     * way resulted in overuse of resources (congestion).  If there is no way   *
     * to route this net, even ignoring congestion, it returns false.  In this  *
     * case the rr_graph is disconnected and you can give up.                   */
    auto& cluster_ctx = g_vpr_ctx.clustering();
    auto& device_ctx = g_vpr_ctx.device();
    const auto& rr_graph = device_ctx.rr_graph;
    auto& route_ctx = g_vpr_ctx.routing();

    printf("Routing net %d\n", size_t(net_id));

    unsigned int num_sinks = cluster_ctx.clb_nlist.net_sinks(net_id).size();

    VTR_LOGV_DEBUG(f_router_debug, "Routing Net %zu (%zu sinks)\n", size_t(net_id), num_sinks);
    
    std::vector<int> factorials = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880};
    int min_detailed_nodes;
    int max_tries = 48;
    int all_permutation_max_fanout = 0;
    int max_sub_iterations = 0; 
    int t_min_incremental_reroute_fanout = router_opts.min_incremental_reroute_fanout;
    int t_num_sinks = num_sinks;


    if (router_opts.shuffle_first_iteration) {
        /*
            (PARSA) Luka, 2025: Shuffle sink order only in the first iteration.
        */
        max_sub_iterations = (router_opts.shuffle_first_iteration) ? router_opts.shuffle1 : 0;
    }
    else {
        /*
            (PARSA) Shashwat, 2024: Original sink order shuffling logic.
        */
        // just the loading of reconvergent nets for detailed router so no need to shuffle
        if (itry == 1) { 
            max_sub_iterations = (router_opts.detailed_router == 1) ? 0 : router_opts.shuffle1;
        }
        else if (itry > 1 && itry < 10){
            if (num_sinks >= factorials.size()){
                t_num_sinks = factorials.size() - 1; 
            }
                max_sub_iterations = std::min(factorials[t_num_sinks], router_opts.shuffle1);
            if (router_opts.shuffle1 > 120){
                    all_permutation_max_fanout = 5;
            }
            else if (router_opts.shuffle1 > 24){
                    all_permutation_max_fanout = 4;
            }
        }
        else if (itry >= 10) {
            if (num_sinks >= factorials.size()){
                t_num_sinks = factorials.size() - 1; 
            }
                max_sub_iterations = std::min(factorials[t_num_sinks], router_opts.shuffle2);
            if (router_opts.shuffle1 > 120){
                    all_permutation_max_fanout = 5;
            }
            else if (router_opts.shuffle1 > 24){
                    all_permutation_max_fanout = 4;
            }
        }
    }


    if (itry > 1) {
        min_detailed_nodes = connections_inf.get_minimum_detailed_nodes();	
    }
    std::vector<int> remaining_targets;
    std::map<std::tuple<int, float>, std::vector<int>> tree_cost_sink_order_map;
    std::map<std::tuple<float, int>, std::vector<int>> tree_cost_first_sink_order_map;
    std::map<std::tuple<int, float, int>, std::vector<int>> tree_ops_sink_order_map;

    const std::set<std::vector<int>>& sink_orders_from_prev_iterations = connections_inf.get_best_sink_orders();
    int additional_sub_iterations;
    if (num_sinks < t_min_incremental_reroute_fanout && router_opts.shuffle1 > 0 && router_opts.shuffle2 > 0) {
    	// complete rip up
	additional_sub_iterations = std::min(static_cast<int>(sink_orders_from_prev_iterations.size()), max_sub_iterations);
    } 
    else {
	additional_sub_iterations = 0;
    }
    max_sub_iterations -= additional_sub_iterations;
    // no need to check if max_sub_iterations can go below zero, as the max value of additional_sub_iterations is max_sub_iterations
    /*if (max_sub_iterations < 0){
	max_sub_iterations = 0;
    }*/
    VTR_ASSERT(max_sub_iterations >= 0);
    VTR_ASSERT(max_sub_iterations + additional_sub_iterations <= router_opts.shuffle1 || max_sub_iterations + additional_sub_iterations <= router_opts.shuffle2);
    t_rt_node* base_rt_root;
    std::vector<size_t> net_heap_pushes;
    std::vector<size_t> net_heap_pops;

    size_t min_heap_pushes;
    size_t min_heap_pops;
    for (int sink_order_itr = 0; sink_order_itr < max_sub_iterations + additional_sub_iterations + 1; sink_order_itr++) {
        // assert to check only one sub iteration when no shuffling
        VTR_ASSERT(!(router_opts.shuffle1 == 0 && sink_order_itr > 0));
        size_t heap_pushes_before = router_stats.heap_pushes;
        size_t heap_pops_before = router_stats.heap_pops;

        t_rt_node* rt_root;
        //VTR_LOG("itry: %d sub itr: %d net_id: %zu\n", itry, sink_order_itr, size_t(net_id));
        if (((itry > 2 && router_opts.detailed_router == 1) || (itry > 1 && router_opts.detailed_router == 0)) && sink_order_itr == 0 && num_sinks >= t_min_incremental_reroute_fanout) { // for standard router: in itry 1, the tree is just OPIN. For detailed router: itry 2 is the first iteration
        base_rt_root = traceback_to_route_tree(net_id);
        }
        else if (((itry > 2 && router_opts.detailed_router == 1) || (itry > 1 && router_opts.detailed_router == 0)) && num_sinks >= t_min_incremental_reroute_fanout){
            pathfinder_update_path_occupancy(route_ctx.trace[net_id].head, -1);
            traceback_from_route_tree(net_id, base_rt_root, num_sinks);
            pathfinder_update_path_occupancy(route_ctx.trace[net_id].head, 1);
            connections_inf.prepare_routing_for_net(net_id);
            /*if (itry >=10 && itry <13) { 
        VTR_LOG("done this sub itr: %d\n", sink_order_itr);
        }*/
        }
        rt_root = setup_routing_resources_incr_route(filename_opts,
            itry,
            net_id,
            num_sinks,
            //router_opts.min_incremental_reroute_fanout,
            t_min_incremental_reroute_fanout,
            connections_inf,
            rt_node_of_sink,
            router_opts,
            check_hold(router_opts, worst_neg_slack));

        bool high_fanout = is_high_fanout(num_sinks, router_opts.high_fanout_threshold);

        SpatialRouteTreeLookup spatial_route_tree_lookup;
        if (high_fanout) {
            spatial_route_tree_lookup = build_route_tree_spatial_lookup(net_id, rt_root);
        }

        // after this point the route tree is correct
        // remaining_targets from this point on are the **pin indices** that have yet to be routed
        //if (sink_order_itr == 0){
        auto& ref_remaining_targets = connections_inf.get_remaining_targets();
        

        if (router_opts.detailed_router == 0 && sink_order_itr == 0) {
            net_order_file << (size_t)net_id << " " << ref_remaining_targets.size() << "/" << num_sinks << " ";
        }
        //}
        if (sink_order_itr == 0){// need to fix remaining targets so that a new permutation  is produced in subsequent iterations
            remaining_targets = ref_remaining_targets;
            // VTR_LOG("Net Id %d", size_t(net_id));
            // for (auto target : remaining_targets) {
            //     VTR_LOG("%d ", target);
            // }
            // VTR_LOG("\n");
        }
        // calculate criticality of remaining target pins
        for (int ipin : remaining_targets) {
            if (timing_info) {
                auto clb_pin = cluster_ctx.clb_nlist.net_pin(net_id, ipin);
                if (!route_ctx.is_clock_net[net_id]) {
                    pin_criticality[ipin] = calculate_clb_net_pin_criticality(*timing_info, netlist_pin_lookup, clb_pin);
                } else {
                    // Use max_criticality for clock nets.
                    // calculate_clb_net_pin_criticality likely doesn't generate
                    // good values for clock nets.
                    //
                    // This will cause them to use min delay paths rather than
                    // avoid congestion. As a future enchancement, the clock nets
                    // should likely route for min slew, but that is a larger
                    // change.
                    pin_criticality[ipin] = router_opts.max_criticality;
                }

                /* Pin criticality is between 0 and 1.
                * Shift it downwards by 1 - max_criticality (max_criticality is 0.99 by default,
                * so shift down by 0.01) and cut off at 0.  This means that all pins with small
                * criticalities (<0.01) get criticality 0 and are ignored entirely, and everything
                * else becomes a bit less critical. This effect becomes more pronounced if
                * max_criticality is set lower. */
                // VTR_ASSERT(pin_criticality[ipin] > -0.01 && pin_criticality[ipin] < 1.01);
                pin_criticality[ipin] = std::max(pin_criticality[ipin] - (1.0 - router_opts.max_criticality), 0.0);

                /* Take pin criticality to some power (1 by default). */
                pin_criticality[ipin] = std::pow(pin_criticality[ipin], router_opts.criticality_exp);

                /* Cut off pin criticality at max_criticality. */
                pin_criticality[ipin] = std::min(pin_criticality[ipin], router_opts.max_criticality);
            } else {
                //No timing info, implies we want a min delay routing, so use criticality of 1.
                pin_criticality[ipin] = 1.;
            }
            //pin_criticality[ipin] = 1.;
        }

        int min_total_detailed_nodes;
        float min_total_cong_cost;

        if (router_opts.steiner_constraints) {
            /*
                (PARSA) Luka, 2025: When routing along RSMT constrained regions, in order to minimize wirecount
                and router effort (since lookahead is not helpful in this case, we minimize heap pops by making sure
                to route to the sink closes to the current partial route tree).
            */
            std::unordered_map<int, int> rsmt_sink_order = g_vpr_ctx.mutable_steiner().steiner_sink_orders[size_t(net_id)];
            std::sort(remaining_targets.begin(), remaining_targets.end(),
                    [&](int a, int b) {
                        return rsmt_sink_order[a] < rsmt_sink_order[b];
                    });
        }

        if (router_opts.shuffle_first_iteration) {
            /*
                (PARSA) Luka, 2025: Only shuffle the sink order in the first iteration of routing, then revert back to the default one.
            */
            if (itry == 1) {
                if (router_opts.detailed_router == 0) {
                    bool is_net_in_target_bracket = false;
                    if (router_opts.target_bracket == 1) {
                        if (cluster_ctx.clb_nlist.net_sinks(net_id).size() <= 5)
                            is_net_in_target_bracket = true;
                    }
                    else if (router_opts.target_bracket == 2) {
                        if (cluster_ctx.clb_nlist.net_sinks(net_id).size() <= 10)
                            is_net_in_target_bracket = true;
                    }
                    else if (router_opts.target_bracket == 3) {
                        if (cluster_ctx.clb_nlist.net_sinks(net_id).size() >= 100)
                            is_net_in_target_bracket = true;
                    }
                    else if (router_opts.target_bracket == 4) {
                        // All nets
                        is_net_in_target_bracket = true;
                    }            

                    if (is_net_in_target_bracket) {
                        std::cout << "Shuffling sink order for net " << size_t(net_id) << " with " << sink_order_itr << std::endl;
                        if (sink_order_itr == max_sub_iterations + additional_sub_iterations) {
                            if (router_opts.tree_type == MIN_TOTAL_NODES) { 
                                auto it = tree_cost_sink_order_map.begin();
                                std::tuple<int, float> minKey = it->first;
                                min_total_detailed_nodes = std::get<0>(minKey);
                                min_total_cong_cost = std::get<1>(minKey);
                                remaining_targets = it->second;
                            }
                            else if (router_opts.tree_type == MIN_CONG_COST) {
                                auto it = tree_cost_first_sink_order_map.begin();
                                std::tuple<float, int> minKey = it->first;
                                min_total_detailed_nodes = std::get<1>(minKey);
                                min_total_cong_cost = std::get<0>(minKey);
                                remaining_targets = it->second;
                            }
                            else if (router_opts.tree_type == MIN_HEAP_PUSHES) {
                                auto it = tree_ops_sink_order_map.begin();
                                std::tuple<int, float, int> minKey = it->first;
                                int total_detailed_nodes_min_ops = std::get<2>(minKey);
                                remaining_targets = it->second;
                            }
                            else if (router_opts.tree_type == MAX_TOTAL_NODES) {
                                auto it = tree_cost_sink_order_map.rbegin();
                                std::tuple<int, float> maxKey = it->first;
                                remaining_targets = it->second;
                            }
                            else {
                                VTR_ASSERT("Wrong flag for argument tree\n");
                            }
                            connections_inf.add_best_sink_order(remaining_targets);	
                        }
                        else if (sink_order_itr == 0) {
                            sort(begin(remaining_targets), end(remaining_targets), [&](int a, int b) {return a < b;});
                        }
                        else if (sink_order_itr > 0) {
                            if (remaining_targets.size() <= all_permutation_max_fanout) {
                                std::next_permutation(remaining_targets.begin(), remaining_targets.end()); 
                            }
                            else {
                                std::random_device rd;  // Seed for random number generator
                                std::mt19937 g(rd());   // Standard Mersenne Twister engine
                                std::shuffle(begin(remaining_targets), end(remaining_targets), g);
                            }
                        }
                    }
                    else {
                        sink_order_itr = max_sub_iterations + additional_sub_iterations;
                        sort(begin(remaining_targets), end(remaining_targets), [&](int a, int b) {return a < b;});
                        connections_inf.add_best_sink_order(remaining_targets);
                    }
                }
            }
            else {
                sort(begin(remaining_targets), end(remaining_targets), [&](int a, int b) {return a < b;});
            }
        }
        else if (router_opts.dependency_graph_sink_order) {
            /*
                (PARSA) Luka, 2025: Use the sink order derived from the net's rectilinear Steiner
                minimal tree.
            */
            std::unordered_map<int, int> dependency_graph_sink_order = g_vpr_ctx.mutable_steiner().steiner_sink_orders[size_t(net_id)];
            std::sort(remaining_targets.begin(), remaining_targets.end(),
                [&](int a, int b) {
                    return dependency_graph_sink_order[a] < dependency_graph_sink_order[b];
                });
        }
        else {
            /* 
                (PARSA) Shashwat 2024: Original sink order shuffling logic.
                [Luka, July 1st 2025]: This should be all of the original code; in case of errors, look into uncommenting lines of code.
            */
            // compare the criticality of different sink nodes
            sort(begin(remaining_targets), end(remaining_targets), Criticality_comp{pin_criticality});
            if (router_opts.detailed_router == 1 && router_opts.preorder_sink_order == 1){
                // // Change sink order only for the first iteration
                // if (itry == 1) {
                //     sort(begin(remaining_targets), end(remaining_targets), [&](int a, int b) {return sink_order_index[a] < sink_order_index[b];});
                // }
            	// else {
                //     sort(begin(remaining_targets), end(remaining_targets), [&](int a, int b) {return a < b;});
                // }
                std::random_device rd;  // Seed for random number generator
                std::mt19937 g(rd());   // Standard Mersenne Twister engine
                std::shuffle(begin(remaining_targets), end(remaining_targets), g);
            }
            else if (sink_order_itr == 0) {
                sort(begin(remaining_targets), end(remaining_targets), Criticality_comp{pin_criticality});
            /*if (size_t(net_id) == 6366) {
                    VTR_LOG("  order:"); 
                    for (unsigned itarget = 0; itarget < remaining_targets.size(); ++itarget) {
                        int target_pin = remaining_targets[itarget];
                        VTR_LOG(" %d", target_pin);
                    }
                    VTR_LOG("\n"); 
            }*/
            }
            else if (sink_order_itr < max_sub_iterations){
                    if (remaining_targets.size() <= all_permutation_max_fanout){
                        std::next_permutation(remaining_targets.begin(), remaining_targets.end()); 
                    }
                    else {
                    // nondeterministic shuffling
                        //std::random_device rd;  // Seed for random number generator
                        //std::mt19937 g(rd());   // Standard Mersenne Twister engine
                        // deterministic shuffling across different instances of VTR
                        std::mt19937 g(itry * 100 + size_t(net_id) + 1000 * sink_order_itr); 
                        std::shuffle(begin(remaining_targets), end(remaining_targets), g);
                    }
                }
            else if (sink_order_itr == max_sub_iterations + additional_sub_iterations){
            if (router_opts.tree_type == MIN_TOTAL_NODES){ 
                    auto it = tree_cost_sink_order_map.begin();
                std::tuple<int, float> minKey = it->first;
                    min_total_detailed_nodes = std::get<0>(minKey);
                    min_total_cong_cost = std::get<1>(minKey);
                    remaining_targets = it->second;
            }
            else if (router_opts.tree_type == MIN_CONG_COST){
                    auto it = tree_cost_first_sink_order_map.begin();
                std::tuple<float, int> minKey = it->first;
                    min_total_detailed_nodes = std::get<1>(minKey);
                    min_total_cong_cost = std::get<0>(minKey);
                    remaining_targets = it->second;
            }
            else if (router_opts.tree_type == MIN_HEAP_PUSHES){
                    auto it = tree_ops_sink_order_map.begin();
                std::tuple<int, float, int> minKey = it->first;
                    int total_detailed_nodes_min_ops = std::get<2>(minKey);
                    remaining_targets = it->second;
                    
                //auto it2 = tree_cost_sink_order_map.begin();
                //std::tuple<int, float> minKey_t = it2->first;
                    //min_total_detailed_nodes = std::get<0>(minKey_t);
                /*if(min_total_detailed_nodes < total_detailed_nodes_min_ops) {
                    VTR_LOG("[FOUND] A tree with minimum heap operations but not minimum possible total nodes that could be achieved with sampled sink orders. Net: %zu total_nodes: %d (min_possible: %d). This implies we cannot do early stopping to find the tree with minimum nodes\n", size_t(net_id), total_detailed_nodes_min_ops, min_total_detailed_nodes);
                }*/
            }
            else {
                VTR_ASSERT("Wrong flag for argument tree\n");
            }
                connections_inf.add_best_sink_order(remaining_targets);	
        
            }

            else if ((sink_order_itr >= max_sub_iterations && sink_order_itr < max_sub_iterations + additional_sub_iterations) && additional_sub_iterations != 0){
            auto it = std::next(sink_orders_from_prev_iterations.begin(), sink_order_itr - max_sub_iterations);
            remaining_targets = *it;
            //remaining_targets = sink_orders_from_prev_iterations[sink_order_itr - max_sub_iterations];
            }
        }

        if (size_t(net_id) == 38908) {
            VTR_LOG("Net 38908 sink order: ");
            for (auto target : remaining_targets) {
                VTR_LOG("%d, ", target);
            }
            VTR_LOG("\n");
        }
        
        /* Update base costs according to fanout and criticality rules */
        update_rr_base_costs(num_sinks);

        t_conn_delay_budget conn_delay_budget;
        t_conn_cost_params cost_params;
        cost_params.astar_fac = router_opts.astar_fac;
        cost_params.bend_cost = router_opts.bend_cost;
        cost_params.pres_fac = pres_fac;
        cost_params.delay_budget = ((budgeting_inf.if_set()) ? &conn_delay_budget : nullptr);
        cost_params.bias = router_opts.sbNode_lookahead_factor;
        cost_params.offpath_penalty = router_opts.offpath_penalty;
        cost_params.detailed_router = router_opts.detailed_router;
        cost_params.relax_hop_order = router_opts.relax_hop_order;
        cost_params.global_occ_factor = router_opts.global_occ_factor;

        // Pre-route to clock source for clock nets (marked as global nets)
        if (cluster_ctx.clb_nlist.net_is_global(net_id) && router_opts.two_stage_clock_routing) {
            //VTR_ASSERT(router_opts.clock_modeling == DEDICATED_NETWORK);
            int sink_node = device_ctx.virtual_clock_network_root_idx;

            enable_router_debug(router_opts, net_id, sink_node, itry, &router);

            VTR_LOGV_DEBUG(f_router_debug, "Pre-routing global net %zu\n", size_t(net_id));

            // Set to the max timing criticality which should intern minimize clock insertion
            // delay by selecting a direct route from the clock source to the virtual sink
            cost_params.criticality = router_opts.max_criticality;
            cost_params.bias = router_opts.sbNode_lookahead_factor;
            if (!timing_driven_pre_route_to_clock_root(
                    router,
                    net_id,
                    sink_node,
                    cost_params,
                    router_opts.high_fanout_threshold,
                    rt_root,
                    spatial_route_tree_lookup,
                    router_stats,
                    is_flat)) {
                return false;
            }
        }

        if (budgeting_inf.if_set()) {
            budgeting_inf.set_should_reroute(net_id, false);
        }

    // explore in order of decreasing criticality (no longer need sink_order array)
        for (unsigned itarget = 0; itarget < remaining_targets.size(); ++itarget) {
            int target_pin = remaining_targets[itarget];

            std::set<int> branch_nodes;
            if (router_opts.detailed_router == 1){
                    branch_nodes = branch_node_map[net_id][target_pin]; // all nodesare part of same global node
            }
            /*VTR_LOG("target_pin: %d\n", target_pin);        		
            for (auto node : branch_nodes){
                VTR_LOG("Branch nodes: %d\n", node);        		
            
            }*/
            int sink_rr = route_ctx.net_rr_terminals[net_id][target_pin];

                enable_router_debug(router_opts, net_id, sink_rr, itry, &router);

                VTR_LOGV_DEBUG(f_router_debug, "Routing Net %zu (%zu sinks) sub_iter: %d\n", size_t(net_id), num_sinks, sink_order_itr);

                cost_params.criticality = pin_criticality[target_pin];
                cost_params.bias = router_opts.sbNode_lookahead_factor;

                if (budgeting_inf.if_set()) {
                    conn_delay_budget.max_delay = budgeting_inf.get_max_delay_budget(net_id, target_pin);
                    conn_delay_budget.target_delay = budgeting_inf.get_delay_target(net_id, target_pin);
                    conn_delay_budget.min_delay = budgeting_inf.get_min_delay_budget(net_id, target_pin);
                    conn_delay_budget.short_path_criticality = budgeting_inf.get_crit_short_path(net_id, target_pin);
                    conn_delay_budget.routing_budgets_algorithm = router_opts.routing_budgets_algorithm;
                }

                profiling::conn_start();

                // build a branch in the route tree to the target
                if (!timing_driven_route_sink(router,
                                            net_id,
                                            itarget,
                                            target_pin,
                                            cost_params,
                                            router_opts,
                                            rt_root, rt_node_of_sink,
                                            spatial_route_tree_lookup,
                                            router_stats,
                                            budgeting_inf,
                                            routing_predictor,
                                            is_flat,
                                            branch_nodes,
                                            itry))
                    return false;

                profiling::conn_finish(route_ctx.net_rr_terminals[net_id][0],
                                    sink_rr,
                                    pin_criticality[target_pin]);

            if (sink_order_itr == 0){
                ++router_stats.connections_routed;
            }
        } // finished all sinks

    std::pair<int, float> cost = get_tree_cost(route_ctx.trace[net_id].head);
    int total_detailed_nodes = cost.first;
    float total_cong_cost = cost.second;
    tree_cost_sink_order_map[std::make_tuple(total_detailed_nodes, total_cong_cost)] = remaining_targets;
    tree_cost_first_sink_order_map[std::make_tuple(total_cong_cost, total_detailed_nodes)] = remaining_targets;


    if (itry == 1 && sink_order_itr == 0) {
        connections_inf.set_minimum_detailed_nodes(total_detailed_nodes);	
    }
    if (sink_order_itr == 0){
        ++router_stats.nets_routed;
    }
    profiling::net_finish();

    /* For later timing analysis. */

    // may have to update timing delay of the previously legally reached sinks since downstream capacitance could be changed
    update_net_delays_from_route_tree(net_delay, rt_node_of_sink, net_id, timing_info.get(), pin_timing_invalidator);

    if (router_opts.update_lower_bound_delays) {
        for (int ipin : remaining_targets) {
            connections_inf.update_lower_bound_connection_delay(net_id, ipin, net_delay[ipin]);
        }
    }

    if (!cluster_ctx.clb_nlist.net_is_ignored(net_id)) {
        for (unsigned ipin = 1; ipin < cluster_ctx.clb_nlist.net_pins(net_id).size(); ++ipin) {
            if (net_delay[ipin] == 0) { // should be SOURCE->OPIN->IPIN->SINK
                VTR_ASSERT(rr_graph.node_type(RRNodeId(rt_node_of_sink[ipin]->parent_node->parent_node->inode)) == OPIN);
            }
        }
    }
    VTR_ASSERT_MSG(route_ctx.rr_node_route_inf[rt_root->inode].occ() <= rr_graph.node_capacity(RRNodeId(rt_root->inode)), "SOURCE should never be congested");

    // route tree is not kept persistent since building it from the traceback the next iteration takes almost 0 time
    VTR_LOGV_DEBUG(f_router_debug, "Routed Net %zu (%zu sinks)\n", size_t(net_id), num_sinks);

    free_route_tree(rt_root);
    router.empty_rcv_route_tree_set();
    /*if ((itry == 10 || itry == 11 || itry == 12) && remaining_targets.size() >= all_permutation_max_fanout){
       VTR_LOG("itry: %d (%d) net_id: %zu total_nodes: %d (min: %d) cost: %f pool_size: %d\n", itry, sink_order_itr, size_t(net_id), total_detailed_nodes, min_detailed_nodes, total_cong_cost, additional_sub_iterations); 
       VTR_LOG("  order:"); 
       for (unsigned itarget = 0; itarget < remaining_targets.size(); ++itarget) {
           int target_pin = remaining_targets[itarget];
       	   VTR_LOG(" %d", target_pin);
       }
       VTR_LOG("\n"); 
    }*/
    /*else if (itry == 2){
       VTR_LOG(" -- itry: %d (%d) net_id: %zu total_nodes: %d (min: %d) cost: %f\n", itry, sink_order_itr, size_t(net_id), total_detailed_nodes, min_detailed_nodes, total_cong_cost); 
    }*/
    /*if (total_detailed_nodes == min_detailed_nodes && total_cong_cost == min_detailed_nodes){
    	break;
    }*/
    // net_id <heap pushes for each permutation> <min> <max>
    // min_heap_pushes += 
    // max_heap_pushes += 

    	

    size_t heap_pushes_after = router_stats.heap_pushes;
    size_t heap_pops_after = router_stats.heap_pops;

    size_t sub_iter_heap_pushes = heap_pushes_after - heap_pushes_before;
    size_t sub_iter_heap_pops = heap_pops_after - heap_pops_before;
    net_heap_pushes.push_back(sub_iter_heap_pushes);
    net_heap_pops.push_back(sub_iter_heap_pops);

    tree_ops_sink_order_map[std::make_tuple(sub_iter_heap_pushes, total_cong_cost, total_detailed_nodes)] = remaining_targets;

    if (sink_order_itr == max_sub_iterations + additional_sub_iterations){// whichever tree_type is used, operations corresponding to that will be reported
    	min_heap_pushes = sub_iter_heap_pushes;
	min_heap_pops = sub_iter_heap_pops;	
    }
    }

    auto [min_heap_pushes_it, max_heap_pushes_it] = std::minmax_element(net_heap_pushes.begin(), net_heap_pushes.end());
    auto [min_heap_pops_it, max_heap_pops_it] = std::minmax_element(net_heap_pops.begin(), net_heap_pops.end());
    int max_index = std::distance(net_heap_pushes.begin(), max_heap_pushes_it);
    size_t corresponding_heap_pops = net_heap_pops[max_index];	
    
    //VTR_ASSERT (*min_heap_pushes_it == min_heap_pushes);
    //ASSERT (*min_heap_pops_it == min_heap_pops);
    router_stats.min_heap_pushes += min_heap_pushes; 
    router_stats.max_heap_pushes += *max_heap_pushes_it; 
    

    router_stats.min_heap_pops += min_heap_pops; 
    //router_stats.max_heap_pops += *max_heap_pops_it; 
    router_stats.max_heap_pops += corresponding_heap_pops; 
    return (true);
}
static t_rt_node* setup_routing_resources_incr_route(const t_file_name_opts& filename_opts,
                                          int itry,
                                          ClusterNetId net_id,
                                          unsigned num_sinks,
                                          int min_incremental_reroute_fanout,
                                          CBRR& connections_inf,
                                          t_rt_node** rt_node_of_sink,
                                          const t_router_opts& router_opts,
                                          bool ripup_high_fanout_nets) {
    /* Build and return a partial route tree from the legal connections from last iteration.
     * along the way do:
     * 	update pathfinder costs to be accurate to the partial route tree
     *	update the net's traceback to be accurate to the partial route tree
     * 	find and store the pins that still need to be reached in incremental_rerouting_resources.remaining_targets
     * 	find and store the rt nodes that have been reached in incremental_rerouting_resources.reached_rt_sinks
     *	mark the rr_node sinks as targets to be reached */

    auto& route_ctx = g_vpr_ctx.routing();
    t_rt_node* rt_root;
    // for nets below a certain size (min_incremental_reroute_fanout), rip up any old routing
    // otherwise, we incrementally reroute by reusing legal parts of the previous iteration
    // convert the previous iteration's traceback into the starting route tree for this iteration
    // For locking and loading branch changing nets
    //if ((int)num_sinks < min_incremental_reroute_fanout || itry == 1 || ripup_high_fanout_nets) {
    if ((int)num_sinks < min_incremental_reroute_fanout || itry == 1 || ripup_high_fanout_nets || (itry == 2 && router_opts.detailed_router == 1)) {
        /*for (int illegal_net_i = 0; illegal_net_i < total_illegal_nets; illegal_net_i++){
            ClusterNetId debug_net = (ClusterNetId)illegal_net_list[illegal_net_i];
	        if (net_id == debug_net){
	    	    VTR_LOG("[COMPLETE RIP UP] NET: %d SINK: %d\n", net_id, num_sinks);
	        }
        }*/    
	//VTR_LOG("[COMPLETE RIP UP] NET: %zu SINK: %d\n", size_t(net_id), num_sinks);
        complete_rip_up_net_count++;
        profiling::net_rerouted();

        // rip up the whole net
        
        //read_route_incr_route(net_id, filename_opts.RouteFile.c_str(), router_opts, filename_opts.verify_file_digests);
        pathfinder_update_path_occupancy(route_ctx.trace[net_id].head, -1);
        free_traceback(net_id);

        rt_root = init_route_tree_to_source(net_id);
        for (unsigned int sink_pin = 1; sink_pin <= num_sinks; ++sink_pin)
            connections_inf.toreach_rr_sink(sink_pin);
        // since all connections will be rerouted for this net, clear all of net's forced reroute flags
        connections_inf.clear_force_reroute_for_net();

        // when we don't prune the tree, we also don't know the sink node indices
        // thus we'll use functions that act on pin indices like mark_ends instead
        // of their versions that act on node indices directly like mark_remaining_ends
        mark_ends(net_id);
        
    } else {
        /*for (int illegal_net_i = 0; illegal_net_i < total_illegal_nets; illegal_net_i++){
            ClusterNetId debug_net = (ClusterNetId)illegal_net_list[illegal_net_i];
	        if (net_id == debug_net){
	    	    VTR_LOG("[PARTIAL RIP UP] NET: %d SINK: %d\n", net_id, num_sinks);
	        }
        }*/
	    //if (net_id == debug_net){
        //    VTR_LOG("[PARTIAL RIP UP] NET: %d SINK: %d\n", net_id, num_sinks);
        //}
        partial_rip_up_net_count++;    
        auto& reached_rt_sinks = connections_inf.get_reached_rt_sinks();
        auto& remaining_targets = connections_inf.get_remaining_targets();

        profiling::net_rebuild_start();

        // convert the previous iteration's traceback into a route tree
        rt_root = traceback_to_route_tree(net_id);

        //Sanity check that route tree and traceback are equivalent before pruning
        VTR_ASSERT_DEBUG(verify_traceback_route_tree_equivalent(route_ctx.trace[net_id].head, rt_root));

        // check for edge correctness
        VTR_ASSERT_SAFE(is_valid_skeleton_tree(rt_root));

        // Skip this check if RCV is enabled, as RCV can use another method to cause reroutes
        VTR_ASSERT_SAFE(should_route_net(net_id, connections_inf, true) || router_opts.routing_budgets_algorithm == YOYO);

        //Prune the branches of the tree that don't legally lead to sinks
        rt_root = prune_route_tree(rt_root, connections_inf);

        //Now that the tree has been pruned, we can free the old traceback
        // NOTE: this must happen *after* pruning since it changes the
        //       recorded congestion
        pathfinder_update_path_occupancy(route_ctx.trace[net_id].head, -1);
        free_traceback(net_id);

        if (rt_root) { //Partially pruned
            profiling::route_tree_preserved();

            //Since we have a valid partial routing (to at least one SINK)
            //we need to make sure the traceback is synchronized to the route tree
            traceback_from_route_tree(net_id, rt_root, reached_rt_sinks.size());

            //Sanity check the traceback for self-consistency
            VTR_ASSERT_DEBUG(validate_traceback(route_ctx.trace[net_id].head));

            //Sanity check that route tree and traceback are equivalent after pruning
            VTR_ASSERT_DEBUG(verify_traceback_route_tree_equivalent(route_ctx.trace[net_id].head, rt_root));

            // put the updated occupancies of the route tree nodes back into pathfinder
            pathfinder_update_path_occupancy(route_ctx.trace[net_id].head, 1);

        } else { //Fully destroyed
            profiling::route_tree_pruned();

            //Initialize only to source
            rt_root = init_route_tree_to_source(net_id);
            //pathfinder_update_path_occupancy(route_ctx.trace[net_id].head, 1);

            //NOTE: We leave the traceback uninitialized, so update_traceback()
            //      will correctly add the SOURCE node when the branch to
            //      the first SINK is found.
            VTR_ASSERT(route_ctx.trace[net_id].head == nullptr);
            VTR_ASSERT(route_ctx.trace[net_id].tail == nullptr);
            VTR_ASSERT(route_ctx.trace_nodes[net_id].empty());
        }

        //Update R/C
        load_new_subtree_R_upstream(rt_root);
        load_new_subtree_C_downstream(rt_root);

        VTR_ASSERT(reached_rt_sinks.size() + remaining_targets.size() == num_sinks);

        //Record current routing
        add_route_tree_to_rr_node_lookup(rt_root);

        // give lookup on the reached sinks
        for (t_rt_node* sink_node : reached_rt_sinks) {
            rt_node_of_sink[sink_node->net_pin_index] = sink_node;
        }

        profiling::net_rebuild_end(num_sinks, remaining_targets.size());

        // check for R_upstream C_downstream and edge correctness
        VTR_ASSERT_SAFE(is_valid_route_tree(rt_root));
        // congestion should've been pruned away
        VTR_ASSERT_SAFE(is_uncongested_route_tree(rt_root));

        // mark remaining ends
        mark_remaining_ends(net_id, remaining_targets);

        // still need to calculate the tree's time delay (0 Tarrival means from SOURCE)
        load_route_tree_Tdel(rt_root, 0);

        // mark the lookup (rr_node_route_inf) for existing tree elements as NO_PREVIOUS so add_to_path stops when it reaches one of them
        load_route_tree_rr_route_inf(rt_root);
    }

    // completed constructing the partial route tree and updated all other data structures to match
    return rt_root;
}
//================================
#ifndef NO_GRAPHICS
//updates router iteration information and checks for router iteration and net id breakpoints
//stops after the specified router iteration or net id is encountered
void update_router_info_and_check_bp(bp_router_type type, int net_id) {
    t_draw_state* draw_state = get_draw_state_vars();
    if (draw_state->list_of_breakpoints.size() != 0) {
        if (type == BP_ROUTE_ITER)
            get_bp_state_globals()->get_glob_breakpoint_state()->router_iter++;
        else if (type == BP_NET_ID)
            get_bp_state_globals()->get_glob_breakpoint_state()->route_net_id = net_id;
        f_router_debug = check_for_breakpoints(false);
        if (f_router_debug) {
            breakpoint_info_window(get_bp_state_globals()->get_glob_breakpoint_state()->bp_description, *get_bp_state_globals()->get_glob_breakpoint_state(), false);
            update_screen(ScreenUpdatePriority::MAJOR, "Breakpoint Encountered", ROUTING, nullptr);
        }
    }
}
#endif
