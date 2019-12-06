#pragma once

#include <vector>
#include <cstring>
#include <algorithm>

#include "rr_graph.h"
#include "rr_graph_uxsdcxx_interface.h"
#include "rr_node.h"
#include "vpr_error.h"
#include "rr_metadata.h"
#include "read_xml_arch_file.h"
#include "vtr_log.h"
#include "vtr_version.h"
#include "vpr_utils.h"
#include "check_rr_graph.h"
#include "rr_graph2.h"
#include "rr_graph_indexed_data.h"

template<typename T>
void make_room_in_vector(T* vec, size_t elem_position) {
    if (elem_position < vec->size()) {
        return;
    }

    size_t capacity = std::max(vec->capacity(), size_t(16));
    while (elem_position >= capacity) {
        capacity *= 2;
    }

    if (capacity >= vec->capacity()) {
        vec->reserve(capacity);
    }

    vec->resize(elem_position + 1);
}

class RrGraphImpl : public uxsd::RrGraphBase {
  public:
    RrGraphImpl(
        bool is_global_graph,
        bool read_edge_metadata,
        t_chan_width* chan_width,
        std::vector<t_rr_node>* rr_nodes,
        std::vector<t_rr_switch_inf>* rr_switch_inf,
        std::vector<t_rr_indexed_data>* rr_indexed_data,
        t_rr_node_indices* rr_node_indices,
        const size_t num_arch_switches,
        const t_arch_switch_inf* arch_switch_inf,
        const std::vector<t_segment_inf>& segment_inf,
        const std::vector<t_physical_tile_type>& physical_tile_types,
        const DeviceGrid& grid,
        const std::unordered_map<int, t_metadata_dict>& rr_node_metadata,
        const std::unordered_map<std::tuple<int, int, short>,
                                 t_metadata_dict>& rr_edge_metadata)
        : is_global_graph_(is_global_graph)
        , read_edge_metadata_(read_edge_metadata)
        , chan_width_(chan_width)
        , rr_nodes_(rr_nodes)
        , rr_switch_inf_(rr_switch_inf)
        , rr_indexed_data_(rr_indexed_data)
        , rr_node_indices_(rr_node_indices)
        , num_arch_switches_(num_arch_switches)
        , arch_switch_inf_(arch_switch_inf)
        , segment_inf_(segment_inf)
        , physical_tile_types_(physical_tile_types)
        , grid_(grid)
        , rr_node_metadata_(rr_node_metadata)
        , rr_edge_metadata_(rr_edge_metadata) {}

    /** Generated for complex type "timing":
     * <xs:complexType name="timing">
     *   <xs:attribute name="R" type="xs:float" />
     *   <xs:attribute name="Cin" type="xs:float" />
     *   <xs:attribute name="Cinternal" type="xs:float" />
     *   <xs:attribute name="Cout" type="xs:float" />
     *   <xs:attribute name="Tdel" type="xs:float" />
     * </xs:complexType>
     */
    inline void set_timing_Cin(float Cin, const void* /*data*/, void* iter) override {
        static_cast<t_rr_switch_inf*>(iter)->Cin = Cin;
    }
    inline float get_timing_Cin(const void* data, void* /*iter*/) override {
        return static_cast<const t_rr_switch_inf*>(data)->Cin;
    }
    inline void set_timing_Cinternal(float Cinternal, const void* /*data*/, void* iter) override {
        static_cast<t_rr_switch_inf*>(iter)->Cinternal = Cinternal;
    }
    inline float get_timing_Cinternal(const void* data, void* /*iter*/) override {
        return static_cast<const t_rr_switch_inf*>(data)->Cinternal;
    }
    inline void set_timing_Cout(float Cout, const void* /*data*/, void* iter) override {
        static_cast<t_rr_switch_inf*>(iter)->Cout = Cout;
    }
    inline float get_timing_Cout(const void* data, void* /*iter*/) override {
        return static_cast<const t_rr_switch_inf*>(data)->Cout;
    }
    inline void set_timing_R(float R, const void* /*data*/, void* iter) override {
        static_cast<t_rr_switch_inf*>(iter)->R = R;
    }
    inline float get_timing_R(const void* data, void* /*iter*/) override {
        return static_cast<const t_rr_switch_inf*>(data)->R;
    }
    inline void set_timing_Tdel(float Tdel, const void* /*data*/, void* iter) override {
        static_cast<t_rr_switch_inf*>(iter)->Tdel = Tdel;
    }
    inline float get_timing_Tdel(const void* data, void* /*iter*/) override {
        return static_cast<const t_rr_switch_inf*>(data)->Tdel;
    }

    /** Generated for complex type "switch":
     * <xs:complexType name="switch">
     *   <xs:all>
     *     <xs:element minOccurs="0" name="timing" type="timing" />
     *     <xs:element name="sizing" type="sizing" />
     *   </xs:all>
     *   <xs:attribute name="id" type="xs:int" use="required" />
     *   <xs:attribute name="name" type="xs:string" use="required" />
     *   <xs:attribute name="type" type="switch_type" />
     * </xs:complexType>
     */
    inline void set_switch_name(const char* name, const void* /*data*/, void* iter) override {
        bool found_arch_name = false;
        for (size_t i = 0; i < num_arch_switches_; ++i) {
            if (strcmp(name, arch_switch_inf_[i].name) == 0) {
                name = arch_switch_inf_[i].name;
                found_arch_name = true;
                break;
            }
        }
        if (!found_arch_name) {
            VPR_FATAL_ERROR(VPR_ERROR_ROUTE,
                            "Switch name '%s' not found in architecture\n", name);
        }

        static_cast<t_rr_switch_inf*>(iter)->name = name;
    }
    inline void set_switch_type(uxsd::enum_switch_type type, const void* /*data*/, void* iter) override {
        SwitchType switch_type = SwitchType::INVALID;
        switch (type) {
            case uxsd::enum_switch_type::TRISTATE:
                switch_type = SwitchType::TRISTATE;
                break;
            case uxsd::enum_switch_type::MUX:
                switch_type = SwitchType::MUX;
                break;
            case uxsd::enum_switch_type::PASS_GATE:
                switch_type = SwitchType::PASS_GATE;
                break;
            case uxsd::enum_switch_type::SHORT:
                switch_type = SwitchType::SHORT;
                break;
            case uxsd::enum_switch_type::BUFFER:
                switch_type = SwitchType::BUFFER;
                break;
            default:
                VPR_FATAL_ERROR(VPR_ERROR_ROUTE,
                                "Invalid switch type '%d'\n", type);
        }

        static_cast<t_rr_switch_inf*>(iter)->set_type(switch_type);
    }
    inline uxsd::enum_switch_type get_switch_type(const void* data, void* /*iter*/) override {
        auto* switch_obj = static_cast<const t_rr_switch_inf*>(data);
        switch (switch_obj->type()) {
            case SwitchType::TRISTATE:
                return uxsd::enum_switch_type::TRISTATE;
            case SwitchType::MUX:
                return uxsd::enum_switch_type::MUX;
            case SwitchType::PASS_GATE:
                return uxsd::enum_switch_type::PASS_GATE;
            case SwitchType::SHORT:
                return uxsd::enum_switch_type::SHORT;
            case SwitchType::BUFFER:
                return uxsd::enum_switch_type::BUFFER;
            default:
                VPR_FATAL_ERROR(VPR_ERROR_ROUTE,
                                "Invalid switch type '%d'\n",
                                switch_obj->type());
        }

        return uxsd::enum_switch_type::UXSD_INVALID;
    }

    inline std::pair<const void*, void*> init_switch_timing(const void* data, void* iter) override {
        return std::make_pair(data, iter);
    }
    inline void finish_switch_timing(const void* /*data*/, void* /*iter*/) override {}
    inline std::pair<const void*, void*> get_switch_timing(const void* data, void* iter) override {
        return std::make_pair(data, iter);
    }
    inline bool has_switch_timing(const void* /*data*/, void* /*iter*/) override {
        return true;
    }

    inline std::pair<const void*, void*> init_switch_sizing(const void* /*data*/, void* iter, float buf_size, float mux_trans_size) override {
        auto* switch_obj = static_cast<t_rr_switch_inf*>(iter);
        switch_obj->buf_size = buf_size;
        switch_obj->mux_trans_size = mux_trans_size;
        return std::make_pair(nullptr, nullptr);
    }
    inline void finish_switch_sizing(const void* /*data*/, void* /*iter*/) override {}
    inline std::pair<const void*, void*> get_switch_sizing(const void* data, void* iter) override {
        return std::make_pair(data, iter);
    }
    inline float get_sizing_buf_size(const void* data, void* /*iter*/) override {
        return static_cast<const t_rr_switch_inf*>(data)->buf_size;
    }
    inline float get_sizing_mux_trans_size(const void* data, void* /*iter*/) override {
        return static_cast<const t_rr_switch_inf*>(data)->mux_trans_size;
    }

    /** Generated for complex type "switches":
     * <xs:complexType name="switches">
     *   <xs:sequence>
     *     <xs:element maxOccurs="unbounded" name="switch" type="switch" />
     *   </xs:sequence>
     * </xs:complexType>
     */
    inline std::pair<const void*, void*> add_switches_switch(const void* /*data*/, void* /*iter*/, int id) override {
        make_room_in_vector(rr_switch_inf_, id);

        (*rr_switch_inf_)[id].R = 0;
        (*rr_switch_inf_)[id].Cin = 0;
        (*rr_switch_inf_)[id].Cout = 0;
        (*rr_switch_inf_)[id].Cinternal = 0;
        (*rr_switch_inf_)[id].Tdel = 0;

        return std::make_pair(nullptr, &(*rr_switch_inf_)[id]);
    }
    inline void finish_switches_switch(const void* /*data*/, void* /*iter*/) override {}
    inline size_t num_switches_switch(const void* /*data*/, void* /*iter*/) override {
        return rr_switch_inf_->size();
    }
    inline std::pair<const void*, void*> get_switches_switch(int n, const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(&(*rr_switch_inf_)[n], nullptr);
    }

    inline int get_switch_id(const void* data, void* /*iter*/) override {
        auto* switch_obj = static_cast<const t_rr_switch_inf*>(data);
        return switch_obj - &(*rr_switch_inf_)[0];
    }
    inline const char* get_switch_name(const void* data, void* /*iter*/) override {
        return static_cast<const t_rr_switch_inf*>(data)->name;
    }

    inline std::pair<const void*, void*> init_rr_graph_switches(const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(nullptr, nullptr);
    }
    inline void finish_rr_graph_switches(const void* /*data*/, void* /*iter*/) override {
        rr_switch_inf_->shrink_to_fit();
    }

  private:
    class MetadataBind {
      public:
        MetadataBind()
            : is_node_(false)
            , is_edge_(false) {}

        void set_name(const char* name) {
            name_.assign(name);
        }
        void set_value(const char* value) {
            value_.assign(value);
        }
        void set_node_target(int inode) {
            VTR_ASSERT(!is_node_);
            VTR_ASSERT(!is_edge_);

            is_node_ = true;
            inode_ = inode;
        }

        void set_edge_target(int source_node, int sink_node, int switch_id) {
            VTR_ASSERT(!is_node_);
            VTR_ASSERT(!is_edge_);

            is_edge_ = true;
            inode_ = source_node;
            sink_node_ = sink_node;
            switch_id_ = switch_id;
        }

        void bind() {
            if (is_node_) {
                vpr::add_rr_node_metadata(inode_,
                                          std::move(name_), std::move(value_));
            } else if (is_edge_) {
                vpr::add_rr_edge_metadata(inode_, sink_node_, switch_id_,
                                          std::move(name_), std::move(value_));
            } else {
                VTR_ASSERT(is_node_ || is_edge_);
            }
        }

        void assert_clear() {
            VTR_ASSERT(name_.empty());
            VTR_ASSERT(value_.empty());
        }

        void finish_meta() {
            is_node_ = false;
            is_edge_ = false;
        }

      private:
        bool is_node_;
        bool is_edge_;
        int inode_;
        int sink_node_;
        int switch_id_;
        std::string name_;
        std::string value_;
    };

    MetadataBind metadata_bind_;

  public:
    /** Generated for complex type "meta":
     * <xs:complexType name="meta">
     *   <xs:simpleContent>
     *     <xs:extension base="xs:string">
     *       <xs:attribute name="name" type="xs:string" use="required" />
     *     </xs:extension>
     *   </xs:simpleContent>
     * </xs:complexType>
     */
    inline void set_meta_name(const char* name, const void* /*data*/, void* iter) override {
        if (iter != nullptr) {
            auto* bind = static_cast<MetadataBind*>(iter);
            bind->set_name(name);
        }
    }
    inline void set_meta_value(const char* value, const void* /*data*/, void* iter) override {
        if (iter != nullptr) {
            auto* bind = static_cast<MetadataBind*>(iter);
            bind->set_value(value);
        }
    }

    /** Generated for complex type "metadata":
     * <xs:complexType name="metadata">
     *   <xs:sequence>
     *     <xs:element maxOccurs="unbounded" name="meta" type="meta" />
     *   </xs:sequence>
     * </xs:complexType>
     */
    inline std::pair<const void*, void*> add_metadata_meta(const void* /*data*/, void* iter) override {
        if (iter != nullptr) {
            auto* bind = static_cast<MetadataBind*>(iter);
            bind->assert_clear();
            return std::make_pair(nullptr, iter);
        } else {
            return std::make_pair(nullptr, nullptr);
        }
    }
    inline void finish_metadata_meta(const void* /*data*/, void* iter) override {
        if (iter != nullptr) {
            auto* bind = static_cast<MetadataBind*>(iter);
            bind->bind();
        }
    }

    /** Generated for complex type "node_loc":
     * <xs:complexType name="node_loc">
     *   <xs:attribute name="xlow" type="xs:int" use="required" />
     *   <xs:attribute name="ylow" type="xs:int" use="required" />
     *   <xs:attribute name="xhigh" type="xs:int" use="required" />
     *   <xs:attribute name="yhigh" type="xs:int" use="required" />
     *   <xs:attribute name="side" type="loc_side" />
     *   <xs:attribute name="ptc" type="xs:int" use="required" />
     * </xs:complexType>
     */
    inline void set_node_loc_side(uxsd::enum_loc_side side, const void* data, void* /*iter*/) override {
        const int* inode = static_cast<const int*>(data);
        auto& node = (*rr_nodes_)[*inode];

        switch (side) {
            case uxsd::enum_loc_side::UXSD_INVALID:
                // node_loc.side is only expected on IPIN/OPIN
                VTR_ASSERT(!(node.type() == IPIN || node.type() == OPIN));
                break;
            case uxsd::enum_loc_side::LEFT:
                node.set_side(LEFT);
                break;
            case uxsd::enum_loc_side::RIGHT:
                node.set_side(RIGHT);
                break;
            case uxsd::enum_loc_side::TOP:
                node.set_side(TOP);
                break;
            case uxsd::enum_loc_side::BOTTOM:
                node.set_side(BOTTOM);
                break;
            default:
                VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                                "Invalid side %d", side);
        }
    }

    inline const char* get_meta_name(const void* data, void* /*iter*/) override {
        auto* meta_value = static_cast<
            const std::pair<
                std::string,
                std::vector<t_metadata_value>>*>(data);
        return meta_value->first.c_str();
    }
    inline const char* get_meta_value(const void* data, void* /*iter*/) override {
        auto* meta_value = static_cast<
            const std::pair<
                std::string,
                std::vector<t_metadata_value>>*>(data);
        VTR_ASSERT(meta_value->second.size() == 1);
        return meta_value->second[0].as_string().c_str();
    }
    inline size_t num_metadata_meta(const void* data, void* /*iter*/) override {
        auto* meta = static_cast<const struct t_metadata_dict*>(data);
        return meta->size();
    }

  private:
    t_metadata_dict::const_iterator iter_;

  public:
    inline std::pair<const void*, void*> get_metadata_meta(int n, const void* data, void* /*iter*/) override {
        auto* meta = static_cast<const struct t_metadata_dict*>(data);
        if (n == 0) {
            iter_ = meta->begin();
        } else {
            ++iter_;
        }
        VTR_ASSERT(iter_ != meta->end());
        return std::make_pair(&*iter_, nullptr);
    }

    inline int get_node_loc_ptc(const void* data, void* /*iter*/) override {
        auto* node = static_cast<const t_rr_node*>(data);
        return node->ptc_num();
    }
    inline uxsd::enum_loc_side get_node_loc_side(const void* data, void* /*iter*/) override {
        auto* node = static_cast<const t_rr_node*>(data);
        if (node->type() == IPIN || node->type() == OPIN) {
            switch (node->side()) {
                case LEFT:
                    return uxsd::enum_loc_side::LEFT;
                case RIGHT:
                    return uxsd::enum_loc_side::RIGHT;
                case TOP:
                    return uxsd::enum_loc_side::TOP;
                case BOTTOM:
                    return uxsd::enum_loc_side::BOTTOM;
                default:
                    VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                                    "Invalid side %d", node->side());
            }
        } else {
            return uxsd::enum_loc_side::UXSD_INVALID;
        }
    }
    inline int get_node_loc_xhigh(const void* data, void* /*iter*/) override {
        auto* node = static_cast<const t_rr_node*>(data);
        return node->xhigh();
    }
    inline int get_node_loc_xlow(const void* data, void* /*iter*/) override {
        auto* node = static_cast<const t_rr_node*>(data);
        return node->xlow();
    }
    inline int get_node_loc_yhigh(const void* data, void* /*iter*/) override {
        auto* node = static_cast<const t_rr_node*>(data);
        return node->yhigh();
    }
    inline int get_node_loc_ylow(const void* data, void* /*iter*/) override {
        auto* node = static_cast<const t_rr_node*>(data);
        return node->ylow();
    }
    inline float get_node_timing_C(const void* data, void* /*iter*/) override {
        auto* node = static_cast<const t_rr_node*>(data);
        return node->C();
    }
    inline float get_node_timing_R(const void* data, void* /*iter*/) override {
        auto* node = static_cast<const t_rr_node*>(data);
        return node->R();
    }
    inline int get_node_segment_segment_id(const void* data, void* /*iter*/) override {
        auto* node = static_cast<const t_rr_node*>(data);
        return (*rr_indexed_data_)[node->cost_index()].seg_index;
    }
    inline unsigned int get_node_capacity(const void* data, void* /*iter*/) override {
        auto* node = static_cast<const t_rr_node*>(data);
        return node->capacity();
    }
    inline uxsd::enum_node_direction get_node_direction(const void* data, void* /*iter*/) override {
        auto* node = static_cast<const t_rr_node*>(data);
        if (node->type() == CHANX || node->type() == CHANY) {
            switch (node->direction()) {
                case INC_DIRECTION:
                    return uxsd::enum_node_direction::INC_DIR;
                case DEC_DIRECTION:
                    return uxsd::enum_node_direction::DEC_DIR;
                case BI_DIRECTION:
                    return uxsd::enum_node_direction::BI_DIR;
                default:
                    VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                                    "Invalid direction %d", node->direction());
            }
        } else {
            return uxsd::enum_node_direction::UXSD_INVALID;
        }
    }
    inline unsigned int get_node_id(const void* data, void* /*iter*/) override {
        auto* node = static_cast<const t_rr_node*>(data);
        return node - &(*rr_nodes_)[0];
    }
    inline uxsd::enum_node_type get_node_type(const void* data, void* /*iter*/) override {
        auto* node = static_cast<const t_rr_node*>(data);
        switch (node->type()) {
            case CHANX:
                return uxsd::enum_node_type::CHANX;
            case CHANY:
                return uxsd::enum_node_type::CHANY;
            case SOURCE:
                return uxsd::enum_node_type::SOURCE;
            case SINK:
                return uxsd::enum_node_type::SINK;
            case OPIN:
                return uxsd::enum_node_type::OPIN;
            case IPIN:
                return uxsd::enum_node_type::IPIN;
            default:
                VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                                "Invalid type %d", node->type());
        }

        return uxsd::enum_node_type::UXSD_INVALID;
    }
    inline std::pair<const void*, void*> get_node_loc(const void* data, void* iter) override {
        return std::make_pair(data, iter);
    }
    inline std::pair<const void*, void*> get_node_timing(const void* data, void* iter) override {
        return std::make_pair(data, iter);
    }
    inline bool has_node_timing(const void* /*data*/, void* /*iter*/) override {
        return true;
    }
    inline std::pair<const void*, void*> get_node_segment(const void* data, void* iter) override {
        return std::make_pair(data, iter);
    }
    inline bool has_node_segment(const void* data, void* /*iter*/) override {
        auto* node = static_cast<const t_rr_node*>(data);
        return (*rr_indexed_data_)[node->cost_index()].seg_index != -1;
    }
    inline std::pair<const void*, void*> get_node_metadata(const void* data, void* iter) override {
        int node = get_node_id(data, iter);
        const auto itr = rr_node_metadata_.find(node);

        return std::make_pair(static_cast<const void*>(&itr->second), nullptr);
    }
    inline bool has_node_metadata(const void* data, void* iter) override {
        int node = get_node_id(data, iter);
        const auto itr = rr_node_metadata_.find(node);
        return itr != rr_node_metadata_.end();
    }
    inline size_t num_rr_nodes_node(const void* /*data*/, void* /*iter*/) override {
        return rr_nodes_->size();
    }
    inline std::pair<const void*, void*> get_rr_nodes_node(int n, const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(&(*rr_nodes_)[n], nullptr);
    }

  private:
    class EdgeWalker {
      public:
        size_t initialize(const std::vector<t_rr_node>* nodes) {
            nodes_ = nodes;
            num_edges_ = 0;
            current_src_inode_ = 0;
            current_edge_ = 0;
            current_idx_ = 0;

            for (const auto& node : *nodes) {
                num_edges_ += node.num_edges();
            }

            return num_edges_;
        }

        int current_src_node() const {
            return current_src_inode_;
        }
        int current_sink_node() const {
            VTR_ASSERT(current_src_inode_ < nodes_->size());
            return (*nodes_)[current_src_inode_].edge_sink_node(current_edge_);
        }
        int current_switch_id_node() const {
            VTR_ASSERT(current_src_inode_ < nodes_->size());
            return (*nodes_)[current_src_inode_].edge_switch(current_edge_);
        }

        size_t advance(int n) {
            VTR_ASSERT(current_idx_ < num_edges_);
            VTR_ASSERT(current_src_inode_ < nodes_->size());

            if (n > 0) {
                current_edge_ += 1;
            }

            if (current_edge_ >= (*nodes_)[current_src_inode_].num_edges()) {
                // Done with current_src_inode_, advance to the end of the
                // node list, or the next node with at least 1 edge.
                current_edge_ = 0;

                do {
                    current_src_inode_ += 1;

                    // This is the last edge, return now.
                    if (current_src_inode_ == nodes_->size()) {
                        VTR_ASSERT(current_idx_ + 1 == num_edges_);
                        return current_idx_++;
                    }
                } while ((*nodes_)[current_src_inode_].num_edges() < 1);
            }

            VTR_ASSERT(current_src_inode_ < nodes_->size());

            return current_idx_++;
        }

      private:
        const std::vector<t_rr_node>* nodes_;
        size_t num_edges_;
        size_t current_src_inode_;
        size_t current_edge_;
        size_t current_idx_;
    };

    EdgeWalker walker_;

  public:
    inline unsigned int get_edge_sink_node(const void* /*data*/, void* /*iter*/) override {
        return walker_.current_sink_node();
    }
    inline unsigned int get_edge_src_node(const void* /*data*/, void* /*iter*/) override {
        return walker_.current_src_node();
    }
    inline unsigned int get_edge_switch_id(const void* /*data*/, void* /*iter*/) override {
        return walker_.current_switch_id_node();
    }
    inline std::pair<const void*, void*> get_edge_metadata(const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(static_cast<const void*>(&rr_edge_metadata_.find(
                                                                             std::make_tuple(
                                                                                 walker_.current_src_node(),
                                                                                 walker_.current_sink_node(),
                                                                                 walker_.current_switch_id_node()))
                                                            ->second),
                              nullptr);
    }
    inline bool has_edge_metadata(const void* /*data*/, void* /*iter*/) override {
        return rr_edge_metadata_.find(
                   std::make_tuple(
                       walker_.current_src_node(),
                       walker_.current_sink_node(),
                       walker_.current_switch_id_node()))
               != rr_edge_metadata_.end();
    }
    inline size_t num_rr_edges_edge(const void* /*data*/, void* /*iter*/) override {
        return walker_.initialize(rr_nodes_);
    }
    inline std::pair<const void*, void*> get_rr_edges_edge(int n, const void* /*data*/, void* /*iter*/) override {
        size_t cur = walker_.advance(n);
        if ((ssize_t)cur != n) {
            VPR_FATAL_ERROR(VPR_ERROR_OTHER, "Incorrect edge index %zu != %d", cur, n);
        }
        return std::make_pair(nullptr, nullptr);
    }

    /** Generated for complex type "node_timing":
     * <xs:complexType name="node_timing">
     *   <xs:attribute name="R" type="xs:float" use="required" />
     *   <xs:attribute name="C" type="xs:float" use="required" />
     * </xs:complexType>
     */

    /** Generated for complex type "node_segment":
     * <xs:complexType name="node_segment">
     *   <xs:attribute name="segment_id" type="xs:int" use="required" />
     * </xs:complexType>
     */

    /** Generated for complex type "node":
     * <xs:complexType name="node">
     *   <xs:all>
     *     <xs:element name="loc" type="node_loc" />
     *     <xs:element minOccurs="0" name="timing" type="node_timing" />
     *     <xs:element minOccurs="0" name="segment" type="node_segment" />
     *     <xs:element minOccurs="0" name="metadata" type="metadata" />
     *   </xs:all>
     *   <xs:attribute name="id" type="xs:unsignedInt" use="required" />
     *   <xs:attribute name="type" type="node_type" use="required" />
     *   <xs:attribute name="direction" type="node_direction" />
     *   <xs:attribute name="capacity" type="xs:unsignedInt" use="required" />
     * </xs:complexType>
     */
    inline void set_node_direction(uxsd::enum_node_direction direction, const void* data, void* /*iter*/) override {
        const int* inode = static_cast<const int*>(data);
        auto& node = (*rr_nodes_)[*inode];
        switch (direction) {
            case uxsd::enum_node_direction::UXSD_INVALID:
                VTR_ASSERT(!(node.type() == CHANX || node.type() == CHANY));
                return;
            case uxsd::enum_node_direction::INC_DIR:
                node.set_direction(INC_DIRECTION);
                break;
            case uxsd::enum_node_direction::DEC_DIR:
                node.set_direction(DEC_DIRECTION);
                break;
            case uxsd::enum_node_direction::BI_DIR:
                node.set_direction(BI_DIRECTION);
                break;
            default:
                VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                                "Invalid node direction %d", direction);
        }
    }

    inline std::pair<const void*, void*> init_node_loc(const void* data, void* iter, int ptc, int xhigh, int xlow, int yhigh, int ylow) override {
        const int* inode = static_cast<const int*>(data);
        auto& node = (*rr_nodes_)[*inode];

        node.set_coordinates(xlow, ylow, xhigh, yhigh);
        node.set_ptc_num(ptc);
        return std::make_pair(data, iter);
    }
    inline void finish_node_loc(const void* /*data*/, void* /*iter*/) override {}

    inline std::pair<const void*, void*> init_node_timing(const void* data, void* /*iter*/, float C, float R) override {
        const int* inode = static_cast<const int*>(data);
        auto& node = (*rr_nodes_)[*inode];
        node.set_rc_index(find_create_rr_rc_data(R, C));
        return std::make_pair(nullptr, nullptr);
    }
    inline void finish_node_timing(const void* /*data*/, void* /*iter*/) override {}
    inline std::pair<const void*, void*> init_node_segment(const void* data, void* /*iter*/, int segment_id) override {
        if (segment_id > (ssize_t)segment_inf_.size()) {
            VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                            "Specified segment %d is larger than number of known segments %zu",
                            segment_inf_.size());
        }

        const int* inode = static_cast<const int*>(data);
        auto& node = (*rr_nodes_)[*inode];
        if (is_global_graph_) {
            node.set_cost_index(0);
        } else if (node.type() == CHANX) {
            node.set_cost_index(CHANX_COST_INDEX_START + segment_id);
            seg_index_[node.cost_index()] = segment_id;
        } else if (node.type() == CHANY) {
            node.set_cost_index(CHANX_COST_INDEX_START + segment_inf_.size() + segment_id);
            seg_index_[node.cost_index()] = segment_id;
        }
        return std::make_pair(nullptr, nullptr);
    }
    inline void finish_node_segment(const void* /*data*/, void* /*iter*/) override {}

    inline std::pair<const void*, void*> init_node_metadata(const void* /*data*/, void* /*iter*/) override {
        metadata_bind_.assert_clear();
        return std::make_pair(nullptr, &metadata_bind_);
    }
    inline void finish_node_metadata(const void* /*data*/, void* /*iter*/) override {
    }

    /** Generated for complex type "rr_nodes":
     * <xs:complexType name="rr_nodes">
     *   <xs:choice maxOccurs="unbounded">
     *     <xs:element name="node" type="node" />
     *   </xs:choice>
     * </xs:complexType>
     */
    inline std::pair<const void*, void*> add_rr_nodes_node(const void* /*data*/, void* /*iter*/, unsigned int capacity, unsigned int id, uxsd::enum_node_type type) override {
        make_room_in_vector(rr_nodes_, id);
        auto& node = (*rr_nodes_)[id];
        metadata_bind_.set_node_target(id);

        node.set_capacity(capacity);

        switch (type) {
            case uxsd::enum_node_type::CHANX:
                node.set_type(CHANX);
                break;
            case uxsd::enum_node_type::CHANY:
                node.set_type(CHANY);
                break;
            case uxsd::enum_node_type::SOURCE:
                node.set_type(SOURCE);
                node.set_cost_index(SOURCE_COST_INDEX);
                break;
            case uxsd::enum_node_type::SINK:
                node.set_type(SINK);
                node.set_cost_index(SINK_COST_INDEX);
                break;
            case uxsd::enum_node_type::OPIN:
                node.set_type(OPIN);
                node.set_cost_index(OPIN_COST_INDEX);
                break;
            case uxsd::enum_node_type::IPIN:
                node.set_type(IPIN);
                node.set_cost_index(IPIN_COST_INDEX);
                break;
            default:
                VPR_FATAL_ERROR(
                    VPR_ERROR_OTHER,
                    "Invalid node type %d",
                    type);
        }

        node.set_rc_index(find_create_rr_rc_data(0, 0));

        id_ = id;
        return std::make_pair(&id_, nullptr);
    }
    inline void finish_rr_nodes_node(const void* data, void* /*iter*/) override {
        const int* inode = static_cast<const int*>(data);
        auto& node = (*rr_nodes_)[*inode];
        metadata_bind_.finish_meta();
        node.set_num_edges(0);
    }

    inline std::pair<const void*, void*> init_rr_graph_rr_nodes(const void* /*data*/, void* /*iter*/) override {
        rr_nodes_->clear();
        seg_index_.resize(CHANX_COST_INDEX_START + 2 * segment_inf_.size(), -1);
        return std::make_pair(nullptr, nullptr);
    }
    inline void finish_rr_graph_rr_nodes(const void* /*data*/, void* /*iter*/) override {
        rr_nodes_->shrink_to_fit();
    }

    /** Generated for complex type "edge":
     * <xs:complexType name="edge">
     *   <xs:all>
     *     <xs:element minOccurs="0" name="metadata" type="metadata" />
     *   </xs:all>
     *   <xs:attribute name="src_node" type="xs:unsignedInt" use="required" />
     *   <xs:attribute name="sink_node" type="xs:unsignedInt" use="required" />
     *   <xs:attribute name="switch_id" type="xs:unsignedInt" use="required" />
     * </xs:complexType>
     */
    inline std::pair<const void*, void*> init_edge_metadata(const void* data, void* iter) override {
        if (read_edge_metadata_) {
            return std::make_pair(data, iter);
        } else {
            return std::make_pair(nullptr, nullptr);
        }
    }
    inline void finish_edge_metadata(const void* /*data*/, void* /*iter*/) override {
    }

    /** Generated for complex type "channels":
     * <xs:complexType name="channels">
     *   <xs:sequence>
     *     <xs:element name="channel" type="channel" />
     *     <xs:element maxOccurs="unbounded" name="x_list" type="x_list" />
     *     <xs:element maxOccurs="unbounded" name="y_list" type="y_list" />
     *   </xs:sequence>
     * </xs:complexType>
     */
    inline int get_channel_chan_width_max(const void* /*data*/, void* /*iter*/) override {
        return chan_width_->max;
    }
    inline int get_channel_x_max(const void* /*data*/, void* /*iter*/) override {
        return chan_width_->x_max;
    }
    inline int get_channel_x_min(const void* /*data*/, void* /*iter*/) override {
        return chan_width_->x_min;
    }
    inline int get_channel_y_max(const void* /*data*/, void* /*iter*/) override {
        return chan_width_->y_max;
    }
    inline int get_channel_y_min(const void* /*data*/, void* /*iter*/) override {
        return chan_width_->y_min;
    }
    inline std::pair<const void*, void*> init_channels_channel(const void* /*data*/, void* /*iter*/, int chan_width_max, int x_max, int x_min, int y_max, int y_min) override {
        chan_width_->max = chan_width_max;
        chan_width_->x_min = x_min;
        chan_width_->y_min = y_min;
        chan_width_->x_max = x_max;
        chan_width_->y_max = y_max;
        chan_width_->x_list.resize(grid_.height());
        chan_width_->y_list.resize(grid_.width());
        return std::make_pair(nullptr, nullptr);
    }
    inline void finish_channels_channel(const void* /*data*/, void* /*iter*/) override {
    }
    inline std::pair<const void*, void*> get_channels_channel(const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(nullptr, nullptr);
    }

    inline std::pair<const void*, void*> add_channels_x_list(const void* /*data*/, void* /*iter*/, unsigned int index, int info) override {
        if (index >= chan_width_->x_list.size()) {
            VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                            "index %d on x_list exceeds x_list size %u",
                            index, chan_width_->x_list.size());
        }
        chan_width_->x_list[index] = info;
        return std::make_pair(nullptr, nullptr);
    }
    inline void finish_channels_x_list(const void* /*data*/, void* /*iter*/) override {}

    inline unsigned int get_x_list_index(const void* data, void* /*iter*/) override {
        return static_cast<const int*>(data) - &chan_width_->x_list[0];
    }
    inline int get_x_list_info(const void* data, void* /*iter*/) override {
        return *static_cast<const int*>(data);
    }
    inline size_t num_channels_x_list(const void* /*data*/, void* /*iter*/) override {
        return chan_width_->x_list.size();
    }
    inline std::pair<const void*, void*> get_channels_x_list(int n, const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(&chan_width_->x_list[n], nullptr);
    }

    inline std::pair<const void*, void*> add_channels_y_list(const void* /*data*/, void* /*iter*/, unsigned int index, int info) override {
        if (index >= chan_width_->y_list.size()) {
            VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                            "index %d on y_list exceeds y_list size %u",
                            index, chan_width_->y_list.size());
        }
        chan_width_->y_list[index] = info;
        return std::make_pair(nullptr, nullptr);
    }
    inline void finish_channels_y_list(const void* /*data*/, void* /*iter*/) override {}

    inline unsigned int get_y_list_index(const void* data, void* /*iter*/) override {
        return static_cast<const int*>(data) - &chan_width_->y_list[0];
    }
    inline int get_y_list_info(const void* data, void* /*iter*/) override {
        return *static_cast<const int*>(data);
    }
    inline size_t num_channels_y_list(const void* /*data*/, void* /*iter*/) override {
        return chan_width_->y_list.size();
    }
    inline std::pair<const void*, void*> get_channels_y_list(int n, const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(&chan_width_->y_list[n], nullptr);
    }

    inline std::pair<const void*, void*> init_rr_graph_channels(const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(nullptr, nullptr);
    }
    inline void finish_rr_graph_channels(const void* /*data*/, void* /*iter*/) override {
    }

  private:
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> edges_;

  public:
    /** Generated for complex type "rr_edges":
     * <xs:complexType name="rr_edges">
     *   <xs:choice maxOccurs="unbounded">
     *     <xs:element name="edge" type="edge" />
     *   </xs:choice>
     * </xs:complexType>
     */
    inline std::pair<const void*, void*> add_rr_edges_edge(const void* /*data*/, void* /*iter*/, unsigned int sink_node, unsigned int src_node, unsigned int switch_id) override {
        if (src_node >= rr_nodes_->size()) {
            VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                            "source_node %d is larger than rr_nodes.size() %d",
                            src_node, rr_nodes_->size());
        }

        metadata_bind_.set_edge_target(src_node, sink_node, switch_id);
        edges_.push_back(std::make_tuple(src_node, sink_node, switch_id));
        return std::make_pair(nullptr, &metadata_bind_);
    }
    inline void finish_rr_edges_edge(const void* /*data*/, void* /*iter*/) override {
        metadata_bind_.finish_meta();
    }
    inline std::pair<const void*, void*> init_rr_graph_rr_edges(const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(nullptr, nullptr);
    }
    inline void finish_rr_graph_rr_edges(const void* /*data*/, void* /*iter*/) override {
        std::vector<size_t> num_edges_for_node(rr_nodes_->size());
        for (const auto& edge : edges_) {
            num_edges_for_node[std::get<0>(edge)]++;
        }

        for (size_t inode = 0; inode < rr_nodes_->size(); inode++) {
            if (num_edges_for_node[inode] > std::numeric_limits<t_edge_size>::max()) {
                VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                                "source node %d edge count %d is too high",
                                inode, num_edges_for_node[inode]);
            }
            (*rr_nodes_)[inode].set_num_edges(num_edges_for_node[inode]);
            num_edges_for_node[inode] = 0;
        }

        /*initialize a vector that keeps track of the number of wire to ipin switches
         * There should be only one wire to ipin switch. In case there are more, make sure to
         * store the most frequent switch */
        std::vector<int> count_for_wire_to_ipin_switches;
        count_for_wire_to_ipin_switches.resize(rr_switch_inf_->size(), 0);
        //first is index, second is count
        std::pair<int, int> most_frequent_switch(-1, 0);

        for (const auto& edge : edges_) {
            auto source_node = std::get<0>(edge);
            auto sink_node = std::get<1>(edge);
            auto switch_id = std::get<2>(edge);

            if (sink_node >= rr_nodes_->size()) {
                VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                                "sink_node %u is larger than rr_nodes.size() %zu",
                                sink_node, rr_nodes_->size());
            }

            if (switch_id >= rr_switch_inf_->size()) {
                VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                                "switch_id %u is larger than num_rr_switches %zu",
                                switch_id, rr_switch_inf_->size());
            }

            auto& node = (*rr_nodes_)[source_node];

            /*Keeps track of the number of the specific type of switch that connects a wire to an ipin
             * use the pair data structure to keep the maximum*/
            if (node.type() == CHANX || node.type() == CHANY) {
                if ((*rr_nodes_)[sink_node].type() == IPIN) {
                    count_for_wire_to_ipin_switches[switch_id]++;
                    if (count_for_wire_to_ipin_switches[switch_id] > most_frequent_switch.second) {
                        most_frequent_switch.first = switch_id;
                        most_frequent_switch.second = count_for_wire_to_ipin_switches[switch_id];
                    }
                }
            }

            //set edge in correct rr_node data structure
            node.set_edge_sink_node(num_edges_for_node[source_node], sink_node);
            node.set_edge_switch(num_edges_for_node[source_node], switch_id);
            num_edges_for_node[source_node]++;
        }

        edges_.clear();
        edges_.shrink_to_fit();
        wire_to_rr_ipin_switch_ = most_frequent_switch.first;
    }

    /** Generated for complex type "rr_graph":
     * <xs:complexType xmlns:xs="http://www.w3.org/2001/XMLSchema">
     *     <xs:all>
     *       <xs:element name="channels" type="channels" />
     *       <xs:element name="switches" type="switches" />
     *       <xs:element name="segments" type="segments" />
     *       <xs:element name="block_types" type="block_types" />
     *       <xs:element name="grid" type="grid_locs" />
     *       <xs:element name="rr_nodes" type="rr_nodes" />
     *       <xs:element name="rr_edges" type="rr_edges" />
     *     </xs:all>
     *     <xs:attribute name="tool_name" type="xs:string" />
     *     <xs:attribute name="tool_version" type="xs:string" />
     *     <xs:attribute name="tool_comment" type="xs:string" />
     *   </xs:complexType>
     */
    inline void set_rr_graph_tool_comment(const char* tool_comment, const void* /*data*/, void* /*iter*/) override {
        std::string correct_string = "Generated from arch file ";
        correct_string += get_arch_file_name();
        if (correct_string != tool_comment) {
            VTR_LOG("\n");
            VTR_LOG_WARN("This RR graph file is based on %s while your input architecture file is %s compatability issues may arise\n",
                         get_arch_file_name(), tool_comment);
            VTR_LOG("\n");
        }
    }
    inline void set_rr_graph_tool_name(const char* /*tool_name*/, const void* /*data*/, void* /*iter*/) override {
    }
    inline void set_rr_graph_tool_version(const char* tool_version, const void* /*data*/, void* /*iter*/) override {
        if (strcmp(tool_version, vtr::VERSION) != 0) {
            VTR_LOG("\n");
            VTR_LOG_WARN("This architecture version is for VPR %s while your current VPR version is %s compatability issues may arise\n",
                         vtr::VERSION, tool_version);
            VTR_LOG("\n");
        }
    }

    /** Generated for complex type "segment_timing":
     * <xs:complexType name="segment_timing">
     *   <xs:attribute name="R_per_meter" type="xs:float" />
     *   <xs:attribute name="C_per_meter" type="xs:float" />
     * </xs:complexType>
     */
    inline float get_segment_timing_C_per_meter(const void* data, void* /*iter*/) override {
        const auto* segment = static_cast<const t_segment_inf*>(data);
        return segment->Cmetal;
    }
    inline void set_segment_timing_C_per_meter(float C_per_meter, const void* data, void* /*iter*/) override {
        const auto* segment = static_cast<const t_segment_inf*>(data);
        if (segment->Cmetal != C_per_meter) {
            VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                            "Architecture file does not match RR graph's segment C_per_meter");
        }
    }
    inline float get_segment_timing_R_per_meter(const void* data, void* /*iter*/) override {
        const auto* segment = static_cast<const t_segment_inf*>(data);
        return segment->Rmetal;
    }
    inline void set_segment_timing_R_per_meter(float R_per_meter, const void* data, void* /*iter*/) override {
        const auto* segment = static_cast<const t_segment_inf*>(data);
        if (segment->Rmetal != R_per_meter) {
            VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                            "Architecture file does not match RR graph's segment R_per_meter");
        }
    }

    /** Generated for complex type "segment":
     * <xs:complexType name="segment">
     *   <xs:all>
     *     <xs:element minOccurs="0" name="timing" type="segment_timing" />
     *   </xs:all>
     *   <xs:attribute name="id" type="xs:int" use="required" />
     *   <xs:attribute name="name" type="xs:string" use="required" />
     * </xs:complexType>
     */
    inline int get_segment_id(const void* data, void* /*iter*/) override {
        const auto* segment = static_cast<const t_segment_inf*>(data);
        return segment - &segment_inf_.at(0);
    }
    inline const char* get_segment_name(const void* data, void* /*iter*/) override {
        const auto* segment = static_cast<const t_segment_inf*>(data);
        return segment->name.c_str();
    }
    inline void set_segment_name(const char* name, const void* data, void* /*iter*/) override {
        const auto* segment = static_cast<const t_segment_inf*>(data);
        if (segment->name != name) {
            VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                            "Architecture file does not match RR graph's segment name: arch uses %s, RR graph uses %s",
                            segment->name.c_str(), name);
        }
    }
    inline std::pair<const void*, void*> init_segment_timing(const void* data, void* iter) override {
        return std::make_pair(data, iter);
    }
    inline void finish_segment_timing(const void* /*data*/, void* /*iter*/) override {}

    inline std::pair<const void*, void*> get_segment_timing(const void* data, void* iter) override {
        return std::make_pair(data, iter);
    }
    inline bool has_segment_timing(const void* /*data*/, void* /*iter*/) override {
        return true;
    }

    /** Generated for complex type "segments":
     * <xs:complexType name="segments">
     *   <xs:sequence>
     *     <xs:element maxOccurs="unbounded" name="segment" type="segment" />
     *   </xs:sequence>
     * </xs:complexType>
     */
    inline std::pair<const void*, void*> add_segments_segment(const void* /*data*/, void* /*iter*/, int id) override {
        return std::make_pair(static_cast<const void*>(&segment_inf_.at(id)), nullptr);
    }
    inline void finish_segments_segment(const void* /*data*/, void* /*iter*/) override {}
    inline size_t num_segments_segment(const void* /*data*/, void* /*iter*/) override {
        return segment_inf_.size();
    }
    inline std::pair<const void*, void*> get_segments_segment(int n, const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(static_cast<const void*>(&segment_inf_.at(n)), nullptr);
    }

    inline std::pair<const void*, void*> init_rr_graph_segments(const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(nullptr, nullptr);
    }
    inline void finish_rr_graph_segments(const void* /*data*/, void* /*iter*/) override {
    }

  private:
    class PinClassChecker {
      public:
        PinClassChecker()
            : tile_type_(nullptr)
            , pin_count_(0) {}

        void set_tile(const t_physical_tile_type* tile_type) {
            VTR_ASSERT(tile_type_ == nullptr);
            tile_type_ = tile_type;
            pin_class_index_ = 0;
        }

        void next_pin_class() {
            VTR_ASSERT(tile_type_ != nullptr);
            if (get_current_pin_class().num_pins != (ssize_t)pin_count_) {
                VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                                "Incorrect number of pins (%zu != %u) in %zu pin_class in block %s",
                                pin_count_, get_current_pin_class().num_pins,
                                pin_class_index_, tile_type_->name);
            }
            pin_class_index_ += 1;
            current_pin_name_.clear();
            pin_count_ = 0;
        }

        bool done_with_pin_classes() const {
            VTR_ASSERT(tile_type_ != nullptr);
            return (ssize_t)pin_class_index_ == tile_type_->num_class;
        }

        const t_class& get_current_pin_class() const {
            VTR_ASSERT(tile_type_ != nullptr);
            VTR_ASSERT((ssize_t)pin_class_index_ < tile_type_->num_class);
            return tile_type_->class_inf[pin_class_index_];
        }

        void clear() {
            tile_type_ = nullptr;
        }

        void set_pin_ptc(int ptc) {
            int index = pin_index_by_num(get_current_pin_class(), ptc);
            if (index < 0) {
                VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                                "Architecture file does not match RR graph's block pin list: invalid ptc for pin class");
            }

            pin_count_ += 1;
            current_pin_name_ = block_type_pin_index_to_name(
                tile_type_,
                get_current_pin_class().pinlist[index]);
        }
        const std::string get_pin_name() const {
            VTR_ASSERT(tile_type_ != nullptr);
            VTR_ASSERT(!current_pin_name_.empty());
            return current_pin_name_;
        }

      private:
        static int pin_index_by_num(const t_class& class_inf, int num) {
            for (int index = 0; index < class_inf.num_pins; ++index) {
                if (num == class_inf.pinlist[index]) {
                    return index;
                }
            }
            return -1;
        }
        const t_physical_tile_type* tile_type_;
        size_t pin_class_index_;
        std::string current_pin_name_;
        size_t pin_count_;
    };

    PinClassChecker pin_class_checker_;

  public:
    /** Generated for complex type "pin":
     * <xs:complexType name="pin">
     *   <xs:simpleContent>
     *     <xs:extension base="xs:string">
     *       <xs:attribute name="ptc" type="xs:int" use="required" />
     *     </xs:extension>
     *   </xs:simpleContent>
     * </xs:complexType>
     */
    inline void set_pin_value(const char* value, const void* /*data*/, void* /*iter*/) override {
        if (pin_class_checker_.get_pin_name() != value) {
            VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                            "Architecture file does not match RR graph's block pin list");
        }
    }

  private:
    const t_physical_tile_type* current_tile_;
    std::string pin_value_;

  public:
    /** Generated for complex type "pin_class":
     * <xs:complexType name="pin_class">
     *   <xs:sequence>
     *     <xs:element maxOccurs="unbounded" name="pin" type="pin" />
     *   </xs:sequence>
     *   <xs:attribute name="type" type="pin_type" use="required" />
     * </xs:complexType>
     */
    inline std::pair<const void*, void*> add_pin_class_pin(const void* /*data*/, void* /*iter*/, int ptc) override {
        pin_class_checker_.set_pin_ptc(ptc);
        return std::make_pair(nullptr, nullptr);
    }
    inline void finish_pin_class_pin(const void* /*data*/, void* /*iter*/) override {
    }

    inline int get_pin_ptc(const void* data, void* /*iter*/) override {
        return *static_cast<const int*>(data);
    }
    inline const char* get_pin_value(const void* data, void* /*iter*/) override {
        int pin = *static_cast<const int*>(data);
        pin_value_ = block_type_pin_index_to_name(current_tile_, pin);
        return pin_value_.c_str();
    }
    inline uxsd::enum_pin_type get_pin_class_type(const void* data, void* /*iter*/) override {
        auto type = static_cast<const t_class*>(data)->type;
        switch (type) {
            case OPEN:
                return uxsd::enum_pin_type::OPEN;
            case DRIVER:
                return uxsd::enum_pin_type::OUTPUT;
            case RECEIVER:
                return uxsd::enum_pin_type::INPUT;
            default:
                VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                                "Unknown pin class type %d", type);
        }
        return uxsd::enum_pin_type::UXSD_INVALID;
    }
    inline size_t num_pin_class_pin(const void* data, void* /*iter*/) override {
        return static_cast<const t_class*>(data)->num_pins;
    }
    inline std::pair<const void*, void*> get_pin_class_pin(int n, const void* data, void* /*iter*/) override {
        auto* class_inf = static_cast<const t_class*>(data);
        return std::make_pair(&class_inf->pinlist[n], nullptr);
    }
    inline int get_block_type_height(const void* data, void* /*iter*/) override {
        return static_cast<const t_physical_tile_type*>(data)->height;
    }
    inline int get_block_type_id(const void* data, void* /*iter*/) override {
        return static_cast<const t_physical_tile_type*>(data)->index;
    }
    inline const char* get_block_type_name(const void* data, void* /*iter*/) override {
        return static_cast<const t_physical_tile_type*>(data)->name;
    }
    inline int get_block_type_width(const void* data, void* /*iter*/) override {
        return static_cast<const t_physical_tile_type*>(data)->width;
    }
    inline size_t num_block_type_pin_class(const void* data, void* /*iter*/) override {
        return static_cast<const t_physical_tile_type*>(data)->num_class;
    }
    inline std::pair<const void*, void*> get_block_type_pin_class(int n, const void* data, void* /*iter*/) override {
        return std::make_pair(&static_cast<const t_physical_tile_type*>(data)->class_inf[n], nullptr);
    }
    inline size_t num_block_types_block_type(const void* /*data*/, void* /*iter*/) override {
        return physical_tile_types_.size();
    }
    inline std::pair<const void*, void*> get_block_types_block_type(int n, const void* /*data*/, void* /*iter*/) override {
        current_tile_ = &physical_tile_types_[n];
        return std::make_pair(static_cast<const void*>(&physical_tile_types_[n]), nullptr);
    }

    /** Generated for complex type "block_type":
     * <xs:complexType name="block_type">
     *   <xs:sequence>
     *     <xs:element maxOccurs="unbounded" minOccurs="0" name="pin_class" type="pin_class" />
     *   </xs:sequence>
     *   <xs:attribute name="id" type="xs:int" use="required" />
     *   <xs:attribute name="name" type="xs:string" use="required" />
     *   <xs:attribute name="width" type="xs:int" use="required" />
     *   <xs:attribute name="height" type="xs:int" use="required" />
     * </xs:complexType>
     */
    inline void set_block_type_name(const char* name, const void* data, void* /*iter*/) override {
        const auto& block_info = *static_cast<const t_physical_tile_type*>(data);
        if (strcmp(block_info.name, name) != 0) {
            VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                            "Architecture file does not match RR graph's block name: arch uses name %s, RR graph uses name %s",
                            block_info.name, name);
        }
    }

    inline std::pair<const void*, void*> add_block_type_pin_class(const void* /*data*/, void* /*iter*/, uxsd::enum_pin_type type) override {
        e_pin_type pin_type;

        switch (type) {
            case uxsd::enum_pin_type::OPEN:
                pin_type = OPEN;
                break;
            case uxsd::enum_pin_type::OUTPUT:
                pin_type = DRIVER;
                break;
            case uxsd::enum_pin_type::INPUT:
                pin_type = RECEIVER;
                break;
            default:
                VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                                "Unknown pin class type %d", type);
        }

        if (pin_class_checker_.get_current_pin_class().type != pin_type) {
            VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                            "Architecture file does not match RR graph's block type");
        }

        return std::make_pair(nullptr, nullptr);
    }

    inline void finish_block_type_pin_class(const void* /*data*/, void* /*iter*/) override {
        pin_class_checker_.next_pin_class();
    }

    /** Generated for complex type "block_types":
     * <xs:complexType name="block_types">
     *   <xs:sequence>
     *     <xs:element maxOccurs="unbounded" name="block_type" type="block_type" />
     *   </xs:sequence>
     * </xs:complexType>
     */
    inline std::pair<const void*, void*> add_block_types_block_type(const void* /*data*/, void* /*iter*/, int height, int id, int width) override {
        const auto& block_info = physical_tile_types_.at(id);
        if (block_info.width != width) {
            VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                            "Architecture file does not match RR graph's block width");
        }
        if (block_info.height != height) {
            VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                            "Architecture file does not match RR graph's block height");
        }

        pin_class_checker_.set_tile(&block_info);
        return std::make_pair(static_cast<const void*>(&block_info), nullptr);
    }
    inline void finish_block_types_block_type(const void* /*data*/, void* /*iter*/) override {
        VTR_ASSERT(pin_class_checker_.done_with_pin_classes());
        pin_class_checker_.clear();
    }
    inline std::pair<const void*, void*> init_rr_graph_block_types(const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(nullptr, nullptr);
    }
    inline void finish_rr_graph_block_types(const void* /*data*/, void* /*iter*/) override {
    }

    /** Generated for complex type "grid_loc":
     * <xs:complexType name="grid_loc">
     *   <xs:attribute name="x" type="xs:int" use="required" />
     *   <xs:attribute name="y" type="xs:int" use="required" />
     *   <xs:attribute name="block_type_id" type="xs:int" use="required" />
     *   <xs:attribute name="width_offset" type="xs:int" use="required" />
     *   <xs:attribute name="height_offset" type="xs:int" use="required" />
     * </xs:complexType>
     */
    /** Generated for complex type "grid_locs":
     * <xs:complexType name="grid_locs">
     *   <xs:sequence>
     *     <xs:element maxOccurs="unbounded" name="grid_loc" type="grid_loc" />
     *   </xs:sequence>
     * </xs:complexType>
     */
    inline std::pair<const void*, void*> add_grid_locs_grid_loc(const void* /*data*/, void* /*iter*/, int block_type_id, int height_offset, int width_offset, int x, int y) override {
        const t_grid_tile& grid_tile = grid_[x][y];

        if (grid_tile.type->index != block_type_id) {
            VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                            "Architecture file does not match RR graph's block_type_id at (%d, %d): arch used ID %d, RR graph used ID %d.", x, y,
                            (grid_tile.type->index), block_type_id);
        }
        if (grid_tile.width_offset != width_offset) {
            VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                            "Architecture file does not match RR graph's width_offset at (%d, %d)", x, y);
        }

        if (grid_tile.height_offset != height_offset) {
            VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                            "Architecture file does not match RR graph's height_offset at (%d, %d)", x, y);
        }
        return std::make_pair(nullptr, nullptr);
    }
    inline void finish_grid_locs_grid_loc(const void* /*data*/, void* /*iter*/) override {}

    inline std::pair<const void*, void*> init_rr_graph_grid(const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(nullptr, nullptr);
    }
    inline void finish_rr_graph_grid(const void* /*data*/, void* /*iter*/) override {
    }

    inline int get_grid_loc_block_type_id(const void* data, void* /*iter*/) override {
        auto* grid_loc = static_cast<const t_grid_tile*>(data);
        return grid_loc->type->index;
    }
    inline int get_grid_loc_height_offset(const void* data, void* /*iter*/) override {
        auto* grid_loc = static_cast<const t_grid_tile*>(data);
        return grid_loc->height_offset;
    }
    inline int get_grid_loc_width_offset(const void* data, void* /*iter*/) override {
        auto* grid_loc = static_cast<const t_grid_tile*>(data);
        return grid_loc->width_offset;
    }
    inline int get_grid_loc_x(const void* data, void* /*iter*/) override {
        auto* grid_loc = static_cast<const t_grid_tile*>(data);
        auto diff = grid_loc - &grid_.matrix().get(0);

        return diff / grid_.matrix().dim_size(1);
    }
    inline int get_grid_loc_y(const void* data, void* /*iter*/) override {
        auto* grid_loc = static_cast<const t_grid_tile*>(data);
        auto diff = grid_loc - &grid_.matrix().get(0);

        return diff % grid_.matrix().dim_size(1);
    }
    inline size_t num_grid_locs_grid_loc(const void* /*data*/, void* /*iter*/) override {
        return grid_.matrix().size();
    }
    inline std::pair<const void*, void*> get_grid_locs_grid_loc(int n, const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(static_cast<const void*>(&grid_.matrix().get(n)), nullptr);
    }

  private:
    std::string tool_comment_;

  public:
    inline const char* get_rr_graph_tool_comment(const void* /*data*/, void* /*iter*/) override {
        tool_comment_.assign("Generated from arch file ");
        tool_comment_ += get_arch_file_name();
        return tool_comment_.c_str();
    }
    inline const char* get_rr_graph_tool_name(const void* /*data*/, void* /*iter*/) override {
        return "vpr";
    }
    inline const char* get_rr_graph_tool_version(const void* /*data*/, void* /*iter*/) override {
        return vtr::VERSION;
    }
    inline std::pair<const void*, void*> get_rr_graph_channels(const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(nullptr, nullptr);
    }
    inline std::pair<const void*, void*> get_rr_graph_switches(const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(nullptr, nullptr);
    }
    inline std::pair<const void*, void*> get_rr_graph_segments(const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(nullptr, nullptr);
    }
    inline std::pair<const void*, void*> get_rr_graph_block_types(const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(nullptr, nullptr);
    }
    inline std::pair<const void*, void*> get_rr_graph_grid(const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(nullptr, nullptr);
    }
    inline std::pair<const void*, void*> get_rr_graph_rr_nodes(const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(nullptr, nullptr);
    }
    inline std::pair<const void*, void*> get_rr_graph_rr_edges(const void* /*data*/, void* /*iter*/) override {
        return std::make_pair(nullptr, nullptr);
    }

    void finish_rr_graph_read(
        const t_graph_type graph_type,
        const enum e_base_cost_type base_cost_type,
        int* wire_to_rr_ipin_switch,
        bool do_check_rr_graph,
        const char* read_rr_graph_name,
        std::string* read_rr_graph_filename) {
        *wire_to_rr_ipin_switch = wire_to_rr_ipin_switch_;

        //Partition the rr graph edges for efficient access to configurable/non-configurable
        //edge subsets. Must be done after RR switches have been allocated
        partition_rr_graph_edges(rr_nodes_);

        process_rr_node_indices();

        init_fan_in(*rr_nodes_, rr_nodes_->size());

        int max_chan_width = (is_global_graph_ ? 1 : chan_width_->max);
        alloc_and_load_rr_indexed_data(
            segment_inf_,
            *rr_node_indices_,
            max_chan_width,
            *wire_to_rr_ipin_switch,
            base_cost_type);

        VTR_ASSERT(rr_indexed_data_->size() == seg_index_.size());
        for (size_t i = 0; i < seg_index_.size(); ++i) {
            (*rr_indexed_data_)[i].seg_index = seg_index_[i];
        }

        read_rr_graph_filename->assign(read_rr_graph_name);

        if (do_check_rr_graph) {
            check_rr_graph(graph_type, grid_, physical_tile_types_);
        }
    }

  private:
    /*Allocates and load the rr_node look up table. SINK and SOURCE, IPIN and OPIN
     *share the same look up table. CHANX and CHANY have individual look ups */
    void process_rr_node_indices() {
        /* Alloc the lookup table */
        auto& indices = *rr_node_indices_;

        indices.resize(NUM_RR_TYPES);

        typedef struct max_ptc {
            short chanx_max_ptc = 0;
            short chany_max_ptc = 0;
        } t_max_ptc;

        /*
         * Local multi-dimensional vector to hold max_ptc for every coordinate.
         * It has same height and width as CHANY and CHANX are inverted
         */
        vtr::Matrix<t_max_ptc> coordinates_max_ptc; /* [x][y] */
        size_t max_coord_size = std::max(grid_.width(), grid_.height());
        coordinates_max_ptc.resize({max_coord_size, max_coord_size}, t_max_ptc());

        /* Alloc the lookup table */
        for (t_rr_type rr_type : RR_TYPES) {
            if (rr_type == CHANX) {
                indices[rr_type].resize(grid_.height());
                for (size_t y = 0; y < grid_.height(); ++y) {
                    indices[rr_type][y].resize(grid_.width());
                    for (size_t x = 0; x < grid_.width(); ++x) {
                        indices[rr_type][y][x].resize(NUM_SIDES);
                    }
                }
            } else {
                indices[rr_type].resize(grid_.width());
                for (size_t x = 0; x < grid_.width(); ++x) {
                    indices[rr_type][x].resize(grid_.height());
                    for (size_t y = 0; y < grid_.height(); ++y) {
                        indices[rr_type][x][y].resize(NUM_SIDES);
                    }
                }
            }
        }

        /*
         * Add the correct node into the vector
         * For CHANX and CHANY no node is added yet, but the maximum ptc is counted for each
         * x/y location. This is needed later to add the correct node corresponding to CHANX
         * and CHANY.
         *
         * Note that CHANX and CHANY 's x and y are swapped due to the chan and seg convention.
         */
        for (size_t inode = 0; inode < rr_nodes_->size(); inode++) {
            auto& node = (*rr_nodes_)[inode];
            if (node.type() == SOURCE || node.type() == SINK) {
                for (int ix = node.xlow(); ix <= node.xhigh(); ix++) {
                    for (int iy = node.ylow(); iy <= node.yhigh(); iy++) {
                        if (node.type() == SOURCE) {
                            indices[SOURCE][ix][iy][0].push_back(inode);
                            indices[SINK][ix][iy][0].push_back(OPEN);
                        } else {
                            VTR_ASSERT(node.type() == SINK);
                            indices[SINK][ix][iy][0].push_back(inode);
                            indices[SOURCE][ix][iy][0].push_back(OPEN);
                        }
                    }
                }
            } else if (node.type() == IPIN || node.type() == OPIN) {
                for (int ix = node.xlow(); ix <= node.xhigh(); ix++) {
                    for (int iy = node.ylow(); iy <= node.yhigh(); iy++) {
                        if (node.type() == OPIN) {
                            indices[OPIN][ix][iy][node.side()].push_back(inode);
                            indices[IPIN][ix][iy][node.side()].push_back(OPEN);
                        } else {
                            VTR_ASSERT(node.type() == IPIN);
                            indices[IPIN][ix][iy][node.side()].push_back(inode);
                            indices[OPIN][ix][iy][node.side()].push_back(OPEN);
                        }
                    }
                }
            } else if (node.type() == CHANX) {
                for (int ix = node.xlow(); ix <= node.xhigh(); ix++) {
                    for (int iy = node.ylow(); iy <= node.yhigh(); iy++) {
                        coordinates_max_ptc[iy][ix].chanx_max_ptc = std::max(coordinates_max_ptc[iy][ix].chanx_max_ptc, node.ptc_num());
                    }
                }
            } else if (node.type() == CHANY) {
                for (int ix = node.xlow(); ix <= node.xhigh(); ix++) {
                    for (int iy = node.ylow(); iy <= node.yhigh(); iy++) {
                        coordinates_max_ptc[ix][iy].chany_max_ptc = std::max(coordinates_max_ptc[ix][iy].chany_max_ptc, node.ptc_num());
                    }
                }
            }
        }

        /* Alloc the lookup table */
        for (t_rr_type rr_type : RR_TYPES) {
            if (rr_type == CHANX) {
                for (size_t y = 0; y < grid_.height(); ++y) {
                    for (size_t x = 0; x < grid_.width(); ++x) {
                        indices[CHANX][y][x][0].resize(coordinates_max_ptc[y][x].chanx_max_ptc + 1, OPEN);
                    }
                }
            } else if (rr_type == CHANY) {
                for (size_t x = 0; x < grid_.width(); ++x) {
                    for (size_t y = 0; y < grid_.height(); ++y) {
                        indices[CHANY][x][y][0].resize(coordinates_max_ptc[x][y].chany_max_ptc + 1, OPEN);
                    }
                }
            }
        }

        int count;
        /* CHANX and CHANY need to reevaluated with its ptc num as the correct index*/
        for (size_t inode = 0; inode < rr_nodes_->size(); inode++) {
            auto& node = (*rr_nodes_)[inode];
            if (node.type() == CHANX) {
                for (int iy = node.ylow(); iy <= node.yhigh(); iy++) {
                    for (int ix = node.xlow(); ix <= node.xhigh(); ix++) {
                        count = node.ptc_num();
                        if (count >= int(indices[CHANX][iy][ix][0].size())) {
                            VPR_FATAL_ERROR(VPR_ERROR_ROUTE,
                                            "Ptc index %d for CHANX (%d, %d) is out of bounds, size = %zu",
                                            count, ix, iy, indices[CHANX][iy][ix][0].size());
                        }
                        indices[CHANX][iy][ix][0][count] = inode;
                    }
                }
            } else if (node.type() == CHANY) {
                for (int ix = node.xlow(); ix <= node.xhigh(); ix++) {
                    for (int iy = node.ylow(); iy <= node.yhigh(); iy++) {
                        count = node.ptc_num();
                        if (count >= int(indices[CHANY][ix][iy][0].size())) {
                            VPR_FATAL_ERROR(VPR_ERROR_ROUTE,
                                            "Ptc index %d for CHANY (%d, %d) is out of bounds, size = %zu",
                                            count, ix, iy, indices[CHANY][ix][iy][0].size());
                        }
                        indices[CHANY][ix][iy][0][count] = inode;
                    }
                }
            }
        }

        //Copy the SOURCE/SINK nodes to all offset positions for blocks with width > 1 and/or height > 1
        // This ensures that look-ups on non-root locations will still find the correct SOURCE/SINK
        for (size_t x = 0; x < grid_.width(); x++) {
            for (size_t y = 0; y < grid_.height(); y++) {
                int width_offset = grid_[x][y].width_offset;
                int height_offset = grid_[x][y].height_offset;
                if (width_offset != 0 || height_offset != 0) {
                    int root_x = x - width_offset;
                    int root_y = y - height_offset;

                    indices[SOURCE][x][y] = indices[SOURCE][root_x][root_y];
                    indices[SINK][x][y] = indices[SINK][root_x][root_y];
                }
            }
        }
    }

    const bool is_global_graph_;
    const bool read_edge_metadata_;
    t_chan_width* chan_width_;
    std::vector<t_rr_node>* rr_nodes_;
    std::vector<t_rr_switch_inf>* rr_switch_inf_;
    std::vector<t_rr_indexed_data>* rr_indexed_data_;
    std::vector<int> seg_index_;
    t_rr_node_indices* rr_node_indices_;
    int id_;
    int wire_to_rr_ipin_switch_;

    const size_t num_arch_switches_;
    const t_arch_switch_inf* arch_switch_inf_;
    const std::vector<t_segment_inf>& segment_inf_;
    const std::vector<t_physical_tile_type>& physical_tile_types_;
    const DeviceGrid& grid_;
    const std::unordered_map<int, t_metadata_dict>& rr_node_metadata_;
    const std::unordered_map<std::tuple<int, int, short>,
                             t_metadata_dict>& rr_edge_metadata_;
};