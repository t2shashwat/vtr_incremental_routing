/**
 * @file
 * @brief Read a .route file and load the route tree and other associated data structure
 *        with the correct values.
 *
 * This is used to perform --analysis only
 */

#ifndef READ_ROUTE_H
#define READ_ROUTE_H

bool read_route(const char* route_file, const t_router_opts& RouterOpts, bool verify_file_digests);
bool read_route_incr_route(ClusterNetId inet_target, const char* route_file, const t_router_opts& RouterOpts, bool verify_file_digests);

#endif /* READ_ROUTE_H */
