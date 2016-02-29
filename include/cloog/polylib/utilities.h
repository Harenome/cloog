/**
 * \file utilities.h
 * \brief Various utilities to ease the use of the PolyLib in CLooG.
 *
 * This header declares various utilities that shall help using the PolyLib in
 * CLooG.
 */

#ifndef CLOOG_POLYLIB_UTILITIES_H
#define CLOOG_POLYLIB_UTILITIES_H

#if defined(__cplusplus)
extern "C"
{
#endif

#include <stdlib.h>

#include <polylib/polylibgmp.h>
#include <cloog/cloog.h>

/**
 * \brief Compute the disjoint union of a polyhedra union and free the input.
 *
 * \param[in] domain      Input domain.
 * \param[in] nb_max_rays Maximum number of rays.
 *
 * \return The disjoint union of the input union of polyhedra.
 *
 * \warning The geometrical dimension must be the same for all polyhedra.
 *          Otherwise, duplicates may appear.
 * \warning This function will free \a domain. Always use the returned value.
 *
 * \details This function computes the result of
 *          cloog_polylib_domain_soft_disjoint_const(), frees the input domain
 *          and returns the aforementioned result.
 *
 * \see cloog_polylib_domain_soft_disjoint_const().
 */
Polyhedron* cloog_polylib_domain_soft_disjoin(Polyhedron* domain,
                                              unsigned nb_max_rays);

/**
 * \brief Compute the disjoint union of an union of polyhedra.
 *
 * \param[in] p           The union of polyhedra.
 * \param[in] nb_max_rays Maximum number of rays.
 *
 * \return The disjoint union of the input union of polyhedra.
 *
 * \warning The geometrical dimension must be the same for all polyhedra.
 *          Otherwise, duplicates may appear.
 *
 * \details This function computes the disjoint union of an union of polyhedra.
 *          The resulting polyhedra shall be completely disjoint intersect (as
 *          opposed to PolyLib's Disjoint_Domain() which may allow – depending
 *          on its \a flag parameter – the resulting polyhedra to have their
 *          facets in common).
 *
 * \note This is a modified version of Disjoint_Domain() from the PolyLib.
 *       Whenever two domains (in the union of polyhedra) intersect, the PolyLib
 *       implementation produces three domains (the intersection and both
 *       domains without the intersection) whereas this implementation produces
 *       only two domains (one of the domains remains intact and the
 *       intersection is substracted from the other domain).
 */
Polyhedron* cloog_polylib_domain_soft_disjoint_const(Polyhedron* p,
                                                     unsigned nb_max_rays);

/**
 * \brief Compute the convex of a polyhedra union if, and only if it is convex.
 *
 * \param[in] domain      The input union of polyhedra.
 * \param[in] nb_max_rays Maximum number of rays.
 *
 * \return The convex version of the union if the union is convex. Otherwise,
 *         the input union remains intact and is returned.
 *
 * \warning The input union may be freed. Always use the returned value or use
 *          cloog_polylib_domain_convex_attempt_const().
 *
 * \see cloog_polylib_domain_convex_attempt_const().
 */
Polyhedron* cloog_polylib_domain_convex_attempt(Polyhedron* domain,
                                                unsigned nb_max_rays);

/**
 * \brief Compute the convex of a polyhedra union if, and only if it is convex.
 *
 * \param[in] domain      The input union of polyhedra.
 * \param[in] nb_max_rays Maximum number of rays.
 *
 * \return The convex version of the union if the union is convex. Otherwise,
 *         the input union remains intact and is returned.
 *
 * \warning The input union may be freed. Always use the returned value.
 */
Polyhedron* cloog_polylib_domain_convex_attempt_const(Polyhedron* domain,
                                                      unsigned nb_max_rays);
/**
 * \brief Count the constraints of a domain.
 *
 * \param[in] domain Input domain.
 *
 * \return The number of constraints.
 */
size_t cloog_polylib_domain_count(const Polyhedron* domain);

/**
 * \brief Get constraints from a domain.
 *
 * \param[in]  domain        Input domain.
 * \param[in] nb_constraints The number of constraints.
 *
 * \return The constraints.
 */
Matrix** cloog_polylib_domain_get_constraints(Polyhedron* domain,
                                              size_t nb_constraints);

/**
 * \brief Get the last domain in a Polyhedron.
 *
 * \param[in] domain Input domain.
 *
 * \return The last domain in the union of domains.
 */
Polyhedron* cloog_polylib_domain_last(Polyhedron* domain);

/**
 * \brief Free a domain if, and only if, it is empty.
 *
 * \param[in] domain Input domain.
 *
 * \return (Polyhedron*) 0
 */
Polyhedron* cloog_polylib_domain_free_empty(Polyhedron* domain);

/**
 * \brief Convert a Polyhedron domain to a CloogLoop.
 *
 * \param[in] domain Input domain.
 * \param[in] state  Input state.
 * \param[in] source Original loop.
 * \param[in] nb_par Parameter number.
 *
 * \return The resulting CloogLoop.
 */
CloogLoop* cloog_loop_from_polyhedron(Polyhedron* domain, CloogState* state,
                                      CloogLoop* source, const size_t nb_par);

#if defined(__cplusplus)
}
#endif

#endif /* CLOOG_POLYLIB_UTILITIES_H */
