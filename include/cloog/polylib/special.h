/**
 * \file special.h
 * \brief Special features thanks to PolyLib and thus specific to CLooG-PolyLib.
 *
 * This header declares functions for special features only available when using
 * the PolyLib.
 */

#ifndef CLOOG_POLYLIB_SPECIAL_H
#define CLOOG_POLYLIB_SPECIAL_H

#if defined(__cplusplus)
extern "C"
{
#endif

#include <polylib/polylibgmp.h>
#include <cloog/cloog.h>

/**
 * \defgroup cloog_polylib_special CLooG-PolyLib special features
 */

/**
 * \brief Split loops.
 *
 * \param[in] state               CLooG's state.
 * \param[in] context             Input context.
 * \param[in] loop                Input loops.
 * \param[in] max_depth           Maximum depth of the split algorithm.
 * \param[in] nb_max_splits       Maximum number of splits.
 * \param[in] nb_max_constraints  Maximum constraint number in split candidates.
 * \param[in] nb_max_dependencies Maximum dependency level in split candidates.
 * \param[in] nb_max_rays         Maximum number of rays.
 *
 * \return The split loops, if any.
 *
 * \warning This function behaves correctly only with CloogLoop loops that have
 *          not been processed yet.
 * \warning The input loop may be freed. Always use the returned value.
 * \details This function attempts to split loops. It is intended to be used at
 *          the earliest stages of CLooG (right after the input is read) or
 *          the first level of cloog_loop_generate() because it expects a
 *          CloogLoop that has not been separated yet.
 *
 *          The arguments \a state and \a context will not be modified. On the
 *          other hand, \a loop will be freed.
 *
 *          Various parameters can be used to tune the splitting operation:
 *          - \a max_depth controls the maximum depth of the algorithm, i.e.
 *            the level of most inner loop to be split.
 *          - \a nb_max_splits defines the maximum number of splits to be
 *            returned at the end of the algorithm.
 *          - \a nb_max_constraints specifies the maximum number of constraints
 *            a domain must not exceed in order to be elligible for the split
 *            algorithm.
 *          - \a nb_max_dependencies limits the dependency level for iterators
 *            in a split candidate (i.e. each iterator of the domain can depend
 *            on, at most, \a nb_max_dependencies other iterators).
 *
 * \ingroup cloog_polylib_special
*/
CloogLoop* cloog_loop_polylib_split(
    CloogState* state, CloogDomain* context, CloogLoop* loop, size_t level,
    size_t max_depth, size_t nb_max_splits, size_t nb_max_constraints,
    size_t nb_max_dependencies, unsigned nb_max_rays);

/**
 * \brief If possible, split the current loop.
 *
 * \param[in] state     CLooG state.
 * \param[in] context   Input context.
 * \param[in] loop      Input loops.
 * \param[in] max_depth Maximum depth of the split algorithm.
 * \param[in] nb_max_splits       Maximum number of splits.
 * \param[in] nb_max_constraints  Maximum constraint number in split candidates.
 * \param[in] nb_max_dependencies Maximum dependency level in split candidates.
 * \param[in] nb_max_rays         Maximum number of rays.
 *
 * \return The unmodified loop or the splits, if any.
 *
 * \note The input loop is not modified.
 * \details This function attempts to split the current CloogLoop (i.e. only the
 *          struct directly pointed by \a loop, not all loops on the same level
 *          pointed by \a loop->next). See cloog_loop_polylib_split() for an
 *          explanation on \a max_depth, \a nb_max_splits,
 *          \a nb_max_constraints, \a nb_max_dependencies and \a nb_max_rays.
 * \see cloog_loop_polylib_split
 *
 *
 * \ingroup cloog_polylib_special
 */
CloogLoop* cloog_loop_polylib_split_current_loop(
    CloogState* state, CloogDomain* context, CloogLoop* loop, size_t level,
    size_t max_depth, size_t nb_max_splits, size_t nb_max_constraints,
    size_t nb_max_dependencies, unsigned nb_max_rays);

/**
 * \brief Split a PolyLib domain.
 *
 * \param[in]  context             Input context.
 * \param[in]  domain              Input domain.
 * \param[in]  nb_parameters       The number of parameters.
 * \param[in]  max_depth           Max depth of the split algorithm.
 * \param[in]  nb_max_splits       Max number of splits.
 * \param[in]  nb_max_constraints  Max constraint number in split candidates.
 * \param[in]  nb_max_dependencies Max dependency level in split candidates.
 * \param[in]  nb_max_rays         Max number of rays.
 * \param[out] nb_splits           The number of splits.
 *
 * \return Splits of the input domain.
 *
 * \note The input parameters \a context and \a domain are not modified.
 * \details If \a max_depth = 0, then the depth is infinite. See
 *          cloog_loop_polylib_split() for an explanation on \a max_depth,
 *          \a nb_max_splits, \a nb_max_constraints, \a nb_max_dependencies and
 *          \a nb_max_rays.
 * \see cloog_loop_polylib_split_current_loop
 *
 * \ingroup cloog_polylib_special
 */
Polyhedron** cloog_polylib_split_polyhedron(
    Polyhedron* context, Polyhedron* domain, size_t nb_parameters, size_t level,
    size_t max_depth, size_t nb_max_split, size_t nb_max_constraints,
    size_t nb_max_dependencies, unsigned nb_max_rays, size_t* nb_splits);

#if defined(__cplusplus)
}
#endif

#endif /* CLOOG_POLYLIB_SPECIAL_H */
