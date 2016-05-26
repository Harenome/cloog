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
CloogLoop* cloog_loop_generate_split(
    CloogOptions* options, size_t level, CloogLoop* loop);

extern CloogLoop * cloog_loop_copy(CloogLoop * source);
extern CloogLoop* cloog_loop_last(CloogLoop* loop);

#if defined(__cplusplus)
}
#endif

#endif /* CLOOG_POLYLIB_SPECIAL_H */
