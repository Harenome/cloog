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
 * \param[in] state     The CLooG state.
 * \param[in] context   The input context.
 * \param[in] loop      The input loops.
 * \param[in] max_depth The maximum depth of the split algorithm.
 *
 * \return The split loops, if any.
 *
 * \warning This function behaves correctly only with CloogLoop loops that have
 *          not been processed yet.
 * \warning The input loop may be freed. Always use the returned value.
 * \details This function attempts to split loops. It is intended to be used at
 *          the earliest stages of CLooG (right after the input is read) because
 *          it expects a CloogLoop that has not been processed yet.
 *
 *          The arguments \a state and \a context will not be modified. On the
 *          other hand, \a loop will be freed.
 *
 * \ingroup cloog_polylib_special
 */
CloogLoop* cloog_loop_polylib_split(CloogState* state, CloogDomain* context,
                                    CloogLoop* loop, size_t max_depth,
                                    unsigned nb_max_rays);

/**
 * \brief If possible, split the current loop.
 *
 * \param[in] state     CLooG state.
 * \param[in] context   Input context.
 * \param[in] loop      Input loops.
 * \param[in] max_depth Maximum depth of the split algorithm.
 *
 * \return The unmodified loop or the splits, if any.
 *
 * \note The input loop is not modified.
 * \details This function attempts to split the current CloogLoop (i.e. only the
 *          struct directly pointed by \a loop, not all loops on the same level
 *          pointed by \a loop->next).
 *
 * \ingroup cloog_polylib_special
 */
CloogLoop* cloog_loop_polylib_split_current_loop(CloogState* state,
                                                 CloogDomain* context,
                                                 CloogLoop* loop,
                                                 size_t max_depth,
                                                 unsigned nb_max_rays);

/**
 * \brief Split a PolyLib domain.
 *
 * \param[in]  context       Input context.
 * \param[in]  domain        Input domain.
 * \param[in]  nb_parameters The number of parameters.
 * \param[in]  max_depth     Maximum depth of the split algorithm.
 * \param[in]  nb_max_rays   Maximum number of rays.
 * \param[out] nb_splits     The number of splits.
 *
 * \return Splits of the input domain.
 *
 * \note The input parameters \a context and \a domain are not modified.
 * \details If \a max_depth = 0, then the depth is infinite.
 *
 * \ingroup cloog_polylib_special
 */
Polyhedron** cloog_polylib_split_polyhedron(Polyhedron* context,
                                            Polyhedron* domain,
                                            size_t nb_parameters,
                                            size_t max_depth,
                                            unsigned nb_max_rays,
                                            size_t* nb_splits);

#if defined(__cplusplus)
}
#endif

#endif /* CLOOG_POLYLIB_SPECIAL_H */
