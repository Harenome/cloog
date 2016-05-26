/**
 * \file special.c
 * \brief Special features thanks to PolyLib and thus specific to CLooG-PolyLib.
 *
 * This file defines functions for special features only available when using
 * the PolyLib.
 */

#include "cloog/polylib/special.h"

#include "cloog/polylib/cloog.h"
#include "cloog/polylib/utilities.h"

#define WS 0

typedef struct polylib_split_data {
  size_t count;
  Polyhedron** domains;
} polylib_split_data;

static void polylib_split_data_init(polylib_split_data*, size_t);
static size_t polylib_split_data_count_all_polyhedra(const polylib_split_data*);
static void polylib_split_data_flatten_domains(polylib_split_data*, const unsigned);
static void polylib_split_data_clean(polylib_split_data*);
static void polylib_split_data_swap(polylib_split_data*, polylib_split_data*);
static void polylib_split_data_grow(polylib_split_data*, size_t);

/*******************************************************************************
 * Static functions prototypes                                                 *
 ******************************************************************************/

/* Utilities for CloogLoop. */
static inline CloogLoop* cloog_loop_copy_current_loop(CloogLoop*);

/* Various utilities. */

/* Utilities for cloog_polylib_split_domain(). */
static int check_domain_dependencies(const Polyhedron*, size_t, size_t);
static inline int domain_is_splittable(const Polyhedron*, size_t, size_t,
                                       size_t);
static void hunt_min_max(Polyhedron*, size_t, size_t, unsigned,
                         polylib_split_data*);

/* Utilities for hunt_min_max(). */
static inline int iterator_is_alone(const Value*, size_t, size_t);
static inline int has_param_or_constant(const Value*, size_t, size_t);

/* Parameterized domain. */
static Matrix* matrix_move_iterators_last(Matrix*, size_t, size_t);
static Polyhedron* domain_move_iterators_last(Polyhedron*, size_t, size_t,
                                              unsigned);
static Param_Polyhedron* compute_param_domain(Polyhedron*, Polyhedron*, size_t,
                                              size_t, unsigned);

/* New splits. */
static size_t param_domain_count_domains(const Param_Domain*);
static void extract_domains(const Param_Polyhedron*, size_t, unsigned,
                            Polyhedron**);
static Polyhedron* fix_boundaries(Polyhedron*, Polyhedron*, unsigned);
static Matrix* matrix_fix_iterators(Matrix*, size_t, size_t);
static Polyhedron* fix_iterators_and_mix(Polyhedron*, Polyhedron*, size_t,
                                         size_t, unsigned);
static void get_new_domains(Polyhedron*, const Param_Polyhedron*, size_t,
                            size_t, size_t, unsigned, Polyhedron**);

static CloogLoop* cloog_loop_split(
    CloogOptions* options, size_t level, CloogLoop* loop);

static CloogLoop* cloog_loop_split_current(
    CloogOptions* options, size_t nb_par, size_t level, CloogLoop* loop);

static Polyhedron** cloog_polylib_split_domain(
    CloogOptions* options, size_t nb_parameters, size_t level,
    Polyhedron* domain, size_t* nb_splits);
/*******************************************************************************
 * cloog/polylib/special.h functions                                           *
 ******************************************************************************/


CloogLoop* cloog_loop_generate_split(
    CloogOptions* options, const size_t level, CloogLoop* loop) {
  CloogLoop* current;

  for (current = loop; current; current = current->next)
    current->inner = cloog_loop_split(options, level, current->inner);

  return loop;
}

CloogLoop* cloog_loop_split(
    CloogOptions* options, const size_t level, CloogLoop* loop) {
   /* CloogLoop is a chained struct. Some structs will be changed, some will
    * remain intact. cloog_loop_split_current() shall never modify
    * nor return any of the structs from `loop`. It must either return a copy
    * of the current struct (in the case where it is not possible to split it)
    * or the splits of the current struct. This way, `loop` can safely be freed
    * as a whole with cloog_loop_free() at the end of this function.
    *
    * (cloog_loop_free() frees the current struct, then moves on to
    * the struct's `next` field. However, at each step, we don't know whether
    * the current struct has been split.)
    */
  const size_t nb_par = loop->domain->nb_par;
  CloogLoop* next, *last, *current, *result;

  result = loop;
  next = last = current = NULL;
  if (!level && loop && loop->domain) {
    /* cloog_loop_split_current() may free loop. */
    next = loop->next;
    result = cloog_loop_split_current(options, nb_par, level, loop);
    last = cloog_loop_last(result);
    for (current = next; current; current = next) {
      /* cloog_loop_split_current() may free current. */
      next = current->next;
      last->next = cloog_loop_split_current(options, nb_par, level, current);
      last = cloog_loop_last(last);
    }
  }

  return result;
}

/**
 * \brief If possible, split the current loop.
 *
 * \param[in] options CLooG options.
 * \param[in] nb_par  Number of parameters.
 * \param[in] level   Current recursion level.
 * \param[in] loop    Input loop.
 *
 * \return The input loop or the splits, if any.
 *
 * \warning The input loop *may* be freed.
 * \details This function attempts to split the current CloogLoop (i.e. only the
 *          struct directly pointed by \a loop, not all loops on the same level
 *          pointed by \a loop->next).
 * \see cloog_loop_split
 *
 * \ingroup cloog_polylib_special
 */
CloogLoop* cloog_loop_split_current(
    CloogOptions* const options, const size_t nb_par, const size_t level,
    CloogLoop* const loop) {
  /* The first step is to copy the input `loop` because it is important that the
   * input `loop` remains unmodified.
   * If the loop can be split, this copy will be freed and the splits will be
   * returned. Otherwise, the copy will be returned.
   */
  size_t count = 0, i = 0;
  Polyhedron** splits = NULL;
  CloogLoop* last = NULL, *result = loop;
  CloogState* const state = options->state;
  Polyhedron* const domain = cloog_domain_polyhedron(loop->domain);

  splits = cloog_polylib_split_domain(options, nb_par, level, domain, &count);
  if (splits) {
    if (count) {
      result = cloog_loop_from_polyhedron(splits[0], state, loop, nb_par);
      last = cloog_loop_last(result);
      for (i = 1; i < count; ++i) {
        last->next = cloog_loop_from_polyhedron(splits[i], state, loop, nb_par);
        last = cloog_loop_last(last);
      }
      /* Only free the current loop! */
      loop->next = NULL;
      cloog_loop_free(loop);
    }
    free(splits);
  }

  return result;
}

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
 *          cloog_loop_split() for an explanation on \a max_depth,
 *          \a nb_max_splits, \a nb_max_constraints, \a nb_max_dependencies and
 *          \a nb_max_rays.
 * \see cloog_loop_split_current
 *
 * \ingroup cloog_polylib_special
 */
Polyhedron** cloog_polylib_split_domain(
    CloogOptions* options, const size_t nb_parameters, const size_t level,
    Polyhedron* const domain, size_t* const result_size) {
  /* Summary:
   * At each level, the first step is to "flatten" the domains to prepare them
   * for compute_param_domain(). Then, each domain is split using
   * compute_param_domain() which returns a Param_Polyhedron. get_new_domains()
   * computes the new domains using this Param_Polyhedron and the algorithm
   * moves on to the next level.
   */
  /* A domain can be an union of polyhedra. However, compute_param_domain()
   * expects single polyhedron domains. flatten_domains() breaks such unions
   * of polyhedra into an array of single polyhedra.
   *
   * compute_param_domain() modifies (copies of) the input context and domain,
   * so the domains in the resulting Param_Polyhedron can not be directly used.
   * Hence the call to get_new_domains() which extracts the new domains, fixes
   * them (i.e. reverts the modifications due to compute_param_domain()) and
   * intersects them with their parent (i.e. the domain passed to
   * compute_param_domain()).
   *
   * As the final number of splits is unknown at the beginning, they are stored
   * using two Polyhedron* arrays (the pointers 'splits' and 'new_splits') that
   * are realloc'd. Domains that must be split are stored in splits. If a domain
   * can not be split, it is copied into new_splits, otherwise its splits are
   * stored into new_splits. The two arrays are swapped/freed each time the end
   * of a given step of the algorithm is reached.
   */
  const size_t max_rays = options->state->backend->MAX_RAYS;
  const size_t max_depth = options->split_depth;
  const size_t nb_max_splits = options->split_max;
  const size_t max_constraints = options->split_constraints;
  const size_t max_dependencies = options->split_dependencies;

  const size_t nb_iterators = domain->Dimension - nb_parameters;
  const size_t max_level =
    (max_depth && max_depth < nb_iterators) ? (max_depth + 1) : nb_iterators;

  /* Param_Polyhedron */
  Param_Polyhedron* PA = NULL;
  Matrix* context_matrix = NULL;
  Polyhedron* context = NULL, *current = NULL;

  polylib_split_data splits, new_splits;
  size_t depth = 1, i = 0, j = 0, nb_temp_domains = 0;

  /* To store both the remaining domains to process and the final results. */
  polylib_split_data_init(&splits, 1);
  polylib_split_data_init(&new_splits, 0);
  /* Clean up the input, if needed. */
  splits.domains[0] = Domain_Copy(domain);
  splits.domains[0] = cloog_polylib_domain_convex_attempt(splits.domains[0],
                                                          max_rays);

  for (depth = 1; depth < max_level; ++depth) {
    /* Flattening will produce, at most, count_all_polyhedra() domains. */
    if (polylib_split_data_count_all_polyhedra(&splits) > nb_max_splits)
      goto cloog_polylib_split_polyhedron_end;
    /* Flatten and split the domains. */
    polylib_split_data_flatten_domains(&splits, max_rays);
    /* Empty context. */
    context_matrix = Matrix_Alloc(1, nb_parameters + depth + 2);
    context = Constraints2Polyhedron(context_matrix, max_rays);
    Matrix_Free(context_matrix);
    /* Split min/max on iterators. */
    for (j = 0; j < splits.count; ++j) {
      if (new_splits.count + 1 > nb_max_splits)
        goto cloog_polylib_split_polyhedron_end;
      current = splits.domains[j];
      if (domain_is_splittable(current, nb_iterators, max_constraints,
                               max_dependencies))
        PA = compute_param_domain(context, current, nb_iterators, depth,
                                  max_rays);
      else
        PA = NULL;
      nb_temp_domains = PA ? param_domain_count_domains(PA->D) : 1;
      if (nb_temp_domains + new_splits.count > nb_max_splits)
        goto cloog_polylib_split_polyhedron_end;
      polylib_split_data_grow(&new_splits, nb_temp_domains);
      get_new_domains(current, PA, nb_temp_domains, nb_iterators, depth,
                      max_rays, &new_splits.domains[new_splits.count]);
      new_splits.count += nb_temp_domains;
      if (PA) {
        Param_Polyhedron_Free(PA);
        PA = NULL;
      }
    }
    polylib_split_data_swap(&splits, &new_splits);
    Domain_Free(context);
    context = NULL;
  }

  /* Split min/max on parameters. */
  polylib_split_data_flatten_domains(&splits, max_rays);
  for (i = 0; i < splits.count; ++i) {
    hunt_min_max(splits.domains[i], nb_iterators, nb_parameters, max_rays,
                 &new_splits);
    if ((splits.count - i + new_splits.count) > nb_max_splits)
      goto cloog_polylib_split_polyhedron_end;
  }
  polylib_split_data_swap(&splits, &new_splits);

  cloog_polylib_split_polyhedron_end:
    if (PA)
      Param_Polyhedron_Free(PA);
    if (context)
      Domain_Free(context);
    polylib_split_data_clean(&new_splits);

    *result_size = splits.count;
    return splits.domains;
}

/*******************************************************************************
 * Static functions definitions                                                *
 ******************************************************************************/

void polylib_split_data_init(polylib_split_data* data, size_t size) {
  if (data) {
    data->count = size;
    data->domains = size ? malloc(size * (sizeof *data->domains)) : NULL;
  }
}

size_t polylib_split_data_count_all_polyhedra(const polylib_split_data* data) {
  size_t nb_domains = 0, count = 0, i = 0;
  const Polyhedron* p = NULL;
  Polyhedron** domains = NULL;

  if (!data || !data->count)
    return 0;

  nb_domains = data->count;
  domains = data->domains;
  count = 0;
  for (i = 0; i < nb_domains; ++i)
    for (p = domains[i]; p; p = p->next)
      ++count;
  return count;
}

/**
 * \brief Break unions of polyhedra into an array of single polyhedra.
 *
 * \param[in]     domains     Array of domains to flatten.
 * \param[in]     nb_max_rays Maximum number of rays.
 * \param[in,out] nb_domains  Number of domains.
 *
 * \return Array of "flattened" domains.
 *
 * \warning The input array \a domain will be freed, always use the returned
 *          value.
 *
 * \details The current implementation of Polyhedron in PolyLib relies on
 *          chained structs to represent unions of polyhedra. If any of the
 *          domains supplied to this function happens to be an union, the
 *          domains of this union are broken into separate Polyhedron structs
 *          and placed in the returned array. Otherwise, the domain is merely
 *          transferred to the returned array.
 */
void polylib_split_data_flatten_domains(polylib_split_data* data,
                                        const unsigned max_rays) {
  /* Flattening: Separate unions of domains into single domains.
   * For each union of domains, we first attempt to find the convex domain and
   * then compute the equivalent disjoint union to ensure they can safely be
   * separated.
   */
  Polyhedron* next = NULL, *p = NULL;
  const size_t nb_domains_local = data->count;
  Polyhedron** domains = data->domains;

  size_t position = 0, i = 0;
  size_t nb_flattened = data->count;
  Polyhedron** result = malloc (nb_flattened * (sizeof *result));
  for (i = 0; i < nb_domains_local; ++i) {
    if (domains[i]->next) {
      domains[i] = cloog_polylib_domain_convex_attempt(domains[i], max_rays);
      domains[i] = cloog_polylib_domain_soft_disjoin(domains[i], max_rays);
      next = NULL;
      nb_flattened += cloog_polylib_domain_count(domains[i]) - 1;
      result = realloc (result, nb_flattened * (sizeof *result));
      for (p = domains[i]; p; p = next) {
        next = p->next;
        p->next = NULL;
        result[position] = p;
        ++position;
      }
    } else {
      result[position] = domains[i];
      ++position;
    }
  }
  free(domains);

  data->count = nb_flattened;
  data->domains = result;
}

void polylib_split_data_clean(polylib_split_data* data) {
  size_t nb_domains = 0, i = 0;
  if (data && data->domains) {
    nb_domains = data->count;
    for (i = 0; i < nb_domains; ++i)
      Domain_Free(data->domains[i]);
    free(data->domains);
    data->count = 0;
    data->domains = NULL;
  }
}

void polylib_split_data_swap(polylib_split_data* old, polylib_split_data* new) {
  if (old && new) {
    polylib_split_data_clean(old);
    old->domains = new->domains;
    old->count = new->count;
    new->domains = NULL;
    new->count = 0;
  }
}

static void polylib_split_data_grow(polylib_split_data* data, size_t growth) {
  const size_t new_size = data->count + growth;
  data->domains = realloc(data->domains, new_size * (sizeof *data->domains));
}

/**
 * \brief Ensure iterators do not depend on too many other iterators.
 *
 * \param[in] p                   Input domain.
 * \param[in] nb_iterators        Total number of iterators.
 * \param[in] nb_max_dependencies Maximum dependency level.
 *
 * \return The result of the test.
 * \retval 0 if at least one iterator depends on \a nb_max_dependencies
 *           (or more) iterators.
 * \retval 1 if the iterators of the domain do not depend on too many iterators.
 */
int check_domain_dependencies(const Polyhedron* p, size_t nb_iterators,
                              size_t nb_max_dependencies) {
  int result = 1;
  size_t row = 0, inner = 0, count = 0, column = 1;
  for (row = 0; result && row < p->NbConstraints; ++row) {
    Value* constraint = p->Constraint[row];
    /* Find the most inner iterator. */
    inner = nb_iterators;
    for ( ; !value_notzero_p(constraint[inner]) && inner >= 1; --inner)
      ;
    /* Count the dependencies. */
    count = 0;
    for (column = 1; result && column < inner; ++column)
      if (value_notzero_p(constraint[column])) {
        ++count;
        if (count >= nb_max_dependencies)
          result = 0;
      }
  }

  return result;
}

/**
 * \brief Decide a domain's splittableness according to various characteristics.
 *
 * \param[in] domain             The split candidate.
 * \param[in] nb_iterators       Total number of iterators.
 * \param[in] nb_max_constraints Maximum number of constraints.
 * \param[in] nb_max_deps        Maximum dependency level.
 *
 * \return The splittableness of the domain.
 * \retval 0 if it is a bad idea to attempt to split the domain.
 * \retval 1 if the domain *seems* to be a good split candidate.
 */
int domain_is_splittable(const Polyhedron* domain, size_t nb_iterators,
                         size_t nb_max_constraints, size_t nb_max_deps) {
  int is_ok = domain->NbConstraints < nb_max_constraints &&
              check_domain_dependencies(domain, nb_iterators, nb_max_deps);
  return is_ok;
}

/**
 * \brief Determine whether an iterator is the sole iterator in a constraint.
 *
 * \param[in] constraint   Input constraint.
 * \param[in] it           Iterator index.
 * \param[in] nb_iterators Total number of iterators.
 *
 * \return The result of the test.
 * \retval 0 if more than one iterator is used in the constraint.
 * \retval 1 if the iterator is the only one used in the constraint.
 */
int iterator_is_alone(const Value* constraint, const size_t it,
                      const size_t nb_iterators) {
  int is_alone = 1;
  size_t column = 1;
  for (column = 1; is_alone && column < it; ++column)
    if (value_notzero_p(constraint[column]))
      is_alone = 0;
  for (column = it + 1; is_alone && column <= nb_iterators; ++column)
    if (value_notzero_p(constraint[column]))
      is_alone = 0;
  return is_alone;
}

/**
 * \brief Determine whether at least a parameter or a constant is used.
 *
 * \param[in] constraint    Input constraint.
 * \param[in] nb_iterators  Total number of iterators.
 * \param[in] nb_parameters Total number of parameters.
 *
 * \return The result of the test.
 * \retval 0 if no parameter or constant is used.
 * \retval 1 if at least one parameter or a constant is used in the constraint.
 */
int has_param_or_constant(const Value* constraint, const size_t nb_iterators,
                          const size_t nb_parameters) {
  int has_one = 0;
  const size_t start = nb_iterators + 1;
  const size_t end = nb_iterators + nb_parameters + 1;
  size_t column = 0;
  for (column = start; !has_one && column <= end; ++column)
    if (value_notzero_p(constraint[column]))
      has_one = 1;
  return has_one;
}

/**
 * \brief Hunt down min/max on parameters.
 *
 * \param[in]     domain        Input domain.
 * \param[in]     nb_iterators  Total number of iterators.
 * \param[in]     nb_parameters Total number of parameters.
 * \param[in]     max_rays      Maximum number of rays.
 * \param[in,out] result        Array of the resulting splits.
 * \param[in,out] result_size   Size of the resulting array.
 */
void hunt_min_max(Polyhedron* domain, const size_t nb_iterators,
                  const size_t nb_parameters, const unsigned max_rays,
                  polylib_split_data* result) {
  const size_t dimension = domain->Dimension + 2;
  int splitted = 0;
  Polyhedron* p = NULL;
  polylib_split_data splits, new_splits;
  Value new_constraint[dimension], * minimums[2], *maximums[2];
  size_t i, j, it, nb_min, nb_max, nb_temp_domains, nb_rows, col, row;

  /* Initializations. */
  i = j = it = nb_min = nb_max = nb_temp_domains = nb_rows = col = row = 0;
  for (i = 0; i < dimension; ++i)
    value_init(new_constraint[i]);
  polylib_split_data_init(&splits, 1);
  polylib_split_data_init(&new_splits, 0);
  splits.domains[0] = Domain_Copy(domain);

  for (it = 1; it <= nb_iterators; ++it) {
    splitted = 1;
    while (splitted) {
      splitted = 0;
      for (i = 0; i < splits.count; ++i) {
        p = splits.domains[i];
        nb_rows = p->NbConstraints;
        /* Find the min/max constraints. */
        nb_min = nb_max = 0;
        for (row = 0; nb_min < 2 && nb_max < 2 && row < nb_rows; ++row) {
          Value* constraint = p->Constraint[row];
          if (!value_zero_p(constraint[0]) &&
              value_notzero_p(constraint[it]) &&
              iterator_is_alone(constraint, it, nb_iterators) &&
              has_param_or_constant(constraint, nb_iterators, nb_parameters)) {
            if (value_pos_p(constraint[it]))
              maximums[nb_max++] = constraint;
            else
              minimums[nb_min++] = constraint;
          }
        }
        nb_temp_domains = 1;
        /* Enlarge the domain list. */
        polylib_split_data_grow(&new_splits, 2);
        /* Split, if possible. */
        if (nb_min >= 2 || nb_max >= 2) {
          splitted = 1;
          ++nb_temp_domains;
          /* Prepare the first new constraint. */
          Value** involved = nb_max >= 2 ? maximums : minimums;
          Vector_Copy(involved[0], new_constraint, dimension);
          for (col = 1; col < dimension; ++col)
            value_subtract(new_constraint[col], new_constraint[col],
                           involved[1][col]);
          /* Use the first constraint for the first split. */
          new_splits.domains[new_splits.count] =
              AddConstraints(new_constraint, 1, p, max_rays);
          /* Prepare the second new constraint. */
          for (col = 1; col < dimension; ++col)
            value_oppose(new_constraint[col], new_constraint[col]);
          value_sub_int(new_constraint[dimension - 1],
                        new_constraint[dimension - 1], 1);
          /* Use the second constraint for the second split. */
          new_splits.domains[new_splits.count + 1] =
              AddConstraints(new_constraint, 1, p, max_rays);
        } else {
          new_splits.domains[new_splits.count] = Domain_Copy(p);
        }
        new_splits.count += nb_temp_domains;
      }
      polylib_split_data_swap(&splits, &new_splits);
    }
  }

  /* Copy the new splits. */
  polylib_split_data_grow(result, splits.count);
  for (j = 0; j < splits.count; ++j)
    result->domains[result->count + j] = splits.domains[j];
  result->count += splits.count;

  if (splits.domains)
    free(splits.domains);
  for (i = 0; i < dimension; ++i)
    value_clear(new_constraint[i]);
}

/**
 * \brief Copy the current CloogLoop.
 *
 * \param[in] loop Input loop.
 *
 * \return The copy.
 */
CloogLoop* cloog_loop_copy_current_loop(CloogLoop* loop) {
  /* At this step, most fields are NULL. */
  CloogLoop* result = cloog_loop_malloc(loop->state);
  result->domain = cloog_domain_copy(loop->domain);
  result->unsimplified = cloog_domain_copy(loop->unsimplified);
  result->otl = loop->otl;
  result->stride = cloog_stride_copy(loop->stride);
  result->block = cloog_block_copy(loop->block);

  return result;
}

/**
 * \brief Move some iterators from the start to the end.
 *
 * \param[in,out] m            Input matrix
 * \param[in]     nb_iterators The total number of iterators.
 * \param[in]     n            The number of iterators to move.
 *
 * \return The resulting matrix (\a m).
 *
 * \warning This function modifies \a m.
 */
Matrix* matrix_move_iterators_last(Matrix* const m, const size_t nb_iterators,
                                   const size_t n) {
  size_t i = 0;
  int j = 1;
  for (i = 0; i < n; ++i)
    for (j = 1; j < (int) (nb_iterators - i); ++j)
      ExchangeColumns(m, j, j + 1);

  return m;
}

/**
 * \brief Move some iterators from the start to the end.
 *
 * \param[in] domain       Input domain
 * \param[in] nb_iterators The total number of iterators.
 * \param[in] n            The number of iterators to move.
 *
 * \return The resulting domain (\a m).
 *
 * \warning Do not use this function if \a domain is an union of polyhedra.
 */
Polyhedron* domain_move_iterators_last(
    Polyhedron* const domain, const size_t nb_iterators, const size_t n,
    unsigned nb_max_rays) {
  Polyhedron* temp = Domain_Copy(domain);
  Matrix* constraints = matrix_move_iterators_last(Polyhedron2Constraints(temp),
                                                   nb_iterators, n);
  Polyhedron* result = Constraints2Polyhedron(constraints, nb_max_rays);
  Matrix_Free(constraints);
  Domain_Free(temp);

  return result;
}

/**
 * \brief Compute the parameterized domain.
 *
 * \param[in] original_context Input context.
 * \param[in] original_domain  Input domain.
 * \param[in] nb_iterators     The number of iterators.
 * \param[in] depth            The current depth in the algorithm.
 *
 * \return The parameterized domain.
 *
 * \note This function does not modify its input.
 */
Param_Polyhedron* compute_param_domain(
    Polyhedron* const original_context, Polyhedron* const original_domain,
    const size_t nb_iterators, const size_t depth, const unsigned nb_max_rays) {
  /* This step temporarily "transforms" iterators into parameters:
   * - move their columns after the remaining iterators in the domain.
   * - add empty columns to the context.
   * Finally, the Param_Polyhedron is computed.
   */
  Polyhedron* context = Domain_Copy(original_context);
  Polyhedron* domain = domain_move_iterators_last(original_domain, nb_iterators,
                                                  depth, nb_max_rays);
  Param_Polyhedron* PA = Polyhedron2Param_Domain(domain, context, WS);
  Domain_Free(context);
  Domain_Free(domain);

  return PA;
}

/**
 * \brief Count the domains of a Param_Domain.
 *
 * \param[in] domain Input Param_Domain.
 *
 * \return The number of domains.
 */
size_t param_domain_count_domains(const Param_Domain* domain) {
  size_t count = 0;
  const Param_Domain* P = domain;
  for (P = domain; P; P = P->next)
    ++count;

  return count;
}

/**
 * \brief Extract the domains of a Param_Polyhedron.
 *
 * \param[in] PA          The input Param_Polyhedron.
 * \param[in] nb_domains  The number of domains to extract.
 * \param[in] nb_max_rays Maximum number of rays
 * \param[in,out]
 *
 * \return The domains.
 */
void extract_domains(const Param_Polyhedron* PA, const size_t nb_domains,
                     const unsigned nb_max_rays, Polyhedron** const domains) {
  size_t count = 0, i = 0;
  const Param_Domain* P = NULL;
  for (P = PA->D; P; P = P->next) {
    domains[count] = cloog_polylib_domain_soft_disjoint_const(P->Domain,
                                                              nb_max_rays);
    for (i = 0; i < count; ++i) {
      domains[i] = fix_boundaries(domains[i], domains[count], nb_max_rays);
      if (domains[i]->next)
        domains[i] = cloog_polylib_domain_soft_disjoin(domains[i], nb_max_rays);
    }
    ++count;
  }
}

/**
 * \brief Fix the boundaries of the left input domain.
 *
 * \param[in] left        Input domain to fix.
 * \param[in] right       Input domain.
 * \param[in] nb_max_rays Max number of rays.
 *
 * \return The fixed domain.
 *
 * \warning \a left might be unusable afterwards. Always use the return value.
 * \details This function fixes the boundaries of the \a left input domain by
 *          removing the right domain from it, when necessary.
 */
Polyhedron* fix_boundaries(Polyhedron* left, Polyhedron* const right,
                           const unsigned nb_max_rays) {
  Polyhedron* result = left;
  Polyhedron* intersection = DomainIntersection(left, right, nb_max_rays);
  if (!emptyQ(intersection)) {
    result = DomainDifference(left, right, nb_max_rays);
    Domain_Free(left);
  }
  Domain_Free(intersection);

  return result;
}

/**
 * \brief Fix iterators.
 *
 * This function adds as many empty columns as there are missing iterators.
 *
 * \param[in] m            Input matrix.
 * \param[in] nb_iterators The number of iterators.
 * \param[in] depth        The current depth in the algorithm.
 *
 * \return The resulting constraints matrix.
 *
 * \warning \a m will be unusable afterwards. Always use the returned value.
 */
Matrix* matrix_fix_iterators(Matrix* m, const size_t nb_iterators,
                             const size_t depth) {
  Matrix* result = m, *temp = NULL;
  size_t i = 0;
  /* Temporarily move the (in)equality column and existing iterators to the
   * end. */
  PutColumnLast (result, 0);
  for (i = 0; i < depth; ++i)
    PutColumnLast (result, (int) (depth - i) - 1);
  /* Add as many iterators as the current depth. */
  for (i = 0; i < nb_iterators - depth; ++i) {
    temp = AddANullColumn (result);
    PutColumnFirst (temp, (int) temp->NbColumns - 1);
    Matrix_Free (result);
    result = temp;
  }
  /* Move back all temporarily moved columns. */
  for (i = 0; i < depth; ++i)
    PutColumnFirst (result, (int) result->NbColumns - 1);
  PutColumnFirst (result, (int) result->NbColumns - 1);

  return result;
}

/**
 * \brief Fix iterators in the input domain and mix it with the original domain.
 *
 * \param[in] original     Original domain.
 * \param[in] domain       Input domain.
 * \param[in] nb_iterators Iterator number.
 * \param[in] depth        Current depth in the algorithm.
 * \param[in] nb_max_rays  Maximum number of rays.
 *
 * \return The resulting domain.
 *
 * \warning The input domain is freed.
 */
Polyhedron* fix_iterators_and_mix(Polyhedron* original, Polyhedron* domain,
                                  const size_t nb_iterators, const size_t depth,
                                  const unsigned nb_max_rays) {
  Matrix* constraints = NULL;
  Polyhedron* fixed, *old_fixed, *current, *old_result, *result;

  fixed = old_fixed = current = old_result = NULL;
  result = Domain_Copy(original);

  if (domain) {
    constraints = Polyhedron2Constraints(domain);
    constraints = matrix_fix_iterators(constraints, nb_iterators, depth);
    fixed = Constraints2Polyhedron(constraints, nb_max_rays);
    Matrix_Free(constraints);

    for (Polyhedron* p = domain->next; p; p = p->next) {
      constraints = Polyhedron2Constraints(p);
      constraints = matrix_fix_iterators(constraints, nb_iterators, depth);
      current = Constraints2Polyhedron(constraints, nb_max_rays);
      Matrix_Free(constraints);

      old_fixed = fixed;
      fixed = DomainUnion(fixed, current, nb_max_rays);
      Domain_Free(old_fixed);
      Domain_Free(current);
    }

    old_result = result;
    result = DomainIntersection(result, fixed, nb_max_rays);
    Domain_Free(old_result);
    Domain_Free(fixed);
    Domain_Free(domain);
  }

  return result;
}

/**
 * \brief Get new domains from a parameterized polyhedron.
 *
 * \param[in] original     Original domain (which was used to compute \a PA).
 * \param[in] PA           Parameterized polyhedron.
 * \param[in] nb_domains   The number of new domains.
 * \param[in] nb_iterators The number of iterators
 * \param[in] depth        The current depth.
 * \param[in] nb_max_rays  The maximum number of rays.
 * \param[in,out] domains  The new domains.
 *
 * \note The input parameter \a original is not modified.
 * \warning \a domains must have been allocated prior the call to this funcion.
 */
void get_new_domains (
    Polyhedron* original, const Param_Polyhedron* PA, const size_t nb_domains,
    const size_t nb_iterators, const size_t depth, const unsigned nb_max_rays,
    Polyhedron** const domains) {
  /* The domains we get via Polyhedron2Param_Domain() may have overlapping
   * boundaries which must be fixed. Then, the iterators must be rearranged
   * (because their order had been modified in order to use
   * Polyhedron2Param_Domain()). Finally, each new domain's constraints are
   * added to the original domain.
   */
  size_t i = 0;
  if (PA) {
    extract_domains(PA, nb_domains, nb_max_rays, domains);
    for (i = 0; i < nb_domains; ++i)
      domains[i] = fix_iterators_and_mix(original, domains[i], nb_iterators,
                                         depth, nb_max_rays);
  } else {
    domains[0] = Domain_Copy(original);
  }
}
