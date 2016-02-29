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
#define MAX_RAYS 2048

/*******************************************************************************
 * Static functions prototypes                                                 *
 ******************************************************************************/

/* Utilities for CloogLoop. */
static inline CloogLoop* cloog_loop_last(CloogLoop*);
static inline CloogLoop* cloog_loop_copy_current_loop(CloogLoop*);

/* Parameterized domain. */
static Polyhedron** flatten_domains(Polyhedron**, unsigned, size_t*);
static Polyhedron* domain_from_constraints(size_t, Matrix*[], unsigned);
static Matrix* matrix_add_parameters(Matrix*, size_t);
static Polyhedron* domain_add_parameters(Polyhedron*, size_t, unsigned);
static Matrix* matrix_move_iterators_last(Matrix*, size_t, size_t);
static Polyhedron* domain_move_iterators_last(Polyhedron*, size_t, size_t,
                                              unsigned);
static Param_Polyhedron* compute_param_domain(Polyhedron*, Polyhedron*, size_t,
                                              size_t, unsigned);

/* New splits. */
static size_t param_domain_count_domains(const Param_Domain*);
static Polyhedron** extract_domains(const Param_Polyhedron*, size_t, unsigned);
static Polyhedron* fix_boundaries(Polyhedron*, Polyhedron*, unsigned);
static Matrix* matrix_fix_iterators(Matrix*, size_t, size_t);
static Polyhedron* fix_iterators_and_mix(Polyhedron*, Polyhedron*, size_t,
                                         size_t, unsigned);
static Polyhedron** get_new_domains(Polyhedron*, const Param_Polyhedron*,
                                    size_t, size_t, size_t, unsigned);


/*******************************************************************************
 * cloog/polylib/special.h functions                                           *
 ******************************************************************************/

CloogLoop* cloog_loop_polylib_split(CloogState* state, CloogDomain* context,
                                    CloogLoop* loop, const size_t max_depth,
                                    const unsigned nb_max_rays) {
   /* CloogLoop is a chained struct. Some structs will be changed, some will
    * remain intact. cloog_loop_polylib_split_current_loop() shall never modify
    * nor return any of the structs from `loop`. It must either return a copy
    * of the current struct (in the case where it is not possible to split it)
    * or the splits of the current struct. This way, `loop` can safely be freed
    * as a whole with cloog_loop_free() at the end of this function.
    *
    * (cloog_loop_free() frees the current struct, then moves on to
    * the struct's `next` field. However, at each step, we don't know whether
    * the current struct has been split.)
    */

  CloogLoop* result = loop;

  if (loop && loop->domain) {
    result = cloog_loop_polylib_split_current_loop(state, context, loop,
                                                   max_depth, nb_max_rays);
    CloogLoop* last = cloog_loop_last(result);
    for (CloogLoop* current = loop->next; current; current = current->next) {
      last->next = cloog_loop_polylib_split_current_loop(state, context,
                                                         current, max_depth,
                                                         nb_max_rays);
      last = cloog_loop_last(last);
    }
    cloog_loop_free(loop);
  }

  return result;
}

CloogLoop* cloog_loop_polylib_split_current_loop(CloogState* state,
                                                 CloogDomain* context,
                                                 CloogLoop* loop,
                                                 const size_t max_depth,
                                                 const unsigned nb_max_rays) {
  /* The first step is to copy the input `loop` because it is important that the
   * input `loop` remains unmodified.
   * If the loop can be split, this copy will be freed and the splits will be
   * returned. Otherwise, the copy will be returned.
   */

  CloogLoop* result = cloog_loop_copy_current_loop(loop);
  const size_t nb_par = context->nb_par;
  Polyhedron* const context_p = cloog_domain_polyhedron(context);
  Polyhedron* const domain_p = cloog_domain_polyhedron(loop->domain);

  size_t nb_splits = 0;
  Polyhedron** splits = cloog_polylib_split_polyhedron(context_p, domain_p, nb_par,
                                               max_depth, nb_max_rays,
                                               &nb_splits);
  if (splits && nb_splits) {
    CloogLoop* old = result;
    result = cloog_loop_from_polyhedron(splits[0], state, loop, nb_par);
    CloogLoop* last = cloog_loop_last(result);
    for (size_t i = 1; i < nb_splits; ++i) {
      last->next = cloog_loop_from_polyhedron(splits[i], state, loop, nb_par);
      last = cloog_loop_last(last);
    }
    free(splits);
    cloog_loop_free(old);
  }

  return result;
}

Polyhedron** cloog_polylib_split_polyhedron(Polyhedron* const context,
                                            Polyhedron* const domain,
                                            const size_t nb_parameters,
                                            const size_t max_depth,
                                            const unsigned nb_max_rays,
                                            size_t* const nb_new_splits) {
  /* Summary:
   * ========
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
   * using two Polyhedron* arrays (the pointers 'splits' and 'new_domains') that
   * are realloc'd. Domains that must be split are stored in splits. If a domain
   * can not be split, it is moved into new_domains, otherwise it is freed and
   * its splits are stored into new_domains. At the end of each iteration of the
   * loop on the depth, the two arrays are swapped/freed.
   */

  const size_t nb_iterators = domain->Dimension - nb_parameters;
  const size_t max_level =
    (max_depth && max_depth < nb_iterators) ? (max_depth + 1) : nb_iterators;

  /* To both store the remaining domains to process and the final results. */
  size_t nb_splits = 1;
  Polyhedron** splits = malloc(nb_splits * (sizeof *splits));
  /* To temporarily store the new domains (or intact domains at a given step. */
  size_t nb_new_domains = 0;
  Polyhedron** new_domains = NULL;
  /* Clean up the input, if needed. */
  splits[0] = Domain_Copy(domain);
  splits[0] = cloog_polylib_domain_convex_attempt(splits[0], nb_max_rays);

  for (size_t depth = 1; depth < max_level; ++depth) {
    splits = flatten_domains(splits, nb_max_rays, &nb_splits);
    for (size_t j = 0; j < nb_splits; ++j) {
      /* Compute the Param_Polyhedron. */
      Param_Polyhedron* PA = compute_param_domain(context, splits[j],
                                                  nb_iterators, depth,
                                                  nb_max_rays);
      if (PA) {
        size_t nb_temp_domains = param_domain_count_domains(PA->D);
        Polyhedron** temp_domains = get_new_domains(splits[j], PA,
                                                    nb_temp_domains,
                                                    nb_iterators, depth,
                                                    nb_max_rays);
        new_domains = realloc(new_domains, (nb_new_domains + nb_temp_domains) *
                              (sizeof *new_domains));
        for (size_t k = 0; k < nb_temp_domains; ++k)
          new_domains[nb_new_domains + k] = temp_domains[k];
        nb_new_domains += nb_temp_domains;

        free(temp_domains);
        Domain_Free(splits[j]);
        Param_Polyhedron_Free(PA);
      } else {
        new_domains = realloc(new_domains,
                              (nb_new_domains + 1) * (sizeof *new_domains));
        new_domains[nb_new_domains] = splits[j];
        ++nb_new_domains;
      }
    }

    if (new_domains && nb_new_domains > 0) {
      free(splits);
      splits = new_domains;
      nb_splits = nb_new_domains;
      new_domains = NULL;
      nb_new_domains = 0;
    }
  }

  *nb_new_splits = nb_splits;
  return splits;
}

/*******************************************************************************
 * Static functions definitions                                                *
 ******************************************************************************/

/* Utilities for CloogLoop. */

/**
 * \brief Get the last CloogLoop.
 *
 * \param[in] loop The input CloogLoop.
 *
 * \return The last CloogLoop.
 */
CloogLoop* cloog_loop_last(CloogLoop* loop) {
  CloogLoop* last = loop;
  while (last->next)
    last = last->next;

  return last;
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
  result->otl = loop->otl;
  result->stride = cloog_stride_copy(loop->stride);
  result->block = cloog_block_copy(loop->block);

  return result;
}

/* Parameterized domain. */

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
Polyhedron** flatten_domains(Polyhedron** domains, const unsigned nb_max_rays,
                             size_t* const nb_domains) {
  /* Flattening: Separate unions of domains into single domains.
   * For each union of domains, we first attempt to find the convex domain and
   * then compute the equivalent disjoint union to ensure they can safely be
   * separated.
   */
  const size_t nb_domains_local = *nb_domains;

  size_t position = 0;
  size_t nb_flattened = *nb_domains;
  Polyhedron** result = malloc (nb_flattened * (sizeof *result));
  for (size_t i = 0; i < nb_domains_local; ++i) {
    if (domains[i]->next) {
      domains[i] = cloog_polylib_domain_convex_attempt(domains[i], nb_max_rays);
      domains[i] = cloog_polylib_domain_soft_disjoin(domains[i], nb_max_rays);
      Polyhedron* next = NULL;
      nb_flattened += cloog_polylib_domain_count(domains[i]) - 1;
      result = realloc (result, nb_flattened * (sizeof *result));
      for (Polyhedron* p = domains[i]; p; p = next) {
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

  *nb_domains = nb_flattened;
  return result;
}

/**
 * \brief Convert an array of constraints into a domain.
 *
 * \param[in] nb_constraints Size of the input array.
 * \param[in] constraints    The input array of constraints.
 *
 * \return The resulting domain.
 *
 * \warning The input array will be freed. Do not use it afterwards.
 */
Polyhedron* domain_from_constraints(const size_t nb_constraints,
                                    Matrix* constraints[nb_constraints],
                                    const unsigned nb_max_rays) {
  Polyhedron* result = Constraints2Polyhedron(constraints[0], nb_max_rays);
  Matrix_Free(constraints[0]);
  for (size_t i = 1; i < nb_constraints; ++i) {
    Polyhedron* old_result = result;
    result = DomainAddConstraints(result, constraints[i], nb_max_rays);
    Domain_Free(old_result);
    Matrix_Free(constraints[i]);
  }
  free(constraints);
  if (result->next)
    result = cloog_polylib_domain_soft_disjoin(result, nb_max_rays);

  return result;
}

/**
 * \brief Add parameters to a constraints matrix.
 *
 * \param[in] m Input matrix.
 * \param[in] n The number of parameters to add.
 *
 * \return The resulting matrix.
 *
 * \warning Do not use \a m aftewards.
 * \details This function adds \a n empty columns to the input matrix right
 *          after the (in)equality column.
 */
Matrix* matrix_add_parameters(Matrix* m, const size_t n) {
  /* Move the (in)equality column to the end. */
  PutColumnLast(m, 0);
  /* Add new columns and move them to the start. */
  for (size_t i = 0; i < n; ++i) {
    Matrix* const temp = AddANullColumn(m);
    PutColumnFirst(temp, (int) temp->NbColumns - 1);
    Matrix_Free(m);
    m = temp;
  }
  /* Move the (in)equality column back to the start. */
  PutColumnFirst(m, (int) m->NbColumns - 1);

  return m;
}

/**
 * \brief Add parameters to a domain.
 *
 * \param[in] domain Input domain.
 * \param[in] n      The number of parameters to add.
 *
 * \return The resulting domain.
 */
Polyhedron* domain_add_parameters(Polyhedron* domain, const size_t n,
                                  const unsigned nb_max_rays) {
  size_t nb_constraints = cloog_polylib_domain_count(domain);
  Matrix** constraints = cloog_polylib_domain_get_constraints(domain,
                                                              nb_constraints);
  for (size_t i = 0; i < nb_constraints; ++i)
    constraints[i] = matrix_add_parameters(constraints[i], n);
  Polyhedron* result = domain_from_constraints(nb_constraints, constraints,
                                               nb_max_rays);

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
  for (size_t i = 0; i < n; ++i)
    for (int j = 1; j < (int) (nb_iterators - i); ++j)
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
Polyhedron* domain_move_iterators_last(Polyhedron* const domain,
                                       const size_t nb_iterators,
                                       const size_t n,
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
Param_Polyhedron* compute_param_domain(Polyhedron* const original_context,
                                       Polyhedron* const original_domain,
                                       const size_t nb_iterators,
                                       const size_t depth,
                                       const unsigned nb_max_rays) {
  /* This step temporarily "transforms" iterators into parameters:
   * - move their columns after the remaining iterators in the domain.
   * - add empty columns to the context.
   * Finally, the Param_Polyhedron is computed.
   */
  Polyhedron* context = domain_add_parameters(original_context, depth,
                                              nb_max_rays);
  Polyhedron* domain = domain_move_iterators_last(original_domain, nb_iterators,
                                                  depth, nb_max_rays);
  Param_Polyhedron* PA = Polyhedron2Param_Domain(domain, context, WS);
  Domain_Free(context);
  Domain_Free(domain);

  return PA;
}

/* New splits. */

/**
 * \brief Count the domains of a Param_Domain.
 *
 * \param[in] domain Input Param_Domain.
 *
 * \return The number of domains.
 */
size_t param_domain_count_domains(const Param_Domain* domain) {
  size_t count = 0;
  for (const Param_Domain* P = domain; P; P = P->next)
    ++count;

  return count;
}

/**
 * \brief Extract the domains of a Param_Polyhedron.
 *
 * \param[in] PA          The input Param_Polyhedron.
 * \param[in] nb_domains  The number of domains to extract.
 * \param[in] nb_max_rays Maximum number of rays
 *
 * \return The domains.
 */
Polyhedron** extract_domains(const Param_Polyhedron* PA,
                             const size_t nb_domains,
                             const unsigned nb_max_rays) {
  size_t count = 0;
  Polyhedron** const domains = malloc(nb_domains * (sizeof *domains));
  for (const Param_Domain* P = PA->D; P; P = P->next) {
    domains[count] = cloog_polylib_domain_soft_disjoint_const(P->Domain,
                                                              nb_max_rays);
    for (size_t i = 0; i < count; ++i) {
      domains[i] = fix_boundaries(domains[i], domains[count], nb_max_rays);
      if (domains[i]->next)
        domains[i] = cloog_polylib_domain_soft_disjoin(domains[i], nb_max_rays);
    }
    ++count;
  }

  return domains;
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
  if (! emptyQ (intersection)) {
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
  Matrix* result = m;
  /* Temporarily move the (in)equality column and existing iterators to the
   * end. */
  PutColumnLast (result, 0);
  for (size_t i = 0; i < depth; ++i)
    PutColumnLast (result, (int) (depth - i) - 1);
  /* Add as many iterators as the current depth. */
  for (size_t i = 0; i < nb_iterators - depth; ++i) {
    Matrix* temp = AddANullColumn (result);
    PutColumnFirst (temp, (int) temp->NbColumns - 1);
    Matrix_Free (result);
    result = temp;
  }
  /* Move back all temporarily moved columns. */
  for (size_t i = 0; i < depth; ++i)
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
  Polyhedron* result = Domain_Copy(original);
  if (domain) {
    Matrix* constraints = Polyhedron2Constraints(domain);
    constraints = matrix_fix_iterators(constraints, nb_iterators, depth);
    Polyhedron* fixed = Constraints2Polyhedron(constraints, nb_max_rays);
    Matrix_Free(constraints);

    for (Polyhedron* p = domain->next; p; p = p->next) {
      constraints = Polyhedron2Constraints(p);
      constraints = matrix_fix_iterators(constraints, nb_iterators, depth);
      Polyhedron* current = Constraints2Polyhedron(constraints, nb_max_rays);
      Matrix_Free(constraints);

      Polyhedron* old_fixed = fixed;
      fixed = DomainUnion(fixed, current, nb_max_rays);
      Domain_Free(old_fixed);
      Domain_Free(current);
    }

    Polyhedron* old_result = result;
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
 *-
 * \return The new domains.
 *
 * \note The input parameter \a original is not modified.
 */
Polyhedron** get_new_domains (Polyhedron* original, const Param_Polyhedron* PA,
                              const size_t nb_domains,
                              const size_t nb_iterators, const size_t depth,
                              const unsigned nb_max_rays) {
  /* The domains we get via Polyhedron2Param_Domain() may have overlapping
   * boundaries which must be fixed. Then, the iterators must be rearranged
   * (because their order had been modified in order to use
   * Polyhedron2Param_Domain()). Finally, each new domain's constraints are
   * added to the original domain.
   */
  Polyhedron** domains = extract_domains(PA, nb_domains, nb_max_rays);
  for (size_t i = 0; i < nb_domains; ++i)
    domains[i] = fix_iterators_and_mix(original, domains[i], nb_iterators,
                                       depth, nb_max_rays);

  return domains;
}
