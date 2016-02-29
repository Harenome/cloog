/**
 * \file utilities.c
 * \brief Various utilities to ease the use of the PolyLib in CLooG.
 *
 * This file defines various utilities that shall help using the PolyLib in
 * CLooG.
 */

#include "cloog/polylib/utilities.h"

#include "cloog/polylib/cloog.h"

/*******************************************************************************
 * Static functions prototypes                                                 *
 ******************************************************************************/

/**
 * \brief Compute the difference 'p1 \ p2', free the input domain.
 *
 * \param[in] p1          Input domain.
 * \param[in] p2          Input polyhedron.
 * \param[in] nb_max_rays Maximum number of rays.
 *
 * \return The difference 'p1 \ p2'.
 *
 * \details This function computes the result of
 *          cloog_polylib_compute_difference(), frees \a p1 and returns the
 *          aforementioned result.
 *
 * \see cloog_polylib_compute_difference_const().
 */
static Polyhedron* cloog_polylib_compute_difference(Polyhedron* p1,
                                                    Polyhedron* p2,
                                                    unsigned nb_max_rays);
/**
 * \brief Compute the difference 'p1 \ p2'.
 *
 * \param[in] p1          Input domain.
 * \param[in] p2          Input polyhedron.
 * \param[in] nb_max_rays Maximum number of rays.
 *
 * \return The difference 'p1 \ p2'.
 *
 * \details This function computes the difference 'p1 \ p2'. It expects \a p1
 *          to be a domain, and \a p2 to be a single polyhedron. Neither \a p1,
 *          nor \a p2 *should* be modified.
 *
 * \see cloog_polylib_compute_difference_const().
 */
static Polyhedron* cloog_polylib_compute_difference_const(Polyhedron* p1,
                                                          Polyhedron* p2,
                                                          unsigned nb_max_rays);

/*******************************************************************************
 * cloog/polylib/utilities.h functions                                         *
 ******************************************************************************/

Polyhedron* cloog_polylib_domain_soft_disjoin(Polyhedron* domain,
                                              const unsigned rays) {
  Polyhedron* result = cloog_polylib_domain_soft_disjoint_const(domain, rays);
  Domain_Free(domain);

  return result;
}

Polyhedron* cloog_polylib_domain_soft_disjoint_const(Polyhedron* P,
                                                     unsigned nb_max_rays) {
  if (!P)
    return (Polyhedron*) 0;
  if (!P->next)
    return Polyhedron_Copy(P);

  Polyhedron* result = (Polyhedron *) 0;
  for (Polyhedron* p1 = P; p1; p1 = p1->next) {
    Polyhedron* remainder = Polyhedron_Copy(p1);
    for (Polyhedron* p2 = p1->next; p2; p2 = p2->next)
      remainder = cloog_polylib_compute_difference(remainder, p2, nb_max_rays);
    remainder = cloog_polylib_domain_free_empty(remainder);
    if (!remainder)
      continue;
    for (Polyhedron* next, * current = remainder; current; current = next) {
      next = current->next;
      current->next = NULL;
      result = AddPolyToDomain(current, result);
    }
  }

  return result;
}

Polyhedron* cloog_polylib_domain_convex_attempt(Polyhedron* domain,
                                                unsigned nb_max_rays) {
  Polyhedron* result = NULL, * discard = NULL;
  Polyhedron* convex = DomainConvex(domain, nb_max_rays);
  Polyhedron* difference = DomainDifference(convex, domain, nb_max_rays);

  if (!emptyQ(difference)) {
    result = domain;
    discard = convex;
  } else {
    result = convex;
    discard = domain;
  }

  Domain_Free(difference);
  Domain_Free(discard);

  return result;
}

Polyhedron* cloog_polylib_domain_convex_attempt_const(Polyhedron* domain,
                                                      unsigned nb_max_rays) {
  Polyhedron* result = domain;
  Polyhedron* convex = DomainConvex(domain, nb_max_rays);
  Polyhedron* difference = DomainDifference(convex, domain, nb_max_rays);

  if (emptyQ(difference)) {
    result = convex;
  } else {
    Domain_Free(convex);
    result = convex;
  }
  Domain_Free(difference);

  return result;
}

size_t cloog_polylib_domain_count(const Polyhedron* domain) {
  size_t count = 0;
  for (const Polyhedron* p = domain; p; p = p->next)
    ++count;

  return count;
}

Matrix** cloog_polylib_domain_get_constraints(Polyhedron* domain,
                                              const size_t nb_constraints) {
  Polyhedron* copy = Domain_Copy(domain);
  Polyhedron* p = copy;
  Matrix** constraints = malloc(nb_constraints * (sizeof *constraints));
  for (size_t i = 0; p && i < nb_constraints; ++i) {
    constraints[i] = Polyhedron2Constraints(p);
    p = p->next;
  }
  Domain_Free(copy);

  return constraints;
}

Polyhedron* cloog_polylib_domain_last(Polyhedron* domain) {
  if (domain)
    while (domain->next)
      domain = domain->next;

  return domain;
}

Polyhedron* cloog_polylib_domain_free_empty(Polyhedron* domain) {
  Polyhedron* result = domain;
  if (domain && emptyQ(domain)) {
    Domain_Free(domain);
    result = (Polyhedron*) 0;
  }

  return result;
}

CloogLoop* cloog_loop_from_polyhedron(Polyhedron* domain, CloogState* state,
                                      CloogLoop* source, const size_t nb_par) {
  CloogLoop* result = cloog_loop_malloc(state);
  result->domain = cloog_domain_from_polylib_polyhedron(state, domain, nb_par);
  result->otl = source->otl;
  result->stride = cloog_stride_copy(source->stride);
  result->block = cloog_block_copy(source->block);

  return result;
}

/*******************************************************************************
 * Static functions definitions                                                *
 ******************************************************************************/

Polyhedron* cloog_polylib_compute_difference(Polyhedron* p1, Polyhedron* p2,
                                             const unsigned nb_max_rays) {
  Polyhedron* res = cloog_polylib_compute_difference_const(p1, p2, nb_max_rays);
  Domain_Free(p1);

  return res;
}

Polyhedron* cloog_polylib_compute_difference_const(Polyhedron* p1,
                                                   Polyhedron* p2,
                                                   const unsigned nb_max_rays) {
  /* Most of this implementation is inspired from the PolyLib function
   * Disjoint_Domain() (see the parts where d1 and d2 are computed, as of
   * PolyLib version 5.22.5).
   */

  Polyhedron *result, *new, *current, *temp;
  result = current = new = temp = (Polyhedron*) 0;

  for (Polyhedron* p = p1; p; p = p->next) {
    current = p;
    for (int i = 0; i < p2->NbConstraints && current; ++i) {
      /* Add '-p2->Constraint[i] >= 1'. */
      temp = SubConstraint(p2->Constraint[i], current, nb_max_rays, 0);
      result = AddPolyToDomain(temp, result);
      if (value_zero_p(p2->Constraint[i][0])) {
        /* Add '+p2->Constraint[i] >= 1'. */
        temp = SubConstraint(p2->Constraint[i], current, nb_max_rays, 1);
        result = AddPolyToDomain(temp, result);
        /* Add constraint 'p2->constraint[i] == 0'. */
        new = AddConstraints(p2->Constraint[i], 1, current, nb_max_rays);
      } else {
        /* Add constraint '+p2->constraint[i] >= 0'. */
        new = SubConstraint(p2->Constraint[i], current, nb_max_rays, 3);
      }
      new = cloog_polylib_domain_free_empty(new);
      if (current != p)
        Domain_Free(current);
      current = new;
    }
    if (current != p)
      Domain_Free(current);
  }

  return result;
}

