#ifndef __GSL_INTEGRATION_H__
#define __GSL_INTEGRATION_H__
#include <config.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/* Workspace for adaptive integrators */

typedef struct
  {
    size_t limit;
    size_t size;
    size_t nrmax;
    size_t i;
    size_t maximum_level;
    double *alist;
    double *blist;
    double *rlist;
    double *elist;
    size_t *order;
    size_t *level;
  }
gsl_integration_workspace;

gsl_integration_workspace *
  gsl_integration_workspace_alloc (const size_t n);

void
  gsl_integration_workspace_free (gsl_integration_workspace * w);

/* Definition of an integration rule */

typedef void gsl_integration_rule (const gsl_function * f,
                                   double a, double b,
                                   double *result, double *abserr,
                                   double *defabs, double *resabs);

void gsl_integration_qk15 (const gsl_function * f, double a, double b,
                           double *result, double *abserr,
                           double *resabs, double *resasc);

enum
  {
    GSL_INTEG_GAUSS15 = 1,      /* 15 point Gauss-Kronrod rule */
  };

void 
gsl_integration_qk (const int n, const double xgk[], 
                    const double wg[], const double wgk[],
                    double fv1[], double fv2[],
                    const gsl_function *f, double a, double b,
                    double * result, double * abserr, 
                    double * resabs, double * resasc);


int gsl_integration_qag (const gsl_function * f,
                         double a, double b,
                         double epsabs, double epsrel, size_t limit,
                         int key,
                         gsl_integration_workspace * workspace,
                         double *result, double *abserr);

__END_DECLS

//CONTENTS OF QAG.C
static int
qag (const gsl_function *f,
     const double a, const double b,
     const double epsabs, const double epsrel,
     const size_t limit,
     gsl_integration_workspace * workspace,
     double * result, double * abserr,
     gsl_integration_rule * q) ;

int
gsl_integration_qag (const gsl_function *f,
                     double a, double b,
                     double epsabs, double epsrel, size_t limit,
                     int key,
                     gsl_integration_workspace * workspace,
                     double * result, double * abserr)
{
  int status ;
  gsl_integration_rule * integration_rule = gsl_integration_qk15 ;

  status = qag (f, a, b, epsabs, epsrel, limit,
                workspace, 
                result, abserr, 
                integration_rule) ;
  
  return status ;
}

static int
qag (const gsl_function * f,
     const double a, const double b,
     const double epsabs, const double epsrel,
     const size_t limit,
     gsl_integration_workspace * workspace,
     double *result, double *abserr,
     gsl_integration_rule * q)
{
  double area, errsum;
  double result0, abserr0, resabs0, resasc0;
  double tolerance;
  size_t iteration = 0;
  int roundoff_type1 = 0, roundoff_type2 = 0, error_type = 0;

  double round_off;     

  /* Initialize results */

  initialise (workspace, a, b);

  *result = 0;
  *abserr = 0;

  if (limit > workspace->limit)
    {
      GSL_ERROR ("iteration limit exceeds available workspace", GSL_EINVAL) ;
    }

  if (epsabs <= 0 && (epsrel < 50 * GSL_DBL_EPSILON || epsrel < 0.5e-28))
    {
      GSL_ERROR ("tolerance cannot be achieved with given epsabs and epsrel",
                 GSL_EBADTOL);
    }

  /* perform the first integration */

  q (f, a, b, &result0, &abserr0, &resabs0, &resasc0);

  set_initial_result (workspace, result0, abserr0);

  /* Test on accuracy */

  tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (result0));

  /* need IEEE rounding here to match original quadpack behavior */

  round_off = GSL_COERCE_DBL (50 * GSL_DBL_EPSILON * resabs0);

  if (abserr0 <= round_off && abserr0 > tolerance)
    {
      *result = result0;
      *abserr = abserr0;

      GSL_ERROR ("cannot reach tolerance because of roundoff error "
                 "on first attempt", GSL_EROUND);
    }
  else if ((abserr0 <= tolerance && abserr0 != resasc0) || abserr0 == 0.0)
    {
      *result = result0;
      *abserr = abserr0;

      return GSL_SUCCESS;
    }
  else if (limit == 1)
    {
      *result = result0;
      *abserr = abserr0;

      GSL_ERROR ("a maximum of one iteration was insufficient", GSL_EMAXITER);
    }

  area = result0;
  errsum = abserr0;

  iteration = 1;

  do
    {
      double a1, b1, a2, b2;
      double a_i, b_i, r_i, e_i;
      double area1 = 0, area2 = 0, area12 = 0;
      double error1 = 0, error2 = 0, error12 = 0;
      double resasc1, resasc2;
      double resabs1, resabs2;

      /* Bisect the subinterval with the largest error estimate */

      retrieve (workspace, &a_i, &b_i, &r_i, &e_i);

      a1 = a_i; 
      b1 = 0.5 * (a_i + b_i);
      a2 = b1;
      b2 = b_i;

      q (f, a1, b1, &area1, &error1, &resabs1, &resasc1);
      q (f, a2, b2, &area2, &error2, &resabs2, &resasc2);

      area12 = area1 + area2;
      error12 = error1 + error2;

      errsum += (error12 - e_i);
      area += area12 - r_i;

      if (resasc1 != error1 && resasc2 != error2)
        {
          double delta = r_i - area12;

          if (fabs (delta) <= 1.0e-5 * fabs (area12) && error12 >= 0.99 * e_i)
            {
              roundoff_type1++;
            }
          if (iteration >= 10 && error12 > e_i)
            {
              roundoff_type2++;
            }
        }

      tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (area));

      if (errsum > tolerance)
        {
          if (roundoff_type1 >= 6 || roundoff_type2 >= 20)
            {
              error_type = 2;   /* round off error */
            }

          /* set error flag in the case of bad integrand behaviour at
             a point of the integration range */

          if (subinterval_too_small (a1, a2, b2))
            {
              error_type = 3;
            }
        }

      update (workspace, a1, b1, area1, error1, a2, b2, area2, error2);

      retrieve (workspace, &a_i, &b_i, &r_i, &e_i);

      iteration++;

    }
  while (iteration < limit && !error_type && errsum > tolerance);

  *result = sum_results (workspace);
  *abserr = errsum;

  if (errsum <= tolerance)
    {
      return GSL_SUCCESS;
    }
  else if (error_type == 2)
    {
      GSL_ERROR ("roundoff error prevents tolerance from being achieved",
                 GSL_EROUND);
    }
  else if (error_type == 3)
    {
      GSL_ERROR ("bad integrand behavior found in the integration interval",
                 GSL_ESING);
    }
  else if (iteration == limit)
    {
      GSL_ERROR ("maximum number of subdivisions reached", GSL_EMAXITER);
    }
  else
    {
      GSL_ERROR ("could not integrate function", GSL_EFAILED);
    }
}

//CONTENTS OF WORKSPACE.C
gsl_integration_workspace *
gsl_integration_workspace_alloc (const size_t n) 
{
  gsl_integration_workspace * w ;
  
  if (n == 0)
    {
      GSL_ERROR_VAL ("workspace length n must be positive integer",
                        GSL_EDOM, 0);
    }

  w = (gsl_integration_workspace *) 
    malloc (sizeof (gsl_integration_workspace));

  if (w == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for workspace struct",
                        GSL_ENOMEM, 0);
    }

  w->alist = (double *) malloc (n * sizeof (double));

  if (w->alist == 0)
    {
      free (w);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for alist ranges",
                        GSL_ENOMEM, 0);
    }

  w->blist = (double *) malloc (n * sizeof (double));

  if (w->blist == 0)
    {
      free (w->alist);
      free (w);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for blist ranges",
                        GSL_ENOMEM, 0);
    }

  w->rlist = (double *) malloc (n * sizeof (double));

  if (w->rlist == 0)
    {
      free (w->blist);
      free (w->alist);
      free (w);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for rlist ranges",
                        GSL_ENOMEM, 0);
    }


  w->elist = (double *) malloc (n * sizeof (double));

  if (w->elist == 0)
    {
      free (w->rlist);
      free (w->blist);
      free (w->alist);
      free (w);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for elist ranges",
                        GSL_ENOMEM, 0);
    }

  w->order = (size_t *) malloc (n * sizeof (size_t));

  if (w->order == 0)
    {
      free (w->elist);
      free (w->rlist);
      free (w->blist);
      free (w->alist);
      free (w);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for order ranges",
                        GSL_ENOMEM, 0);
    }

  w->level = (size_t *) malloc (n * sizeof (size_t));

  if (w->level == 0)
    {
      free (w->order);
      free (w->elist);
      free (w->rlist);
      free (w->blist);
      free (w->alist);
      free (w);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for order ranges",
                        GSL_ENOMEM, 0);
    }

  w->size = 0 ;
  w->limit = n ;
  w->maximum_level = 0 ;
  
  return w ;
}

void
gsl_integration_workspace_free (gsl_integration_workspace * w)
{
  RETURN_IF_NULL (w);
  free (w->level) ;
  free (w->order) ;
  free (w->elist) ;
  free (w->rlist) ;
  free (w->blist) ;
  free (w->alist) ;
  free (w) ;
}

//CONTENTS OF QK15.C
/* Gauss quadrature weights and kronrod quadrature abscissae and
   weights as evaluated with 80 decimal digit arithmetic by
   L. W. Fullerton, Bell Labs, Nov. 1981. */

static const double xgk[8] =    /* abscissae of the 15-point kronrod rule */
{
  0.991455371120812639206854697526329,
  0.949107912342758524526189684047851,
  0.864864423359769072789712788640926,
  0.741531185599394439863864773280788,
  0.586087235467691130294144838258730,
  0.405845151377397166906606412076961,
  0.207784955007898467600689403773245,
  0.000000000000000000000000000000000
};

/* xgk[1], xgk[3], ... abscissae of the 7-point gauss rule. 
   xgk[0], xgk[2], ... abscissae to optimally extend the 7-point gauss rule */

static const double wg[4] =     /* weights of the 7-point gauss rule */
{
  0.129484966168869693270611432679082,
  0.279705391489276667901467771423780,
  0.381830050505118944950369775488975,
  0.417959183673469387755102040816327
};

static const double wgk[8] =    /* weights of the 15-point kronrod rule */
{
  0.022935322010529224963732008058970,
  0.063092092629978553290700663189204,
  0.104790010322250183839876322541518,
  0.140653259715525918745189590510238,
  0.169004726639267902826583426598550,
  0.190350578064785409913256402421014,
  0.204432940075298892414161999234649,
  0.209482141084727828012999174891714
};

void
gsl_integration_qk15 (const gsl_function * f, double a, double b,
      double *result, double *abserr,
      double *resabs, double *resasc)
{
  double fv1[8], fv2[8];
  gsl_integration_qk (8, xgk, wg, wgk, fv1, fv2, f, a, b, result, abserr, resabs, resasc);
}

//CONTENTS OF QK.C
void
gsl_integration_qk (const int n, 
                    const double xgk[], const double wg[], const double wgk[],
                    double fv1[], double fv2[],
                    const gsl_function * f, double a, double b,
                    double *result, double *abserr,
                    double *resabs, double *resasc)
{

  const double center = 0.5 * (a + b);
  const double half_length = 0.5 * (b - a);
  const double abs_half_length = fabs (half_length);
  const double f_center = GSL_FN_EVAL (f, center);

  double result_gauss = 0;
  double result_kronrod = f_center * wgk[n - 1];

  double result_abs = fabs (result_kronrod);
  double result_asc = 0;
  double mean = 0, err = 0;

  int j;

  if (n % 2 == 0)
    {
      result_gauss = f_center * wg[n / 2 - 1];
    }

  for (j = 0; j < (n - 1) / 2; j++)
    {
      const int jtw = j * 2 + 1;  /* in original fortran j=1,2,3 jtw=2,4,6 */
      const double abscissa = half_length * xgk[jtw];
      const double fval1 = GSL_FN_EVAL (f, center - abscissa);
      const double fval2 = GSL_FN_EVAL (f, center + abscissa);
      const double fsum = fval1 + fval2;
      fv1[jtw] = fval1;
      fv2[jtw] = fval2;
      result_gauss += wg[j] * fsum;
      result_kronrod += wgk[jtw] * fsum;
      result_abs += wgk[jtw] * (fabs (fval1) + fabs (fval2));
    }

  for (j = 0; j < n / 2; j++)
    {
      int jtwm1 = j * 2;
      const double abscissa = half_length * xgk[jtwm1];
      const double fval1 = GSL_FN_EVAL (f, center - abscissa);
      const double fval2 = GSL_FN_EVAL (f, center + abscissa);
      fv1[jtwm1] = fval1;
      fv2[jtwm1] = fval2;
      result_kronrod += wgk[jtwm1] * (fval1 + fval2);
      result_abs += wgk[jtwm1] * (fabs (fval1) + fabs (fval2));
    };

  mean = result_kronrod * 0.5;

  result_asc = wgk[n - 1] * fabs (f_center - mean);

  for (j = 0; j < n - 1; j++)
    {
      result_asc += wgk[j] * (fabs (fv1[j] - mean) + fabs (fv2[j] - mean));
    }

  /* scale by the width of the integration region */

  err = (result_kronrod - result_gauss) * half_length;

  result_kronrod *= half_length;
  result_abs *= abs_half_length;
  result_asc *= abs_half_length;

  *result = result_kronrod;
  *resabs = result_abs;
  *resasc = result_asc;
  *abserr = rescale_error (err, result_abs, result_asc);

}

//CONTENTS OF ERR.C
static double rescale_error (double err, const double result_abs, const double result_asc) ;

static double
rescale_error (double err, const double result_abs, const double result_asc)
{
  err = fabs(err) ;

  if (result_asc != 0 && err != 0)
      {
        double scale = pow((200 * err / result_asc), 1.5) ;
        
        if (scale < 1)
          {
            err = result_asc * scale ;
          }
        else 
          {
            err = result_asc ;
          }
      }
  if (result_abs > GSL_DBL_MIN / (50 * GSL_DBL_EPSILON))
    {
      double min_err = 50 * GSL_DBL_EPSILON * result_abs ;

      if (min_err > err) 
        {
          err = min_err ;
        }
    }
  
  return err ;
}

//CONTENTS OF SET_INTIAL.C
static inline
void set_initial_result (gsl_integration_workspace * workspace, 
                         double result, double error);

static inline
void set_initial_result (gsl_integration_workspace * workspace, 
                         double result, double error)
{
  workspace->size = 1;
  workspace->rlist[0] = result;
  workspace->elist[0] = error;
}

//CONTENTS OF INITIALISE.C
static inline
void initialise (gsl_integration_workspace * workspace, double a, double b);

static inline
void initialise (gsl_integration_workspace * workspace, double a, double b)
{
  workspace->size = 0;
  workspace->nrmax = 0;
  workspace->i = 0;
  workspace->alist[0] = a;
  workspace->blist[0] = b;
  workspace->rlist[0] = 0.0;
  workspace->elist[0] = 0.0;
  workspace->order[0] = 0;
  workspace->level[0] = 0;

  workspace->maximum_level = 0;
}

//CONTENTS OF UTIL.C
static inline
void update (gsl_integration_workspace * workspace,
                 double a1, double b1, double area1, double error1,
                 double a2, double b2, double area2, double error2);

static inline void
retrieve (const gsl_integration_workspace * workspace, 
          double * a, double * b, double * r, double * e);



static inline
void update (gsl_integration_workspace * workspace,
             double a1, double b1, double area1, double error1,
             double a2, double b2, double area2, double error2)
{
  double * alist = workspace->alist ;
  double * blist = workspace->blist ;
  double * rlist = workspace->rlist ;
  double * elist = workspace->elist ;
  size_t * level = workspace->level ;

  const size_t i_max = workspace->i ;
  const size_t i_new = workspace->size ;

  const size_t new_level = workspace->level[i_max] + 1;

  /* append the newly-created intervals to the list */
  
  if (error2 > error1)
    {
      alist[i_max] = a2;        /* blist[maxerr] is already == b2 */
      rlist[i_max] = area2;
      elist[i_max] = error2;
      level[i_max] = new_level;
      
      alist[i_new] = a1;
      blist[i_new] = b1;
      rlist[i_new] = area1;
      elist[i_new] = error1;
      level[i_new] = new_level;
    }
  else
    {
      blist[i_max] = b1;        /* alist[maxerr] is already == a1 */
      rlist[i_max] = area1;
      elist[i_max] = error1;
      level[i_max] = new_level;
      
      alist[i_new] = a2;
      blist[i_new] = b2;
      rlist[i_new] = area2;
      elist[i_new] = error2;
      level[i_new] = new_level;
    }
  
  workspace->size++;

  if (new_level > workspace->maximum_level)
    {
      workspace->maximum_level = new_level;
    }

  qpsrt (workspace) ;
}

static inline void
retrieve (const gsl_integration_workspace * workspace, 
          double * a, double * b, double * r, double * e)
{
  const size_t i = workspace->i;
  double * alist = workspace->alist;
  double * blist = workspace->blist;
  double * rlist = workspace->rlist;
  double * elist = workspace->elist;

  *a = alist[i] ;
  *b = blist[i] ;
  *r = rlist[i] ;
  *e = elist[i] ;
}

static inline double
sum_results (const gsl_integration_workspace * workspace);

static inline double
sum_results (const gsl_integration_workspace * workspace)
{
  const double * const rlist = workspace->rlist ;
  const size_t n = workspace->size;

  size_t k;
  double result_sum = 0;

  for (k = 0; k < n; k++)
    {
      result_sum += rlist[k];
    }
  
  return result_sum;
}

static inline int
subinterval_too_small (double a1, double a2, double b2);

static inline int
subinterval_too_small (double a1, double a2, double b2)
{
  const double e = GSL_DBL_EPSILON;
  const double u = GSL_DBL_MIN;

  double tmp = (1 + 100 * e) * (fabs (a2) + 1000 * u);

  int status = fabs (a1) <= tmp && fabs (b2) <= tmp;

  return status;
}

//CONTENTS OF QPSRT.C
static inline void 
qpsrt (gsl_integration_workspace * workspace);

static inline
void qpsrt (gsl_integration_workspace * workspace)
{
  const size_t last = workspace->size - 1;
  const size_t limit = workspace->limit;

  double * elist = workspace->elist;
  size_t * order = workspace->order;

  double errmax ;
  double errmin ;
  int i, k, top;

  size_t i_nrmax = workspace->nrmax;
  size_t i_maxerr = order[i_nrmax] ;
  
  /* Check whether the list contains more than two error estimates */

  if (last < 2) 
    {
      order[0] = 0 ;
      order[1] = 1 ;
      workspace->i = i_maxerr ;
      return ;
    }

  errmax = elist[i_maxerr] ;

  /* This part of the routine is only executed if, due to a difficult
     integrand, subdivision increased the error estimate. In the normal
     case the insert procedure should start after the nrmax-th largest
     error estimate. */

  while (i_nrmax > 0 && errmax > elist[order[i_nrmax - 1]]) 
    {
      order[i_nrmax] = order[i_nrmax - 1] ;
      i_nrmax-- ;
    } 

  /* Compute the number of elements in the list to be maintained in
     descending order. This number depends on the number of
     subdivisions still allowed. */
  
  if(last < (limit/2 + 2)) 
    {
      top = last ;
    }
  else
    {
      top = limit - last + 1;
    }
  
  /* Insert errmax by traversing the list top-down, starting
     comparison from the element elist(order(i_nrmax+1)). */
  
  i = i_nrmax + 1 ;
  
  /* The order of the tests in the following line is important to
     prevent a segmentation fault */

  while (i < top && errmax < elist[order[i]])
    {
      order[i-1] = order[i] ;
      i++ ;
    }
  
  order[i-1] = i_maxerr ;
  
  /* Insert errmin by traversing the list bottom-up */
  
  errmin = elist[last] ;
  
  k = top - 1 ;
  
  while (k > i - 2 && errmin >= elist[order[k]])
    {
      order[k+1] = order[k] ;
      k-- ;
    }
  
  order[k+1] = last ;

  /* Set i_max and e_max */

  i_maxerr = order[i_nrmax] ;
  
  workspace->i = i_maxerr ;
  workspace->nrmax = i_nrmax ;
}
#endif /* __GSL_INTEGRATION_H__ */
