#include <R.h>
#include <Rmath.h>
#include <math.h>

// Compute residual error differences
double get_rikjl(double *X, double *sigma,
		 int *N, int *p, int ik_idx, int jl_idx) {
  double *xdif = Calloc(*p, double);
  double rikjl = 0.0;
  int m = 0, q = 0;
  for (m = 0; m < *p; m++) {
    xdif[m] = 0.0;
    xdif[m] = X[ik_idx + m * *N] - X[jl_idx + m * *N];
  }
  for (m = 0; m < *p; m++) {
    for (q = 0; q < *p; q++) {
      rikjl += xdif[m] * sigma[m * *p + q] * xdif[q];
    }
  }
  rikjl = sqrt(rikjl);
  Free(xdif);
  return(rikjl);
}

// Gehan type objective function (smooth); old name = lfun
void gehan_s_obj(double *beta, double *Y, double *X, double *delta,
		 int *clsize, double *sigma, int *n,
	  int *p, int *N, double *weights, double *gehanweights, double *ln) {
  int i, j, k, l, ik_idx = 0, jl_idx;
  double rikjl, edif, z, H, h, sqrtn = sqrt(*n);
  double *e = Calloc(*N, double);
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }
  *ln = 0.0;
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      if (delta[ik_idx] != 0) {
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    rikjl = get_rikjl(X, sigma, N, p, ik_idx, jl_idx);
	    if (rikjl != 0) {
	      edif = e[jl_idx] - e[ik_idx];
	      z = sqrtn * edif / rikjl;
	      H = pnorm(z, 0.0, 1.0, 1, 0);
	      h = dnorm(z, 0.0, 1.0, 0);
	      *ln += gehanweights[ik_idx] * weights[ik_idx] * weights[jl_idx] *
		(edif * H + rikjl * h / sqrtn);
	    }
	    jl_idx++;
	  }
	}
      }
      ik_idx++;
    }
  }
  Free(e);
  /* *ln /= (*n * (*n - 1)); */
}

// Gehan type estimating function (non-smooth) old name = unsfun
void gehan_ns_est(double *beta, double *Y, double *X, double *delta, int *clsize,
	    int *n, int *p, int *N, double *weights, double *gehanweights, double *sn) {
  int i, j, k, l, ik_idx = 0, jl_idx, r;
  double *e = Calloc(*N, double), *xdif = Calloc(*p, double);
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      if (delta[ik_idx] != 0) {
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    if (e[ik_idx] - e[jl_idx] <= 0) {
	      for (r = 0; r < *p; r++) {
		xdif[r] = 0.0;
		xdif[r] = X[ik_idx + r * *N] - X[jl_idx + r * *N];
		sn[r] += gehanweights[ik_idx] * weights[ik_idx] * weights[jl_idx] * xdif[r];
	      } // end e[ik] - e[jl]
	    }
	    jl_idx++;
	  }
	}
      }
      ik_idx++;
    }
  }
  Free(xdif);
  Free(e);
}

// Gehan type estimating function (smooth); old name = ufun
void gehan_s_est(double *beta, double *Y, double *X, double *delta, int *clsize,
	  double *sigma, int *n, int *p, int *N, double *weights, double *gehanweights, double *sn) {
  int i, j, k, l, ik_idx = 0, jl_idx, r;
  double *e = Calloc(*N, double), *xdif = Calloc(*p, double);
  double rikjl, z, H, sqrtn = sqrt(*n), edif;
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      if (delta[ik_idx] != 0) {
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    rikjl = get_rikjl(X, sigma, N, p, ik_idx, jl_idx);
	    if (rikjl != 0) {
	      edif = e[jl_idx] - e[ik_idx];
	      z = sqrtn * edif / rikjl;
	      H = pnorm(z, 0.0, 1.0, 1, 0);
	      for (r = 0; r < *p; r++) {
		xdif[r] = 0.0;
		xdif[r] = X[ik_idx + r * *N] - X[jl_idx + r * *N];
		sn[r] += gehanweights[ik_idx] * weights[ik_idx] * weights[jl_idx] * xdif[r] * H;
	      } // end for r
	    }
	    jl_idx++;
	  }
	}
      }
      ik_idx++;
    }
  }
  Free(xdif);
  Free(e);
}

// log-rank type estimating function (non-smooth); old name = ulognsfun
void log_ns_est(double *beta, double *Y, double *X, double *delta, int *clsize,
		int *n, int *p, int *N, double *weights, double *gw, double *sn) {
  int i, j, k, l, ik_idx = 0, jl_idx, r;
  double *e = Calloc(*N, double), *nu = Calloc(*p, double);
  double de = 0.0;
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      if (delta[ik_idx] != 0) {
	/* reset nu */
	for ( r = 0; r < *p; r++) {
	  nu[r] = 0.0;
	}
	de = 0.0;
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    if (e[ik_idx] - e[jl_idx] <= 0) {
	      for ( r = 0; r < *p; r++) {
		nu[r] += X[jl_idx + r * *N] * weights[jl_idx];
	      }
	      de += weights[jl_idx];
	    }
	    jl_idx++;
	  }
	}  // end jl
	for (r =  0; r < *p; r++) {
	  sn[r] += weights[ik_idx] * gw[ik_idx] * (X[ik_idx + r * *N] - nu[r] / de);
	}
      }
      ik_idx++;
    }
  }
  Free(nu);
  Free(e);
}

// log-rank type estimating function (smooth-equivent form); old name = ulogfun
void log_s_est(double *beta, double *Y, double *X, double *delta, int *clsize,
	     double *sigma, int *n, int *p, int *N, double *weights, double *gw, double *sn) {
  int i, j, k, l, ik_idx = 0, jl_idx, r;
  double *e = Calloc(*N, double), *nu = Calloc(*p, double);
  double rikjl, z, H, sqrtn = sqrt(*n), edif, de;
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      if (delta[ik_idx] != 0) {
	/* reset nu */
	for ( r = 0; r < *p; r++) {
	  nu[r] = 0.0;
	}
	de = 0.0;
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    rikjl = get_rikjl(X, sigma, N, p, ik_idx, jl_idx);
	    if (rikjl != 0) {
	      edif = e[jl_idx] - e[ik_idx];
	      z = sqrtn * edif / rikjl;
	      H = pnorm(z, 0.0, 1.0, 1, 0);
	      for (r = 0; r < *p; r++) {
		nu[r] += weights[jl_idx] * X[jl_idx + r * *N] * H;
	      }
	      de += weights[jl_idx] * H;
	    }
	    jl_idx++;
	  }
	}
	for ( r = 0; r < *p; r ++) {
	  sn[r] += weights[ik_idx] * gw[ik_idx] * (X[ik_idx + r * *N] - nu[r] / de);
	}
      }
      ik_idx++;
    }
  }
  //  Free(edif);
  Free(nu);
  //  Free(de);
  Free(e);
  /* for (r = 0; r < *p; r++) { */
  /*   sn[r] /= (*n * (*n - 1)); */
  /* } */
}

// Compute non-smooth gehan weight; used to prepare for method #3 and #4; old name = getnsgehan
void gehan_ns_wt(double *beta, double *Y, double *X, int *clsize,
	     int *n, int *p, int *N, double *weights, double *de) {
  int i, j, k, l, ik_idx = 0, jl_idx = 0;
  double *e = Calloc(*N, double);
  double edif;
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    edif = e[jl_idx] - e[ik_idx];
	    if (edif >= 0) {
	      de[ik_idx] += weights[jl_idx];
	    }
	    jl_idx++;
	  } // end for l
	} // end for j
	ik_idx++;
    } // end for k
  } // end for i
  Free(e);
}

// Compute smooth gehan weight; used to prepare for method #3 and #4; old name = getgehan
void gehan_s_wt(double *beta, double *Y, double *X, int *clsize,
	       double *sigma, int *n, int *p, int *N, double *weights, double *de) {
  int i, j, k, l, ik_idx = 0, jl_idx;
  double *e = Calloc(*N, double);
  double rikjl, edif, H, z, sqrtn = sqrt(*n);
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    rikjl = get_rikjl(X, sigma, N, p, ik_idx, jl_idx);
	    if (rikjl != 0) {
	      edif = e[jl_idx] - e[ik_idx];
	      z = sqrtn * edif / rikjl;
	      H = pnorm(z, 0.0, 1.0, 1, 0);
	      de[ik_idx] += weights[jl_idx] * H;
	    } // end if rikjl != 0
	    jl_idx++;
	  } // end for l
	} // end for j
      ik_idx++;
    } // end for k
  } // end for i
  Free(e);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// March 15
////////////////////////////////////////////////////////////////////////////////////////////////////

void abargehanfunC(double *beta, double *Y, double *X, 
		  double *delta, int *clsize, double *sigma,
		  int *n, int *p, int *N,
		  double *weights,
		  double *gehanWeights,
		  //output
		  double *abar) {
  int i, j, k, l, ik_idx = 0, jl_idx, r, s, xdif_idx;
  double *e = Calloc(*N, double), *xdif = Calloc(*p, double),
    *nu1 = Calloc(*p, double), *nu2 = Calloc(*p, double), *nu3 = Calloc(*p, double);
  double rikjl, edif, h, z, sqrtn = sqrt(*n), coef;

  /* compute the e vector */
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }

  /*The main abar part */
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      if (delta[ik_idx] != 0) {
	/* reset */
	for ( r = 0; r < *p; r++) {
	  nu1[r] = 0.0;
	  nu2[r] = 0.0;
	  nu3[r] = 0.0;
	}
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    rikjl = get_rikjl(X, sigma, N, p, ik_idx, jl_idx);
	    if (rikjl != 0) {
	      edif = e[jl_idx] - e[ik_idx];
	      z = sqrtn * edif / rikjl;
	      h = dnorm(z, 0.0, 1.0, 0);
	      coef = gehanWeights[ik_idx] * weights[ik_idx] * weights[jl_idx] * h * sqrtn / rikjl;
	      for (r = 0; r < *p; r++) {
		xdif[r] = 0.0;
		xdif[r] = X[ik_idx + r * *N] - X[jl_idx + r * *N];
	      }
	      xdif_idx = 0;
	      for (r = 0; r < *p; r++) {
		for (s = 0; s < *p; s++) {
		  abar[xdif_idx] += coef * xdif[r] * xdif[s];  // p by p matrix
		  xdif_idx++;
		} // end for s
	      } // end for r
	    } // end if rikjl != 0
	    jl_idx++;
	  } // end for l
	} // end for j
      } // end if delta[ik_idx] != 0
      ik_idx++;
    } // end for k
  } // end for i

  Free(xdif);
  Free(e);
  Free(nu1);
  Free(nu2);
  Free(nu3);
  /* for (r = 0; r < (*p * *p); r++) { */
  /*   abar[r] /= (*n * (*n - 1)); */
  /* } */
}


void abarlogfunC(double *beta, double *Y, double *X, double *delta, int *clsize, double *pw,
	     double *sigma,
	     int *n,
	     int *p,
	     int *N,
	     double *weights,
		/* output */
	     double *abar) {
  int i, j, k, l, ik_idx = 0, jl_idx, r, s, xdif_idx, xdif_idx2;
  double *e = Calloc(*N, double), *xdif = Calloc(*p, double),
    *nu1 = Calloc(*p * *p, double), *nu2 = Calloc(*p, double), *nu3 = Calloc(*p, double);
  double rikjl, edif, H, h, z, sqrtn = sqrt(*n), coef, de;

  /* compute the e vector */
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }

  /*The main abar part */
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      if (delta[ik_idx] != 0) {
	/* reset */
	xdif_idx = 0;
	for ( r = 0; r < *p; r++) {
	  for ( s = 0; s < *p; s++) {
	    nu1[xdif_idx] = 0.0;
	    xdif_idx++;
	  }
	  nu2[r] = 0.0;
	  nu3[r] = 0.0;
	}
	de = 0.0;
	coef = 0.0;
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    rikjl = get_rikjl(X, sigma, N, p, ik_idx, jl_idx);
	    if (rikjl != 0) {
	      edif = e[jl_idx] - e[ik_idx];
	      z = sqrtn * edif / rikjl;
	      h = dnorm(z, 0.0, 1.0, 0);
	      H = pnorm(z, 0.0, 1.0, 1, 0);
	      coef = weights[ik_idx] * weights[jl_idx] * h; // * sqrtn / rikjl;
	      for (r = 0; r < *p; r++) {
		xdif[r] = 0.0;
		xdif[r] = X[ik_idx + r * *N] - X[jl_idx + r * *N];
	      }
	      xdif_idx = 0;
	      for (r = 0; r < *p; r++) {
		for (s = 0; s < *p; s++) {
		  nu1[xdif_idx] += coef * X[jl_idx + s * *N] * xdif[r] * sqrtn / rikjl;  // s, r
		  xdif_idx++;
		}
		nu2[r] += weights[jl_idx] * H * X[jl_idx + r * *N];
		nu3[r] += coef * xdif[r] * sqrtn / rikjl;
	      }
	      de += weights[jl_idx] * H;
	    } // end if rikjl != 0
	    jl_idx++;
	  } // end for l
	} // end for j
	xdif_idx2 = 0;
	for (r = 0; r < *p; r++) {
	  for (s = 0; s < *p; s++) {
	    abar[xdif_idx2] += -1 * pw[ik_idx] * (de * nu1[xdif_idx2] - nu2[s] * nu3[r]) / (de * de);
	    xdif_idx2++;
	  }
	}
      } // end if delta[ik_idx] != 0
      ik_idx++;
    } // end for k
  } // end for i

  Free(xdif);
  Free(e);
  Free(nu1);
  Free(nu2);
  Free(nu3);
  /* for (r = 0; r < (*p * *p); r++) { */
  /*   abar[r] /= (*n * (*n - 1)); */
  /* } */
}


void omegafun(double *beta,
	      double *Y,
	      double *X,
	      double *delta,
	      int *clsize,
	      int *n,
	      int *p,
	      int *N,
	      double *weights,
	      //output
	      double *omega) {
  int i, k, j, l, m, r, s, ik_idx = 0, jl_idx, rs_idx, emk = 0, ind = 0, omega_idx;
  double *xdif = Calloc(*p, double), *e = Calloc(*N, double), *ksi = Calloc(*p * *N, double);

  // compute e
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
   }

  // compute ksi
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      jl_idx = 0;
      for (j = 0; j < *n; j++) {
	for (l = 0; l < clsize[j]; l++) {

	  if (delta[ik_idx] != 0) {
	    if (e[ik_idx] - e[jl_idx] < 0) {
	      for (m = 0; m < *p; m++) {
		ksi[ik_idx + m * *N] += delta[ik_idx] * (X[ik_idx + m * *N] - X[jl_idx + m * *N]) * weights[jl_idx] / *n;  //changed
	      } // end for m
	    } // end if e - e <0
	  } // end if delta[ik_idx] != 0

	  if (delta[jl_idx] != 0) {
	    if (e[ik_idx] - e[jl_idx] >= 0) {
	      rs_idx = 0;
	      emk = 0;
	      for (r = 0; r < *n; r++) {
		for (s = 0; s < clsize[r]; s++) {
		  if (e[rs_idx] - e[jl_idx] >= 0) {
		    for (m = 0; m < *p; m++) {
		      xdif[m] += weights[rs_idx] * ( X[ik_idx + m * *N] - X[rs_idx + m * *N] );  //changed
		    } // end for m
		    emk++;
		  } // end if
		  rs_idx++;
		} // end for s
	      } // end for r

	      for (m = 0; m < *p; m++) {
		ksi[ik_idx + m * *N] -= xdif[m] / (*n * emk);
		xdif[m] = 0;
	      } // end for m
	    } // end if e -e >= 0
	  }// end if delta jl != 0
	  jl_idx++;
	} // end for l
      } // end for j
      ik_idx++;
    } // end for k1
  } // end for i

  for (i = 0; i < *n; i++) {
    ind += clsize[i];
    for (k = 0; k < clsize[i]; k++) {
      for (l = 0; l < clsize[i]; l++) {
	omega_idx = 0;
	for (r = 0; r < *p; r++) {
	  for (s = 0; s < *p; s++) {
	    omega[omega_idx] += ksi[ind + k - clsize[0] + r * *N] * ksi[ind + l - clsize[0] + s * *N];
	    omega_idx++;
	  } // end for r
	} // end for s
      } // end for l
    } // end for k
  } // end for i

  /* for (r = 0; r < *p * *p; r++) { */
  /* omega[r] /= *n; */
  /* } */
  Free(e);
  Free(xdif);
  Free(ksi);
}



void abarpwfunC(double *beta, double *Y, double *X, double *delta, 
	       int *clsize, double *pw, double *sigma,
	       int *n, int *p, int *N,
	       double *weights,
	       //output
	       double *abar) {
  int i, j, k, l, ik_idx = 0, jl_idx, r, s, xdif_idx, xdif_idx2;
  double *e = Calloc(*N, double), *nu = Calloc(*N * *p, double);
  double rikjl, edif, H, z, sqrtn = sqrt(*n), de;

  /* compute the e vector */
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }

  /*The main abar part */
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      if (delta[ik_idx] != 0) {
	/* reset */
	xdif_idx = 0;
	for ( r = 0; r < *p; r++) {
	  nu[xdif_idx] = 0.0;
	  xdif_idx++;
	}
	de = 0.0;
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    rikjl = get_rikjl(X, sigma, N, p, ik_idx, jl_idx);
	    if (rikjl != 0) {
	      edif = e[jl_idx] - e[ik_idx];
	      z = sqrtn * edif / rikjl;
	      H = pnorm(z, 0.0, 1.0, 1, 0);
	      xdif_idx = 0;
	      for (r = 0; r < *p; r++) {
		// nu[jl_idx + r * *N] += weights[jl_idx] * H * X[jl_idx + r * *N];
		nu[r] += weights[jl_idx] * H * X[jl_idx + r * *N];
	      }
	      de += weights[jl_idx] * H;
	    } // end if rikjl != 0
	    jl_idx++;
	  } // end for l
	} // end for j
	xdif_idx2 = 0;
	for (r = 0; r < *p; r++) {
	  for (s = 0; s < *p; s++) {
	    abar[xdif_idx2] += -1 * pw[ik_idx] * X[ik_idx + r * *N] * (X[ik_idx + s * *N] - nu[s] / (de));
	    // abar[xdif_idx2] += pw[ik_idx] * (X[ik_idx + r * *N] - nu[r] / de);
	  xdif_idx2++;
	  } // end s
	}
      } // end if delta[ik_idx] != 0
      ik_idx++;
    } // end for k
  } // end for i
  Free(e);
  Free(nu);
  /* for (r = 0; r < (*p * *p); r++) { */
  /*   abar[r] /= (*n * (*n - 1)); */
  /* } */
}


void getnsgehan (double *beta, double *Y, double *X, int *clsize,
	     double *sigma,
	     int *n,
	     int *p,
	     int *N,
	     double *weights,
		/* output */
	     double *de) {
  int i, j, k, l, ik_idx = 0, jl_idx;
  double *e = Calloc(*N, double);
  double edif;

  /* compute the e vector */
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }

  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    //  rikjl = get_rikjl(X, sigma, N, p, ik_idx, jl_idx);
	    //  if (rikjl != 0) {
	      edif = e[jl_idx] - e[ik_idx];
	      if (edif >= 0) {
		de[ik_idx] += weights[jl_idx];
	      }
	      //  } // end if rikjl != 0
	    jl_idx++;
	  } // end for l
	} // end for j
      ik_idx++;
    } // end for k
  } // end for i
  Free(e);
  /* for (r = 0; r < (*p * *p); r++) { */
  /*   abar[r] /= (*n * (*n - 1)); */
  /* } */
}

void getgehan (double *beta, double *Y, double *X, int *clsize,
	     double *sigma,
	     int *n,
	     int *p,
	     int *N,
	     double *weights,
		/* output */
	     double *de) {
  int i, j, k, l, ik_idx = 0, jl_idx;
  double *e = Calloc(*N, double);
  double rikjl, edif, H, z, sqrtn = sqrt(*n);
  /* compute the e vector */
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }

  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    rikjl = get_rikjl(X, sigma, N, p, ik_idx, jl_idx);
	    if (rikjl != 0) {
	      edif = e[jl_idx] - e[ik_idx];
	      z = sqrtn * edif / rikjl;
	      H = pnorm(z, 0.0, 1.0, 1, 0);
	      de[ik_idx] += weights[jl_idx] * H;
	    } // end if rikjl != 0
	    jl_idx++;
	  } // end for l
	} // end for j
      ik_idx++;
    } // end for k
  } // end for i
  Free(e);
  /* for (r = 0; r < (*p * *p); r++) { */
  /*   abar[r] /= (*n * (*n - 1)); */
  /* } */
}


void matAfun(double *X, double *sigma, double *delta,
	  int *clsize, int *n, int *p, int *N, 
	  double *matA) {
  int i, j, k, ikId = 0, jlId, l, r, m, q;
  double *xdif = Calloc(*p, double);
  double rikjl = 0.0;
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      if (delta[ikId] != 0) {
	jlId = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {

	    for (r = 0; r < *p; r++) {
	      xdif[r] = 0.0;
	      xdif[r] = X[ikId + r * *N] - X[jlId + r * *N];
	    }

	    rikjl = get_rikjl(X, sigma, N, p, ikId, jlId);
	    if (rikjl != 0) {
	      for (m = 0; m < *p; m++) {
		for (q = 0; q < *p; q++) {
		  matA[q + m * *p] += xdif[m] * xdif[q] * sqrt(*n)/ (rikjl);
		}
	      }
	    }
	    jlId++;
	  } // end l
	} // end j
      } // end delta[ikId] != 0
      ikId++;
    } // end k 
  } // end i
  Free(xdif);
}


// Length biase

void ulblk(double*beta, double *Y, double *X, double *delta, double *a, double *lbw,
	   int *n, int *p, 
	   //output 
	   double *out) {
  int i, j, r;
  double *e = Calloc(*n, double); 
  double *ei = Calloc(*n, double); 
  double *ea = Calloc(*n, double); 
  double *nu = Calloc(*p, double); 
  double de;

  for (i = 0; i < *n; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *n + i] * beta[j];
    }
    ei[i] = Y[i] - e[i];
    ea[i] = a[i] - e[i];
  }

  for (i = 0; i < *n; i++) {
    if (delta[i] != 0) {
      for ( r = 0; r < *p; r++) {
	nu[r] = 0.0;
      }
      de = 0.0;
      for (j = 0; j < *n; j++) {
	if (ei[i] <= ei[j] && ei[i] > ea[j]) {
	  for (r = 0; r < *p; r++) {
	    nu[r] += lbw[j + i * *n] * X[j + r * *n];
	}
	  de += lbw[j + i * *n];
	}
      } // end j
      for (r = 0; r < *p; r++) {
	out[r] += X[i + r * *n] - (nu[r] / de);
      }
    } // end delta
  } // end i
  Free(e);
  Free(ea);
  Free(ei);
  Free(nu);
}


void ulblk2(double*beta, double *Y, double *X, double *delta, double *lbw,
	   int *n, int *p, 
	   //output 
	   double *out) {
  int i, j, r;
  double *e = Calloc(*n, double); 
  double *ei = Calloc(*n, double); 
  double *nu = Calloc(*p, double); 
  double de;
  for (i = 0; i < *n; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *n + i] * beta[j];
    }
    ei[i] = Y[i] - e[i];
  }

  for (i = 0; i < *n; i++) {
    if (delta[i] != 0) {
      for ( r = 0; r < *p; r++) {
	nu[r] = 0.0;
      }
      de = 0.0;
      for (j = 0; j < *n; j++) {
	if (ei[i] <= ei[j]) {
	  for (r = 0; r < *p; r++) {
	    nu[r] += lbw[j + i * *n] * X[j + r * *n];
	}
	  de += lbw[j + i * *n];
	}
      } // end j
      for (r = 0; r < *p; r++) {
	out[r] += X[i + r * *n] - (nu[r] / de);
      }
    } // end delta
  } // end i
  Free(e);
  Free(ei);
  Free(nu);
}

void ulbge(double*beta, double *Y, double *X, double *delta, double *a, int *n, int *p, 
	   //output 
	   double *out) {
  int i, j, r;
  double *e = Calloc(*n, double); 
  double *ei = Calloc(*n, double); 
  double *ea = Calloc(*n, double); 

  for (i = 0; i < *n; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *n + i] * beta[j];
    }
    ei[i] = Y[i] - e[i];
    ea[i] = a[i] - e[i];
  }

  for (i = 0; i < *n; i++) {
    if (delta[i] != 0) {    
      for (j = 0; j < *n; j++) {
    	if (ei[i] <= ei[j] && ei[i] > ea[j]) {
	  for (r = 0; r < *p; r++) {
	    out[r] += X[i + r * *n] - X[j + r * *n];
	  }
	}
      }
    }
  }
  Free(e);
  Free(ea);
  Free(ei);
}

void ulbge2(double*beta, double *Y, double *X, double *delta, double *lbw, 
	    int *n, int *p,
	   //output 
	   double *out) {
  int i, j, r;
  double *e = Calloc(*n, double); 
  double *ei = Calloc(*n, double); 

  for (i = 0; i < *n; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *n + i] * beta[j];
    }
    ei[i] = Y[i] - e[i];
  }

  for (i = 0; i < *n; i++) {
    if (delta[i] != 0) {    
      for (j = 0; j < *n; j++) {
    	if (ei[i] <= ei[j]) {
	  for (r = 0; r < *p; r++) {
	    out[r] += (X[i + r * *n] - X[j + r * *n]) * lbw[j + i * *n];
	  }
	}
      }
    }
  }
  Free(e);
  Free(ei);
}

void ulbsge(double*beta, double *Y, double *X, 
	    double *delta, double *a, double *sigma, int *n, int *p, 
	   //output 
	   double *out) {
  int i, j, r;
  double *e = Calloc(*n, double); 
  double *ei = Calloc(*n, double); 
  double *ea = Calloc(*n, double); 
  double z1, z2, H1, H2, rikjl;
  for (i = 0; i < *n; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *n + i] * beta[j];
    }
    ei[i] = Y[i] - e[i];
    ea[i] = a[i] - e[i];
  }

  for (i = 0; i < *n; i++) {
    if (delta[i] != 0) {    
      for (j = 0; j < *n; j++) {
	rikjl = get_rikjl(X, sigma, n, p, i, j);
	if (rikjl != 0) {
	  z1 = sqrt(*n) * (ea[j] - ei[i]) / rikjl;
	  z2 = sqrt(*n) * (ei[j] - ei[i]) / rikjl;
	  H1 = pnorm(z1, 0.0, 1.0, 1, 0);
	  H2 = pnorm(z2, 0.0, 1.0, 1, 0);
	  for (r = 0; r < *p; r++) {
	    out[r] += (X[i + r * *n] - X[j + r * *n]) * (H2 - H1);
	  }
	}
      }
    }
  }
  Free(e);
  Free(ea);
  Free(ei);
}


void uuge(double*beta, double *Y, double *X, double *delta, double *lbw, double *sigma,
	  int *n, int *p,
	   //output 
	   double *out) {
  int i, j, r;
  double z, H, rikjl;
  double *e = Calloc(*n, double); 
  double *ei = Calloc(*n, double); 

  for (i = 0; i < *n; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *n + i] * beta[j];
    }
    ei[i] = Y[i] - e[i];
  }

  for (i = 0; i < *n; i++) {
    if (delta[i] != 0) {    
      for (j = 0; j < *n; j++) {
	rikjl = get_rikjl(X, sigma, n, p, i, j);
	if (rikjl != 0) {
	  z = sqrt(*n) * (ei[j] - ei[i]) / rikjl;
	  H = pnorm(z, 0.0, 1.0, 1, 0);
	  for (r = 0; r < *p; r++) {
	    out[r] += (X[i + r * *n] - X[j + r * *n]) * H * lbw[j + i * *n];
	  }
	}
      }
    }
    
  }
  Free(e);
  Free(ei);
}
