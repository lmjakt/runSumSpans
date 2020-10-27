#include <R.h>
#include <Rinternals.h>
#include <vector>

// This code is a generalised rewrite of functions within cpg_distribute
// written by Gos Micklem + Tim Cutts (2004)
// but taking a sequence of scores and their associated positions,
// and a separation penalty.

// This code uses the .Call interface.

// Originally written to search for differentially methylated regions in RRBS data
// which is not continuous requiring that locations are specified.

// values are taken at face value.. 

// The positions should be sorted within regions; that is:
//  1, 12, ..., 234, 23, 45, ...
// is allowed. Discontinuities are taken to indicate borders
// which spans cannot cross.

extern "C" SEXP find_spans(SEXP, SEXP, SEXP);

// Functions to find spans using a running sum..
struct span {
  span(unsigned int start, unsigned int end, double pos_start, double pos_end, double score) : 
    start(start), end(end), pos_start(pos_start), pos_end(pos_end), score(score){}
  span() :
    start(0), end(0), pos_start(0), pos_end(0), score(0) {}
  unsigned int start;
  unsigned int end;
  double pos_start;
  double pos_end;
  double score;
};


// The separation penalty is multiplied by the distance between adjacent
// values.
void find_span(int start, int end, double *values, double *positions,
	       int l, double separation_penalty, 
	       std::vector<span>& spans){
  if(start > end)
    return;
  if(!l || end >= l)
    return;

  //  Rprintf("find_span: %d -> %d (%d), %f\n", start, end, l, separation_penalty);

  double score = 0;
  double last_score = 0;   // this keeps track of whether we have entered a span region
  double maxScore = 0;
  int maxI = start;
  int minI = start;
  int i = start;
  
  while(i < end){
    last_score = score;
    // we do not allow a score to go across a chromosomal boundary.
    // when the positions are appropriately sorted this will show as a negative distances
    if(i > 0 && positions[i] < positions[i-1]){
      // we are crossing a chromosomal boundary. We don't want to do that.
      score = 0;
    }else{
      if(!score){
	score += values[i];
      }else{
	score += values[i] - separation_penalty * (positions[i] - positions[i-1]);
      }
    }
    score = score < 0 ? 0 : score; // negative score not allowed.
    if(last_score == 0 && score > 0)
      minI = i;
    // if last_score is positive and score is negative, then we are at the
    // end of a span and we should do something useful.
    if(score == 0 && last_score > 0){
      spans.push_back( span(minI, maxI, positions[minI], positions[maxI], maxScore) );
      //      Rprintf("pushed back: %d -> %d %f  : %f -> %f\n", minI, maxI, maxScore, positions[minI], positions[maxI]);
      find_span(maxI+1, i, values, positions, l, separation_penalty, spans);
      score = last_score = maxScore = 0;
      maxI = minI = i;
    }
    maxI = (score > maxScore) ? i : maxI;
    maxScore = (score >= maxScore) ? score : maxScore;
    minI = (score == 0) ? i : minI;
    ++i;
  }
  // Unfortunately the above code doesn't handle the last span
  if(score > 0){
    spans.push_back( span(minI, maxI, positions[minI], positions[maxI], maxScore) );
    //    Rprintf("pushed back: %d -> %d %f  : %f -> %f\n", minI, maxI, maxScore, positions[minI], positions[maxI]);
    find_span(maxI+1, end, values, positions, l, separation_penalty, spans);
  }
}

// and an interface function



SEXP find_spans( SEXP scores_r, SEXP pos_r, SEXP sep_p_r ){
  if( TYPEOF(scores_r) != REALSXP || length( scores_r ) < 2 )
    return(R_NilValue);
  if(TYPEOF(pos_r) != REALSXP || length( scores_r ) != length( pos_r ))
    error("positions should be real and of the same length as the scores");
  if(TYPEOF(sep_p_r) != REALSXP || length( sep_p_r ) != 1 )
    error("sep_p_r should be real and of length 1");
  
  int l = length(scores_r);
  double *scores = REAL( scores_r );
  double *pos = REAL( pos_r );
  double sep_p = REAL( sep_p_r )[0];
  std::vector<span> spans;
  find_span( 0, l-1, scores, pos, l, sep_p, spans );
  if(spans.size() == 0)
    return(R_NilValue);
  // make a data frame with 6 columns to which we copy the values
  // note that for this it would be better to change the span struct
  // to be a set of vectors and to have a push() function that grows
  // each vector separately. (Well maybe, that increases the
  // number of memory allocations, so may not be better).
  SEXP ret_value = PROTECT( allocVector( VECSXP, 5));
  SEXP class_r = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT( class_r, 0, mkChar("data.frame"));
  
  const char *colnames[] = {"i1", "i2", "start", "end", "score"};
  SEXP colnames_r = PROTECT(allocVector(STRSXP, 5));
  for(int i=0; i < 5; ++i)
    SET_STRING_ELT(colnames_r, i, mkChar( colnames[i] ));
  namesgets( ret_value, colnames_r);
  
  SET_VECTOR_ELT(ret_value, 0, allocVector(INTSXP, spans.size()));
  SET_VECTOR_ELT(ret_value, 1, allocVector(INTSXP, spans.size()));
  SET_VECTOR_ELT(ret_value, 2, allocVector(REALSXP, spans.size()));
  SET_VECTOR_ELT(ret_value, 3, allocVector(REALSXP, spans.size()));
  SET_VECTOR_ELT(ret_value, 4, allocVector(REALSXP, spans.size()));
  
  // And set up some pointers:
  int *start = INTEGER( VECTOR_ELT(ret_value, 0));
  int *end = INTEGER( VECTOR_ELT(ret_value, 1));
  double *pos_start = REAL( VECTOR_ELT(ret_value, 2));
  double *pos_end = REAL( VECTOR_ELT(ret_value, 3));
  double *score = REAL( VECTOR_ELT(ret_value, 4));
  for(size_t i=0; i < spans.size(); ++i){
    start[i] = (int)spans[i].start + 1;
    end[i] = (int)spans[i].end + 1;
    pos_start[i] = spans[i].pos_start;
    pos_end[i] = spans[i].pos_end;
    score[i] = spans[i].score;
  }
  // Setting the class attribute to 
  //    classgets( ret_value, class_r ); // shortcut for setAttrib
  // setAttrib( ret_value, R_ClassSymbol, class_r);
  // both of these result in an empty data.frame. This really is easier
  // in Rcpp. Easy to wrap in R code though.
  
  UNPROTECT(3);
  return( ret_value );
}

// NOTE: R_init_runSumSpans does not appear to be called when
// loading the library. This makes using this rather dangerous.

static const R_CallMethodDef callMethods[] = {
  {"rs_spans", (DL_FUNC)&find_spans, 3},
  {NULL, NULL, 0}
};

void R_init_runSumSpans(DllInfo *info)
{
  Rprintf("Trying to register dllinfo \n");
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}



// [[Rcpp::export]]
// List find_spans(NumericVector scores, NumericVector positions, float seperation_penalty){
//   std::vector<float> s = as<std::vector<float> >(scores);
//   std::vector<int> p = as<std::vector<int> >(positions);
//   std::vector<span> spans;
  
//   find_span(0, s.size() - 1, s, p, spans, seperation_penalty);
  
//   std::vector<unsigned int> starts(spans.size());
//   std::vector<unsigned int> ends(spans.size());
//   std::vector<int> pos_start(spans.size());
//   std::vector<int> pos_end(spans.size());
//   std::vector<float> span_scores(spans.size());

//   for(unsigned int i=0; i < spans.size(); ++i){
//     // starts and ends are indices. In R these count from 0, not from 1.
//     starts[i] = spans[i].start + 1;
//     ends[i] = spans[i].end + 1;
//     pos_start[i] = spans[i].pos_start;
//     pos_end[i] = spans[i].pos_end;
//     span_scores[i] = spans[i].score;
//   }
//   return(List::create(Named("start", starts),
// 		      Named("end", ends),
// 		      Named("pos_start", pos_start),
// 		      Named("pos_end", pos_end),
// 		      Named("score", span_scores)
// 		      ));
 // }
